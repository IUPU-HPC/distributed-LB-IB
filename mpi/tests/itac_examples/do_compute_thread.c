#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>

#include "semaphore.h"
#include "lv.h"
#include "tasknode.h"
#include "do_compute_thread.h"
#include "bins.h"
#include "tblas.h"

//#include <VT.h>

//#define V_T

double ONE = 1.0, NONE = -1.0;
char Lower='L', Transpose='T', NoTranspose='N', Forward='F', Columnwise='C', Left='L', NonUnit='N', Right='R', Up='U', NormM='M';


void* do_compute_thread(void *v) {
  GV               gv;
  LV               lv;
  int              thrd_id, pid;
  Matrix*          matrices;
  int              NB, nblks;
  ReadyNode        task;
  double           t0, t1, t2;
  double           *tau, *work; //temporary working space
  
  t1         = 0;
  t2         = 0;
  lv         = (LV) v;
  gv         = (GV) lv->gv;
  pid        = gv->proc_id;
  thrd_id    = lv->thrd_id;
  matrices   = gv->matrices;
  NB         = gv->NB;
  nblks      = gv->nblks;

  //pin_me_now(thrd_id);

  //printf("(P%d, t%d) computing thread starts ...\n", pid, thrd_id);

  tau   = (double*) malloc(sizeof(double)*NB);
  work  = (double*) malloc(sizeof(double)*NB*NB);

  //hist = 0;
  while(!gv->done) {
    if((task = get_task(lv)) != NULL) {//task is already removed from ready_q
      t0 = get_cur_time();
      solve_task(task, matrices, nblks, NB, gv, tau, work); //compute it
      t1 += get_cur_time() - t0;
      /*if(task->task->type == TASK_TYPE_TTQRT) {//progress report
	printf("(P%d,t%d) computed a task: ", pid, thrd_id);
	print_task(task->task);
      }*/
      t0 = get_cur_time();
      finish_task(lv, task);//update different queues and fire new tasks
      t2 += get_cur_time() - t0;
      lv->solved_tasks++;
      //if(pid == 0)
      //printf("(P%d,t%d):R[0,0]=%p\n", pid, thrd_id, matrices[1]->blk_ptrs[0]);
    } 
  } //end of while loop

  //added because there are a few leftover StoreR tasks!
  while( (task=get_task(lv)) != NULL ) {
    //printf("trailing tasks after Exit:\n");
    //print_task(task->task);
    t0 = get_cur_time();
    solve_task(task, matrices, nblks, NB, gv, tau, work); //compute it
    t1 += get_cur_time() - t0;
#ifdef PLASMA_DEBUG
    printf("solved a trailing storeR task\n");
#endif
  }
  
  lv->comp_time = t1;
  lv->fire_time = t2;

  free(tau);
  free(work);

  return NULL;
}

//called by router when you need router to do Computations!! not for now
void do_compute_one_task(void *v, double* tau, double* work) {
  GV               gv;
  LV               lv;
  Matrix*          matrices;
  int              NB, nblks;
  ReadyNode        task;
  double           t0;
  
  lv         = (LV) v;
  gv         = (GV) lv->gv;
  matrices   = gv->matrices;
  NB         = gv->NB;
  nblks      = gv->nblks;

  if((task = get_task(lv)) != NULL) {        //task is already removed from ready_q
    t0 = get_cur_time();
    solve_task(task, matrices, nblks, NB, gv, tau, work); //compute it
    lv->comp_time += get_cur_time() - t0;
    /*if(task->task->type == TASK_TYPE_GEQRT) {//progress report
      printf("(P%d,t%d) computed a task: ", gv->proc_id, lv->thrd_id);
      print_task(task->task);
    }*/
    t0 = get_cur_time();
    finish_task(lv, task);              //update different queues and fire new tasks
    lv->fire_time += get_cur_time() - t0;
    lv->solved_tasks++;
  }

}

//Only called by router: convert the alien output to output-msg! not for now
void do_compute_one_task2(void *v, double* tau, double* work) {
  GV               gv;
  LV               lv;
  ReadyNode        task;
  
  lv         = (LV) v;
  gv         = (GV) lv->gv;

  if((task = get_task2(lv)) != NULL) {//fetch the alien messages -> critic_q
    //if(task->task->type == TASK_TYPE_NOOP) {
    finish_task(lv, task);          //critic_q -> out-mbox
      //}
  }

}

//if found a ready2go task, return 1
ReadyNode get_task(LV lv) {
  ReadyNode               ready;
  GV                      gv;
  MailBox                 inbox; //all about owner tasks
  MailBox                 emergent_inbox;
  GeneratedTaskQueue      generated_q;
  ReadyTaskQueue          critic_q;
  ReadyTaskQueue          regular_q;
  ReadyTaskQueue          emergent_q;
  Msg*                    buf, msg;
  int                     num_removed, i;
  double                  t0;

  gv             = lv->gv;
  buf            = lv->msg_scratch_buf;
  inbox          = gv->in_box;
  critic_q       = gv->critic_ready_q;
  regular_q      = gv->regular_ready_q;
  emergent_q     = gv->emergent_ready_q;
  generated_q    = gv->gen_task_q;
  emergent_inbox = gv->emergent_in_box;

  t0 = get_cur_time();

  num_removed = rm_msgs_inbox(emergent_inbox, buf);
  if(num_removed > 0) {
#ifdef PLASMA_DEBUG_DEEP
    printf("(P%d,t%d) fetches %d msg from in-box\n", gv->proc_id, lv->thrd_id, num_removed);
    print_in_msg(gv->proc_id, buf[0]);
#endif
    if(num_removed > lv->msg_scratch_buf_size) {fprintf(stderr, "Error: msg_scratch_buf overflow!\n"); exit(21);}
    for(i = 0; i < num_removed; i++) {
      msg = buf[i];
      process_recv_msg(msg, lv); 
      free_msg(gv->msg_pool, msg);//buf[i] stores a pointer to msg which was created by router.
    }
  
  } else {  
    //read in-mailbox, move the msgs out of the in-box. That is, an array of pointers
    num_removed = rm_msgs_inbox(inbox, buf);
    if(num_removed > 0) {
#ifdef PLASMA_DEBUG
      printf("(P%d,t%d) fetches %d msg from in-box\n", gv->proc_id, lv->thrd_id, num_removed);
      //print_in_msg(gv->proc_id, buf[0]);
#endif
      if(num_removed > lv->msg_scratch_buf_size) {fprintf(stderr, "Error: msg_scratch_buf overflow!\n"); exit(22);}
      for(i = 0; i < num_removed; i++) {
	msg = buf[i];
	process_recv_msg(msg, lv); //either inputs(for owners) OR alien complementary outpu
	free_msg(gv->msg_pool, msg);//buf[i] stores a pointer to msg which was created by router.
      }
    }
  }
  
  //check ready queues
  ready = popfront_ready_q(emergent_q);
  if(ready == NULL) {
    ready = popfront_ready_q(critic_q);
    if(ready == NULL) {
      ready = popfront_ready_q(regular_q);    
    } 
  }
  
#ifdef PLASMA_DEBUG
  if(ready != NULL) {
    printf("(P%d, t%d) got a task: ", gv->proc_id, lv->thrd_id);
    if(ready->task == NULL) {
      printf("Error: (P%d, t%d) just got a ready->NULL task\n", gv->proc_id, lv->thrd_id);
      exit(23);
    }    
    print_task(ready->task);
  }
#endif  

  lv->getask_time += get_cur_time() - t0;

  //printf("P%d: hello, %p\n", gv->proc_id, ready);
  return ready;
}


//only search Emergent In-Box for No-OP tasks
ReadyNode get_task2(LV lv) {
  ReadyNode               node;
  GV                      gv;
  MailBox                 inbox; //all about owner tasks
  //GeneratedTaskQueue      generated_q;
  ReadyTaskQueue          critic_q;
  //ReadyTaskQueue          regular_q;
  Msg*                    buf, msg;
  int                     num_removed, i;
  double                  t0;

  gv          = lv->gv;
  buf         = lv->msg_scratch_buf;
  inbox       = gv->emergent_in_box;
  critic_q    = gv->critic_ready_q;
  //regular_q   = gv->regular_ready_q;
  //generated_q = gv->gen_task_q;


  t0 = get_cur_time();
  
  //Q: should it be emergent in-box?
  num_removed = rm_msgs_inbox(inbox, buf);
  if(num_removed > 0) {
    //printf("(P%d,t%d) fetches %d msg from in-box\n", gv->proc_id, lv->thrd_id, num_removed);
    //print_in_msg(gv->proc_id, buf[0]);
    if(num_removed > lv->msg_scratch_buf_size) {fprintf(stderr, "Error: msg_scratch_buf overflow!\n"); exit(24);}
    for(i = 0; i < num_removed; i++) {
      msg = buf[i];
      process_recv_msg(msg, lv); 
      free_msg(gv->msg_pool, msg);//buf[i] stores a pointer to msg which was created by router.
    }
  }
  
  lv->getask_time += get_cur_time() - t0;

  if(critic_q->counter <= 0) return NULL;

  pthread_mutex_lock(&critic_q->lock);

  if(critic_q->counter <= 0) {
    pthread_mutex_unlock(&critic_q->lock);
    return NULL;
  }

  for(node = critic_q->senti.next; node != &critic_q->senti; node = node->next) {
    if(node->task->type == TASK_TYPE_NOOP) {
      critic_q->counter--;
      node->next->pre = node->pre;
      node->pre->next = node->next;  
      break;
    } 
  }
  if(node == &critic_q->senti) {//not found
    node = NULL;
  }

  pthread_mutex_unlock(&critic_q->lock);

  return node;
}

//if found a ready2go task, return 1. Not Used Any More!!
void check_emergent_inbox_quick(LV lv) {
  //ReadyNode               ready;
  GV                      gv;
  MailBox                 inbox; //all about owner tasks
  MailBox                 emergent_inbox;
  GeneratedTaskQueue      generated_q;
  ReadyTaskQueue          critic_q;
  ReadyTaskQueue          regular_q;
  Msg*                    buf, msg;
  int                     num_removed, i;
  double                  t0;

  gv          = lv->gv;
  buf         = lv->msg_scratch_buf;
  inbox       = gv->in_box;
  critic_q    = gv->critic_ready_q;
  regular_q   = gv->regular_ready_q;
  generated_q = gv->gen_task_q;

  emergent_inbox = gv->emergent_in_box;

  t0 = get_cur_time();

  num_removed = rm_msgs_inbox(emergent_inbox, buf);
  if(num_removed > 0) {
    //printf("(P%d,t%d) fetches %d msg from in-box\n", gv->proc_id, lv->thrd_id, num_removed);
    //print_in_msg(gv->proc_id, buf[0]);
    if(num_removed > lv->msg_scratch_buf_size) {fprintf(stderr, "Error: msg_scratch_buf overflow!\n"); exit(25);}
    for(i = 0; i < num_removed; i++) {
      msg = buf[i];
      process_recv_msg(msg, lv); 
      free_msg(gv->msg_pool, msg);//buf[i] stores a pointer to msg which was created by router.
    }
  
  } /*else {  
    //read in-mailbox, move the msgs out of the in-box. That is, an array of pointers
    num_removed = rm_msgs_inbox(inbox, buf);
    if(num_removed > 0) {
      //printf("(P%d,t%d) fetches %d msg from in-box\n", gv->proc_id, lv->thrd_id, num_removed);
      //print_in_msg(gv->proc_id, buf[0]);
      if(num_removed > lv->msg_scratch_buf_size) {fprintf(stderr, "Error: msg_scratch_buf overflow!\n"); exit(1);}
      for(i = 0; i < num_removed; i++) {
	msg = buf[i];
	process_recv_msg(msg, lv); 
	free_msg(gv->msg_pool, msg);//buf[i] stores a pointer to msg which was created by router.
      }
    }
    }*/
  
  lv->getask_time += get_cur_time() - t0;

  //return ready;
}


void deallocate_block(GV gv, Matrix mat, int idx);

void solve_task(ReadyNode ready, Matrix* matrices, int nblks, int NB, GV gv, double* tau, double* work) {
  TaskNode    task;
  int         idx, idx2, idx3, r;
  //int         i, j;
  int         INFO;
  BlockNode   *A_blks, *R_blks, *T_blks, *tmp_blks, *newT_blks; //for A, R, T
  //debug:
  double      *data_A, *data_A2, *data_R, *data_T, *data_V, *data_R1, *data_R2;
  static int  IB;
  
#ifdef V_T
  int class_id;
  int qr1_id, up1_id, qr2_id, up2_id, merge_id, merge_up_id;
  
  VT_classdef( "Computation", &class_id );
  VT_funcdef("QR1", class_id, &qr1_id);
  VT_funcdef("UP1", class_id, &up1_id);
  VT_funcdef("QR2", class_id, &qr2_id);
  VT_funcdef("UP2", class_id, &up2_id);
  VT_funcdef("MERGE", class_id, &merge_id);
  VT_funcdef("MUPDATE", class_id, &merge_up_id);
#endif

  IB = NB / 5;
  if(IB == 0) IB = NB;

  //IB = NB; //debug purpose

  A_blks   = matrices[0]->blk_ptrs;
  R_blks   = matrices[1]->blk_ptrs;
  T_blks   = matrices[2]->blk_ptrs;
  newT_blks= matrices[3]->blk_ptrs;

  task     = ready->task;

  switch(task->type) {
    case TASK_TYPE_GEQRT:
      //printf("b4 qr1()\n");
      idx    = task->output[1]*nblks + task->output[2];
      data_A = get_data_ptr(A_blks[idx]->data, NB);
      //task->private_ptr = malloc_blocknode(gv->block_pool);
      task->private_ptr = R_blks[idx]; //reuse preallocated R_kk blocks
      pthread_mutex_lock(&gv->blk_ptrs_lock);
      if(T_blks[idx] == NULL) T_blks[idx] = malloc_blocknode(gv->block_pool);
      pthread_mutex_unlock(&gv->blk_ptrs_lock);
      data_T = get_data_ptr(T_blks[idx]->data, NB);
      T_blks[idx]->nconsumers = 0;
      *(get_loc_ptr(data_T, NB))   = 2;
      *(get_loc_ptr(data_T, NB)+1) = task->output[1];
      *(get_loc_ptr(data_T, NB)+2) = task->output[2];
      //printf("P%d: QR1 wrote T %d %d\n", gv->proc_id, task->output[1], task->output[2]);
      print_one_block(data_A, NB);
#ifdef V_T
      VT_begin(qr1_id);
#endif
      core_dgeqrt_(&NB, &NB, &NB, data_A, &NB, data_T, &NB, tau, work, &INFO);
#ifdef V_T
      VT_end(qr1_id);
#endif
      //printf("P%d: computed QR1\n", gv->proc_id);
      //copy to private R buffer
      memcpy(get_data_ptr(task->private_ptr->data, NB), data_A, sizeof(double)*NB*NB);
      print_one_block(data_A, NB);
      //printf("after qr1()\n");
      if(task->private_ptr == NULL) {
	printf("Error: QR1 still allocates a NULL ptr to R block after trying\n");
	exit(26);
      }
      break;
    case  TASK_TYPE_TSQRT:
      //printf("b4 qr2()\n");
      idx    = task->inputs[1]*nblks + task->inputs[2];//[k,k]
      idx2   = task->output[1]*nblks + task->output[2];//[i,k]
      //R_blks = matrices[task->compl_outputs[0]]->blk_ptrs; //either R_kk or A_kk
      //data_R = get_data_ptr(R_blks[idx]->data, NB);
      if(task->private_ptr == NULL) {
	r = task->real_loc[0]*nblks + task->real_loc[1];
	task->private_ptr = *(R_blks+r);
	*(R_blks+r) = NULL;
	printf("QR2 has a NULL private_ptr, allocate a block for it\n");
      }
      //printf("task->private_ptr = %p\n", task->private_ptr);
      data_R = get_data_ptr(task->private_ptr->data, NB);
      data_A = get_data_ptr(A_blks[idx2]->data, NB);
      pthread_mutex_lock(&gv->blk_ptrs_lock);
      if(T_blks[idx2] == NULL) T_blks[idx2] = malloc_blocknode(gv->block_pool);
      pthread_mutex_unlock(&gv->blk_ptrs_lock);
      data_T = get_data_ptr(T_blks[idx2]->data, NB);
      T_blks[idx2]->nconsumers = 0;
      *(get_loc_ptr(data_T, NB))   = 2;
      *(get_loc_ptr(data_T, NB)+1) = task->output[1];
      *(get_loc_ptr(data_T, NB)+2) = task->output[2];
      //printf("P%d: QR1 wrote T %d %d\n", gv->proc_id, task->output[1], task->output[2]);
      //if(/*task->output[1] == 3 &&*/ task->output[2] == 2) {
      //print_one_block(data_R, NB);
	//}
#ifdef V_T
      VT_begin(qr2_id);
#endif
      core_dtsqrt_(&NB, &NB, &IB, data_R, &NB, data_A, &NB, data_T, &NB, tau, work, &INFO);
#ifdef V_T
      VT_end(qr2_id);
#endif
      //if(/*task->output[1] == 3 &&*/ task->output[2] == 2) {
      //print_one_block(data_R, NB);
	//}
      //printf("after qr2()\n");
      //application specific free
      //deallocate_block(gv, matrices[1], idx); //R_kk
      //relay R blocks
      idx = task->real_loc[0]*nblks + task->real_loc[1];
      *(R_blks+idx) = task->private_ptr;
      if(task->private_ptr == NULL) {
	printf("Error: QR2 still got a NULL ptr to R block after trying\n");
	exit(27);
      }

      break;
    case TASK_TYPE_ORMQR:
      //printf("b4 up1()\n");
      idx     = task->inputs[1]*nblks  + task->inputs[2];//[k,k]
      idx2    = task->output[1]*nblks  + task->output[2];//[k,j]
      data_A  = get_data_ptr(A_blks[idx]->data, NB);
      if(T_blks[idx] == NULL) {
        fprintf(stderr, "\n\n\nError P%d: UP1 null data T %d\n\n\n", gv->proc_id, idx); 
        print_task(task); exit(28);
      }
      data_T  = get_data_ptr(T_blks[idx]->data, NB);
      //printf("P%d: input T M%d [%d %d]\n", gv->proc_id, *(get_loc_ptr(data_T, NB)), 
      //     *(get_loc_ptr(data_T, NB)+1), *(get_loc_ptr(data_T, NB)+2));
      data_A2 = get_data_ptr(A_blks[idx2]->data, NB); //original matrix A
      task->private_ptr = malloc_blocknode(gv->block_pool);//prepare for output matrix R
      //printf("UP1: R[%d,%d] = %p\n", task->real_loc[0], task->real_loc[1], task->private_ptr);
#ifdef V_T
      VT_begin(up1_id);
#endif
      core_dormqr_(&Left, &Transpose, &NB, &NB, &NB, &NB,data_A/*V*/, &NB, data_T, &NB, data_A2/*A*/, &NB, work, &NB, &INFO);
#ifdef V_T
      VT_end(up1_id);
#endif
      memcpy(get_data_ptr(task->private_ptr->data,NB), data_A2, sizeof(double)*NB*NB);
      deallocate_block(gv, matrices[0], idx); //V_kk
      deallocate_block(gv, matrices[2], idx); //T_kk
      //printf("after up1()\n");
      //relay R blocks
      idx = task->real_loc[0]*nblks + task->real_loc[1];
      *(R_blks+idx) = task->private_ptr;
      if(task->private_ptr == NULL) {
	printf("Error: UP1 still allocates a NULL ptr to R block after trying\n");
	exit(29);
      }
      //printf("private_ptr to forward: %p\n", task->private_ptr);
      //print_one_block(get_data_ptr(R_blks[idx]->data, NB), NB);
      //print_one_block(get_data_ptr(task->private_ptr->data, NB), NB);
      break;
    case TASK_TYPE_TSSSMQR:
      //printf("b4 up2()\n");
      idx      = task->inputs[1]*nblks  + task->inputs[2];//[i,k]
      idx2     = task->inputs[7]*nblks  + task->inputs[8];//[k,j]
      idx3     = task->output[1]*nblks  + task->output[2];//[i,j]
      //data_R   = get_data_ptr(A_blks[idx2]->data, NB);
      //printf("private_ptr received: %p\n", task->private_ptr);
      if(task->private_ptr == NULL) {
#ifdef PLASMA_DEBUG
	printf("Warning: UP2 resorts to relay an R block at [%d,%d]\n", task->real_loc[0], task->real_loc[1]);
	print_task(task);
#endif
	r = task->real_loc[0]*nblks + task->real_loc[1];
	task->private_ptr = *(R_blks+r);
	*(R_blks+r) = NULL;
      }
      if(task->private_ptr == NULL) {
	printf("Error: UP2 still got a NULL ptr to R block after trying\n");
	exit(30);
      }
      data_R   = get_data_ptr(task->private_ptr->data, NB);
      //printf("UP2: R[%d,%d] = %p\n", task->real_loc[0], task->real_loc[1], task->private_ptr);
      data_A   = get_data_ptr(A_blks[idx3]->data, NB);
      data_A2  = get_data_ptr(A_blks[idx]->data, NB);
      if(T_blks[idx] == NULL) {
        fprintf(stderr, "\n\n\nError P%d: UP2 null data T %d\n\n\n", gv->proc_id, idx); 
        print_task(task); exit(31);
        //T_blks[idx] = malloc_blocknode(gv->block_pool);
      }
      data_T   = get_data_ptr(T_blks[idx]->data, NB);
      //printf("P%d: input T M%d [%d %d]\n", gv->proc_id, *(get_loc_ptr(data_T, NB)), 
      //     *(get_loc_ptr(data_T, NB)+1), *(get_loc_ptr(data_T, NB)+2));
      //if(task->output[1] == 3 && task->output[2] == 2) {
      //print_one_block(data_R, NB);
	//}
#ifdef V_T
      VT_begin(up2_id);
#endif
      core_dtsssmqr_(&Left, &Transpose, &NB, &NB, &NB, &IB/*IB*/, &NB,data_R, &NB, data_A, &NB, data_A2/*V*/, &NB, data_T, &NB, work, &NB, &INFO); 
#ifdef V_T
       VT_end(up2_id);
#endif
      //if(task->output[1] == 3 && task->output[2] == 2) {
      //print_one_block(data_R, NB);
	//}
      //printf("after up2() (1)\n");
      //application specific free
      deallocate_block(gv, matrices[0], idx);   //V_ik
      deallocate_block(gv, matrices[2], idx);   //T_ik
      //deallocate_block(gv, matrices[0], idx2);//A_kj
      //relay R blocks
      idx = task->real_loc[0]*nblks + task->real_loc[1];
      *(R_blks+idx) = task->private_ptr;
      //printf("after up2() (2)\n");
      break;
      
  case TASK_TYPE_TTQRT:
    //printf("b4 merge\n");
    idx     = task->inputs[1]*nblks + task->inputs[2]; //[i1,k];
    idx2    = task->inputs[4]*nblks + task->inputs[5]; //[i2,k];
    //R blocks should be in place!?
    //printf("R_blks=%p, idx=%d\n", R_blks, idx);
    if(R_blks[idx] == NULL) {
      printf("P%d: MERGE task has a NULL R block\n", gv->proc_id);
      //R_blks[idx] = malloc_blocknode(gv->block_pool);
      exit(1);
    }
    data_R1 = get_data_ptr(R_blks[idx]->data, NB);     //R_i1,k
#ifdef PLASMA_DEBUG
    printf("P%d: wanna compute MERGE: R1[%d,%d] R2[%d, %d]\n", 
	   gv->proc_id, task->inputs[1], task->inputs[2],
	   task->inputs[4], task->inputs[5]);
#endif
    data_R2 = get_data_ptr(R_blks[idx2]->data, NB);    //R_i2,k
    data_A  = get_data_ptr(A_blks[idx2]->data, NB);    //V_i2,k


    //allocate newT block on demand
    pthread_mutex_lock(&gv->blk_ptrs_lock);
    if(newT_blks[idx2] == NULL) newT_blks[idx2] = malloc_blocknode(gv->block_pool);
    pthread_mutex_unlock(&gv->blk_ptrs_lock);

    data_T  = get_data_ptr(newT_blks[idx2]->data, NB); //NewT_i2,k

    //print_one_block(data_R1, NB);
    //print_one_block(data_R2, NB);
    //call the computational subroutine:
#ifdef V_T
    VT_begin(merge_id);
#endif
    core_dttqrt_(&NB, &NB, &IB/*IB*/, data_R1, &NB, data_R2, &NB, data_T, &IB/*IB*/, tau, work, &INFO);
#ifdef V_T
    VT_end(merge_id);
#endif
    //printf("P%d: computed MERGE\n", gv->proc_id);
    //R2 will be overwritten by upper triangular V2, need partial copy?
    //print_one_block(data_R1, NB);
    //print_one_block(data_R2, NB);
    memcpy(data_A, data_R2, sizeof(double)*NB*NB);    
    task->private_ptr = R_blks[idx];
    //printf("after merge\n");
    break;
    
    
  case TASK_TYPE_TTSSMQR:
    idx     = task->inputs[1]*nblks  + task->inputs[2]; //[i2,k]
    idx2    = task->inputs[7]*nblks  + task->inputs[8]; //[i1,j]
    idx3    = task->inputs[10]*nblks + task->inputs[11];//[i2,j]
    data_V  = get_data_ptr(A_blks[idx]->data, NB);      //V_i2,k
    data_A  = get_data_ptr(A_blks[idx2]->data, NB);     //A_i1,j
    data_A2 = get_data_ptr(A_blks[idx3]->data, NB);     //A_i2,j
    data_T  = get_data_ptr(newT_blks[idx]->data, NB);   //T_i2,k
#ifdef V_T
    VT_begin(merge_up_id);
#endif
    //calll the computational subroutine:
    core_dttssmqr_(&Left, &Transpose, &NB, &NB, &NB, &IB/*IB*/, &NB, data_A, &NB, data_A2, &NB, data_V, &NB, data_T, &IB/*IB*/, work, &NB, &INFO);
#ifdef V_T
    VT_end(merge_up_id);
#endif
    //added on March 15, 2010.
    deallocate_block(gv, matrices[3], idx);   //new T_i2,k
    break;
    
  case TASK_TYPE_NOOP:
    //printf("Warning: NOOP operation is found\n");
    break;
    
  case TASK_TYPE_STORE_R:
      //source:
      if(task->private_ptr == NULL) {
	r = task->real_loc[0]*nblks + task->real_loc[1];
#ifdef PLASMA_DEBUG
	printf("Warning: StoreR resorts to relay an R block\n");
#endif
	if(*(R_blks+r) != NULL) {
	  task->private_ptr = *(R_blks+r);
	} else {
	  fprintf(stderr, "Error: R[%d,%d] is NULL\n", task->real_loc[0], task->real_loc[1]);
	}
      }
      //destination:
      tmp_blks   = matrices[task->output[0]]->blk_ptrs;//could be Rkk or A
      idx        = task->output[1]*nblks + task->output[2];
      if(tmp_blks[idx] == NULL) {//R_ii could be freed already
	if(task->output[1] != task->output[2]) {//impossible case
	  fprintf(stderr, "Error: A_kj block is NULL\n"); exit(32);
	}
	tmp_blks[idx] = malloc_blocknode(gv->block_pool);
      }
      data_R   = get_data_ptr(tmp_blks[idx]->data, NB); //Rkk or Akj

      //print_one_block(get_data_ptr(task->private_ptr,NB), NB);
      //printf("memcpy: &R[%d,%d] = %p\n", task->output[1], task->output[2], task->private_ptr);
      memcpy(data_R, get_data_ptr(task->private_ptr->data, NB), sizeof(double)*NB*NB);
      //printf("task Store-R[%d,%d]\n", task->real_loc[0], task->real_loc[1]);
      //print_one_block(get_data_ptr(task->private_ptr->data, NB),  NB);
      break;


  }


}


void finish_task(LV lv, ReadyNode ready) {
  TaskNode ptr;

  ptr = ready->task;

  fire_tasks(ptr, lv);//scan through tasks in ptr+1, ptr+2, ...

  free_readynode(lv->gv->ready_pool, ready);
#ifdef PLASMA_DEBUG_DEEP
  printf("(P%d, t%d) Ready Q: ", lv->gv->proc_id, lv->thrd_id);
  print_ready_q(lv->gv->regular_ready_q);
#endif
}


void deallocate_block(GV gv, Matrix mat, int idx) {
  BlockNode   del;
  double      t0;

  del = NULL;


  t0 = get_cur_time();

  //printf("enter deallocate\n");
  pthread_mutex_lock(&gv->blk_ptrs_lock);
  if(mat->blk_ptrs[idx] != NULL) {
    (mat->blk_ptrs[idx]->nconsumers)++; //critical section
    if(mat->shadow_blks_nchildren_status[idx] == 1 /*correct counter*/&& 
       mat->blk_ptrs[idx]->nconsumers == mat->shadow_blks_nchildren[idx] /*use the counter*/) {
      //printf("del %d\n", idx);
      del                = mat->blk_ptrs[idx];
      mat->blk_ptrs[idx] = NULL;
    }
  } else {
    printf("Panic: block %d becomes NULL immediately after using it\n", idx);
    //exit(33);
  }
  pthread_mutex_unlock(&gv->blk_ptrs_lock);
  
  if(del != NULL) free_blocknode(gv->block_pool, del);

  gv->dealloc_time += get_cur_time() - t0;

  //printf("exit deallocate\n");
}
