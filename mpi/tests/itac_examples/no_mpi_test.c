#include <stdio.h>
#include <VT.h>

int func1(){
    printf("this is func1");
}

int func2(){
    printf("this is func1");
}


int main(int argc, char *argv[]){

      int class_id;
      int func1_id, func2_id;


      VT_initialize(&argc, &argv);
      
      VT_classdef( "Computation", &class_id );
      VT_funcdef("ADVSTEP", class_id, &func1_id);
      VT_funcdef("GETBUF", class_id, &func2_id);

      VT_begin(func1_id);
      func1();
      VT_end(func2_id);


      VT_begin(func2_id);
      func2();
      VT_end(func2_id);

      VT_finalize();

    return 0;
}
