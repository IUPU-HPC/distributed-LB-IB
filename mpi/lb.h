#ifndef LB_H
#define LB_H

extern const int c[19][3];
extern const double t[19];

// compute local equilibrium from rho and u
inline double _computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uz, double uSqr)
{
    double c_u = c[iPop][0]*ux + c[iPop][1]*uy + c[iPop][2]*uz;
    return rho * t[iPop] * (
               1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr
           );
}

// initialize a node to its local equilibrium term
inline void computeEquilibrium(Fluidnode* node) {
  int iPop;
  double rho = node->rho;
  double ux = node->vel_x;
  double uy = node->vel_y;
  double uz = node->vel_z;

  double uSqr = ux*ux + uy*uy + uz*uz;
  for (iPop = 0; iPop < 19; ++iPop) {
      node->dfeq[iPop] =
          _computeEquilibrium(iPop, rho, ux, uy, uz, uSqr);
  }
}

#endif //LB_H