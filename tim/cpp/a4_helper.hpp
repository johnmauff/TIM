#ifndef A4_HELPER_H_
#define A4_HELPER_H_
#include <AMReX_Array4.H>

namespace A4 {

struct A4Box
{
    amrex::Box bx;
    amrex::Real* data = nullptr;
    amrex::Array4<amrex::Real> arr;

    int nx, ny, nz, ncomp;
};

   A4Box make_a4(int nx, int ny, int nz, int ncomp=1);
   void free_a4(A4Box& a4);
   void copy_fh_to_a4(const double* f, A4Box& a4);
   void copy_a4_to_fh(const A4Box& a4, double* f);
   }
#endif
