#ifndef A4_HELPER_H_
#define A4_HELPER_H_
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>
#include <string>

namespace turbotmp {

struct A4Box
{
    amrex::Box bx;
    amrex::Real* data = nullptr;     // Raw C-pointer on device C order
    amrex::Array4<amrex::Real> arr;  //  AMReX view C order

    amrex::Real* data_f = nullptr;  // Raw C-pointer on device Fortran order

    int nx, ny, nz, ncomp;
};

   A4Box make_array4(int nx, int ny, int nz, int ncomp);
   void free_array4(A4Box& a4);
   void copy_FortranHost_to_array4(const double* f, A4Box& a4);
   void copy_array4_to_FortranHost(const A4Box& a4, double* f);
}
#endif
