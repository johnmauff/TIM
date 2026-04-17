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

   A4Box make_a4(int nx, int ny, int nz, int ncomp);
   void free_a4(A4Box& a4);
   void copy_fh_to_a4(const double* f, A4Box& a4);
   void copy_a4_to_fh(const A4Box& a4, double* f);
   amrex::MultiFab make_mf_from_a4(const A4Box& a4);
   void fill_mf_from_a4(amrex::MultiFab & Mf, const A4Box& at);
   void write_a4_vismf(const A4Box & a4, const std::string & name);
}
#endif
