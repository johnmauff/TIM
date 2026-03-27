#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

using namespace amrex;


#include "ppm_limit_K.hpp"
#include "a4_helper.hpp"

extern "C"
void ppm_limit_pos_f (Real* h_in_h,
                      Real* h_L_h,
                      Real* h_R_h,
                      const Real* h_min,
                      const int* lo_i, 
		      const int* hi_i,
                      const int* lo_j,
		      const int* hi_j,
                      const int* i_min,
                      const int* i_max,
                      const int* j_min,
                      const int* j_max)
{
    // Active domain (kernel launch only on real cells)
    Box bx(IntVect(*lo_i, *lo_j, 0),
	       IntVect(*hi_i, *hi_j, 0));

    // Create A4 containers for the Fortran arrays
    auto H_IN = A4::make_a4(*i_max, *j_max, 1);
    auto HL   = A4::make_a4(*i_max, *j_max, 1);
    auto HR   = A4::make_a4(*i_max, *j_max, 1);

    // Copy from Fortran arrays to A4 container
    A4::copy_f_to_a4(h_in_h,H_IN);
    A4::copy_f_to_a4(h_L_h,HL);
    A4::copy_f_to_a4(h_R_h,HR);

    Real hmin = *h_min;

    // Launch kernel
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        ppm_limit_pos_cell(HL.arr(i,j,k), HR.arr(i,j,k), H_IN.arr(i,j,k), hmin);
    });

    // Ensure kernel is done before copying back
    Gpu::synchronize();

    // Copy device → host
    A4::copy_a4_to_f(HL, h_L_h);
    A4::copy_a4_to_f(HR, h_R_h);

    // Free device memory
    A4::free_a4(H_IN);
    A4::free_a4(HR);
    A4::free_a4(HL);
}
