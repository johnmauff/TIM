#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

using namespace amrex;


#include "mom_continuity_ppm.hpp"
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
    // Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(*lo_i, *lo_j, 0),
	       IntVect(*hi_i, *hi_j, 0));

    // Create A4 containers for the Fortran arrays
    auto H_IN = A4::make_a4(*i_max, *j_max, 1);
    auto HL   = A4::make_a4(*i_max, *j_max, 1);
    auto HR   = A4::make_a4(*i_max, *j_max, 1);

    // Copy from Fortran arrays to A4 container
    A4::copy_fh_to_a4(h_in_h,H_IN);
    A4::copy_fh_to_a4(h_L_h,HL);
    A4::copy_fh_to_a4(h_R_h,HR);

    // Launch kernel
    ppm_limit_pos(bx,H_IN.arr, HL.arr, HR.arr, *h_min);

    // Ensure kernel is done before copying back
    Gpu::synchronize();

    // Copy device → host
    A4::copy_a4_to_fh(HL, h_L_h);
    A4::copy_a4_to_fh(HR, h_R_h);

    // Free device memory
    A4::free_a4(H_IN);
    A4::free_a4(HR);
    A4::free_a4(HL);
}

extern "C"
void ppm_limit_cw84_f (Real* h_in_h,
                      Real* h_L_h,
                      Real* h_R_h,
                      const int* lo_i, 
		      const int* hi_i,
                      const int* lo_j,
		      const int* hi_j,
                      const int* i_min,
                      const int* i_max,
                      const int* j_min,
                      const int* j_max)
{
    // Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(*lo_i, *lo_j, 0),
	       IntVect(*hi_i, *hi_j, 0));

    // Create A4 containers for the Fortran arrays
    auto H_IN = A4::make_a4(*i_max, *j_max, 1);
    auto HL   = A4::make_a4(*i_max, *j_max, 1);
    auto HR   = A4::make_a4(*i_max, *j_max, 1);

    // Copy from Fortran arrays to A4 container
    A4::copy_fh_to_a4(h_in_h,H_IN);
    A4::copy_fh_to_a4(h_L_h,HL);
    A4::copy_fh_to_a4(h_R_h,HR);

    // Launch kernel
    ppm_limit_cw84(bx,H_IN.arr, HL.arr, HR.arr);

    // Ensure kernel is done before copying back
    Gpu::synchronize();

    // Copy device → host
    A4::copy_a4_to_fh(HL, h_L_h);
    A4::copy_a4_to_fh(HR, h_R_h);

    // Free device memory
    A4::free_a4(H_IN);
    A4::free_a4(HR);
    A4::free_a4(HL);
}
