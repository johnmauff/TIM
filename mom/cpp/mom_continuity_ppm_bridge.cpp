#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

using namespace amrex;

#include "mom_continuity_ppm.hpp"
#include "tim_helper.hpp"
extern "C"
void ppm_limit_pos_c (Real* h_in_h,
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
                      const int* j_max,
		      const int* mode)
{
    bool debug = std::getenv("PPM_LIMIT_POS_DEBUG") != nullptr;

    std::string tag = "debug_ppm_limit_pos/";
    std::ofstream meta( tag + "meta.txt");

    // Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(*lo_i, *lo_j, 0),
	       IntVect(*hi_i, *hi_j, 0));

    // Create A4 containers for the Fortran arrays
    auto H_IN = timh::make_a4(*i_max, *j_max, 1, 1);
    auto HL   = timh::make_a4(*i_max, *j_max, 1, 1);
    auto HR   = timh::make_a4(*i_max, *j_max, 1, 1);

    // Copy from Fortran arrays to A4 container
    timh::copy_fh_to_a4(h_in_h,H_IN);
    timh::copy_fh_to_a4(h_L_h,HL);
    timh::copy_fh_to_a4(h_R_h,HR);

    switch (*mode) {
      case 101:  // FIXME: elimiate hardcoded constants
         //-------------------------------------------------
         // Capture input
         //-------------------------------------------------
         Gpu::streamSynchronize();
         timh::write_a4_vismf(H_IN, tag + "_H_IN_before");
         timh::write_a4_vismf(HL,   tag + "_HL_before");
         timh::write_a4_vismf(HR,   tag + "_HR_before");
         meta << "h_min = " << *h_min << "\n";
         meta << "lo_i = " << *lo_i << "\n";
         meta << "hi_i = " << *hi_i << "\n";
         meta << "lo_j = " << *lo_j << "\n";
         meta << "hi_j = " << *hi_j << "\n";
         //-------------------------------------------------
      case 103: // FIXME: eliminate hardcoded constants
         //-------------------------------------------------
         //  Execute kernel
         //-------------------------------------------------
         ppm_limit_pos(bx,H_IN.arr, HL.arr, HR.arr, *h_min);

         // Ensure kernel is done before copying back
         Gpu::synchronize();

         // Copy device → host
         timh::copy_a4_to_fh(HL, h_L_h);
         timh::copy_a4_to_fh(HR, h_R_h);
         //-------------------------------------------------
      case 102: // FIXME: eliminate hardcode constants
         //-------------------------------------------------
         // Capture Output
         //-------------------------------------------------
         timh::write_a4_vismf(HL, tag + "_HL_after");
         timh::write_a4_vismf(HR, tag + "_HR_after");
         //-------------------------------------------------
    }

    // Free device memory
    timh::free_a4(H_IN);
    timh::free_a4(HR);
    timh::free_a4(HL);
}

extern "C"
void ppm_limit_cw84_c (Real* h_in_h,
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
    auto H_IN = timh::make_a4(*i_max, *j_max, 1, 1);
    auto HL   = timh::make_a4(*i_max, *j_max, 1, 1);
    auto HR   = timh::make_a4(*i_max, *j_max, 1, 1);

    // Copy from Fortran arrays to A4 container
    timh::copy_fh_to_a4(h_in_h,H_IN);
    timh::copy_fh_to_a4(h_L_h,HL);
    timh::copy_fh_to_a4(h_R_h,HR);

    // Launch kernel
    ppm_limit_cw84(bx,H_IN.arr, HL.arr, HR.arr);

    // Ensure kernel is done before copying back
    Gpu::synchronize();

    // Copy device → host
    timh::copy_a4_to_fh(HL, h_L_h);
    timh::copy_a4_to_fh(HR, h_R_h);

    // Free device memory
    timh::free_a4(H_IN);
    timh::free_a4(HR);
    timh::free_a4(HL);
}
