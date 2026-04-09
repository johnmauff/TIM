#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

using namespace amrex;

#include "mom_continuity_ppm.hpp"
#include "tim_helper.hpp"

struct RealArray_C {
    double* data;
    int* shape;
    int* lb;
    int* ub;
    int rank;
};
struct Box_C {
    int* idxS;
    int* idxE;
};

/**
 * @brief Bridge for the function PPM_limit_pos function
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++ 
 * implementation based on the setting of the @p mode parameter.
 *
 * @param h_in_h Layer thickness [H → m or kg m^-2] 
 * 	on the host in Fortran order
 * @param h_L_h, Left thickness of the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2]
 * @param h_R_h, Right thickness in the reconstruction {host, Fortran order} 
 * 	[H → m or kg m^-2] 
 * @param hmin Minimum thickness allowed by the parabolic fit (host, Fortran order) 
 * 	[H → m or kg m^-2]
 * @param lb  Lower loop bounds
 * @param ub  Upper loop bounds
 * @param mode Execution of the mode of the bridge.
 *
 * @return Modified thickness values @p h_L_h and @p h_R_h
 */
extern "C"
void ppm_limit_pos_bridge(const Box_C* bx_h,
		           const RealArray_C* h_in_h,
			   RealArray_C* h_L_h,
			   RealArray_C* h_R_h,
	                   const Real* h_min,
			   const int* mode)
{ 
    /// Define the output directory for captured data
    std::string dir = "debug/ppm_limit_pos/";

    /// Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(bx_h->idxS[0]-1, bx_h->idxS[1]-1, bx_h->idxS[2]-1),
	   IntVect(bx_h->idxE[0]-1, bx_h->idxE[1]-1, bx_h->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto H_IN = timh::make_a4(h_in_h->shape[0], h_in_h->shape[1], 1, 1);
    auto HL   = timh::make_a4(h_L_h->shape[0],  h_L_h->shape[1],  1, 1);
    auto HR   = timh::make_a4(h_R_h->shape[0],  h_R_h->shape[1],  1, 1);

    /// Copy from Fortran arrays to A4 container
    timh::copy_fh_to_a4(h_in_h->data, H_IN);
    timh::copy_fh_to_a4(h_L_h->data, HL);
    timh::copy_fh_to_a4(h_R_h->data, HR);

    switch (*mode) {
      case 101:  { /// FIXME: eliminate hardcoded constants
         ///-------------------------------------------------
         /// Capture input
         ///-------------------------------------------------
         Gpu::streamSynchronize();
         timh::write_a4_vismf(H_IN, dir + "_H_IN_before");
         timh::write_a4_vismf(HL,   dir + "_HL_before");
         timh::write_a4_vismf(HR,   dir + "_HR_before");

         /// Output the scalar values...
         std::ofstream meta( dir + "meta.txt");
         meta << "h_min = " << *h_min << "\n";
	 meta << "bx = " << bx << "\n";
         break;
         ///-------------------------------------------------
      }
      case 103: { /// FIXME: eliminate hardcoded constants
         ///-------------------------------------------------
         ///  Execute kernel
         ///-------------------------------------------------
         ppm_limit_pos(bx,H_IN.arr, HL.arr, HR.arr, *h_min);

         /// Ensure kernel is done before copying back
         Gpu::synchronize();

         /// Copy device → host
         timh::copy_a4_to_fh(HL, h_L_h->data);
         timh::copy_a4_to_fh(HR, h_R_h->data);
         break;
         ///-------------------------------------------------
      }
      case 102: { /// FIXME: eliminate hardcode constants
         ///-------------------------------------------------
         /// Capture Output
         ///-------------------------------------------------
         timh::write_a4_vismf(HL, dir + "_HL_after");
         timh::write_a4_vismf(HR, dir + "_HR_after");
         break;
         ///-------------------------------------------------
      }
    }

    /// Free a4 container
    timh::free_a4(H_IN);
    timh::free_a4(HR);
    timh::free_a4(HL);
}

/**
 * @brief Bridge for the function PPM_limit_cw84 function
 *
 * This function acts as a bridge between a Fortran interface
 * and an AMReX C++ implementation. It also provides the ability
 * to either capture the input, or output or execute the AMReX C++
 * implementation based on the setting of the @p mode parameter.
 *
 * @param bx_h   Box over which to iterate 
 * @param h_in_h Layer thickness [H → m or kg m^-2]
 *      on the host in Fortran order
 * @param h_L_h, Left thickness of the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 * @param h_R_h, Right thickness in the reconstruction {host, Fortran order}
 *      [H → m or kg m^-2]
 * @param mode Execution of the mode of the bridge.
 *
 * @return Modified thickness values @p h_L_h and @p h_R_h
 */
extern "C"
void ppm_limit_cw84_bridge(const Box_C* bx_h,
		             const RealArray_C* h_in_h,
                             RealArray_C* h_L_h,
                             RealArray_C* h_R_h,
		      	     const int* mode)
{
    /// Define the output directory for captured data
    std::string dir = "debug/ppm_limit_cw84/";
    /// Define Active domain (kernel launch only on real cells)
    Box bx(IntVect(bx_h->idxS[0]-1, bx_h->idxS[1]-1, bx_h->idxS[2]-1),
           IntVect(bx_h->idxE[0]-1, bx_h->idxE[1]-1, bx_h->idxE[2]-1));

    /// Create A4 containers for the Fortran arrays
    auto H_IN = timh::make_a4(h_in_h->shape[0], h_in_h->shape[1], 1, 1);
    auto HL   = timh::make_a4(h_L_h->shape[0],  h_L_h->shape[1],  1, 1);
    auto HR   = timh::make_a4(h_R_h->shape[0],  h_R_h->shape[1],  1, 1);


    /// Copy from Fortran arrays to A4 container
    timh::copy_fh_to_a4(h_in_h->data, H_IN);
    timh::copy_fh_to_a4(h_L_h->data, HL);
    timh::copy_fh_to_a4(h_R_h->data, HR);

    switch(*mode) {
      case 101: { /// FIXME: eliminate hardcoded constants 
         ///-------------------------------------------------
         /// Capture input
         ///-------------------------------------------------
         Gpu::streamSynchronize();
         timh::write_a4_vismf(H_IN, dir + "_H_IN_before");
         timh::write_a4_vismf(HL,   dir + "_HL_before");
         timh::write_a4_vismf(HR,   dir + "_HR_before");

         /// Output the scalar values...
         std::ofstream meta( dir + "meta.txt");
	 meta << "bx = " << bx << "\n";
         break;
      }
      case 103: { /// FIXME: eliminate hardcoded constants
         ///-------------------------------------------------
         ///  Execute kernel
         ///-------------------------------------------------
        ppm_limit_cw84(bx,H_IN.arr, HL.arr, HR.arr);

        /// Ensure kernel is done before copying back
        Gpu::synchronize();

        /// Copy device → host
        timh::copy_a4_to_fh(HL, h_L_h->data);
        timh::copy_a4_to_fh(HR, h_R_h->data);
	break;
      }
      case 102: { /// FIXME: eliminate hardcode constants
   	 ///-------------------------------------------------
         /// Capture Output
         ///-------------------------------------------------
         timh::write_a4_vismf(HL, dir + "_HL_after");
         timh::write_a4_vismf(HR, dir + "_HR_after");
         break;
         ///-------------------------------------------------
      }
    }
    /// Free memory from a4 containers
    timh::free_a4(H_IN);
    timh::free_a4(HR);
    timh::free_a4(HL);
}

