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
    Box active(IntVect(*lo_i, *lo_j, 0),
	       IntVect(*hi_i, *hi_j, 0));

    Box bx_halo(IntVect(*i_min, *j_min, 0),
                IntVect(*i_max, *j_max, 0));

    const Long npts = bx_halo.numPts();

    // This is trying out the A4 helper capability
    auto H_IN = make_a4(*i_max,*j_max,1);
    copy_f_to_a4(h_in_h,H_IN);
    // Do some flops here

    // Allocate device memory
    Real* h_in_d = (Real*) The_Device_Arena()->alloc(npts * sizeof(Real));
    Real* h_L_d  = (Real*) The_Device_Arena()->alloc(npts * sizeof(Real));
    Real* h_R_d  = (Real*) The_Device_Arena()->alloc(npts * sizeof(Real));

    // Copy host → device
    Gpu::copy(Gpu::hostToDevice, h_in_h, h_in_h + npts, h_in_d);
    Gpu::copy(Gpu::hostToDevice, h_L_h,  h_L_h  + npts, h_L_d);
    Gpu::copy(Gpu::hostToDevice, h_R_h,  h_R_h  + npts, h_R_d);

    // Wrap as Array4
    Array4<Real const> h_in(h_in_d, lbound(bx_halo), ubound(bx_halo),1);
    Array4<Real>         hL(h_L_d,  lbound(bx_halo), ubound(bx_halo),1);
    Array4<Real>         hR(h_R_d,  lbound(bx_halo), ubound(bx_halo),1);

    Real hmin = *h_min;


    // Launch kernel
    ParallelFor(active, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        ppm_limit_pos_cell(hL(i,j,k), hR(i,j,k), H_IN.arr(i,j,k), hmin);
    });

    // Ensure kernel is done before copying back
    Gpu::synchronize();

    // Copy device → host
    Gpu::copy(Gpu::deviceToHost, h_L_d, h_L_d + npts, h_L_h);
    Gpu::copy(Gpu::deviceToHost, h_R_d, h_R_d + npts, h_R_h);

    // Free device memory
    The_Device_Arena()->free(h_in_d);
    The_Device_Arena()->free(h_L_d);
    The_Device_Arena()->free(h_R_d);

    free_a4(H_IN);
}
