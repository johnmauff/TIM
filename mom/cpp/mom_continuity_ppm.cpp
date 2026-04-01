// mom_continuity_ppm.cpp

#include <AMReX_Gpu.H>
#include <AMReX_Box.H>
#include "mom_continuity_ppm_kernel.hpp"

using namespace amrex;

void ppm_limit_pos(Box const& bx,
			  Array4<Real> const& H_IN,
			  Array4<Real> & HL,
			  Array4<Real> & HR,
                          Real hmin)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This limiter prevents undershooting minima within the domain with
        //  values less than h_min.
        ppm_limit_pos_cell(HL(i,j,k),
                           HR(i,j,k),
                           H_IN(i,j,k),
                           hmin);
    });
}

void ppm_limit_cw84(Box const& bx,
			  Array4<Real> const& H_IN,
			  Array4<Real> & HL,
			  Array4<Real> & HR)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This limiter monotonizes the parabola following
        // Colella and Woodward, 1984, Eq. 1.10
        ppm_limit_cw84_cell(HL(i,j,k),
                           HR(i,j,k),
                           H_IN(i,j,k));
    });
}
