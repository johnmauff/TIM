#ifndef PPM_LIMIT_K_H_
#define PPM_LIMIT_K_H_
#include <AMReX_REAL.H>
#include <AMReX_Math.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ppm_limit_pos_cell (amrex::Real& hL,
                         amrex::Real& hR,
                         const amrex::Real& hin,
                         amrex::Real  hmin) noexcept
{
    //This limiter prevents undershooting minima within the domain with
    // values less than h_min.
    amrex::Real curv = 3.0 * ((hL + hR) - 2.0 * hin);

    if (curv > 0.0)  // Only minima are limited
    {
        amrex::Real dh = hR - hL;

        if (amrex::Math::abs(dh) < curv)
        { // The parabola's minimum is within the cell
            if (hin <= hmin)
            {
                hL = hin;
                hR = hin;
            }
            else if (12.0 * curv * (hin - hmin)
                     < (curv*curv + 3.0*dh*dh))
            {
                //  The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
                //  be limited in this case.  0 < scale < 1.
                amrex::Real scale = 12.0 * curv * (hin - hmin)
                                  / (curv*curv + 3.0*dh*dh);

                hL = hin + scale * (hL - hin);
                hR = hin + scale * (hR - hin);
            }
        }
    }
}
#endif
