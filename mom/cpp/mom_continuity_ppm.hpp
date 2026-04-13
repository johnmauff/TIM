// mom_continuity_ppm.hpp

#ifndef CONTINUITY_PPM_H_
#define CONTINUITY_PPM_H_

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

/**
 * @brief Piecewise parabolic limiter
 */
void ppm_limit_pos(const amrex::Box & bx,
                   amrex::Array4<const amrex::Real> const& H_IN,
                   amrex::Array4<amrex::Real> const& HL,
                   amrex::Array4<amrex::Real> const& HR,
                   const amrex::Real hmin);

/**
 * @brief Piecewise parabolic limiter of Colella and Woodward, 1984
 */
void ppm_limit_cw84(const amrex::Box & bx,
                    amrex::Array4<const amrex::Real> const& H_IN,
                    amrex::Array4<amrex::Real> const& HL,
                    amrex::Array4<amrex::Real> const& HR);

/**
 * @brief Piecewise reconstruction in the y dimension
 */
void PPM_reconstruction_y(
    const amrex::Box&,
    amrex::Array4<const amrex::Real> const&,
    amrex::Array4<amrex::Real> const&,
    amrex::Array4<amrex::Real> const&,
    amrex::Array4<const amrex::Real> const&,
    amrex::Real,
    bool,
    bool,
    OceanOBC*);

#endif
