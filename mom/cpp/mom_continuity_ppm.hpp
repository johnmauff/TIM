// mom_continuity_ppm.hpp
#pragma once

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

namespace MOM {
using amrex::Box;
using amrex::Array4;
/**
 * @brief Piecewise parabolic limiter
 */
void ppm_limit_pos(const Box &,
                   Array4<const Real> const&,
                   Array4<Real> const&,
                   Array4<Real> const&,
                   const Real);

/**
 * @brief Piecewise parabolic limiter of Colella and Woodward, 1984
 */
void ppm_limit_cw84(const Box&,
                    Array4<const Real> const&,
                    Array4<Real> const&,
                    Array4<Real> const&);

/**
 * @brief Piecewise reconstruction in the y dimension
 */
void PPM_reconstruction_y(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Piecewise reconstruction in the x dimension
 */
void PPM_reconstruction_x(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Zonal edge thickness — upwind copy or x-direction PPM reconstruction
 */
void zonal_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Meridional edge thickness — upwind copy or y-direction PPM reconstruction
 */
void meridional_edge_thickness(
    const Box&,
    Array4<const Real> const&,
    Array4<Real> const&,
    Array4<Real> const&,
    Array4<const Real> const&,
    Real,
    bool,
    bool,
    bool,
    OceanOBC*);

/**
 * @brief Effective thickness at meridional (V) faces for continuity flux
 *
 * Ports meridional_flux_thickness_fortran from MOM_continuity_PPM.F90.
 * The iteration box bxV = growLo(bxC, y, 1) is derived internally.
 * OBC support is not yet implemented; a non-null obc aborts.
 */
void meridional_flux_thickness(
    const Box&,                    // bxC — H-grid continuity box
    Array4<const Real> const&,     // v          — meridional velocity [L T-1]
    Array4<const Real> const&,     // h          — layer thickness [H]
    Array4<const Real> const&,     // h_S        — south edge thickness [H]
    Array4<const Real> const&,     // h_N        — north edge thickness [H]
    Array4<Real> const&,           // h_v (inout)— effective thickness at V-faces [H]
    Real,                          // dt
    bool,                          // vol_CFL
    bool,                          // marginal
    Array4<const Real> const&,     // dx_Cv      — meridional face length [L] (2-D)
    Array4<const Real> const&,     // IareaT     — inverse tracer cell area [L-2] (2-D)
    Array4<const Real> const&,     // IdyT       — inverse tracer cell y-extent [L-1] (2-D)
    OceanOBC*,                     // obc
    Array4<const Real> const&,     // por_face_areaV — fractional open area of V-faces
    Array4<const Real> const&      // visc_rem_v     — viscous remainder fraction
);

/**
 * @brief Effective thickness at zonal (U) faces for continuity flux
 *
 * Ports zonal_flux_thickness_fortran from MOM_continuity_PPM.F90.
 * The iteration box bxU = growLo(bxC, x, 1) is derived internally.
 * OBC support is not yet implemented; a non-null obc aborts.
 */
void zonal_flux_thickness(
    const Box&,                    // bxC — H-grid continuity box
    Array4<const Real> const&,     // u          — zonal velocity [L T-1]
    Array4<const Real> const&,     // h          — layer thickness [H]
    Array4<const Real> const&,     // h_W        — west edge thickness [H]
    Array4<const Real> const&,     // h_E        — east edge thickness [H]
    Array4<Real> const&,           // h_u (inout)— effective thickness at U-faces [H]
    Real,                          // dt
    bool,                          // vol_CFL
    bool,                          // marginal
    Array4<const Real> const&,     // dy_Cu      — zonal face length [L] (2-D)
    Array4<const Real> const&,     // IareaT     — inverse tracer cell area [L-2] (2-D)
    Array4<const Real> const&,     // IdxT       — inverse tracer cell x-extent [L-1] (2-D)
    OceanOBC*,                     // obc
    Array4<const Real> const&,     // por_face_areaU — fractional open area of U-faces
    Array4<const Real> const&      // visc_rem_u     — viscous remainder fraction
);
}
