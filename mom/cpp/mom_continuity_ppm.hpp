// mom_continuity_ppm.hpp
// SKILLS: 0.3.1
#pragma once
/**
 * @file mom_continuity_ppm.hpp
 * @brief Box-level AMReX kernel declarations for MOM6 PPM continuity
 *        (piecewise parabolic reconstruction and edge-thickness routines).
 */

#include "mom_continuity_ppm_kernel.hpp"

struct OceanOBC;    // Undefined at the moment

/// @brief Options controlling the transport adjustment and barotropic-consistency
/// iteration used by the continuity solver. Field-for-field mirror of the Fortran
/// `bind(C)` type `transport_adjust_CS_C` -- order and types must not change.
struct transport_adjust_CS_C {
    double tol_eta;            ///< Tolerance for free-surface height discrepancies.
    double tol_vel;            ///< Tolerance for barotropic velocity discrepancies.
    double CFL_limit_adjust;   ///< Maximum CFL of the adjusted velocities.
    bool   aggress_adjust;     ///< If true, allow a larger relative CFL change.
    bool   vol_CFL;            ///< If true, use the ratio of open face lengths to
                                ///< tracer cell areas when estimating CFL numbers.
    bool   better_iter;        ///< If true, use a velocity-based iteration criterion.
    bool   use_visc_rem_max;   ///< If true, use limiting bounds for viscous columns.
    bool   marginal_faces;     ///< If true, use marginal face areas as barotropic weights.
};

/// @brief AMReX ports of MOM6 numerical kernels.
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
 * @brief Meridional volume/thickness flux — PPM-reconstructed edge
 * thickness advected by the meridional velocity, scaled by viscosity
 * remnant and open-face area
 */
void meridional_flux_thickness(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<const Real> const&,   //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,         //!< Effective thickness at meridional faces [H ~> m or kg m-2]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< The grid cell's 1/dyT [L-1 ~> m-1]
    bool,                        //!< If true, rescale the ratio of face areas to the cell areas when estimating the CFL number
    bool,                        //!< If true, report the marginal face thicknesses; otherwise report transport-averaged thicknesses
    OceanOBC*,                   //!< Open boundaries control structure
    Array4<const Real> const&,   //!< Fractional open area of V-faces [nondim]
    Array4<const Real> const&);  //!< Fraction of momentum remaining after viscosity, and fraction of a
                                  //!< time-step's barotropic acceleration a layer experiences [nondim];
                                  //!< between 0 (bottom) and 1 (far above the bottom)

/**
 * @brief Zonal volume/thickness flux — PPM-reconstructed edge
 * thickness advected by the zonal velocity, scaled by viscosity
 * remnant and open-face area
 */
void zonal_flux_thickness(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<const Real> const&,   //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,         //!< Effective thickness at zonal faces [H ~> m or kg m-2]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< The grid cell's 1/dxT [L-1 ~> m-1]
    bool,                        //!< If true, rescale the ratio of face areas to the cell areas when estimating the CFL number
    bool,                        //!< If true, report the marginal face thicknesses; otherwise report transport-averaged thicknesses
    OceanOBC*,                   //!< Open boundaries control structure
    Array4<const Real> const&,   //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const&);  //!< Fraction of momentum remaining after viscosity, and fraction of a
                                  //!< time-step's barotropic acceleration a layer experiences [nondim];
                                  //!< between 0 (bottom) and 1 (far above the bottom)

/**
 * @brief Zonal continuity update — advances layer thickness by the
 * convergence of the zonal thickness flux
 */
void continuity_zonal_convergence(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<Real> const&,         //!< Final layer thickness [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< Zonal thickness flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< Initial layer thickness [H ~> m or kg m-2]; may be absent
                                  //!< (.p == nullptr), in which case the final thickness is also
                                  //!< used as the initial thickness
    Real);                       //!< The minimum layer thickness [H ~> m or kg m-2]

/**
 * @brief Meridional continuity update — advances layer thickness by the
 * convergence of the meridional thickness flux
 */
void continuity_meridional_convergence(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<Real> const&,         //!< Final layer thickness [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< Meridional thickness flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< Initial layer thickness [H ~> m or kg m-2]; may be absent
                                  //!< (.p == nullptr), in which case the final thickness is also
                                  //!< used as the initial thickness
    Real);                       //!< The minimum layer thickness [H ~> m or kg m-2]
}
