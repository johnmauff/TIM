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

/// @brief Options controlling the edge-value reconstruction scheme used by the
/// continuity solver. Field-for-field mirror of the Fortran `bind(C)` type
/// `reconstruction_CS_C` -- order and types must not change.
struct reconstruction_CS_C {
    bool upwind_1st;  ///< If true, use a first-order upwind scheme.
    bool monotonic;   ///< If true, use the Colella & Woodward monotonic limiter.
    bool simple_2nd;  ///< If true, use a simple second order interpolation.
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

/**
 * @brief Sets the effective open face areas and barotropic-velocity
 * corrections at zonal faces as a function of barotropic flow, for use
 * by the barotropic solver's transport-adjustment iteration
 */
void set_zonal_BT_cont(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<const Real> const&,   //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,         //!< Effective open face area, west, 0 transport
    Array4<Real> const&,         //!< Effective open face area, east, 0 transport
    Array4<Real> const&,         //!< Effective open face area, westerly test velocity
    Array4<Real> const&,         //!< Effective open face area, easterly test velocity
    Array4<Real> const&,         //!< Westerly correction to the barotropic velocity
    Array4<Real> const&,         //!< Easterly correction to the barotropic velocity
    Array4<const Real> const&,   //!< Barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const&,   //!< Partial derivative of du_err with du at 0 adjustment [H L ~> m2 or kg m-1]
    Array4<const Real> const&,   //!< Maximum acceptable value of du [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Minimum acceptable value of du [L T-1 ~> m s-1]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's u-point x-extent [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< The grid cell's 1/dxT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    Array4<const Real> const&,   //!< Fraction of momentum/barotropic acceleration remaining
                                  //!< after viscosity [nondim]
    Array4<const Real> const&,   //!< Maximum allowable viscosity remnant [nondim]
    Array4<const int> const&,    //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const&);  //!< Fractional open area of U-faces [nondim]

/**
 * @brief Sets the effective open face areas and barotropic-velocity
 * corrections at meridional faces as a function of barotropic flow, for
 * use by the barotropic solver's transport-adjustment iteration
 */
void set_merid_BT_cont(
    const Box&,                  //!< Iteration box for continuity solver
    Array4<const Real> const&,   //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,   //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,         //!< Effective open face area, south, 0 transport
    Array4<Real> const&,         //!< Effective open face area, north, 0 transport
    Array4<Real> const&,         //!< Effective open face area, southerly test velocity
    Array4<Real> const&,         //!< Effective open face area, northerly test velocity
    Array4<Real> const&,         //!< Southerly correction to the barotropic velocity
    Array4<Real> const&,         //!< Northerly correction to the barotropic velocity
    Array4<const Real> const&,   //!< Barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const&,   //!< Partial derivative of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1]
    Array4<const Real> const&,   //!< Maximum acceptable value of dv [L T-1 ~> m s-1]
    Array4<const Real> const&,   //!< Minimum acceptable value of dv [L T-1 ~> m s-1]
    Real,                        //!< Time increment [T ~> s]
    Array4<const Real> const&,   //!< The grid cell's v-point y-extent [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,   //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,   //!< The grid cell's 1/dyT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    Array4<const Real> const&,   //!< Fraction of momentum/barotropic acceleration remaining
                                  //!< after viscosity [nondim]
    Array4<const Real> const&,   //!< Maximum allowable viscosity remnant [nondim]
    Array4<const int> const&,    //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const&);  //!< Fractional open area of V-faces [nondim]

/**
 * @brief Newton-iterates a barotropic velocity correction per zonal face so
 * that the vertically-summed zonal mass/volume transport matches the target
 * barotropic transport, to within the transport-adjustment iteration's
 * tolerance
 */
void zonal_flux_adjust(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const&,    //!< Partial derivative of du_err with du at 0 adjustment [H L ~> m2 or kg m-1]
    Array4<Real> const&,          //!< The barotropic velocity adjustment [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Maximum acceptable value of du [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Minimum acceptable value of du [L T-1 ~> m s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dxT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity [nondim]
    Array4<const int> const&,     //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const&,    //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const&,    //!< Summed volume flux through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Volume flux through zonal faces, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    OceanOBC*);                   //!< Open boundary control structure

/**
 * @brief Newton-iterates a barotropic velocity correction per meridional
 * face so that the vertically-summed meridional mass/volume transport
 * matches the target barotropic transport, to within the
 * transport-adjustment iteration's tolerance
 */
void meridional_flux_adjust(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Summed transport with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<const Real> const&,    //!< Partial derivative of dv_err with dv at 0 adjustment [H L ~> m2 or kg m-1]
    Array4<Real> const&,          //!< The barotropic velocity adjustment [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Maximum acceptable value of dv [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Minimum acceptable value of dv [L T-1 ~> m s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dyT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity [nondim]
    Array4<const int> const&,     //!< Logical flag (0/1) indicating which I values to work on
    Array4<const Real> const&,    //!< Fractional open area of V-faces [nondim]
    Array4<const Real> const&,    //!< Summed volume flux through meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Volume flux through meridional faces, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    OceanOBC*);                   //!< Open boundary control structure

/**
 * @brief Zonal mass/volume flux orchestrator -- computes the zonal PPM
 * transport, then (when a barotropic target transport and/or BT_cont
 * output is requested) the transport-adjustment correction via
 * zonal_flux_adjust, set_zonal_BT_cont, and zonal_flux_thickness
 */
void zonal_mass_flux(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< West edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< East edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,          //!< Volume flux through zonal faces, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const&,    //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const&,    //!< The x-extent of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< 0 for land points, 1 for ocean points at u-locations [nondim]
    Array4<const Real> const&,    //!< The grid cell's u-point x-extent [L ~> m]
    Real,                         //!< A negligibly small thickness used to avoid division
                                   //!< by zero [H ~> m or kg m-2]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&,    //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const&,    //!< Summed volume flux through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity [nondim]; may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Zonal velocity with barotropic correction [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, west, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, east, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, westerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, easterly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Westerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Easterly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&);         //!< Zonal velocity increment from u that gives uhbt as the
                                   //!< depth-integrated transport [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)

/**
 * @brief Meridional mass/volume flux orchestrator -- computes the
 * meridional PPM transport, then (when a barotropic target transport
 * and/or BT_cont output is requested) the transport-adjustment
 * correction via meridional_flux_adjust and set_merid_BT_cont
 */
void meridional_mass_flux(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< South edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< North edge thickness in the reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,          //!< Volume flux through meridional faces, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Array4<const Real> const&,    //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const&,    //!< The y-extent of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< 0 for land points, 1 for ocean points at v-locations [nondim]
    Array4<const Real> const&,    //!< The grid cell's v-point y-extent [L ~> m]
    int,                          //!< Start i-index of the data domain (0-based)
    int,                          //!< End i-index of the data domain (0-based)
    Real,                         //!< A negligibly small thickness used to avoid division
                                   //!< by zero [H ~> m or kg m-2]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&,    //!< Fractional open area of V-faces [nondim]
    Array4<const Real> const&,    //!< Summed volume flux through meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity [nondim]; may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Meridional velocity with barotropic correction [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, south, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, north, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, southerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, northerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Southerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Northerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&);         //!< Meridional velocity increment from v that gives vhbt as the
                                   //!< depth-integrated transport [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)

/**
 * @brief Monolithic continuity solver -- reconstructs edge thicknesses, then
 * advects (via zonal_mass_flux/meridional_mass_flux and
 * continuity_zonal_convergence/continuity_meridional_convergence) first in
 * one direction and then the other, in the order set by x_first
 */
void continuity_PPM(
    Array4<const Real> const&,    //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Initial layer thickness [H ~> m or kg m-2]
    Array4<Real> const&,          //!< Final layer thickness [H ~> m or kg m-2]
    Array4<Real> const&,          //!< Zonal volume flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<Real> const&,          //!< Meridional volume flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    const Box&,                   //!< The core (unstencilled) iteration box
    int,                          //!< The continuity solver stencil width
    bool,                         //!< If true, advect zonally before meridionally
    Array4<const Real> const&,    //!< Cell land/ocean mask [nondim]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const&,    //!< The area of the h-cell [L2 ~> m2]
    Array4<const Real> const&,    //!< The x-extent of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< 0 for land points, 1 for ocean points at u-locations [nondim]
    Array4<const Real> const&,    //!< The grid cell's u-point x-extent [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Array4<const Real> const&,    //!< The y-extent of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< 0 for land points, 1 for ocean points at v-locations [nondim]
    Array4<const Real> const&,    //!< The grid cell's v-point y-extent [L ~> m]
    int,                          //!< Start i-index of the data domain (0-based)
    int,                          //!< End i-index of the data domain (0-based)
    Real,                         //!< A one-Angstrom thickness [H ~> m or kg m-2]
    Real,                         //!< A negligibly small thickness used to avoid division
                                   //!< by zero [H ~> m or kg m-2]
    const reconstruction_CS_C&,   //!< Options controlling the edge-value reconstruction scheme
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&,    //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const&,    //!< Fractional open area of V-faces [nondim]
    Array4<const Real> const&,    //!< Summed volume flux through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<const Real> const&,    //!< Summed volume flux through meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity, zonal [nondim]; may be absent (.p == nullptr)
    Array4<const Real> const&,    //!< Fraction of momentum/barotropic acceleration remaining
                                   //!< after viscosity, meridional [nondim]; may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Zonal velocity with barotropic correction [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Meridional velocity with barotropic correction [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, west, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, east, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, westerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, easterly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Westerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Easterly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, south, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, north, 0 transport;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, southerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Effective open face area, northerly test velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Southerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Northerly correction to the barotropic velocity;
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&,          //!< Zonal velocity increment from u that gives uhbt as the
                                   //!< depth-integrated transport [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)
    Array4<Real> const&);         //!< Meridional velocity increment from v that gives vhbt as the
                                   //!< depth-integrated transport [L T-1 ~> m s-1];
                                   //!< may be absent (.p == nullptr)

/**
 * @brief Sums the zonal PPM transport over all layers to give the
 * barotropic (depth-integrated) zonal transport, uhbt
 */
void zonal_BT_mass_flux(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Western edge thickness in the PPM reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Eastern edge thickness in the PPM reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,          //!< The summed volume flux through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dxT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&);   //!< Fractional open area of U-faces [nondim]

/**
 * @brief Sums the meridional PPM transport over all layers to give the
 * barotropic (depth-integrated) meridional transport, vhbt
 */
void meridional_BT_mass_flux(
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness used to calculate fluxes [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Southern edge thickness in the PPM reconstruction [H ~> m or kg m-2]
    Array4<const Real> const&,    //!< Northern edge thickness in the PPM reconstruction [H ~> m or kg m-2]
    Array4<Real> const&,          //!< The summed volume flux through meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dyT [L-1 ~> m-1]
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&);   //!< Fractional open area of V-faces [nondim]

/**
 * @brief Reconstructs zonal and meridional edge thicknesses, then computes
 * the barotropic (depth-integrated) zonal and meridional transports uhbt
 * and vhbt via zonal_BT_mass_flux and meridional_BT_mass_flux
 */
void continuity_PPM_2d_fluxes(
    Array4<const Real> const&,    //!< Zonal velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Meridional velocity [L T-1 ~> m s-1]
    Array4<const Real> const&,    //!< Layer thickness [H ~> m or kg m-2]
    Array4<Real> const&,          //!< Vertically summed thickness flux through zonal
                                   //!< faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Array4<Real> const&,          //!< Vertically summed thickness flux through meridional
                                   //!< faces [H L2 T-1 ~> m3 s-1 or kg s-1]
    Real,                         //!< Time increment [T ~> s]
    const Box&,                   //!< Iteration box for continuity solver
    Array4<const Real> const&,    //!< Cell land/ocean mask [nondim]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the u-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/areaT [L-2 ~> m-2]
    Array4<const Real> const&,    //!< The grid cell's 1/dxT [L-1 ~> m-1]
    Array4<const Real> const&,    //!< The grid cell's unblocked lengths of the v-faces of the h-cell [L ~> m]
    Array4<const Real> const&,    //!< The grid cell's 1/dyT [L-1 ~> m-1]
    Real,                         //!< A one-Angstrom thickness [H ~> m or kg m-2]
    const reconstruction_CS_C&,   //!< Options controlling the edge-value reconstruction scheme
    const transport_adjust_CS_C&, //!< Options controlling the transport adjustment and
                                   //!< barotropic-consistency iteration
    OceanOBC*,                    //!< Open boundary control structure
    Array4<const Real> const&,    //!< Fractional open area of U-faces [nondim]
    Array4<const Real> const&);   //!< Fractional open area of V-faces [nondim]
}
