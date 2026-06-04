// Unit tests for MOM::ppm_limit_pos / PPM_reconstruction_x / PPM_reconstruction_y.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"
#include "mom_continuity_ppm.hpp"

using test_mom::expect_arrays_equal;
using test_mom::to_host_fab;

// -------------------------------------------------------------------------
// ppm_limit_pos
// -------------------------------------------------------------------------
TEST(PpmLimitPos, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_limit_pos");

    const auto   bx        = captured.box("_bx");
    const auto   h_in      = captured.fab_device("_h_in");
    auto         h_L       = captured.fab_device("_h_L_before");
    auto         h_R       = captured.fab_device("_h_R_before");
    const auto   h_L_after = captured.fab_host("_h_L_after");
    const auto   h_R_after = captured.fab_host("_h_R_after");
    const double h_min     = captured.real64("_h_min");

    MOM::ppm_limit_pos(bx,
                       h_in.const_array(),
                       h_L.array(),
                       h_R.array(),
                       h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_L_after, to_host_fab(h_L), "h_L");
    expect_arrays_equal(h_R_after, to_host_fab(h_R), "h_R");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_x
// -------------------------------------------------------------------------
TEST(PpmReconstructionX, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_x");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_W        = captured.fab_device("_h_W_before");
    auto         h_E        = captured.fab_device("_h_E_before");
    const auto   mask2d     = captured.fab_device("_mask2d_t");
    const auto   h_W_after  = captured.fab_host("_h_W_after");
    const auto   h_E_after  = captured.fab_host("_h_E_after");
    const double h_min      = captured.real64("_h_min");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::PPM_reconstruction_x(bxH,
                              h_in.const_array(),
                              h_W.array(),
                              h_E.array(),
                              mask2d.const_array(),
                              h_min,
                              monotonic,
                              simple_2nd,
                              /*OBC=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_W_after, to_host_fab(h_W), "h_W");
    expect_arrays_equal(h_E_after, to_host_fab(h_E), "h_E");
}

// -------------------------------------------------------------------------
// PPM_reconstruction_y
// -------------------------------------------------------------------------
TEST(PpmReconstructionY, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "ppm_reconstruction_y");

    const auto   bxH        = captured.box("_bxH");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_S        = captured.fab_device("_h_S_before");
    auto         h_N        = captured.fab_device("_h_N_before");
    const auto   mask2d     = captured.fab_device("_mask2d_t");
    const auto   h_S_after  = captured.fab_host("_h_S_after");
    const auto   h_N_after  = captured.fab_host("_h_N_after");
    const double h_min      = captured.real64("_h_min");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::PPM_reconstruction_y(bxH,
                              h_in.const_array(),
                              h_S.array(),
                              h_N.array(),
                              mask2d.const_array(),
                              h_min,
                              monotonic,
                              simple_2nd,
                              /*OBC=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_S_after, to_host_fab(h_S), "h_S");
    expect_arrays_equal(h_N_after, to_host_fab(h_N), "h_N");
}

// -------------------------------------------------------------------------
// ppm_limit_cw84 -- no capture available yet
// -------------------------------------------------------------------------
TEST(PpmLimitCw84, MatchesFortranCapture) {
    GTEST_SKIP() << "no captured ppm_limit_cw84.{bin,meta} fixture yet";
}

// -------------------------------------------------------------------------
// meridional_flux_thickness
// -------------------------------------------------------------------------
TEST(MeridionalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_flux_thickness");

    const auto   bxC            = captured.box("_bxC");
    const auto   v              = captured.fab_device("_v");
    const auto   h              = captured.fab_device("_h");
    const auto   h_S            = captured.fab_device("_h_S");
    const auto   h_N            = captured.fab_device("_h_N");
    auto         h_v            = captured.fab_device("_h_v_before");
    const auto   dx_Cv          = captured.fab_device("_dx_Cv");
    const auto   IareaT         = captured.fab_device("_IareaT");
    const auto   IdyT           = captured.fab_device("_IdyT");
    const auto   por_face_areaV = captured.fab_device("_por_face_areaV");
    const auto   visc_rem_v     = captured.fab_device("_visc_rem_v");
    const auto   h_v_after      = captured.fab_host("_h_v_after");
    const double dt             = captured.real64("_dt");
    const bool   vol_CFL        = captured.logical("_vol_CFL");
    const bool   marginal       = captured.logical("_marginal");

    MOM::meridional_flux_thickness(bxC,
                                   v.const_array(),
                                   h.const_array(),
                                   h_S.const_array(),
                                   h_N.const_array(),
                                   h_v.array(),
                                   dt,
                                   vol_CFL,
                                   marginal,
                                   dx_Cv.const_array(),
                                   IareaT.const_array(),
                                   IdyT.const_array(),
                                   /*obc=*/nullptr,
                                   por_face_areaV.const_array(),
                                   visc_rem_v.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_v_after, to_host_fab(h_v), "h_v");
}

// -------------------------------------------------------------------------
// zonal_flux_thickness
// -------------------------------------------------------------------------
TEST(ZonalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_flux_thickness");

    const auto   bxC            = captured.box("_bxC");
    const auto   u              = captured.fab_device("_u");
    const auto   h              = captured.fab_device("_h");
    const auto   h_W            = captured.fab_device("_h_W");
    const auto   h_E            = captured.fab_device("_h_E");
    auto         h_u            = captured.fab_device("_h_u_before");
    const auto   dy_Cu          = captured.fab_device("_dy_Cu");
    const auto   IareaT         = captured.fab_device("_IareaT");
    const auto   IdxT           = captured.fab_device("_IdxT");
    const auto   por_face_areaU = captured.fab_device("_por_face_areaU");
    const auto   visc_rem_u     = captured.fab_device("_visc_rem_u");
    const auto   h_u_after      = captured.fab_host("_h_u_after");
    const double dt             = captured.real64("_dt");
    const bool   vol_CFL        = captured.logical("_vol_CFL");
    const bool   marginal       = captured.logical("_marginal");

    MOM::zonal_flux_thickness(bxC,
                              u.const_array(),
                              h.const_array(),
                              h_W.const_array(),
                              h_E.const_array(),
                              h_u.array(),
                              dt,
                              vol_CFL,
                              marginal,
                              dy_Cu.const_array(),
                              IareaT.const_array(),
                              IdxT.const_array(),
                              /*obc=*/nullptr,
                              por_face_areaU.const_array(),
                              visc_rem_u.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_u_after, to_host_fab(h_u), "h_u");
}
