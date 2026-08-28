// Unit tests for MOM::ppm_limit_pos / PPM_reconstruction_x / PPM_reconstruction_y.
//
// Each test loads a captured Fortran (input, expected-output) pair from
// <data-dir>/<name>.{bin,meta}, runs the C++ kernel over equivalent AMReX
// containers, and compares the result against the captured "after" arrays.

// SKILLS: 0.3.1

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
// meridional_edge_thickness
// -------------------------------------------------------------------------
//
// No capture/meridional_edge_thickness.{bin,meta} fixture exists yet, so
// this test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1870-1966,
// the meridional_edge_thickness shim's TIMH_capture arm (rec%add(...)
// calls at lines 1914-1922 and 1931-1932). Cross-checked against the
// bind(C) interface (lines 136-153) and against
// turbotmp_meridional_edge_thickness_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all three agree.
//
// MOM::meridional_edge_thickness(bxC, h_in, h_S, h_N, mask2dT, h_min,
//                                upwind_1st, monotonic, simple_2nd, obc)
// internally branches on upwind_1st between a 1st-order upwind copy and a
// full PPM_reconstruction_y call; whichever fixture is eventually captured
// will exercise only whichever branch was live at record time. OBC is
// never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
TEST(MeridionalEdgeThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_edge_thickness");

    const auto   bxC        = captured.box("_bxC");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_S        = captured.fab_device("_h_S_before");
    auto         h_N        = captured.fab_device("_h_N_before");
    const auto   mask2dT    = captured.fab_device("_mask2dT");
    const auto   h_S_after  = captured.fab_host("_h_S_after");
    const auto   h_N_after  = captured.fab_host("_h_N_after");
    const double h_min      = captured.real64("_h_min");
    const bool   upwind_1st = captured.logical("_upwind_1st");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::meridional_edge_thickness(bxC,
                                   h_in.const_array(),
                                   h_S.array(),
                                   h_N.array(),
                                   mask2dT.const_array(),
                                   h_min,
                                   upwind_1st,
                                   monotonic,
                                   simple_2nd,
                                   /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_S_after, to_host_fab(h_S), "h_S");
    expect_arrays_equal(h_N_after, to_host_fab(h_N), "h_N");
}

// -------------------------------------------------------------------------
// zonal_edge_thickness
// -------------------------------------------------------------------------
//
// No capture/zonal_edge_thickness.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1737-1833,
// the zonal_edge_thickness shim's TIMH_capture arm (rec%add(...) calls at
// lines 1781-1789 and 1798-1799). Cross-checked against the bind(C)
// interface (lines 113-131) and against
// turbotmp_zonal_edge_thickness_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all three agree.
//
// MOM::zonal_edge_thickness(bxC, h_in, h_W, h_E, mask2dT, h_min,
//                           upwind_1st, monotonic, simple_2nd, obc)
// internally branches on upwind_1st between a 1st-order upwind copy and a
// full PPM_reconstruction_x call; whichever fixture is eventually captured
// will exercise only whichever branch was live at record time. OBC is
// never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
TEST(ZonalEdgeThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_edge_thickness");

    const auto   bxC        = captured.box("_bxC");
    const auto   h_in       = captured.fab_device("_h_in");
    auto         h_W        = captured.fab_device("_h_W_before");
    auto         h_E        = captured.fab_device("_h_E_before");
    const auto   mask2dT    = captured.fab_device("_mask2dT");
    const auto   h_W_after  = captured.fab_host("_h_W_after");
    const auto   h_E_after  = captured.fab_host("_h_E_after");
    const double h_min      = captured.real64("_h_min");
    const bool   upwind_1st = captured.logical("_upwind_1st");
    const bool   monotonic  = captured.logical("_monotonic");
    const bool   simple_2nd = captured.logical("_simple_2nd");

    MOM::zonal_edge_thickness(bxC,
                              h_in.const_array(),
                              h_W.array(),
                              h_E.array(),
                              mask2dT.const_array(),
                              h_min,
                              upwind_1st,
                              monotonic,
                              simple_2nd,
                              /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_W_after, to_host_fab(h_W), "h_W");
    expect_arrays_equal(h_E_after, to_host_fab(h_E), "h_E");
}

// -------------------------------------------------------------------------
// zonal_flux_thickness
// -------------------------------------------------------------------------
//
// No capture/zonal_flux_thickness.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3023-3142,
// the zonal_flux_thickness shim's TIMH_capture arm (rec%add(...) calls at
// lines 3084-3097 and 3105). Cross-checked against the bind(C) interface
// (lines 337-375) and against turbotmp_zonal_flux_thickness_bridge's
// parameter order in mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all
// three agree.
//
// MOM::zonal_flux_thickness(bxC, u, h, h_W, h_E, h_u, dt, dy_Cu, IareaT,
//                           IdxT, vol_CFL, marginal, obc, por_face_areaU,
//                           visc_rem_u)
// h_W/h_E are pure inputs here (only h_u is modified). visc_rem_u is
// recorded by the Fortran shim only when associated, so a fixture
// captured with it unassociated would be missing the _visc_rem_u field.
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
TEST(ZonalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_flux_thickness");

    const auto   bxC            = captured.box("_bxC");
    const auto   u               = captured.fab_device("_u");
    const auto   h               = captured.fab_device("_h");
    const auto   h_W             = captured.fab_device("_h_W");
    const auto   h_E             = captured.fab_device("_h_E");
    auto         h_u             = captured.fab_device("_h_u_before");
    const auto   h_u_after       = captured.fab_host("_h_u_after");
    const double dt              = captured.real64("_dt");
    const auto   dy_Cu           = captured.fab_device("_dy_Cu");
    const auto   IareaT          = captured.fab_device("_IareaT");
    const auto   IdxT            = captured.fab_device("_IdxT");
    const bool   vol_CFL         = captured.logical("_vol_CFL");
    const bool   marginal        = captured.logical("_marginal");
    const auto   por_face_areaU  = captured.fab_device("_por_face_areaU");
    const auto   visc_rem_u      = captured.fab_device("_visc_rem_u");

    MOM::zonal_flux_thickness(bxC,
                              u.const_array(),
                              h.const_array(),
                              h_W.const_array(),
                              h_E.const_array(),
                              h_u.array(),
                              dt,
                              dy_Cu.const_array(),
                              IareaT.const_array(),
                              IdxT.const_array(),
                              vol_CFL,
                              marginal,
                              /*obc=*/nullptr,
                              por_face_areaU.const_array(),
                              visc_rem_u.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_u_after, to_host_fab(h_u), "h_u");
}

// -------------------------------------------------------------------------
// continuity_meridional_convergence
// -------------------------------------------------------------------------
//
// No capture/continuity_meridional_convergence.{bin,meta} fixture exists
// yet, so this test's field mapping is grounded directly in Fortran source
// (not a .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1626-1700,
// the continuity_meridional_convergence shim's TIMH_capture arm (rec%add(...)
// calls at lines 1665-1671 and 1677). Cross-checked against the bind(C)
// interface (lines 176-193) -- both agree.
//
// MOM::continuity_meridional_convergence(bxC, h, vh, dt, IareaT, hin, h_min)
// advances h in place by the convergence of the meridional thickness flux.
// hin is recorded unconditionally by the shim (unlike visc_rem_u/v in the
// flux-thickness kernels), but the Fortran source itself allows hin to be
// unassociated -- if the captured call site had it that way, this fixture's
// _hin field will be empty and CapturedFile::fab_device("_hin") will throw
// rather than silently misbehave.
TEST(ContinuityMeridionalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_meridional_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   vh      = captured.fab_device("_vh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const auto   hin     = captured.fab_device("_hin");
    const double h_min   = captured.real64("_h_min");

    MOM::continuity_meridional_convergence(bxC,
                                           h.array(),
                                           vh.const_array(),
                                           dt,
                                           IareaT.const_array(),
                                           hin.const_array(),
                                           h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after, to_host_fab(h), "h");
}

// -------------------------------------------------------------------------
// continuity_zonal_convergence
// -------------------------------------------------------------------------
//
// No capture/continuity_zonal_convergence.{bin,meta} fixture exists yet, so
// this test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:1504-1578,
// the continuity_zonal_convergence shim's TIMH_capture arm (rec%add(...)
// calls at lines 1543-1549 and 1555). Cross-checked against the bind(C)
// interface (lines 158-172) -- both agree.
//
// MOM::continuity_zonal_convergence(bxC, h, uh, dt, IareaT, hin, h_min)
// advances h in place by the convergence of the zonal thickness flux. hin
// is recorded unconditionally by the shim (unlike visc_rem_u/v in the
// flux-thickness kernels), but the Fortran source itself allows hin to be
// unassociated -- if the captured call site had it that way, this fixture's
// _hin field will be empty and CapturedFile::fab_device("_hin") will throw
// rather than silently misbehave.
TEST(ContinuityZonalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_zonal_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   uh      = captured.fab_device("_uh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const auto   hin     = captured.fab_device("_hin");
    const double h_min   = captured.real64("_h_min");

    MOM::continuity_zonal_convergence(bxC,
                                      h.array(),
                                      uh.const_array(),
                                      dt,
                                      IareaT.const_array(),
                                      hin.const_array(),
                                      h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after, to_host_fab(h), "h");
}
