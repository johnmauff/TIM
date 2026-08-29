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
    const double h_min   = captured.real64("_h_min");
    // _hin is captured unconditionally but may still be null-encoded
    // (Fortran container unassociated at capture time).
    amrex::FArrayBox hin_fab;
    amrex::Array4<const amrex::Real> hin{};
    if (captured.is_associated("_hin")) {
        hin_fab = captured.fab_device("_hin");
        hin = hin_fab.const_array();
    }

    MOM::continuity_meridional_convergence(bxC,
                                           h.array(),
                                           vh.const_array(),
                                           dt,
                                           IareaT.const_array(),
                                           hin,
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
    const double h_min   = captured.real64("_h_min");
    // _hin is captured unconditionally but may still be null-encoded
    // (Fortran container unassociated at capture time).
    amrex::FArrayBox hin_fab;
    amrex::Array4<const amrex::Real> hin{};
    if (captured.is_associated("_hin")) {
        hin_fab = captured.fab_device("_hin");
        hin = hin_fab.const_array();
    }

    MOM::continuity_zonal_convergence(bxC,
                                      h.array(),
                                      uh.const_array(),
                                      dt,
                                      IareaT.const_array(),
                                      hin,
                                      h_min);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after, to_host_fab(h), "h");
}

// -------------------------------------------------------------------------
// set_merid_BT_cont
// -------------------------------------------------------------------------
//
// No capture/set_merid_bt_cont.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5564-5724, the
// set_merid_BT_cont shim's TIMH_capture arm (rec%add(...) calls at lines
// 5637-5662 and 5671-5676). Cross-checked against the bind(C) interface and
// against turbotmp_set_merid_bt_cont_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all three agree.
//
// MOM::set_merid_BT_cont computes the effective open face areas and
// barotropic-velocity corrections at meridional faces as a function of
// barotropic flow, for use by the barotropic solver's transport-adjustment
// iteration. do_I is captured as a LogicalArray_t, read here via
// int_fab_device() into an amrex::IArrayBox. Only CS.vol_CFL is captured --
// the kernel reads no other transport_adjust_CS_C field, so the rest of CS
// is left default-initialized.
TEST(SetMeridBtCont, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "set_merid_bt_cont");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    auto         FA_v_S0          = captured.fab_device("_FA_v_S0_before");
    auto         FA_v_N0          = captured.fab_device("_FA_v_N0_before");
    auto         FA_v_SS          = captured.fab_device("_FA_v_SS_before");
    auto         FA_v_NN          = captured.fab_device("_FA_v_NN_before");
    auto         vBT_SS           = captured.fab_device("_vBT_SS_before");
    auto         vBT_NN           = captured.fab_device("_vBT_NN_before");
    const auto   FA_v_S0_after    = captured.fab_host("_FA_v_S0_after");
    const auto   FA_v_N0_after    = captured.fab_host("_FA_v_N0_after");
    const auto   FA_v_SS_after    = captured.fab_host("_FA_v_SS_after");
    const auto   FA_v_NN_after    = captured.fab_host("_FA_v_NN_after");
    const auto   vBT_SS_after     = captured.fab_host("_vBT_SS_after");
    const auto   vBT_NN_after     = captured.fab_host("_vBT_NN_after");
    const auto   dv0              = captured.fab_device("_dv0");
    const auto   vh_tot_0         = captured.fab_device("_vh_tot_0");
    const auto   dvhdv_tot_0      = captured.fab_device("_dvhdv_tot_0");
    const auto   dv_max_CFL       = captured.fab_device("_dv_max_CFL");
    const auto   dv_min_CFL       = captured.fab_device("_dv_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dyCv             = captured.fab_device("_dyCv");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   visc_rem_max     = captured.fab_device("_visc_rem_max");
    const auto   do_I             = captured.int_fab_device("_do_I");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");

    MOM::set_merid_BT_cont(bxC,
                           v.const_array(),
                           h_in.const_array(),
                           h_S.const_array(),
                           h_N.const_array(),
                           FA_v_S0.array(),
                           FA_v_N0.array(),
                           FA_v_SS.array(),
                           FA_v_NN.array(),
                           vBT_SS.array(),
                           vBT_NN.array(),
                           dv0.const_array(),
                           vh_tot_0.const_array(),
                           dvhdv_tot_0.const_array(),
                           dv_max_CFL.const_array(),
                           dv_min_CFL.const_array(),
                           dt,
                           dyCv.const_array(),
                           dx_Cv.const_array(),
                           IareaT.const_array(),
                           IdyT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           visc_rem_max.const_array(),
                           do_I.const_array(),
                           por_face_areaV.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(FA_v_S0_after, to_host_fab(FA_v_S0), "FA_v_S0");
    expect_arrays_equal(FA_v_N0_after, to_host_fab(FA_v_N0), "FA_v_N0");
    expect_arrays_equal(FA_v_SS_after, to_host_fab(FA_v_SS), "FA_v_SS");
    expect_arrays_equal(FA_v_NN_after, to_host_fab(FA_v_NN), "FA_v_NN");
    expect_arrays_equal(vBT_SS_after,  to_host_fab(vBT_SS),  "vBT_SS");
    expect_arrays_equal(vBT_NN_after,  to_host_fab(vBT_NN),  "vBT_NN");
}

// -------------------------------------------------------------------------
// set_zonal_BT_cont
// -------------------------------------------------------------------------
//
// No capture/set_zonal_bt_cont.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3724-3883, the
// set_zonal_BT_cont shim's TIMH_capture arm (rec%add(...) calls at lines
// 3796-3821 and 3830-3835). Cross-checked against the bind(C) interface and
// against turbotmp_set_zonal_bt_cont_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all three agree.
//
// MOM::set_zonal_BT_cont computes the effective open face areas and
// barotropic-velocity corrections at zonal faces as a function of
// barotropic flow, for use by the barotropic solver's transport-adjustment
// iteration. do_I is captured as a LogicalArray_t, read here via
// int_fab_device() into an amrex::IArrayBox. Only CS.vol_CFL is captured --
// the kernel reads no other transport_adjust_CS_C field, so the rest of CS
// is left default-initialized.
TEST(SetZonalBtCont, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "set_zonal_bt_cont");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    auto         FA_u_W0          = captured.fab_device("_FA_u_W0_before");
    auto         FA_u_E0          = captured.fab_device("_FA_u_E0_before");
    auto         FA_u_WW          = captured.fab_device("_FA_u_WW_before");
    auto         FA_u_EE          = captured.fab_device("_FA_u_EE_before");
    auto         uBT_WW           = captured.fab_device("_uBT_WW_before");
    auto         uBT_EE           = captured.fab_device("_uBT_EE_before");
    const auto   FA_u_W0_after    = captured.fab_host("_FA_u_W0_after");
    const auto   FA_u_E0_after    = captured.fab_host("_FA_u_E0_after");
    const auto   FA_u_WW_after    = captured.fab_host("_FA_u_WW_after");
    const auto   FA_u_EE_after    = captured.fab_host("_FA_u_EE_after");
    const auto   uBT_WW_after     = captured.fab_host("_uBT_WW_after");
    const auto   uBT_EE_after     = captured.fab_host("_uBT_EE_after");
    const auto   du0              = captured.fab_device("_du0");
    const auto   uh_tot_0         = captured.fab_device("_uh_tot_0");
    const auto   duhdu_tot_0      = captured.fab_device("_duhdu_tot_0");
    const auto   du_max_CFL       = captured.fab_device("_du_max_CFL");
    const auto   du_min_CFL       = captured.fab_device("_du_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dxCu             = captured.fab_device("_dxCu");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    transport_adjust_CS_C CS{};
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   visc_rem_max     = captured.fab_device("_visc_rem_max");
    const auto   do_I             = captured.int_fab_device("_do_I");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");

    MOM::set_zonal_BT_cont(bxC,
                           u.const_array(),
                           h_in.const_array(),
                           h_W.const_array(),
                           h_E.const_array(),
                           FA_u_W0.array(),
                           FA_u_E0.array(),
                           FA_u_WW.array(),
                           FA_u_EE.array(),
                           uBT_WW.array(),
                           uBT_EE.array(),
                           du0.const_array(),
                           uh_tot_0.const_array(),
                           duhdu_tot_0.const_array(),
                           du_max_CFL.const_array(),
                           du_min_CFL.const_array(),
                           dt,
                           dxCu.const_array(),
                           dy_Cu.const_array(),
                           IareaT.const_array(),
                           IdxT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           visc_rem_max.const_array(),
                           do_I.const_array(),
                           por_face_areaU.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(FA_u_W0_after, to_host_fab(FA_u_W0), "FA_u_W0");
    expect_arrays_equal(FA_u_E0_after, to_host_fab(FA_u_E0), "FA_u_E0");
    expect_arrays_equal(FA_u_WW_after, to_host_fab(FA_u_WW), "FA_u_WW");
    expect_arrays_equal(FA_u_EE_after, to_host_fab(FA_u_EE), "FA_u_EE");
    expect_arrays_equal(uBT_WW_after,  to_host_fab(uBT_WW),  "uBT_WW");
    expect_arrays_equal(uBT_EE_after,  to_host_fab(uBT_EE),  "uBT_EE");
}

// -------------------------------------------------------------------------
// meridional_flux_adjust
// -------------------------------------------------------------------------
//
// No capture/meridional_flux_adjust.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5210-5362, the
// meridional_flux_adjust shim's TIMH_capture arm (rec%add(...) calls at
// lines 5283-5305 and 5313-5314). Cross-checked against the bind(C)
// interface and against turbotmp_meridional_flux_adjust_bridge's parameter
// order in mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all agree.
//
// MOM::meridional_flux_adjust Newton-iterates a barotropic velocity
// correction per meridional face so that the vertically-summed meridional
// mass/volume transport matches the target barotropic transport, to within
// the transport-adjustment iteration's tolerance. _tol_eta/_tol_vel/
// _better_iter/_vol_CFL are the only transport_adjust_CS_C fields this
// kernel reads, and are the only ones the shim captures; the rest of CS is
// left default-initialized. do_I_in is captured as a LogicalArray_t, read
// via int_fab_device(). vhbt/vh_3d are captured only when associated at
// record time -- a fixture recorded with either unassociated would be
// missing the corresponding field(s). obc is never captured -- pass
// nullptr, matching the existing PPM_reconstruction_x/_y tests
// (OBC-inactive configs only).
TEST(MeridionalFluxAdjust, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_flux_adjust");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    const auto   vh_tot_0         = captured.fab_device("_vh_tot_0");
    const auto   dvhdv_tot_0      = captured.fab_device("_dvhdv_tot_0");
    auto         dv               = captured.fab_device("_dv_before");
    const auto   dv_after         = captured.fab_host("_dv_after");
    const auto   dv_max_CFL       = captured.fab_device("_dv_max_CFL");
    const auto   dv_min_CFL       = captured.fab_device("_dv_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    transport_adjust_CS_C CS{};
    CS.tol_eta                    = captured.real64("_tol_eta");
    CS.tol_vel                    = captured.real64("_tol_vel");
    CS.better_iter                = captured.logical("_better_iter");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   do_I_in          = captured.int_fab_device("_do_I_in");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");
    const auto   vhbt             = captured.fab_device("_vhbt");
    auto         vh_3d            = captured.fab_device("_vh_3d_before");
    const auto   vh_3d_after      = captured.fab_host("_vh_3d_after");

    MOM::meridional_flux_adjust(bxC,
                                v.const_array(),
                                h_in.const_array(),
                                h_S.const_array(),
                                h_N.const_array(),
                                vh_tot_0.const_array(),
                                dvhdv_tot_0.const_array(),
                                dv.array(),
                                dv_max_CFL.const_array(),
                                dv_min_CFL.const_array(),
                                dt,
                                dx_Cv.const_array(),
                                IareaT.const_array(),
                                IdyT.const_array(),
                                CS,
                                visc_rem.const_array(),
                                do_I_in.const_array(),
                                por_face_areaV.const_array(),
                                vhbt.const_array(),
                                vh_3d.array(),
                                /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(dv_after,     to_host_fab(dv),     "dv");
    expect_arrays_equal(vh_3d_after,  to_host_fab(vh_3d),  "vh_3d");
}

// -------------------------------------------------------------------------
// zonal_flux_adjust
// -------------------------------------------------------------------------
//
// No capture/zonal_flux_adjust.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:3371-3526, the
// zonal_flux_adjust shim's TIMH_capture arm (rec%add(...) calls at lines
// 3447-3469 and 3477-3478). Cross-checked against the bind(C) interface and
// against turbotmp_zonal_flux_adjust_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all agree.
//
// MOM::zonal_flux_adjust Newton-iterates a barotropic velocity correction
// per zonal face so that the vertically-summed zonal mass/volume transport
// matches the target barotropic transport, to within the
// transport-adjustment iteration's tolerance. _tol_eta/_tol_vel/
// _better_iter/_vol_CFL are the only transport_adjust_CS_C fields this
// kernel reads, and are the only ones the shim captures; the rest of CS is
// left default-initialized. do_I_in is captured as a LogicalArray_t, read
// via int_fab_device(). uhbt/uh_3d are captured only when associated at
// record time -- a fixture recorded with either unassociated would be
// missing the corresponding field(s). obc is never captured -- pass
// nullptr, matching the existing PPM_reconstruction_x/_y tests
// (OBC-inactive configs only).
TEST(ZonalFluxAdjust, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_flux_adjust");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    const auto   uh_tot_0         = captured.fab_device("_uh_tot_0");
    const auto   duhdu_tot_0      = captured.fab_device("_duhdu_tot_0");
    auto         du               = captured.fab_device("_du_before");
    const auto   du_after         = captured.fab_host("_du_after");
    const auto   du_max_CFL       = captured.fab_device("_du_max_CFL");
    const auto   du_min_CFL       = captured.fab_device("_du_min_CFL");
    const double dt               = captured.real64("_dt");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    transport_adjust_CS_C CS{};
    CS.tol_eta                    = captured.real64("_tol_eta");
    CS.tol_vel                    = captured.real64("_tol_vel");
    CS.better_iter                = captured.logical("_better_iter");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    const auto   visc_rem         = captured.fab_device("_visc_rem");
    const auto   do_I_in          = captured.int_fab_device("_do_I_in");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");
    const auto   uhbt             = captured.fab_device("_uhbt");
    auto         uh_3d            = captured.fab_device("_uh_3d_before");
    const auto   uh_3d_after      = captured.fab_host("_uh_3d_after");

    MOM::zonal_flux_adjust(bxC,
                           u.const_array(),
                           h_in.const_array(),
                           h_W.const_array(),
                           h_E.const_array(),
                           uh_tot_0.const_array(),
                           duhdu_tot_0.const_array(),
                           du.array(),
                           du_max_CFL.const_array(),
                           du_min_CFL.const_array(),
                           dt,
                           dy_Cu.const_array(),
                           IareaT.const_array(),
                           IdxT.const_array(),
                           CS,
                           visc_rem.const_array(),
                           do_I_in.const_array(),
                           por_face_areaU.const_array(),
                           uhbt.const_array(),
                           uh_3d.array(),
                           /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(du_after,    to_host_fab(du),    "du");
    expect_arrays_equal(uh_3d_after, to_host_fab(uh_3d), "uh_3d");
}

// -------------------------------------------------------------------------
// meridional_mass_flux
// -------------------------------------------------------------------------
//
// No capture/meridional_mass_flux.{bin,meta} fixture exists yet, so this
// test's field mapping is grounded directly in Fortran source (not a .meta
// file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:4290-4502, the
// meridional_mass_flux shim's TIMH_capture arm (rec%add(...) calls at lines
// 4383-4417 and 4427-4437). Cross-checked against the bind(C) interface and
// against turbotmp_meridional_mass_flux_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all agree.
//
// MOM::meridional_mass_flux is the orchestrator that computes the
// meridional PPM transport vh, then -- when vhbt/BT_cont output is
// requested -- the transport-adjustment correction via
// meridional_flux_adjust and set_merid_BT_cont. _isd/_ied are Fortran
// 1-based data-domain i-bounds; the bridge normally converts these to
// 0-based before calling the kernel, so this test does the same
// subtraction directly. _vhbt/_visc_rem_v are captured unconditionally by
// the shim (unlike v_cor/dv_cor/the six FA_v_*/vBT_* fields, which are
// captured only when associated). _CFL_limit_adjust/_aggress_adjust/
// _vol_CFL/_use_visc_rem_max/_marginal_faces are the transport_adjust_CS_C
// fields the shim captures; the kernel itself does not read marginal_faces
// (that would only matter for the BT_cont%h_v/meridional_flux_thickness
// tail, which is out of scope -- see the kernel's own source comment). obc
// is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
TEST(MeridionalMassFlux, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_mass_flux");

    const auto   bxC             = captured.box("_bxC");
    const auto   v                = captured.fab_device("_v");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_S              = captured.fab_device("_h_S");
    const auto   h_N              = captured.fab_device("_h_N");
    auto         vh                = captured.fab_device("_vh_before");
    const auto   vh_after          = captured.fab_host("_vh_after");
    const double dt               = captured.real64("_dt");
    const auto   dx_Cv            = captured.fab_device("_dx_Cv");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdyT             = captured.fab_device("_IdyT");
    const auto   areaT            = captured.fab_device("_areaT");
    const auto   dyT              = captured.fab_device("_dyT");
    const auto   mask2dCv         = captured.fab_device("_mask2dCv");
    const auto   dyCv             = captured.fab_device("_dyCv");
    const int    isd              = captured.integer("_isd") - 1;
    const int    ied              = captured.integer("_ied") - 1;
    const double H_subroundoff    = captured.real64("_H_subroundoff");
    transport_adjust_CS_C CS{};
    CS.CFL_limit_adjust           = captured.real64("_CFL_limit_adjust");
    CS.aggress_adjust             = captured.logical("_aggress_adjust");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    CS.use_visc_rem_max           = captured.logical("_use_visc_rem_max");
    CS.marginal_faces             = captured.logical("_marginal_faces");
    const auto   por_face_areaV   = captured.fab_device("_por_face_areaV");
    // _vhbt/_visc_rem_v are captured unconditionally but may still be
    // null-encoded (Fortran container unassociated at capture time).
    amrex::FArrayBox vhbt_fab, visc_rem_v_fab;
    amrex::Array4<const amrex::Real> vhbt{}, visc_rem_v{};
    if (captured.is_associated("_vhbt")) {
        vhbt_fab = captured.fab_device("_vhbt");
        vhbt = vhbt_fab.const_array();
    }
    if (captured.is_associated("_visc_rem_v")) {
        visc_rem_v_fab = captured.fab_device("_visc_rem_v");
        visc_rem_v = visc_rem_v_fab.const_array();
    }
    auto         v_cor             = captured.fab_device("_v_cor_before");
    const auto   v_cor_after       = captured.fab_host("_v_cor_after");
    auto         FA_v_S0           = captured.fab_device("_FA_v_S0_before");
    auto         FA_v_N0           = captured.fab_device("_FA_v_N0_before");
    auto         FA_v_SS           = captured.fab_device("_FA_v_SS_before");
    auto         FA_v_NN           = captured.fab_device("_FA_v_NN_before");
    auto         vBT_SS            = captured.fab_device("_vBT_SS_before");
    auto         vBT_NN            = captured.fab_device("_vBT_NN_before");
    const auto   FA_v_S0_after     = captured.fab_host("_FA_v_S0_after");
    const auto   FA_v_N0_after     = captured.fab_host("_FA_v_N0_after");
    const auto   FA_v_SS_after     = captured.fab_host("_FA_v_SS_after");
    const auto   FA_v_NN_after     = captured.fab_host("_FA_v_NN_after");
    const auto   vBT_SS_after      = captured.fab_host("_vBT_SS_after");
    const auto   vBT_NN_after      = captured.fab_host("_vBT_NN_after");
    auto         dv_cor            = captured.fab_device("_dv_cor_before");
    const auto   dv_cor_after      = captured.fab_host("_dv_cor_after");

    MOM::meridional_mass_flux(bxC,
                              v.const_array(),
                              h_in.const_array(),
                              h_S.const_array(),
                              h_N.const_array(),
                              vh.array(),
                              dt,
                              dx_Cv.const_array(),
                              IareaT.const_array(),
                              IdyT.const_array(),
                              areaT.const_array(),
                              dyT.const_array(),
                              mask2dCv.const_array(),
                              dyCv.const_array(),
                              isd,
                              ied,
                              H_subroundoff,
                              CS,
                              /*obc=*/nullptr,
                              por_face_areaV.const_array(),
                              vhbt,
                              visc_rem_v,
                              v_cor.array(),
                              FA_v_S0.array(),
                              FA_v_N0.array(),
                              FA_v_SS.array(),
                              FA_v_NN.array(),
                              vBT_SS.array(),
                              vBT_NN.array(),
                              dv_cor.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(vh_after,     to_host_fab(vh),     "vh");
    expect_arrays_equal(v_cor_after,  to_host_fab(v_cor),  "v_cor");
    expect_arrays_equal(FA_v_S0_after, to_host_fab(FA_v_S0), "FA_v_S0");
    expect_arrays_equal(FA_v_N0_after, to_host_fab(FA_v_N0), "FA_v_N0");
    expect_arrays_equal(FA_v_SS_after, to_host_fab(FA_v_SS), "FA_v_SS");
    expect_arrays_equal(FA_v_NN_after, to_host_fab(FA_v_NN), "FA_v_NN");
    expect_arrays_equal(vBT_SS_after,  to_host_fab(vBT_SS),  "vBT_SS");
    expect_arrays_equal(vBT_NN_after,  to_host_fab(vBT_NN),  "vBT_NN");
    expect_arrays_equal(dv_cor_after, to_host_fab(dv_cor), "dv_cor");
}

// -------------------------------------------------------------------------
// zonal_mass_flux
// -------------------------------------------------------------------------
//
// No capture/zonal_mass_flux.{bin,meta} fixture exists yet, so this test's
// field mapping is grounded directly in Fortran source (not a .meta file):
// submodules/MOM6/src/core/MOM_continuity_PPM.F90:2371-2570, the
// zonal_mass_flux shim's TIMH_capture arm (rec%add(...) calls at lines
// 2455-2487 and 2496-2506). Cross-checked against the bind(C) interface
// and against turbotmp_zonal_mass_flux_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all agree.
//
// MOM::zonal_mass_flux is the orchestrator that computes the zonal PPM
// transport uh, then -- when uhbt/BT_cont output is requested -- the
// transport-adjustment correction via zonal_flux_adjust and
// set_zonal_BT_cont. _uhbt/_visc_rem_u are captured unconditionally by
// the shim (unlike u_cor/du_cor/the six FA_u_*/uBT_* fields, which are
// captured only when associated). _CFL_limit_adjust/_aggress_adjust/
// _vol_CFL/_use_visc_rem_max/_marginal_faces are the transport_adjust_CS_C
// fields the shim captures; the kernel itself does not read
// marginal_faces (that would only matter for the BT_cont%h_u/
// zonal_flux_thickness tail, which is out of scope -- see the kernel's
// own source comment). obc is never captured -- pass nullptr, matching
// the existing PPM_reconstruction_x/_y tests (OBC-inactive configs only).
TEST(ZonalMassFlux, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "zonal_mass_flux");

    const auto   bxC             = captured.box("_bxC");
    const auto   u                = captured.fab_device("_u");
    const auto   h_in             = captured.fab_device("_h_in");
    const auto   h_W              = captured.fab_device("_h_W");
    const auto   h_E              = captured.fab_device("_h_E");
    auto         uh                = captured.fab_device("_uh_before");
    const auto   uh_after          = captured.fab_host("_uh_after");
    const double dt               = captured.real64("_dt");
    const auto   dy_Cu            = captured.fab_device("_dy_Cu");
    const auto   IareaT           = captured.fab_device("_IareaT");
    const auto   IdxT             = captured.fab_device("_IdxT");
    const auto   areaT            = captured.fab_device("_areaT");
    const auto   dxT              = captured.fab_device("_dxT");
    const auto   mask2dCu         = captured.fab_device("_mask2dCu");
    const auto   dxCu             = captured.fab_device("_dxCu");
    const double H_subroundoff    = captured.real64("_H_subroundoff");
    transport_adjust_CS_C CS{};
    CS.CFL_limit_adjust           = captured.real64("_CFL_limit_adjust");
    CS.aggress_adjust             = captured.logical("_aggress_adjust");
    CS.vol_CFL                    = captured.logical("_vol_CFL");
    CS.use_visc_rem_max           = captured.logical("_use_visc_rem_max");
    CS.marginal_faces             = captured.logical("_marginal_faces");
    const auto   por_face_areaU   = captured.fab_device("_por_face_areaU");
    // _uhbt/_visc_rem_u are captured unconditionally but may still be
    // null-encoded (Fortran container unassociated at capture time).
    amrex::FArrayBox uhbt_fab, visc_rem_u_fab;
    amrex::Array4<const amrex::Real> uhbt{}, visc_rem_u{};
    if (captured.is_associated("_uhbt")) {
        uhbt_fab = captured.fab_device("_uhbt");
        uhbt = uhbt_fab.const_array();
    }
    if (captured.is_associated("_visc_rem_u")) {
        visc_rem_u_fab = captured.fab_device("_visc_rem_u");
        visc_rem_u = visc_rem_u_fab.const_array();
    }
    auto         u_cor             = captured.fab_device("_u_cor_before");
    const auto   u_cor_after       = captured.fab_host("_u_cor_after");
    auto         FA_u_W0           = captured.fab_device("_FA_u_W0_before");
    auto         FA_u_E0           = captured.fab_device("_FA_u_E0_before");
    auto         FA_u_WW           = captured.fab_device("_FA_u_WW_before");
    auto         FA_u_EE           = captured.fab_device("_FA_u_EE_before");
    auto         uBT_WW            = captured.fab_device("_uBT_WW_before");
    auto         uBT_EE            = captured.fab_device("_uBT_EE_before");
    const auto   FA_u_W0_after     = captured.fab_host("_FA_u_W0_after");
    const auto   FA_u_E0_after     = captured.fab_host("_FA_u_E0_after");
    const auto   FA_u_WW_after     = captured.fab_host("_FA_u_WW_after");
    const auto   FA_u_EE_after     = captured.fab_host("_FA_u_EE_after");
    const auto   uBT_WW_after      = captured.fab_host("_uBT_WW_after");
    const auto   uBT_EE_after      = captured.fab_host("_uBT_EE_after");
    auto         du_cor            = captured.fab_device("_du_cor_before");
    const auto   du_cor_after      = captured.fab_host("_du_cor_after");

    MOM::zonal_mass_flux(bxC,
                         u.const_array(),
                         h_in.const_array(),
                         h_W.const_array(),
                         h_E.const_array(),
                         uh.array(),
                         dt,
                         dy_Cu.const_array(),
                         IareaT.const_array(),
                         IdxT.const_array(),
                         areaT.const_array(),
                         dxT.const_array(),
                         mask2dCu.const_array(),
                         dxCu.const_array(),
                         H_subroundoff,
                         CS,
                         /*obc=*/nullptr,
                         por_face_areaU.const_array(),
                         uhbt,
                         visc_rem_u,
                         u_cor.array(),
                         FA_u_W0.array(),
                         FA_u_E0.array(),
                         FA_u_WW.array(),
                         FA_u_EE.array(),
                         uBT_WW.array(),
                         uBT_EE.array(),
                         du_cor.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(uh_after,     to_host_fab(uh),     "uh");
    expect_arrays_equal(u_cor_after,  to_host_fab(u_cor),  "u_cor");
    expect_arrays_equal(FA_u_W0_after, to_host_fab(FA_u_W0), "FA_u_W0");
    expect_arrays_equal(FA_u_E0_after, to_host_fab(FA_u_E0), "FA_u_E0");
    expect_arrays_equal(FA_u_WW_after, to_host_fab(FA_u_WW), "FA_u_WW");
    expect_arrays_equal(FA_u_EE_after, to_host_fab(FA_u_EE), "FA_u_EE");
    expect_arrays_equal(uBT_WW_after,  to_host_fab(uBT_WW),  "uBT_WW");
    expect_arrays_equal(uBT_EE_after,  to_host_fab(uBT_EE),  "uBT_EE");
    expect_arrays_equal(du_cor_after, to_host_fab(du_cor), "du_cor");
}

// -------------------------------------------------------------------------
// continuity_PPM
// -------------------------------------------------------------------------
//
// No capture/continuity_ppm.{bin,meta} fixture exists yet, so this test's
// field mapping is grounded directly in Fortran source (not a .meta file):
// submodules/MOM6/src/core/MOM_continuity_PPM.F90:1148-1467, the
// continuity_PPM shim's TIMH_capture arm (rec%add(...) calls at lines
// 1278-1339 and 1350-1370). Cross-checked against the bind(C) interface
// and against turbotmp_continuity_ppm_bridge's parameter order in
// mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp -- all agree.
//
// MOM::continuity_PPM is the monolithic continuity solver: reconstructs
// edge thicknesses, then advects first in one direction and then the
// other (order set by x_first), via zonal_mass_flux/meridional_mass_flux
// and continuity_zonal_convergence/continuity_meridional_convergence.
// _isd/_ied are Fortran 1-based data-domain i-bounds; the bridge normally
// converts these to 0-based before calling the kernel, so this test does
// the same subtraction directly. Every field of both reconstruction_CS_C
// (upwind_1st/monotonic/simple_2nd) and transport_adjust_CS_C (the
// remaining 8 fields) is captured here, unlike some of the leaf kernels'
// fixtures. uhbt/vhbt/visc_rem_u/visc_rem_v are captured unconditionally;
// u_cor/v_cor/du_cor/dv_cor and the twelve FA_u_*/uBT_*/FA_v_*/vBT_*
// fields are captured only when associated. obc is never captured -- pass
// nullptr, matching the existing PPM_reconstruction_x/_y tests
// (OBC-inactive configs only).
TEST(ContinuityPPM, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_ppm");

    const auto   u                 = captured.fab_device("_u");
    const auto   v                 = captured.fab_device("_v");
    const auto   hin               = captured.fab_device("_hin");
    auto         h                  = captured.fab_device("_h_before");
    const auto   h_after            = captured.fab_host("_h_after");
    auto         uh                 = captured.fab_device("_uh_before");
    const auto   uh_after           = captured.fab_host("_uh_after");
    auto         vh                 = captured.fab_device("_vh_before");
    const auto   vh_after           = captured.fab_host("_vh_after");
    const double dt                = captured.real64("_dt");
    const auto   bx0                = captured.box("_bx0");
    const int    stencil           = captured.integer("_stencil");
    const bool   x_first           = captured.logical("_x_first");
    const auto   mask2dT           = captured.fab_device("_mask2dT");
    const auto   dy_Cu             = captured.fab_device("_dy_Cu");
    const auto   IareaT            = captured.fab_device("_IareaT");
    const auto   IdxT              = captured.fab_device("_IdxT");
    const auto   areaT             = captured.fab_device("_areaT");
    const auto   dxT               = captured.fab_device("_dxT");
    const auto   mask2dCu          = captured.fab_device("_mask2dCu");
    const auto   dxCu              = captured.fab_device("_dxCu");
    const auto   dx_Cv             = captured.fab_device("_dx_Cv");
    const auto   IdyT              = captured.fab_device("_IdyT");
    const auto   dyT               = captured.fab_device("_dyT");
    const auto   mask2dCv          = captured.fab_device("_mask2dCv");
    const auto   dyCv              = captured.fab_device("_dyCv");
    const int    isd               = captured.integer("_isd") - 1;
    const int    ied               = captured.integer("_ied") - 1;
    const double Angstrom_H        = captured.real64("_Angstrom_H");
    const double H_subroundoff     = captured.real64("_H_subroundoff");
    reconstruction_CS_C reconstruction_CS{};
    reconstruction_CS.upwind_1st   = captured.logical("_upwind_1st");
    reconstruction_CS.monotonic    = captured.logical("_monotonic");
    reconstruction_CS.simple_2nd   = captured.logical("_simple_2nd");
    transport_adjust_CS_C transport_adjust_CS{};
    transport_adjust_CS.tol_eta          = captured.real64("_tol_eta");
    transport_adjust_CS.tol_vel          = captured.real64("_tol_vel");
    transport_adjust_CS.CFL_limit_adjust = captured.real64("_CFL_limit_adjust");
    transport_adjust_CS.aggress_adjust   = captured.logical("_aggress_adjust");
    transport_adjust_CS.vol_CFL          = captured.logical("_vol_CFL");
    transport_adjust_CS.better_iter      = captured.logical("_better_iter");
    transport_adjust_CS.use_visc_rem_max = captured.logical("_use_visc_rem_max");
    transport_adjust_CS.marginal_faces   = captured.logical("_marginal_faces");
    const auto   por_face_areaU    = captured.fab_device("_por_face_areaU");
    const auto   por_face_areaV    = captured.fab_device("_por_face_areaV");
    // _uhbt/_vhbt/_visc_rem_u/_visc_rem_v are captured unconditionally but
    // may still be null-encoded (Fortran container unassociated at capture
    // time).
    amrex::FArrayBox uhbt_fab, vhbt_fab, visc_rem_u_fab, visc_rem_v_fab;
    amrex::Array4<const amrex::Real> uhbt{}, vhbt{}, visc_rem_u{}, visc_rem_v{};
    if (captured.is_associated("_uhbt")) {
        uhbt_fab = captured.fab_device("_uhbt");
        uhbt = uhbt_fab.const_array();
    }
    if (captured.is_associated("_vhbt")) {
        vhbt_fab = captured.fab_device("_vhbt");
        vhbt = vhbt_fab.const_array();
    }
    if (captured.is_associated("_visc_rem_u")) {
        visc_rem_u_fab = captured.fab_device("_visc_rem_u");
        visc_rem_u = visc_rem_u_fab.const_array();
    }
    if (captured.is_associated("_visc_rem_v")) {
        visc_rem_v_fab = captured.fab_device("_visc_rem_v");
        visc_rem_v = visc_rem_v_fab.const_array();
    }
    auto         u_cor              = captured.fab_device("_u_cor_before");
    const auto   u_cor_after        = captured.fab_host("_u_cor_after");
    auto         v_cor              = captured.fab_device("_v_cor_before");
    const auto   v_cor_after        = captured.fab_host("_v_cor_after");
    auto         FA_u_W0            = captured.fab_device("_FA_u_W0_before");
    auto         FA_u_E0            = captured.fab_device("_FA_u_E0_before");
    auto         FA_u_WW            = captured.fab_device("_FA_u_WW_before");
    auto         FA_u_EE            = captured.fab_device("_FA_u_EE_before");
    auto         uBT_WW             = captured.fab_device("_uBT_WW_before");
    auto         uBT_EE             = captured.fab_device("_uBT_EE_before");
    const auto   FA_u_W0_after      = captured.fab_host("_FA_u_W0_after");
    const auto   FA_u_E0_after      = captured.fab_host("_FA_u_E0_after");
    const auto   FA_u_WW_after      = captured.fab_host("_FA_u_WW_after");
    const auto   FA_u_EE_after      = captured.fab_host("_FA_u_EE_after");
    const auto   uBT_WW_after       = captured.fab_host("_uBT_WW_after");
    const auto   uBT_EE_after       = captured.fab_host("_uBT_EE_after");
    auto         FA_v_S0            = captured.fab_device("_FA_v_S0_before");
    auto         FA_v_N0            = captured.fab_device("_FA_v_N0_before");
    auto         FA_v_SS            = captured.fab_device("_FA_v_SS_before");
    auto         FA_v_NN            = captured.fab_device("_FA_v_NN_before");
    auto         vBT_SS             = captured.fab_device("_vBT_SS_before");
    auto         vBT_NN             = captured.fab_device("_vBT_NN_before");
    const auto   FA_v_S0_after      = captured.fab_host("_FA_v_S0_after");
    const auto   FA_v_N0_after      = captured.fab_host("_FA_v_N0_after");
    const auto   FA_v_SS_after      = captured.fab_host("_FA_v_SS_after");
    const auto   FA_v_NN_after      = captured.fab_host("_FA_v_NN_after");
    const auto   vBT_SS_after       = captured.fab_host("_vBT_SS_after");
    const auto   vBT_NN_after       = captured.fab_host("_vBT_NN_after");
    auto         du_cor             = captured.fab_device("_du_cor_before");
    const auto   du_cor_after       = captured.fab_host("_du_cor_after");
    auto         dv_cor             = captured.fab_device("_dv_cor_before");
    const auto   dv_cor_after       = captured.fab_host("_dv_cor_after");

    MOM::continuity_PPM(u.const_array(),
                        v.const_array(),
                        hin.const_array(),
                        h.array(),
                        uh.array(),
                        vh.array(),
                        dt,
                        bx0,
                        stencil,
                        x_first,
                        mask2dT.const_array(),
                        dy_Cu.const_array(),
                        IareaT.const_array(),
                        IdxT.const_array(),
                        areaT.const_array(),
                        dxT.const_array(),
                        mask2dCu.const_array(),
                        dxCu.const_array(),
                        dx_Cv.const_array(),
                        IdyT.const_array(),
                        dyT.const_array(),
                        mask2dCv.const_array(),
                        dyCv.const_array(),
                        isd,
                        ied,
                        Angstrom_H,
                        H_subroundoff,
                        reconstruction_CS,
                        transport_adjust_CS,
                        /*obc=*/nullptr,
                        por_face_areaU.const_array(),
                        por_face_areaV.const_array(),
                        uhbt,
                        vhbt,
                        visc_rem_u,
                        visc_rem_v,
                        u_cor.array(),
                        v_cor.array(),
                        FA_u_W0.array(),
                        FA_u_E0.array(),
                        FA_u_WW.array(),
                        FA_u_EE.array(),
                        uBT_WW.array(),
                        uBT_EE.array(),
                        FA_v_S0.array(),
                        FA_v_N0.array(),
                        FA_v_SS.array(),
                        FA_v_NN.array(),
                        vBT_SS.array(),
                        vBT_NN.array(),
                        du_cor.array(),
                        dv_cor.array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after,  to_host_fab(h),  "h");
    expect_arrays_equal(uh_after, to_host_fab(uh), "uh");
    expect_arrays_equal(vh_after, to_host_fab(vh), "vh");
    expect_arrays_equal(u_cor_after, to_host_fab(u_cor), "u_cor");
    expect_arrays_equal(v_cor_after, to_host_fab(v_cor), "v_cor");
    expect_arrays_equal(FA_u_W0_after, to_host_fab(FA_u_W0), "FA_u_W0");
    expect_arrays_equal(FA_u_E0_after, to_host_fab(FA_u_E0), "FA_u_E0");
    expect_arrays_equal(FA_u_WW_after, to_host_fab(FA_u_WW), "FA_u_WW");
    expect_arrays_equal(FA_u_EE_after, to_host_fab(FA_u_EE), "FA_u_EE");
    expect_arrays_equal(uBT_WW_after,  to_host_fab(uBT_WW),  "uBT_WW");
    expect_arrays_equal(uBT_EE_after,  to_host_fab(uBT_EE),  "uBT_EE");
    expect_arrays_equal(FA_v_S0_after, to_host_fab(FA_v_S0), "FA_v_S0");
    expect_arrays_equal(FA_v_N0_after, to_host_fab(FA_v_N0), "FA_v_N0");
    expect_arrays_equal(FA_v_SS_after, to_host_fab(FA_v_SS), "FA_v_SS");
    expect_arrays_equal(FA_v_NN_after, to_host_fab(FA_v_NN), "FA_v_NN");
    expect_arrays_equal(vBT_SS_after,  to_host_fab(vBT_SS),  "vBT_SS");
    expect_arrays_equal(vBT_NN_after,  to_host_fab(vBT_NN),  "vBT_NN");
    expect_arrays_equal(du_cor_after, to_host_fab(du_cor), "du_cor");
    expect_arrays_equal(dv_cor_after, to_host_fab(dv_cor), "dv_cor");
}

// -------------------------------------------------------------------------
// meridional_flux_thickness
// -------------------------------------------------------------------------
//
// No capture/meridional_flux_thickness.{bin,meta} fixture exists yet, so
// this test's field mapping is grounded directly in Fortran source (not a
// .meta file): submodules/MOM6/src/core/MOM_continuity_PPM.F90:5292-5411,
// the meridional_flux_thickness shim's TIMH_capture arm (rec%add(...)
// calls at lines 5352-5365 and 5373). Cross-checked against the bind(C)
// interface and against turbotmp_meridional_flux_thickness_bridge's
// parameter order in mom/cpp/turbotmp_mom_continuity_ppm_bridge.cpp --
// all three agree.
//
// MOM::meridional_flux_thickness computes the effective thickness at
// meridional faces, scaled down to account for the effects of viscosity
// and the fractional open area. visc_rem_v is recorded by the Fortran
// shim only when associated, so a fixture captured with it unassociated
// would be missing the _visc_rem_v field. OBC is never captured -- pass
// nullptr, matching the existing PPM_reconstruction_x/_y tests
// (OBC-inactive configs only).
TEST(MeridionalFluxThickness, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "meridional_flux_thickness");

    const auto   bxC            = captured.box("_bxC");
    const auto   v               = captured.fab_device("_v");
    const auto   h               = captured.fab_device("_h");
    const auto   h_S             = captured.fab_device("_h_S");
    const auto   h_N             = captured.fab_device("_h_N");
    auto         h_v             = captured.fab_device("_h_v_before");
    const auto   h_v_after       = captured.fab_host("_h_v_after");
    const double dt              = captured.real64("_dt");
    const auto   dx_Cv           = captured.fab_device("_dx_Cv");
    const auto   IareaT          = captured.fab_device("_IareaT");
    const auto   IdyT            = captured.fab_device("_IdyT");
    const bool   vol_CFL         = captured.logical("_vol_CFL");
    const bool   marginal        = captured.logical("_marginal");
    const auto   por_face_areaV  = captured.fab_device("_por_face_areaV");
    const auto   visc_rem_v      = captured.fab_device("_visc_rem_v");

    MOM::meridional_flux_thickness(bxC,
                                   v.const_array(),
                                   h.const_array(),
                                   h_S.const_array(),
                                   h_N.const_array(),
                                   h_v.array(),
                                   dt,
                                   dx_Cv.const_array(),
                                   IareaT.const_array(),
                                   IdyT.const_array(),
                                   vol_CFL,
                                   marginal,
                                   /*obc=*/nullptr,
                                   por_face_areaV.const_array(),
                                   visc_rem_v.const_array());
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_v_after, to_host_fab(h_v), "h_v");
}
