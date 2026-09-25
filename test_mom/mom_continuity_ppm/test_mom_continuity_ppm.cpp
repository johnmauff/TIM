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

namespace {

struct OptionalInOutArray {
    amrex::FArrayBox before_fab;
    amrex::FArrayBox after_fab;
    amrex::Array4<amrex::Real> arr{};
    bool present = false;
};

OptionalInOutArray bind_optional_inout(const test_mom::CapturedFile& captured,
                                       const std::string& field) {
    OptionalInOutArray o;
    o.present = captured.is_associated("_" + field + "_before");
    if (o.present) {
        o.before_fab = captured.fab_device("_" + field + "_before");
        o.arr = o.before_fab.array();
        o.after_fab = captured.fab_host("_" + field + "_after");
    }
    return o;
}

} // namespace

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
// OBC is never captured -- pass nullptr, matching the existing
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
// OBC is never captured -- pass nullptr, matching the existing
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
    amrex::FArrayBox visc_rem_u_fab;
    amrex::Array4<const amrex::Real> visc_rem_u{};
    if (captured.is_associated("_visc_rem_u")) {
        visc_rem_u_fab = captured.fab_device("_visc_rem_u");
        visc_rem_u = visc_rem_u_fab.const_array();
    }

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
                              visc_rem_u);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_u_after, to_host_fab(h_u), "h_u");
}

// -------------------------------------------------------------------------
// continuity_meridional_convergence
// -------------------------------------------------------------------------
TEST(ContinuityMeridionalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_meridional_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   vh      = captured.fab_device("_vh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const double h_min   = captured.real64("_h_min");
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
TEST(ContinuityZonalConvergence, MatchesFortranCapture) {
    test_mom::CapturedFile captured(test_mom::data_dir / "continuity_zonal_convergence");

    const auto   bxC     = captured.box("_bxC");
    auto         h       = captured.fab_device("_h_before");
    const auto   h_after = captured.fab_host("_h_after");
    const auto   uh      = captured.fab_device("_uh");
    const double dt      = captured.real64("_dt");
    const auto   IareaT  = captured.fab_device("_IareaT");
    const double h_min   = captured.real64("_h_min");
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
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
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
    amrex::FArrayBox vhbt_fab;
    amrex::Array4<const amrex::Real> vhbt{};
    if (captured.is_associated("_vhbt")) {
        vhbt_fab = captured.fab_device("_vhbt");
        vhbt = vhbt_fab.const_array();
    }
    auto vh_3d = bind_optional_inout(captured, "vh_3d");

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
                                vhbt,
                                vh_3d.arr,
                                /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(dv_after, to_host_fab(dv), "dv");
    if (vh_3d.present) expect_arrays_equal(vh_3d.after_fab, to_host_fab(vh_3d.before_fab), "vh_3d");
}

// -------------------------------------------------------------------------
// zonal_flux_adjust
// -------------------------------------------------------------------------
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
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
    amrex::FArrayBox uhbt_fab;
    amrex::Array4<const amrex::Real> uhbt{};
    if (captured.is_associated("_uhbt")) {
        uhbt_fab = captured.fab_device("_uhbt");
        uhbt = uhbt_fab.const_array();
    }
    auto uh_3d = bind_optional_inout(captured, "uh_3d");

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
                           uhbt,
                           uh_3d.arr,
                           /*obc=*/nullptr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(du_after, to_host_fab(du), "du");
    if (uh_3d.present) expect_arrays_equal(uh_3d.after_fab, to_host_fab(uh_3d.before_fab), "uh_3d");
}

// -------------------------------------------------------------------------
// meridional_mass_flux
// -------------------------------------------------------------------------
// OBC is never captured -- pass nullptr, matching the existing
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
    auto v_cor  = bind_optional_inout(captured, "v_cor");
    auto FA_v_S0 = bind_optional_inout(captured, "FA_v_S0");
    auto FA_v_N0 = bind_optional_inout(captured, "FA_v_N0");
    auto FA_v_SS = bind_optional_inout(captured, "FA_v_SS");
    auto FA_v_NN = bind_optional_inout(captured, "FA_v_NN");
    auto vBT_SS  = bind_optional_inout(captured, "vBT_SS");
    auto vBT_NN  = bind_optional_inout(captured, "vBT_NN");
    auto h_v     = bind_optional_inout(captured, "h_v");
    auto dv_cor  = bind_optional_inout(captured, "dv_cor");

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
                              v_cor.arr,
                              FA_v_S0.arr,
                              FA_v_N0.arr,
                              FA_v_SS.arr,
                              FA_v_NN.arr,
                              vBT_SS.arr,
                              vBT_NN.arr,
                              h_v.arr,
                              dv_cor.arr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(vh_after, to_host_fab(vh), "vh");
    if (v_cor.present)  expect_arrays_equal(v_cor.after_fab,  to_host_fab(v_cor.before_fab),  "v_cor");
    if (FA_v_S0.present) expect_arrays_equal(FA_v_S0.after_fab, to_host_fab(FA_v_S0.before_fab), "FA_v_S0");
    if (FA_v_N0.present) expect_arrays_equal(FA_v_N0.after_fab, to_host_fab(FA_v_N0.before_fab), "FA_v_N0");
    if (FA_v_SS.present) expect_arrays_equal(FA_v_SS.after_fab, to_host_fab(FA_v_SS.before_fab), "FA_v_SS");
    if (FA_v_NN.present) expect_arrays_equal(FA_v_NN.after_fab, to_host_fab(FA_v_NN.before_fab), "FA_v_NN");
    if (vBT_SS.present)  expect_arrays_equal(vBT_SS.after_fab,  to_host_fab(vBT_SS.before_fab),  "vBT_SS");
    if (vBT_NN.present)  expect_arrays_equal(vBT_NN.after_fab,  to_host_fab(vBT_NN.before_fab),  "vBT_NN");
    if (h_v.present)     expect_arrays_equal(h_v.after_fab,     to_host_fab(h_v.before_fab),     "h_v");
    if (dv_cor.present) expect_arrays_equal(dv_cor.after_fab, to_host_fab(dv_cor.before_fab), "dv_cor");
}

// -------------------------------------------------------------------------
// zonal_mass_flux
// -------------------------------------------------------------------------
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
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
    auto u_cor  = bind_optional_inout(captured, "u_cor");
    auto FA_u_W0 = bind_optional_inout(captured, "FA_u_W0");
    auto FA_u_E0 = bind_optional_inout(captured, "FA_u_E0");
    auto FA_u_WW = bind_optional_inout(captured, "FA_u_WW");
    auto FA_u_EE = bind_optional_inout(captured, "FA_u_EE");
    auto uBT_WW  = bind_optional_inout(captured, "uBT_WW");
    auto uBT_EE  = bind_optional_inout(captured, "uBT_EE");
    auto h_u     = bind_optional_inout(captured, "h_u");
    auto du_cor  = bind_optional_inout(captured, "du_cor");

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
                         u_cor.arr,
                         FA_u_W0.arr,
                         FA_u_E0.arr,
                         FA_u_WW.arr,
                         FA_u_EE.arr,
                         uBT_WW.arr,
                         uBT_EE.arr,
                         h_u.arr,
                         du_cor.arr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(uh_after, to_host_fab(uh), "uh");
    if (u_cor.present)  expect_arrays_equal(u_cor.after_fab,  to_host_fab(u_cor.before_fab),  "u_cor");
    if (FA_u_W0.present) expect_arrays_equal(FA_u_W0.after_fab, to_host_fab(FA_u_W0.before_fab), "FA_u_W0");
    if (FA_u_E0.present) expect_arrays_equal(FA_u_E0.after_fab, to_host_fab(FA_u_E0.before_fab), "FA_u_E0");
    if (FA_u_WW.present) expect_arrays_equal(FA_u_WW.after_fab, to_host_fab(FA_u_WW.before_fab), "FA_u_WW");
    if (FA_u_EE.present) expect_arrays_equal(FA_u_EE.after_fab, to_host_fab(FA_u_EE.before_fab), "FA_u_EE");
    if (uBT_WW.present)  expect_arrays_equal(uBT_WW.after_fab,  to_host_fab(uBT_WW.before_fab),  "uBT_WW");
    if (uBT_EE.present)  expect_arrays_equal(uBT_EE.after_fab,  to_host_fab(uBT_EE.before_fab),  "uBT_EE");
    if (h_u.present)     expect_arrays_equal(h_u.after_fab,     to_host_fab(h_u.before_fab),     "h_u");
    if (du_cor.present) expect_arrays_equal(du_cor.after_fab, to_host_fab(du_cor.before_fab), "du_cor");
}

// -------------------------------------------------------------------------
// continuity_PPM
// -------------------------------------------------------------------------
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
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
    auto u_cor  = bind_optional_inout(captured, "u_cor");
    auto v_cor  = bind_optional_inout(captured, "v_cor");
    auto FA_u_W0 = bind_optional_inout(captured, "FA_u_W0");
    auto FA_u_E0 = bind_optional_inout(captured, "FA_u_E0");
    auto FA_u_WW = bind_optional_inout(captured, "FA_u_WW");
    auto FA_u_EE = bind_optional_inout(captured, "FA_u_EE");
    auto uBT_WW  = bind_optional_inout(captured, "uBT_WW");
    auto uBT_EE  = bind_optional_inout(captured, "uBT_EE");
    auto FA_v_S0 = bind_optional_inout(captured, "FA_v_S0");
    auto FA_v_N0 = bind_optional_inout(captured, "FA_v_N0");
    auto FA_v_SS = bind_optional_inout(captured, "FA_v_SS");
    auto FA_v_NN = bind_optional_inout(captured, "FA_v_NN");
    auto vBT_SS  = bind_optional_inout(captured, "vBT_SS");
    auto vBT_NN  = bind_optional_inout(captured, "vBT_NN");
    auto h_u     = bind_optional_inout(captured, "h_u");
    auto h_v     = bind_optional_inout(captured, "h_v");
    auto du_cor  = bind_optional_inout(captured, "du_cor");
    auto dv_cor  = bind_optional_inout(captured, "dv_cor");

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
                        u_cor.arr,
                        v_cor.arr,
                        FA_u_W0.arr,
                        FA_u_E0.arr,
                        FA_u_WW.arr,
                        FA_u_EE.arr,
                        uBT_WW.arr,
                        uBT_EE.arr,
                        FA_v_S0.arr,
                        FA_v_N0.arr,
                        FA_v_SS.arr,
                        FA_v_NN.arr,
                        vBT_SS.arr,
                        vBT_NN.arr,
                        h_u.arr,
                        h_v.arr,
                        du_cor.arr,
                        dv_cor.arr);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_after,  to_host_fab(h),  "h");
    expect_arrays_equal(uh_after, to_host_fab(uh), "uh");
    expect_arrays_equal(vh_after, to_host_fab(vh), "vh");
    if (u_cor.present)  expect_arrays_equal(u_cor.after_fab,  to_host_fab(u_cor.before_fab),  "u_cor");
    if (v_cor.present)  expect_arrays_equal(v_cor.after_fab,  to_host_fab(v_cor.before_fab),  "v_cor");
    if (FA_u_W0.present) expect_arrays_equal(FA_u_W0.after_fab, to_host_fab(FA_u_W0.before_fab), "FA_u_W0");
    if (FA_u_E0.present) expect_arrays_equal(FA_u_E0.after_fab, to_host_fab(FA_u_E0.before_fab), "FA_u_E0");
    if (FA_u_WW.present) expect_arrays_equal(FA_u_WW.after_fab, to_host_fab(FA_u_WW.before_fab), "FA_u_WW");
    if (FA_u_EE.present) expect_arrays_equal(FA_u_EE.after_fab, to_host_fab(FA_u_EE.before_fab), "FA_u_EE");
    if (uBT_WW.present)  expect_arrays_equal(uBT_WW.after_fab,  to_host_fab(uBT_WW.before_fab),  "uBT_WW");
    if (uBT_EE.present)  expect_arrays_equal(uBT_EE.after_fab,  to_host_fab(uBT_EE.before_fab),  "uBT_EE");
    if (FA_v_S0.present) expect_arrays_equal(FA_v_S0.after_fab, to_host_fab(FA_v_S0.before_fab), "FA_v_S0");
    if (FA_v_N0.present) expect_arrays_equal(FA_v_N0.after_fab, to_host_fab(FA_v_N0.before_fab), "FA_v_N0");
    if (FA_v_SS.present) expect_arrays_equal(FA_v_SS.after_fab, to_host_fab(FA_v_SS.before_fab), "FA_v_SS");
    if (FA_v_NN.present) expect_arrays_equal(FA_v_NN.after_fab, to_host_fab(FA_v_NN.before_fab), "FA_v_NN");
    if (vBT_SS.present)  expect_arrays_equal(vBT_SS.after_fab,  to_host_fab(vBT_SS.before_fab),  "vBT_SS");
    if (vBT_NN.present)  expect_arrays_equal(vBT_NN.after_fab,  to_host_fab(vBT_NN.before_fab),  "vBT_NN");
    if (h_u.present)     expect_arrays_equal(h_u.after_fab,     to_host_fab(h_u.before_fab),     "h_u");
    if (h_v.present)     expect_arrays_equal(h_v.after_fab,     to_host_fab(h_v.before_fab),     "h_v");
    if (du_cor.present) expect_arrays_equal(du_cor.after_fab, to_host_fab(du_cor.before_fab), "du_cor");
    if (dv_cor.present) expect_arrays_equal(dv_cor.after_fab, to_host_fab(dv_cor.before_fab), "dv_cor");
}

// -------------------------------------------------------------------------
// meridional_flux_thickness
// -------------------------------------------------------------------------
// OBC is never captured -- pass nullptr, matching the existing
// PPM_reconstruction_x/_y tests (OBC-inactive configs only).
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
    amrex::FArrayBox visc_rem_v_fab;
    amrex::Array4<const amrex::Real> visc_rem_v{};
    if (captured.is_associated("_visc_rem_v")) {
        visc_rem_v_fab = captured.fab_device("_visc_rem_v");
        visc_rem_v = visc_rem_v_fab.const_array();
    }

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
                                   visc_rem_v);
    amrex::Gpu::synchronize();

    expect_arrays_equal(h_v_after, to_host_fab(h_v), "h_v");
}
