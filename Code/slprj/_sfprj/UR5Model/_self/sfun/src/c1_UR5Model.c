/* Include files */

#include <stddef.h>
#include "blas.h"
#include "UR5Model_sfun.h"
#include "c1_UR5Model.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "UR5Model_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c1_debug_family_names[15] = { "q1", "q2", "q3", "q4", "q5",
  "q6", "x", "J", "k", "nargin", "nargout", "x_des", "xd_des", "q", "dq" };

/* Function Declarations */
static void initialize_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void initialize_params_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance);
static void enable_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void disable_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance);
static void set_sim_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_st);
static void finalize_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void sf_gateway_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_chartstep_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void initSimStructsc1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_dq, const char_T *c1_identifier, real_T c1_y[6]);
static void c1_b_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_c_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[36]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_b_info_helper(const mxArray **c1_info);
static void c1_c_info_helper(const mxArray **c1_info);
static void c1_d_info_helper(const mxArray **c1_info);
static real_T c1_eml_div(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x,
  real_T c1_y);
static void c1_eml_switch_helper(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_pinv(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_A[36],
                    real_T c1_X[36]);
static void c1_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_eml_xgesvd(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_A[36], real_T c1_U[36], real_T c1_S[6], real_T c1_V[36]);
static real_T c1_eml_xnrm2(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[36], int32_T c1_ix0);
static void c1_threshold(SFc1_UR5ModelInstanceStruct *chartInstance);
static real_T c1_abs(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x);
static void c1_realmin(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_check_forloop_overflow_error(SFc1_UR5ModelInstanceStruct
  *chartInstance, boolean_T c1_overflow);
static void c1_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_b_x[36]);
static void c1_b_threshold(SFc1_UR5ModelInstanceStruct *chartInstance);
static real_T c1_eml_xdotc(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0);
static void c1_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0, real_T
  c1_b_y[36]);
static void c1_c_threshold(SFc1_UR5ModelInstanceStruct *chartInstance);
static real_T c1_b_eml_xnrm2(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[6], int32_T c1_ix0);
static void c1_b_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_b_x[6]);
static void c1_b_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[6], int32_T
  c1_iy0, real_T c1_b_y[6]);
static void c1_c_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_y[36], int32_T
  c1_iy0, real_T c1_b_y[36]);
static void c1_c_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_b_x[36]);
static void c1_eps(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_b_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_b_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance);
static real_T c1_sqrt(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x);
static void c1_c_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_eml_xrotg(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_a,
  real_T c1_b, real_T *c1_b_a, real_T *c1_b_b, real_T *c1_c, real_T *c1_s);
static void c1_eml_xrot(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x
  [36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_c, real_T c1_s, real_T c1_b_x
  [36]);
static void c1_eml_xswap(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_b_x[36]);
static void c1_d_threshold(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_eml_xgemm(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_k, real_T c1_A[36], real_T c1_B[36], real_T c1_C[36], real_T c1_b_C[36]);
static void c1_e_threshold(SFc1_UR5ModelInstanceStruct *chartInstance);
static void c1_c_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_e_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_f_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_UR5Model, const char_T *c1_identifier);
static uint8_T c1_g_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0);
static void c1_d_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0);
static void c1_e_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0);
static void c1_e_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[6], int32_T
  c1_iy0);
static void c1_f_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_y[36], int32_T
  c1_iy0);
static void c1_f_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_a, real_T c1_x[36], int32_T c1_ix0);
static void c1_b_sqrt(SFc1_UR5ModelInstanceStruct *chartInstance, real_T *c1_x);
static void c1_b_eml_xrotg(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  *c1_a, real_T *c1_b, real_T *c1_c, real_T *c1_s);
static void c1_b_eml_xrot(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_c, real_T c1_s);
static void c1_b_eml_xswap(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0);
static void c1_b_eml_xgemm(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_k, real_T c1_A[36], real_T c1_B[36], real_T c1_C[36]);
static void init_dsm_address_info(SFc1_UR5ModelInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_UR5Model = 0U;
}

static void initialize_params_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[6];
  const mxArray *c1_b_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T (*c1_dq)[6];
  c1_dq = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(2, 1), false);
  for (c1_i0 = 0; c1_i0 < 6; c1_i0++) {
    c1_u[c1_i0] = (*c1_dq)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_UR5Model;
  c1_b_u = c1_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[6];
  int32_T c1_i1;
  real_T (*c1_dq)[6];
  c1_dq = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)), "dq",
                      c1_dv0);
  for (c1_i1 = 0; c1_i1 < 6; c1_i1++) {
    (*c1_dq)[c1_i1] = c1_dv0[c1_i1];
  }

  chartInstance->c1_is_active_c1_UR5Model = c1_f_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 1)), "is_active_c1_UR5Model");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_UR5Model(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  int32_T c1_i2;
  int32_T c1_i3;
  int32_T c1_i4;
  int32_T c1_i5;
  real_T (*c1_q)[6];
  real_T (*c1_xd_des)[6];
  real_T (*c1_dq)[6];
  real_T (*c1_x_des)[6];
  c1_q = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
  c1_xd_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c1_dq = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_x_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i2 = 0; c1_i2 < 6; c1_i2++) {
    _SFD_DATA_RANGE_CHECK((*c1_x_des)[c1_i2], 0U);
  }

  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_UR5Model(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_UR5ModelMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c1_i3 = 0; c1_i3 < 6; c1_i3++) {
    _SFD_DATA_RANGE_CHECK((*c1_dq)[c1_i3], 1U);
  }

  for (c1_i4 = 0; c1_i4 < 6; c1_i4++) {
    _SFD_DATA_RANGE_CHECK((*c1_xd_des)[c1_i4], 2U);
  }

  for (c1_i5 = 0; c1_i5 < 6; c1_i5++) {
    _SFD_DATA_RANGE_CHECK((*c1_q)[c1_i5], 3U);
  }
}

static void c1_chartstep_c1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  int32_T c1_i6;
  real_T c1_x_des[6];
  int32_T c1_i7;
  real_T c1_xd_des[6];
  int32_T c1_i8;
  real_T c1_q[6];
  uint32_T c1_debug_family_var_map[15];
  real_T c1_q1;
  real_T c1_q2;
  real_T c1_q3;
  real_T c1_q4;
  real_T c1_q5;
  real_T c1_q6;
  real_T c1_x[6];
  real_T c1_J[36];
  real_T c1_k[36];
  real_T c1_nargin = 3.0;
  real_T c1_nargout = 1.0;
  real_T c1_dq[6];
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_A;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_y;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_b_A;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_u_x;
  real_T c1_b_y;
  real_T c1_v_x;
  real_T c1_w_x;
  real_T c1_x_x;
  real_T c1_y_x;
  real_T c1_ab_x;
  real_T c1_bb_x;
  real_T c1_c_A;
  real_T c1_cb_x;
  real_T c1_db_x;
  real_T c1_eb_x;
  real_T c1_c_y;
  real_T c1_fb_x;
  real_T c1_gb_x;
  real_T c1_d_A;
  real_T c1_hb_x;
  real_T c1_ib_x;
  real_T c1_jb_x;
  real_T c1_d_y;
  real_T c1_kb_x;
  real_T c1_lb_x;
  real_T c1_mb_x;
  real_T c1_nb_x;
  real_T c1_ob_x;
  real_T c1_pb_x;
  real_T c1_e_A;
  real_T c1_qb_x;
  real_T c1_rb_x;
  real_T c1_sb_x;
  real_T c1_e_y;
  real_T c1_tb_x;
  real_T c1_ub_x;
  real_T c1_vb_x;
  real_T c1_wb_x;
  real_T c1_xb_x;
  real_T c1_yb_x;
  real_T c1_ac_x;
  real_T c1_bc_x;
  real_T c1_cc_x;
  real_T c1_dc_x;
  real_T c1_ec_x;
  real_T c1_fc_x;
  real_T c1_gc_x;
  real_T c1_hc_x;
  real_T c1_f_A;
  real_T c1_ic_x;
  real_T c1_jc_x;
  real_T c1_kc_x;
  real_T c1_f_y;
  real_T c1_lc_x;
  real_T c1_mc_x;
  real_T c1_nc_x;
  real_T c1_oc_x;
  real_T c1_g_A;
  real_T c1_pc_x;
  real_T c1_qc_x;
  real_T c1_rc_x;
  real_T c1_g_y;
  real_T c1_sc_x;
  real_T c1_tc_x;
  real_T c1_uc_x;
  real_T c1_vc_x;
  real_T c1_wc_x;
  real_T c1_xc_x;
  real_T c1_yc_x;
  real_T c1_ad_x;
  real_T c1_bd_x;
  real_T c1_cd_x;
  real_T c1_dd_x;
  real_T c1_ed_x;
  real_T c1_fd_x;
  real_T c1_gd_x;
  real_T c1_hd_x;
  real_T c1_id_x;
  real_T c1_h_A;
  real_T c1_jd_x;
  real_T c1_kd_x;
  real_T c1_ld_x;
  real_T c1_h_y;
  real_T c1_md_x;
  real_T c1_nd_x;
  real_T c1_od_x;
  real_T c1_pd_x;
  real_T c1_qd_x;
  real_T c1_rd_x;
  real_T c1_i_A;
  real_T c1_sd_x;
  real_T c1_td_x;
  real_T c1_ud_x;
  real_T c1_i_y;
  real_T c1_vd_x;
  real_T c1_wd_x;
  real_T c1_xd_x;
  real_T c1_yd_x;
  real_T c1_ae_x;
  real_T c1_be_x;
  real_T c1_j_A;
  real_T c1_ce_x;
  real_T c1_de_x;
  real_T c1_ee_x;
  real_T c1_j_y;
  real_T c1_fe_x;
  real_T c1_ge_x;
  real_T c1_he_x;
  real_T c1_ie_x;
  real_T c1_je_x;
  real_T c1_ke_x;
  real_T c1_le_x;
  real_T c1_me_x;
  real_T c1_ne_x;
  real_T c1_oe_x;
  real_T c1_pe_x;
  real_T c1_qe_x;
  real_T c1_re_x;
  real_T c1_se_x;
  real_T c1_te_x;
  real_T c1_ue_x;
  real_T c1_k_A;
  real_T c1_ve_x;
  real_T c1_we_x;
  real_T c1_xe_x;
  real_T c1_k_y;
  real_T c1_ye_x;
  real_T c1_af_x;
  real_T c1_bf_x;
  real_T c1_cf_x;
  real_T c1_l_A;
  real_T c1_df_x;
  real_T c1_ef_x;
  real_T c1_ff_x;
  real_T c1_l_y;
  real_T c1_gf_x;
  real_T c1_hf_x;
  real_T c1_if_x;
  real_T c1_jf_x;
  real_T c1_m_A;
  real_T c1_kf_x;
  real_T c1_lf_x;
  real_T c1_mf_x;
  real_T c1_m_y;
  real_T c1_nf_x;
  real_T c1_of_x;
  real_T c1_pf_x;
  real_T c1_qf_x;
  real_T c1_n_A;
  real_T c1_rf_x;
  real_T c1_sf_x;
  real_T c1_tf_x;
  real_T c1_n_y;
  real_T c1_uf_x;
  real_T c1_vf_x;
  real_T c1_wf_x;
  real_T c1_xf_x;
  real_T c1_yf_x;
  real_T c1_ag_x;
  real_T c1_bg_x;
  real_T c1_cg_x;
  real_T c1_o_A;
  real_T c1_dg_x;
  real_T c1_eg_x;
  real_T c1_fg_x;
  real_T c1_o_y;
  real_T c1_gg_x;
  real_T c1_hg_x;
  real_T c1_p_A;
  real_T c1_ig_x;
  real_T c1_jg_x;
  real_T c1_kg_x;
  real_T c1_p_y;
  real_T c1_lg_x;
  real_T c1_mg_x;
  real_T c1_ng_x;
  real_T c1_og_x;
  real_T c1_pg_x;
  real_T c1_qg_x;
  real_T c1_rg_x;
  real_T c1_sg_x;
  real_T c1_tg_x;
  real_T c1_ug_x;
  real_T c1_vg_x;
  real_T c1_wg_x;
  real_T c1_q_A;
  real_T c1_xg_x;
  real_T c1_yg_x;
  real_T c1_ah_x;
  real_T c1_q_y;
  real_T c1_bh_x;
  real_T c1_ch_x;
  real_T c1_dh_x;
  real_T c1_eh_x;
  real_T c1_fh_x;
  real_T c1_gh_x;
  real_T c1_r_A;
  real_T c1_hh_x;
  real_T c1_ih_x;
  real_T c1_jh_x;
  real_T c1_r_y;
  real_T c1_kh_x;
  real_T c1_lh_x;
  real_T c1_mh_x;
  real_T c1_nh_x;
  real_T c1_oh_x;
  real_T c1_ph_x;
  real_T c1_s_A;
  real_T c1_qh_x;
  real_T c1_rh_x;
  real_T c1_sh_x;
  real_T c1_s_y;
  real_T c1_th_x;
  real_T c1_uh_x;
  real_T c1_vh_x;
  real_T c1_wh_x;
  real_T c1_xh_x;
  real_T c1_yh_x;
  real_T c1_ai_x;
  real_T c1_bi_x;
  real_T c1_ci_x;
  real_T c1_di_x;
  real_T c1_ei_x;
  real_T c1_fi_x;
  real_T c1_gi_x;
  real_T c1_hi_x;
  real_T c1_ii_x;
  real_T c1_ji_x;
  real_T c1_ki_x;
  real_T c1_li_x;
  real_T c1_t_A;
  real_T c1_mi_x;
  real_T c1_ni_x;
  real_T c1_oi_x;
  real_T c1_t_y;
  real_T c1_pi_x;
  real_T c1_qi_x;
  real_T c1_ri_x;
  real_T c1_si_x;
  real_T c1_u_A;
  real_T c1_ti_x;
  real_T c1_ui_x;
  real_T c1_vi_x;
  real_T c1_u_y;
  real_T c1_wi_x;
  real_T c1_xi_x;
  real_T c1_yi_x;
  real_T c1_aj_x;
  real_T c1_v_A;
  real_T c1_bj_x;
  real_T c1_cj_x;
  real_T c1_dj_x;
  real_T c1_v_y;
  real_T c1_ej_x;
  real_T c1_fj_x;
  real_T c1_gj_x;
  real_T c1_hj_x;
  real_T c1_w_A;
  real_T c1_ij_x;
  real_T c1_jj_x;
  real_T c1_kj_x;
  real_T c1_w_y;
  real_T c1_lj_x;
  real_T c1_mj_x;
  real_T c1_nj_x;
  real_T c1_oj_x;
  real_T c1_x_A;
  real_T c1_pj_x;
  real_T c1_qj_x;
  real_T c1_rj_x;
  real_T c1_x_y;
  real_T c1_sj_x;
  real_T c1_tj_x;
  real_T c1_uj_x;
  real_T c1_vj_x;
  real_T c1_wj_x;
  real_T c1_xj_x;
  real_T c1_y_A;
  real_T c1_yj_x;
  real_T c1_ak_x;
  real_T c1_bk_x;
  real_T c1_y_y;
  real_T c1_ck_x;
  real_T c1_dk_x;
  real_T c1_ek_x;
  real_T c1_fk_x;
  real_T c1_ab_A;
  real_T c1_gk_x;
  real_T c1_hk_x;
  real_T c1_ik_x;
  real_T c1_ab_y;
  real_T c1_jk_x;
  real_T c1_kk_x;
  real_T c1_lk_x;
  real_T c1_mk_x;
  real_T c1_bb_A;
  real_T c1_nk_x;
  real_T c1_ok_x;
  real_T c1_pk_x;
  real_T c1_bb_y;
  real_T c1_qk_x;
  real_T c1_rk_x;
  real_T c1_sk_x;
  real_T c1_tk_x;
  real_T c1_cb_A;
  real_T c1_uk_x;
  real_T c1_vk_x;
  real_T c1_wk_x;
  real_T c1_cb_y;
  real_T c1_xk_x;
  real_T c1_yk_x;
  real_T c1_al_x;
  real_T c1_bl_x;
  real_T c1_cl_x;
  real_T c1_dl_x;
  real_T c1_el_x;
  real_T c1_fl_x;
  real_T c1_gl_x;
  real_T c1_hl_x;
  real_T c1_db_A;
  real_T c1_il_x;
  real_T c1_jl_x;
  real_T c1_kl_x;
  real_T c1_db_y;
  real_T c1_ll_x;
  real_T c1_ml_x;
  real_T c1_nl_x;
  real_T c1_ol_x;
  real_T c1_pl_x;
  real_T c1_ql_x;
  real_T c1_rl_x;
  real_T c1_sl_x;
  real_T c1_tl_x;
  real_T c1_ul_x;
  real_T c1_eb_A;
  real_T c1_vl_x;
  real_T c1_wl_x;
  real_T c1_xl_x;
  real_T c1_eb_y;
  real_T c1_yl_x;
  real_T c1_am_x;
  real_T c1_bm_x;
  real_T c1_cm_x;
  real_T c1_dm_x;
  real_T c1_em_x;
  real_T c1_fm_x;
  real_T c1_gm_x;
  real_T c1_hm_x;
  real_T c1_im_x;
  real_T c1_fb_A;
  real_T c1_jm_x;
  real_T c1_km_x;
  real_T c1_lm_x;
  real_T c1_fb_y;
  real_T c1_mm_x;
  real_T c1_nm_x;
  real_T c1_om_x;
  real_T c1_pm_x;
  real_T c1_qm_x;
  real_T c1_rm_x;
  real_T c1_sm_x;
  real_T c1_tm_x;
  real_T c1_um_x;
  real_T c1_vm_x;
  real_T c1_gb_A;
  real_T c1_wm_x;
  real_T c1_xm_x;
  real_T c1_ym_x;
  real_T c1_gb_y;
  real_T c1_an_x;
  real_T c1_bn_x;
  real_T c1_cn_x;
  real_T c1_dn_x;
  real_T c1_en_x;
  real_T c1_fn_x;
  real_T c1_hb_A;
  real_T c1_gn_x;
  real_T c1_hn_x;
  real_T c1_in_x;
  real_T c1_hb_y;
  real_T c1_jn_x;
  real_T c1_kn_x;
  real_T c1_ln_x;
  real_T c1_mn_x;
  real_T c1_nn_x;
  real_T c1_on_x;
  real_T c1_pn_x;
  real_T c1_qn_x;
  real_T c1_ib_A;
  real_T c1_rn_x;
  real_T c1_sn_x;
  real_T c1_tn_x;
  real_T c1_ib_y;
  real_T c1_un_x;
  real_T c1_vn_x;
  real_T c1_wn_x;
  real_T c1_xn_x;
  real_T c1_yn_x;
  real_T c1_ao_x;
  real_T c1_jb_A;
  real_T c1_bo_x;
  real_T c1_co_x;
  real_T c1_do_x;
  real_T c1_jb_y;
  real_T c1_eo_x;
  real_T c1_fo_x;
  real_T c1_kb_A;
  real_T c1_go_x;
  real_T c1_ho_x;
  real_T c1_io_x;
  real_T c1_kb_y;
  real_T c1_jo_x;
  real_T c1_ko_x;
  real_T c1_lo_x;
  real_T c1_mo_x;
  real_T c1_no_x;
  real_T c1_oo_x;
  real_T c1_lb_A;
  real_T c1_po_x;
  real_T c1_qo_x;
  real_T c1_ro_x;
  real_T c1_lb_y;
  real_T c1_so_x;
  real_T c1_to_x;
  real_T c1_uo_x;
  real_T c1_vo_x;
  real_T c1_wo_x;
  real_T c1_xo_x;
  real_T c1_yo_x;
  real_T c1_ap_x;
  real_T c1_bp_x;
  real_T c1_cp_x;
  real_T c1_dp_x;
  real_T c1_ep_x;
  real_T c1_fp_x;
  real_T c1_gp_x;
  real_T c1_hp_x;
  real_T c1_ip_x;
  real_T c1_jp_x;
  real_T c1_kp_x;
  real_T c1_lp_x;
  real_T c1_mp_x;
  real_T c1_np_x;
  real_T c1_op_x;
  real_T c1_pp_x;
  real_T c1_qp_x;
  real_T c1_mb_A;
  real_T c1_rp_x;
  real_T c1_sp_x;
  real_T c1_tp_x;
  real_T c1_mb_y;
  real_T c1_up_x;
  real_T c1_vp_x;
  real_T c1_wp_x;
  real_T c1_xp_x;
  real_T c1_nb_A;
  real_T c1_yp_x;
  real_T c1_aq_x;
  real_T c1_bq_x;
  real_T c1_nb_y;
  real_T c1_cq_x;
  real_T c1_dq_x;
  real_T c1_eq_x;
  real_T c1_fq_x;
  real_T c1_ob_A;
  real_T c1_gq_x;
  real_T c1_hq_x;
  real_T c1_iq_x;
  real_T c1_ob_y;
  real_T c1_jq_x;
  real_T c1_kq_x;
  real_T c1_lq_x;
  real_T c1_mq_x;
  real_T c1_pb_A;
  real_T c1_nq_x;
  real_T c1_oq_x;
  real_T c1_pq_x;
  real_T c1_pb_y;
  real_T c1_qq_x;
  real_T c1_rq_x;
  real_T c1_sq_x;
  real_T c1_tq_x;
  real_T c1_qb_A;
  real_T c1_uq_x;
  real_T c1_vq_x;
  real_T c1_wq_x;
  real_T c1_qb_y;
  real_T c1_xq_x;
  real_T c1_yq_x;
  real_T c1_ar_x;
  real_T c1_br_x;
  real_T c1_cr_x;
  real_T c1_dr_x;
  real_T c1_rb_A;
  real_T c1_er_x;
  real_T c1_fr_x;
  real_T c1_gr_x;
  real_T c1_rb_y;
  real_T c1_hr_x;
  real_T c1_ir_x;
  real_T c1_jr_x;
  real_T c1_kr_x;
  real_T c1_sb_A;
  real_T c1_lr_x;
  real_T c1_mr_x;
  real_T c1_nr_x;
  real_T c1_sb_y;
  real_T c1_or_x;
  real_T c1_pr_x;
  real_T c1_qr_x;
  real_T c1_rr_x;
  real_T c1_tb_A;
  real_T c1_sr_x;
  real_T c1_tr_x;
  real_T c1_ur_x;
  real_T c1_tb_y;
  real_T c1_vr_x;
  real_T c1_wr_x;
  real_T c1_xr_x;
  real_T c1_yr_x;
  real_T c1_as_x;
  real_T c1_bs_x;
  real_T c1_cs_x;
  real_T c1_ds_x;
  real_T c1_es_x;
  real_T c1_fs_x;
  real_T c1_ub_A;
  real_T c1_gs_x;
  real_T c1_hs_x;
  real_T c1_is_x;
  real_T c1_ub_y;
  real_T c1_js_x;
  real_T c1_ks_x;
  real_T c1_ls_x;
  real_T c1_ms_x;
  real_T c1_ns_x;
  real_T c1_os_x;
  real_T c1_ps_x;
  real_T c1_qs_x;
  real_T c1_rs_x;
  real_T c1_ss_x;
  real_T c1_vb_A;
  real_T c1_ts_x;
  real_T c1_us_x;
  real_T c1_vs_x;
  real_T c1_vb_y;
  real_T c1_ws_x;
  real_T c1_xs_x;
  real_T c1_ys_x;
  real_T c1_at_x;
  real_T c1_wb_A;
  real_T c1_bt_x;
  real_T c1_ct_x;
  real_T c1_dt_x;
  real_T c1_wb_y;
  real_T c1_et_x;
  real_T c1_ft_x;
  real_T c1_gt_x;
  real_T c1_ht_x;
  real_T c1_it_x;
  real_T c1_jt_x;
  real_T c1_kt_x;
  real_T c1_lt_x;
  real_T c1_mt_x;
  real_T c1_nt_x;
  real_T c1_xb_A;
  real_T c1_ot_x;
  real_T c1_pt_x;
  real_T c1_qt_x;
  real_T c1_xb_y;
  real_T c1_rt_x;
  real_T c1_st_x;
  real_T c1_tt_x;
  real_T c1_ut_x;
  real_T c1_vt_x;
  real_T c1_wt_x;
  real_T c1_xt_x;
  real_T c1_yt_x;
  real_T c1_au_x;
  real_T c1_bu_x;
  real_T c1_yb_A;
  real_T c1_cu_x;
  real_T c1_du_x;
  real_T c1_eu_x;
  real_T c1_yb_y;
  real_T c1_fu_x;
  real_T c1_gu_x;
  real_T c1_ac_A;
  real_T c1_hu_x;
  real_T c1_iu_x;
  real_T c1_ju_x;
  real_T c1_ac_y;
  real_T c1_ku_x;
  real_T c1_lu_x;
  real_T c1_bc_A;
  real_T c1_mu_x;
  real_T c1_nu_x;
  real_T c1_ou_x;
  real_T c1_bc_y;
  real_T c1_pu_x;
  real_T c1_qu_x;
  real_T c1_ru_x;
  real_T c1_su_x;
  real_T c1_tu_x;
  real_T c1_uu_x;
  real_T c1_cc_A;
  real_T c1_vu_x;
  real_T c1_wu_x;
  real_T c1_xu_x;
  real_T c1_cc_y;
  real_T c1_yu_x;
  real_T c1_av_x;
  real_T c1_dc_A;
  real_T c1_bv_x;
  real_T c1_cv_x;
  real_T c1_dv_x;
  real_T c1_dc_y;
  real_T c1_ev_x;
  real_T c1_fv_x;
  real_T c1_ec_A;
  real_T c1_gv_x;
  real_T c1_hv_x;
  real_T c1_iv_x;
  real_T c1_ec_y;
  real_T c1_jv_x;
  real_T c1_kv_x;
  real_T c1_lv_x;
  real_T c1_mv_x;
  real_T c1_fc_A;
  real_T c1_nv_x;
  real_T c1_ov_x;
  real_T c1_pv_x;
  real_T c1_fc_y;
  real_T c1_qv_x;
  real_T c1_rv_x;
  real_T c1_gc_A;
  real_T c1_sv_x;
  real_T c1_tv_x;
  real_T c1_uv_x;
  real_T c1_gc_y;
  real_T c1_vv_x;
  real_T c1_wv_x;
  real_T c1_hc_A;
  real_T c1_xv_x;
  real_T c1_yv_x;
  real_T c1_aw_x;
  real_T c1_hc_y;
  real_T c1_bw_x;
  real_T c1_cw_x;
  real_T c1_ic_A;
  real_T c1_dw_x;
  real_T c1_ew_x;
  real_T c1_fw_x;
  real_T c1_ic_y;
  real_T c1_gw_x;
  real_T c1_hw_x;
  real_T c1_jc_A;
  real_T c1_iw_x;
  real_T c1_jw_x;
  real_T c1_kw_x;
  real_T c1_jc_y;
  real_T c1_lw_x;
  real_T c1_mw_x;
  real_T c1_kc_A;
  real_T c1_nw_x;
  real_T c1_ow_x;
  real_T c1_pw_x;
  real_T c1_kc_y;
  real_T c1_qw_x;
  real_T c1_rw_x;
  real_T c1_sw_x;
  real_T c1_tw_x;
  real_T c1_uw_x;
  real_T c1_vw_x;
  real_T c1_ww_x;
  real_T c1_xw_x;
  real_T c1_yw_x;
  real_T c1_ax_x;
  real_T c1_bx_x;
  real_T c1_cx_x;
  real_T c1_dx_x;
  real_T c1_ex_x;
  real_T c1_fx_x;
  real_T c1_gx_x;
  real_T c1_hx_x;
  real_T c1_ix_x;
  real_T c1_jx_x;
  real_T c1_kx_x;
  real_T c1_lx_x;
  real_T c1_mx_x;
  real_T c1_nx_x;
  real_T c1_ox_x;
  real_T c1_px_x;
  real_T c1_qx_x;
  real_T c1_rx_x;
  real_T c1_sx_x;
  real_T c1_tx_x;
  real_T c1_ux_x;
  real_T c1_vx_x;
  real_T c1_wx_x;
  real_T c1_xx_x;
  real_T c1_yx_x;
  real_T c1_ay_x;
  real_T c1_by_x;
  real_T c1_cy_x;
  real_T c1_dy_x;
  real_T c1_ey_x;
  real_T c1_fy_x;
  real_T c1_gy_x;
  real_T c1_hy_x;
  real_T c1_iy_x;
  real_T c1_jy_x;
  real_T c1_ky_x;
  real_T c1_ly_x;
  int32_T c1_i9;
  static real_T c1_a[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c1_i10;
  real_T c1_b[6];
  int32_T c1_i11;
  real_T c1_lc_y[6];
  int32_T c1_i12;
  int32_T c1_i13;
  int32_T c1_i14;
  real_T c1_b_J[36];
  real_T c1_b_a[36];
  int32_T c1_i15;
  int32_T c1_i16;
  int32_T c1_i17;
  int32_T c1_i18;
  int32_T c1_i19;
  int32_T c1_i20;
  int32_T c1_i21;
  int32_T c1_i22;
  int32_T c1_i23;
  int32_T c1_i24;
  int32_T c1_i25;
  real_T (*c1_b_dq)[6];
  real_T (*c1_b_q)[6];
  real_T (*c1_b_xd_des)[6];
  real_T (*c1_b_x_des)[6];
  c1_b_q = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
  c1_b_xd_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_dq = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_b_x_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i6 = 0; c1_i6 < 6; c1_i6++) {
    c1_x_des[c1_i6] = (*c1_b_x_des)[c1_i6];
  }

  for (c1_i7 = 0; c1_i7 < 6; c1_i7++) {
    c1_xd_des[c1_i7] = (*c1_b_xd_des)[c1_i7];
  }

  for (c1_i8 = 0; c1_i8 < 6; c1_i8++) {
    c1_q[c1_i8] = (*c1_b_q)[c1_i8];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q1, 0U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q2, 1U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q3, 2U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q4, 3U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q5, 4U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_q6, 5U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_x, 6U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_J, 7U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_k, 8U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 9U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 10U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_x_des, 11U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_xd_des, 12U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_q, 13U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dq, 14U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 3);
  c1_q1 = c1_q[0];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 4);
  c1_q2 = c1_q[1];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 5);
  c1_q3 = c1_q[2];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 6);
  c1_q4 = c1_q[3];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 7);
  c1_q5 = c1_q[4];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 8);
  c1_q6 = c1_q[5];
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 11);
  c1_b_x = c1_q2 + c1_q3;
  c1_c_x = c1_b_x;
  c1_c_x = muDoubleScalarCos(c1_c_x);
  c1_d_x = c1_q1;
  c1_e_x = c1_d_x;
  c1_e_x = muDoubleScalarCos(c1_e_x);
  c1_f_x = c1_q4;
  c1_g_x = c1_f_x;
  c1_g_x = muDoubleScalarSin(c1_g_x);
  c1_A = 379.0 * c1_c_x * c1_e_x * c1_g_x;
  c1_h_x = c1_A;
  c1_i_x = c1_h_x;
  c1_j_x = c1_i_x;
  c1_y = c1_j_x / 4.0;
  c1_k_x = c1_q1;
  c1_l_x = c1_k_x;
  c1_l_x = muDoubleScalarCos(c1_l_x);
  c1_m_x = c1_q2;
  c1_n_x = c1_m_x;
  c1_n_x = muDoubleScalarCos(c1_n_x);
  c1_o_x = c1_q5;
  c1_p_x = c1_o_x;
  c1_p_x = muDoubleScalarCos(c1_p_x);
  c1_q_x = c1_q1;
  c1_r_x = c1_q_x;
  c1_r_x = muDoubleScalarSin(c1_r_x);
  c1_b_A = 165.0 * c1_p_x * c1_r_x;
  c1_s_x = c1_b_A;
  c1_t_x = c1_s_x;
  c1_u_x = c1_t_x;
  c1_b_y = c1_u_x / 2.0;
  c1_v_x = (c1_q2 + c1_q3) + c1_q4;
  c1_w_x = c1_v_x;
  c1_w_x = muDoubleScalarCos(c1_w_x);
  c1_x_x = c1_q1;
  c1_y_x = c1_x_x;
  c1_y_x = muDoubleScalarCos(c1_y_x);
  c1_ab_x = c1_q5;
  c1_bb_x = c1_ab_x;
  c1_bb_x = muDoubleScalarSin(c1_bb_x);
  c1_c_A = 165.0 * c1_w_x * c1_y_x * c1_bb_x;
  c1_cb_x = c1_c_A;
  c1_db_x = c1_cb_x;
  c1_eb_x = c1_db_x;
  c1_c_y = c1_eb_x / 2.0;
  c1_fb_x = c1_q1;
  c1_gb_x = c1_fb_x;
  c1_gb_x = muDoubleScalarSin(c1_gb_x);
  c1_d_A = 1093.0 * c1_gb_x;
  c1_hb_x = c1_d_A;
  c1_ib_x = c1_hb_x;
  c1_jb_x = c1_ib_x;
  c1_d_y = c1_jb_x / 10.0;
  c1_kb_x = c1_q2 + c1_q3;
  c1_lb_x = c1_kb_x;
  c1_lb_x = muDoubleScalarSin(c1_lb_x);
  c1_mb_x = c1_q1;
  c1_nb_x = c1_mb_x;
  c1_nb_x = muDoubleScalarCos(c1_nb_x);
  c1_ob_x = c1_q4;
  c1_pb_x = c1_ob_x;
  c1_pb_x = muDoubleScalarCos(c1_pb_x);
  c1_e_A = 379.0 * c1_lb_x * c1_nb_x * c1_pb_x;
  c1_qb_x = c1_e_A;
  c1_rb_x = c1_qb_x;
  c1_sb_x = c1_rb_x;
  c1_e_y = c1_sb_x / 4.0;
  c1_tb_x = c1_q1;
  c1_ub_x = c1_tb_x;
  c1_ub_x = muDoubleScalarCos(c1_ub_x);
  c1_vb_x = c1_q2;
  c1_wb_x = c1_vb_x;
  c1_wb_x = muDoubleScalarCos(c1_wb_x);
  c1_xb_x = c1_q3;
  c1_yb_x = c1_xb_x;
  c1_yb_x = muDoubleScalarCos(c1_yb_x);
  c1_ac_x = c1_q1;
  c1_bc_x = c1_ac_x;
  c1_bc_x = muDoubleScalarCos(c1_bc_x);
  c1_cc_x = c1_q2;
  c1_dc_x = c1_cc_x;
  c1_dc_x = muDoubleScalarSin(c1_dc_x);
  c1_ec_x = c1_q3;
  c1_fc_x = c1_ec_x;
  c1_fc_x = muDoubleScalarSin(c1_fc_x);
  c1_gc_x = c1_q1;
  c1_hc_x = c1_gc_x;
  c1_hc_x = muDoubleScalarCos(c1_hc_x);
  c1_f_A = 1093.0 * c1_hc_x;
  c1_ic_x = c1_f_A;
  c1_jc_x = c1_ic_x;
  c1_kc_x = c1_jc_x;
  c1_f_y = c1_kc_x / 10.0;
  c1_lc_x = c1_q1;
  c1_mc_x = c1_lc_x;
  c1_mc_x = muDoubleScalarCos(c1_mc_x);
  c1_nc_x = c1_q5;
  c1_oc_x = c1_nc_x;
  c1_oc_x = muDoubleScalarCos(c1_oc_x);
  c1_g_A = 165.0 * c1_mc_x * c1_oc_x;
  c1_pc_x = c1_g_A;
  c1_qc_x = c1_pc_x;
  c1_rc_x = c1_qc_x;
  c1_g_y = c1_rc_x / 2.0;
  c1_sc_x = c1_q2;
  c1_tc_x = c1_sc_x;
  c1_tc_x = muDoubleScalarCos(c1_tc_x);
  c1_uc_x = c1_q1;
  c1_vc_x = c1_uc_x;
  c1_vc_x = muDoubleScalarSin(c1_vc_x);
  c1_wc_x = c1_q1;
  c1_xc_x = c1_wc_x;
  c1_xc_x = muDoubleScalarSin(c1_xc_x);
  c1_yc_x = c1_q2;
  c1_ad_x = c1_yc_x;
  c1_ad_x = muDoubleScalarSin(c1_ad_x);
  c1_bd_x = c1_q3;
  c1_cd_x = c1_bd_x;
  c1_cd_x = muDoubleScalarSin(c1_cd_x);
  c1_dd_x = (c1_q2 + c1_q3) + c1_q4;
  c1_ed_x = c1_dd_x;
  c1_ed_x = muDoubleScalarCos(c1_ed_x);
  c1_fd_x = c1_q1;
  c1_gd_x = c1_fd_x;
  c1_gd_x = muDoubleScalarSin(c1_gd_x);
  c1_hd_x = c1_q5;
  c1_id_x = c1_hd_x;
  c1_id_x = muDoubleScalarSin(c1_id_x);
  c1_h_A = 165.0 * c1_ed_x * c1_gd_x * c1_id_x;
  c1_jd_x = c1_h_A;
  c1_kd_x = c1_jd_x;
  c1_ld_x = c1_kd_x;
  c1_h_y = c1_ld_x / 2.0;
  c1_md_x = c1_q2 + c1_q3;
  c1_nd_x = c1_md_x;
  c1_nd_x = muDoubleScalarCos(c1_nd_x);
  c1_od_x = c1_q1;
  c1_pd_x = c1_od_x;
  c1_pd_x = muDoubleScalarSin(c1_pd_x);
  c1_qd_x = c1_q4;
  c1_rd_x = c1_qd_x;
  c1_rd_x = muDoubleScalarSin(c1_rd_x);
  c1_i_A = 379.0 * c1_nd_x * c1_pd_x * c1_rd_x;
  c1_sd_x = c1_i_A;
  c1_td_x = c1_sd_x;
  c1_ud_x = c1_td_x;
  c1_i_y = c1_ud_x / 4.0;
  c1_vd_x = c1_q2 + c1_q3;
  c1_wd_x = c1_vd_x;
  c1_wd_x = muDoubleScalarSin(c1_wd_x);
  c1_xd_x = c1_q4;
  c1_yd_x = c1_xd_x;
  c1_yd_x = muDoubleScalarCos(c1_yd_x);
  c1_ae_x = c1_q1;
  c1_be_x = c1_ae_x;
  c1_be_x = muDoubleScalarSin(c1_be_x);
  c1_j_A = 379.0 * c1_wd_x * c1_yd_x * c1_be_x;
  c1_ce_x = c1_j_A;
  c1_de_x = c1_ce_x;
  c1_ee_x = c1_de_x;
  c1_j_y = c1_ee_x / 4.0;
  c1_fe_x = c1_q2;
  c1_ge_x = c1_fe_x;
  c1_ge_x = muDoubleScalarCos(c1_ge_x);
  c1_he_x = c1_q3;
  c1_ie_x = c1_he_x;
  c1_ie_x = muDoubleScalarCos(c1_ie_x);
  c1_je_x = c1_q1;
  c1_ke_x = c1_je_x;
  c1_ke_x = muDoubleScalarSin(c1_ke_x);
  c1_le_x = c1_q2 + c1_q3;
  c1_me_x = c1_le_x;
  c1_me_x = muDoubleScalarSin(c1_me_x);
  c1_ne_x = c1_q2;
  c1_oe_x = c1_ne_x;
  c1_oe_x = muDoubleScalarSin(c1_oe_x);
  c1_pe_x = c1_q5;
  c1_qe_x = c1_pe_x;
  c1_qe_x = muDoubleScalarSin(c1_qe_x);
  c1_re_x = c1_q2 + c1_q3;
  c1_se_x = c1_re_x;
  c1_se_x = muDoubleScalarCos(c1_se_x);
  c1_te_x = c1_q4;
  c1_ue_x = c1_te_x;
  c1_ue_x = muDoubleScalarSin(c1_ue_x);
  c1_k_A = 165.0 * c1_se_x * c1_ue_x;
  c1_ve_x = c1_k_A;
  c1_we_x = c1_ve_x;
  c1_xe_x = c1_we_x;
  c1_k_y = c1_xe_x / 2.0;
  c1_ye_x = c1_q2 + c1_q3;
  c1_af_x = c1_ye_x;
  c1_af_x = muDoubleScalarSin(c1_af_x);
  c1_bf_x = c1_q4;
  c1_cf_x = c1_bf_x;
  c1_cf_x = muDoubleScalarCos(c1_cf_x);
  c1_l_A = 165.0 * c1_af_x * c1_cf_x;
  c1_df_x = c1_l_A;
  c1_ef_x = c1_df_x;
  c1_ff_x = c1_ef_x;
  c1_l_y = c1_ff_x / 2.0;
  c1_gf_x = c1_q2 + c1_q3;
  c1_hf_x = c1_gf_x;
  c1_hf_x = muDoubleScalarCos(c1_hf_x);
  c1_if_x = c1_q4;
  c1_jf_x = c1_if_x;
  c1_jf_x = muDoubleScalarCos(c1_jf_x);
  c1_m_A = 379.0 * c1_hf_x * c1_jf_x;
  c1_kf_x = c1_m_A;
  c1_lf_x = c1_kf_x;
  c1_mf_x = c1_lf_x;
  c1_m_y = c1_mf_x / 4.0;
  c1_nf_x = c1_q2 + c1_q3;
  c1_of_x = c1_nf_x;
  c1_of_x = muDoubleScalarSin(c1_of_x);
  c1_pf_x = c1_q4;
  c1_qf_x = c1_pf_x;
  c1_qf_x = muDoubleScalarSin(c1_qf_x);
  c1_n_A = 379.0 * c1_of_x * c1_qf_x;
  c1_rf_x = c1_n_A;
  c1_sf_x = c1_rf_x;
  c1_tf_x = c1_sf_x;
  c1_n_y = c1_tf_x / 4.0;
  c1_x[0] = ((((((c1_y - 425.0 * c1_l_x * c1_n_x) - c1_b_y) - c1_c_y) - c1_d_y)
              + c1_e_y) - 392.0 * c1_ub_x * c1_wb_x * c1_yb_x) + 392.0 * c1_bc_x
    * c1_dc_x * c1_fc_x;
  c1_x[1] = ((((((c1_f_y + c1_g_y) - 425.0 * c1_tc_x * c1_vc_x) + 392.0 *
                c1_xc_x * c1_ad_x * c1_cd_x) - c1_h_y) + c1_i_y) + c1_j_y) -
    392.0 * c1_ge_x * c1_ie_x * c1_ke_x;
  c1_x[2] = ((((392.0 * c1_me_x + 425.0 * c1_oe_x) + c1_qe_x * (c1_k_y + c1_l_y))
              + c1_m_y) - c1_n_y) + 89.2;
  c1_x[3] = 0.0;
  c1_x[4] = 0.0;
  c1_x[5] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 17);
  c1_uf_x = c1_q2;
  c1_vf_x = c1_uf_x;
  c1_vf_x = muDoubleScalarCos(c1_vf_x);
  c1_wf_x = c1_q1;
  c1_xf_x = c1_wf_x;
  c1_xf_x = muDoubleScalarSin(c1_xf_x);
  c1_yf_x = c1_q1;
  c1_ag_x = c1_yf_x;
  c1_ag_x = muDoubleScalarCos(c1_ag_x);
  c1_bg_x = c1_q5;
  c1_cg_x = c1_bg_x;
  c1_cg_x = muDoubleScalarCos(c1_cg_x);
  c1_o_A = 165.0 * c1_ag_x * c1_cg_x;
  c1_dg_x = c1_o_A;
  c1_eg_x = c1_dg_x;
  c1_fg_x = c1_eg_x;
  c1_o_y = c1_fg_x / 2.0;
  c1_gg_x = c1_q1;
  c1_hg_x = c1_gg_x;
  c1_hg_x = muDoubleScalarCos(c1_hg_x);
  c1_p_A = 1093.0 * c1_hg_x;
  c1_ig_x = c1_p_A;
  c1_jg_x = c1_ig_x;
  c1_kg_x = c1_jg_x;
  c1_p_y = c1_kg_x / 10.0;
  c1_lg_x = c1_q1;
  c1_mg_x = c1_lg_x;
  c1_mg_x = muDoubleScalarSin(c1_mg_x);
  c1_ng_x = c1_q2;
  c1_og_x = c1_ng_x;
  c1_og_x = muDoubleScalarSin(c1_og_x);
  c1_pg_x = c1_q3;
  c1_qg_x = c1_pg_x;
  c1_qg_x = muDoubleScalarSin(c1_qg_x);
  c1_rg_x = (c1_q2 + c1_q3) + c1_q4;
  c1_sg_x = c1_rg_x;
  c1_sg_x = muDoubleScalarCos(c1_sg_x);
  c1_tg_x = c1_q1;
  c1_ug_x = c1_tg_x;
  c1_ug_x = muDoubleScalarSin(c1_ug_x);
  c1_vg_x = c1_q5;
  c1_wg_x = c1_vg_x;
  c1_wg_x = muDoubleScalarSin(c1_wg_x);
  c1_q_A = 165.0 * c1_sg_x * c1_ug_x * c1_wg_x;
  c1_xg_x = c1_q_A;
  c1_yg_x = c1_xg_x;
  c1_ah_x = c1_yg_x;
  c1_q_y = c1_ah_x / 2.0;
  c1_bh_x = c1_q2 + c1_q3;
  c1_ch_x = c1_bh_x;
  c1_ch_x = muDoubleScalarCos(c1_ch_x);
  c1_dh_x = c1_q1;
  c1_eh_x = c1_dh_x;
  c1_eh_x = muDoubleScalarSin(c1_eh_x);
  c1_fh_x = c1_q4;
  c1_gh_x = c1_fh_x;
  c1_gh_x = muDoubleScalarSin(c1_gh_x);
  c1_r_A = 379.0 * c1_ch_x * c1_eh_x * c1_gh_x;
  c1_hh_x = c1_r_A;
  c1_ih_x = c1_hh_x;
  c1_jh_x = c1_ih_x;
  c1_r_y = c1_jh_x / 4.0;
  c1_kh_x = c1_q2 + c1_q3;
  c1_lh_x = c1_kh_x;
  c1_lh_x = muDoubleScalarSin(c1_lh_x);
  c1_mh_x = c1_q4;
  c1_nh_x = c1_mh_x;
  c1_nh_x = muDoubleScalarCos(c1_nh_x);
  c1_oh_x = c1_q1;
  c1_ph_x = c1_oh_x;
  c1_ph_x = muDoubleScalarSin(c1_ph_x);
  c1_s_A = 379.0 * c1_lh_x * c1_nh_x * c1_ph_x;
  c1_qh_x = c1_s_A;
  c1_rh_x = c1_qh_x;
  c1_sh_x = c1_rh_x;
  c1_s_y = c1_sh_x / 4.0;
  c1_th_x = c1_q2;
  c1_uh_x = c1_th_x;
  c1_uh_x = muDoubleScalarCos(c1_uh_x);
  c1_vh_x = c1_q3;
  c1_wh_x = c1_vh_x;
  c1_wh_x = muDoubleScalarCos(c1_wh_x);
  c1_xh_x = c1_q1;
  c1_yh_x = c1_xh_x;
  c1_yh_x = muDoubleScalarSin(c1_yh_x);
  c1_ai_x = c1_q1;
  c1_bi_x = c1_ai_x;
  c1_bi_x = muDoubleScalarCos(c1_bi_x);
  c1_ci_x = c1_q2 + c1_q3;
  c1_di_x = c1_ci_x;
  c1_di_x = muDoubleScalarSin(c1_di_x);
  c1_ei_x = c1_q2;
  c1_fi_x = c1_ei_x;
  c1_fi_x = muDoubleScalarSin(c1_fi_x);
  c1_gi_x = c1_q5;
  c1_hi_x = c1_gi_x;
  c1_hi_x = muDoubleScalarSin(c1_hi_x);
  c1_ii_x = c1_q2 + c1_q3;
  c1_ji_x = c1_ii_x;
  c1_ji_x = muDoubleScalarCos(c1_ji_x);
  c1_ki_x = c1_q4;
  c1_li_x = c1_ki_x;
  c1_li_x = muDoubleScalarSin(c1_li_x);
  c1_t_A = 165.0 * c1_ji_x * c1_li_x;
  c1_mi_x = c1_t_A;
  c1_ni_x = c1_mi_x;
  c1_oi_x = c1_ni_x;
  c1_t_y = c1_oi_x / 2.0;
  c1_pi_x = c1_q2 + c1_q3;
  c1_qi_x = c1_pi_x;
  c1_qi_x = muDoubleScalarSin(c1_qi_x);
  c1_ri_x = c1_q4;
  c1_si_x = c1_ri_x;
  c1_si_x = muDoubleScalarCos(c1_si_x);
  c1_u_A = 165.0 * c1_qi_x * c1_si_x;
  c1_ti_x = c1_u_A;
  c1_ui_x = c1_ti_x;
  c1_vi_x = c1_ui_x;
  c1_u_y = c1_vi_x / 2.0;
  c1_wi_x = c1_q2 + c1_q3;
  c1_xi_x = c1_wi_x;
  c1_xi_x = muDoubleScalarCos(c1_xi_x);
  c1_yi_x = c1_q4;
  c1_aj_x = c1_yi_x;
  c1_aj_x = muDoubleScalarCos(c1_aj_x);
  c1_v_A = 379.0 * c1_xi_x * c1_aj_x;
  c1_bj_x = c1_v_A;
  c1_cj_x = c1_bj_x;
  c1_dj_x = c1_cj_x;
  c1_v_y = c1_dj_x / 4.0;
  c1_ej_x = c1_q2 + c1_q3;
  c1_fj_x = c1_ej_x;
  c1_fj_x = muDoubleScalarSin(c1_fj_x);
  c1_gj_x = c1_q4;
  c1_hj_x = c1_gj_x;
  c1_hj_x = muDoubleScalarSin(c1_hj_x);
  c1_w_A = 379.0 * c1_fj_x * c1_hj_x;
  c1_ij_x = c1_w_A;
  c1_jj_x = c1_ij_x;
  c1_kj_x = c1_jj_x;
  c1_w_y = c1_kj_x / 4.0;
  c1_lj_x = c1_q1;
  c1_mj_x = c1_lj_x;
  c1_mj_x = muDoubleScalarCos(c1_mj_x);
  c1_nj_x = (c1_q2 + c1_q3) + c1_q4;
  c1_oj_x = c1_nj_x;
  c1_oj_x = muDoubleScalarCos(c1_oj_x);
  c1_x_A = 379.0 * c1_oj_x;
  c1_pj_x = c1_x_A;
  c1_qj_x = c1_pj_x;
  c1_rj_x = c1_qj_x;
  c1_x_y = c1_rj_x / 4.0;
  c1_sj_x = c1_q2 + c1_q3;
  c1_tj_x = c1_sj_x;
  c1_tj_x = muDoubleScalarSin(c1_tj_x);
  c1_uj_x = (c1_q2 + c1_q3) + c1_q4;
  c1_vj_x = c1_uj_x;
  c1_vj_x = muDoubleScalarSin(c1_vj_x);
  c1_wj_x = c1_q5;
  c1_xj_x = c1_wj_x;
  c1_xj_x = muDoubleScalarSin(c1_xj_x);
  c1_y_A = 165.0 * c1_vj_x * c1_xj_x;
  c1_yj_x = c1_y_A;
  c1_ak_x = c1_yj_x;
  c1_bk_x = c1_ak_x;
  c1_y_y = c1_bk_x / 2.0;
  c1_ck_x = c1_q1;
  c1_dk_x = c1_ck_x;
  c1_dk_x = muDoubleScalarCos(c1_dk_x);
  c1_ek_x = (c1_q2 + c1_q3) + c1_q4;
  c1_fk_x = c1_ek_x;
  c1_fk_x = muDoubleScalarCos(c1_fk_x);
  c1_ab_A = 379.0 * c1_fk_x;
  c1_gk_x = c1_ab_A;
  c1_hk_x = c1_gk_x;
  c1_ik_x = c1_hk_x;
  c1_ab_y = c1_ik_x / 4.0;
  c1_jk_x = (c1_q2 + c1_q3) + c1_q4;
  c1_kk_x = c1_jk_x;
  c1_kk_x = muDoubleScalarSin(c1_kk_x);
  c1_lk_x = c1_q5;
  c1_mk_x = c1_lk_x;
  c1_mk_x = muDoubleScalarSin(c1_mk_x);
  c1_bb_A = 165.0 * c1_kk_x * c1_mk_x;
  c1_nk_x = c1_bb_A;
  c1_ok_x = c1_nk_x;
  c1_pk_x = c1_ok_x;
  c1_bb_y = c1_pk_x / 2.0;
  c1_qk_x = c1_q1;
  c1_rk_x = c1_qk_x;
  c1_rk_x = muDoubleScalarSin(c1_rk_x);
  c1_sk_x = c1_q5;
  c1_tk_x = c1_sk_x;
  c1_tk_x = muDoubleScalarSin(c1_tk_x);
  c1_cb_A = 165.0 * c1_rk_x * c1_tk_x;
  c1_uk_x = c1_cb_A;
  c1_vk_x = c1_uk_x;
  c1_wk_x = c1_vk_x;
  c1_cb_y = c1_wk_x / 2.0;
  c1_xk_x = c1_q1;
  c1_yk_x = c1_xk_x;
  c1_yk_x = muDoubleScalarCos(c1_yk_x);
  c1_al_x = c1_q2;
  c1_bl_x = c1_al_x;
  c1_bl_x = muDoubleScalarCos(c1_bl_x);
  c1_cl_x = c1_q3;
  c1_dl_x = c1_cl_x;
  c1_dl_x = muDoubleScalarCos(c1_dl_x);
  c1_el_x = c1_q4;
  c1_fl_x = c1_el_x;
  c1_fl_x = muDoubleScalarCos(c1_fl_x);
  c1_gl_x = c1_q5;
  c1_hl_x = c1_gl_x;
  c1_hl_x = muDoubleScalarCos(c1_hl_x);
  c1_db_A = 165.0 * c1_yk_x * c1_bl_x * c1_dl_x * c1_fl_x * c1_hl_x;
  c1_il_x = c1_db_A;
  c1_jl_x = c1_il_x;
  c1_kl_x = c1_jl_x;
  c1_db_y = c1_kl_x / 2.0;
  c1_ll_x = c1_q1;
  c1_ml_x = c1_ll_x;
  c1_ml_x = muDoubleScalarCos(c1_ml_x);
  c1_nl_x = c1_q2;
  c1_ol_x = c1_nl_x;
  c1_ol_x = muDoubleScalarCos(c1_ol_x);
  c1_pl_x = c1_q5;
  c1_ql_x = c1_pl_x;
  c1_ql_x = muDoubleScalarCos(c1_ql_x);
  c1_rl_x = c1_q3;
  c1_sl_x = c1_rl_x;
  c1_sl_x = muDoubleScalarSin(c1_sl_x);
  c1_tl_x = c1_q4;
  c1_ul_x = c1_tl_x;
  c1_ul_x = muDoubleScalarSin(c1_ul_x);
  c1_eb_A = 165.0 * c1_ml_x * c1_ol_x * c1_ql_x * c1_sl_x * c1_ul_x;
  c1_vl_x = c1_eb_A;
  c1_wl_x = c1_vl_x;
  c1_xl_x = c1_wl_x;
  c1_eb_y = c1_xl_x / 2.0;
  c1_yl_x = c1_q1;
  c1_am_x = c1_yl_x;
  c1_am_x = muDoubleScalarCos(c1_am_x);
  c1_bm_x = c1_q3;
  c1_cm_x = c1_bm_x;
  c1_cm_x = muDoubleScalarCos(c1_cm_x);
  c1_dm_x = c1_q5;
  c1_em_x = c1_dm_x;
  c1_em_x = muDoubleScalarCos(c1_em_x);
  c1_fm_x = c1_q2;
  c1_gm_x = c1_fm_x;
  c1_gm_x = muDoubleScalarSin(c1_gm_x);
  c1_hm_x = c1_q4;
  c1_im_x = c1_hm_x;
  c1_im_x = muDoubleScalarSin(c1_im_x);
  c1_fb_A = 165.0 * c1_am_x * c1_cm_x * c1_em_x * c1_gm_x * c1_im_x;
  c1_jm_x = c1_fb_A;
  c1_km_x = c1_jm_x;
  c1_lm_x = c1_km_x;
  c1_fb_y = c1_lm_x / 2.0;
  c1_mm_x = c1_q1;
  c1_nm_x = c1_mm_x;
  c1_nm_x = muDoubleScalarCos(c1_nm_x);
  c1_om_x = c1_q4;
  c1_pm_x = c1_om_x;
  c1_pm_x = muDoubleScalarCos(c1_pm_x);
  c1_qm_x = c1_q5;
  c1_rm_x = c1_qm_x;
  c1_rm_x = muDoubleScalarCos(c1_rm_x);
  c1_sm_x = c1_q2;
  c1_tm_x = c1_sm_x;
  c1_tm_x = muDoubleScalarSin(c1_tm_x);
  c1_um_x = c1_q3;
  c1_vm_x = c1_um_x;
  c1_vm_x = muDoubleScalarSin(c1_vm_x);
  c1_gb_A = 165.0 * c1_nm_x * c1_pm_x * c1_rm_x * c1_tm_x * c1_vm_x;
  c1_wm_x = c1_gb_A;
  c1_xm_x = c1_wm_x;
  c1_ym_x = c1_xm_x;
  c1_gb_y = c1_ym_x / 2.0;
  c1_an_x = c1_q2 + c1_q3;
  c1_bn_x = c1_an_x;
  c1_bn_x = muDoubleScalarCos(c1_bn_x);
  c1_cn_x = c1_q1;
  c1_dn_x = c1_cn_x;
  c1_dn_x = muDoubleScalarCos(c1_dn_x);
  c1_en_x = c1_q4;
  c1_fn_x = c1_en_x;
  c1_fn_x = muDoubleScalarSin(c1_fn_x);
  c1_hb_A = 379.0 * c1_bn_x * c1_dn_x * c1_fn_x;
  c1_gn_x = c1_hb_A;
  c1_hn_x = c1_gn_x;
  c1_in_x = c1_hn_x;
  c1_hb_y = c1_in_x / 4.0;
  c1_jn_x = c1_q1;
  c1_kn_x = c1_jn_x;
  c1_kn_x = muDoubleScalarCos(c1_kn_x);
  c1_ln_x = c1_q2;
  c1_mn_x = c1_ln_x;
  c1_mn_x = muDoubleScalarCos(c1_mn_x);
  c1_nn_x = c1_q5;
  c1_on_x = c1_nn_x;
  c1_on_x = muDoubleScalarCos(c1_on_x);
  c1_pn_x = c1_q1;
  c1_qn_x = c1_pn_x;
  c1_qn_x = muDoubleScalarSin(c1_qn_x);
  c1_ib_A = 165.0 * c1_on_x * c1_qn_x;
  c1_rn_x = c1_ib_A;
  c1_sn_x = c1_rn_x;
  c1_tn_x = c1_sn_x;
  c1_ib_y = c1_tn_x / 2.0;
  c1_un_x = (c1_q2 + c1_q3) + c1_q4;
  c1_vn_x = c1_un_x;
  c1_vn_x = muDoubleScalarCos(c1_vn_x);
  c1_wn_x = c1_q1;
  c1_xn_x = c1_wn_x;
  c1_xn_x = muDoubleScalarCos(c1_xn_x);
  c1_yn_x = c1_q5;
  c1_ao_x = c1_yn_x;
  c1_ao_x = muDoubleScalarSin(c1_ao_x);
  c1_jb_A = 165.0 * c1_vn_x * c1_xn_x * c1_ao_x;
  c1_bo_x = c1_jb_A;
  c1_co_x = c1_bo_x;
  c1_do_x = c1_co_x;
  c1_jb_y = c1_do_x / 2.0;
  c1_eo_x = c1_q1;
  c1_fo_x = c1_eo_x;
  c1_fo_x = muDoubleScalarSin(c1_fo_x);
  c1_kb_A = 1093.0 * c1_fo_x;
  c1_go_x = c1_kb_A;
  c1_ho_x = c1_go_x;
  c1_io_x = c1_ho_x;
  c1_kb_y = c1_io_x / 10.0;
  c1_jo_x = c1_q2 + c1_q3;
  c1_ko_x = c1_jo_x;
  c1_ko_x = muDoubleScalarSin(c1_ko_x);
  c1_lo_x = c1_q1;
  c1_mo_x = c1_lo_x;
  c1_mo_x = muDoubleScalarCos(c1_mo_x);
  c1_no_x = c1_q4;
  c1_oo_x = c1_no_x;
  c1_oo_x = muDoubleScalarCos(c1_oo_x);
  c1_lb_A = 379.0 * c1_ko_x * c1_mo_x * c1_oo_x;
  c1_po_x = c1_lb_A;
  c1_qo_x = c1_po_x;
  c1_ro_x = c1_qo_x;
  c1_lb_y = c1_ro_x / 4.0;
  c1_so_x = c1_q1;
  c1_to_x = c1_so_x;
  c1_to_x = muDoubleScalarCos(c1_to_x);
  c1_uo_x = c1_q2;
  c1_vo_x = c1_uo_x;
  c1_vo_x = muDoubleScalarCos(c1_vo_x);
  c1_wo_x = c1_q3;
  c1_xo_x = c1_wo_x;
  c1_xo_x = muDoubleScalarCos(c1_xo_x);
  c1_yo_x = c1_q1;
  c1_ap_x = c1_yo_x;
  c1_ap_x = muDoubleScalarCos(c1_ap_x);
  c1_bp_x = c1_q2;
  c1_cp_x = c1_bp_x;
  c1_cp_x = muDoubleScalarSin(c1_cp_x);
  c1_dp_x = c1_q3;
  c1_ep_x = c1_dp_x;
  c1_ep_x = muDoubleScalarSin(c1_ep_x);
  c1_fp_x = c1_q1;
  c1_gp_x = c1_fp_x;
  c1_gp_x = muDoubleScalarSin(c1_gp_x);
  c1_hp_x = c1_q2 + c1_q3;
  c1_ip_x = c1_hp_x;
  c1_ip_x = muDoubleScalarSin(c1_ip_x);
  c1_jp_x = c1_q2;
  c1_kp_x = c1_jp_x;
  c1_kp_x = muDoubleScalarSin(c1_kp_x);
  c1_lp_x = c1_q5;
  c1_mp_x = c1_lp_x;
  c1_mp_x = muDoubleScalarSin(c1_mp_x);
  c1_np_x = c1_q2 + c1_q3;
  c1_op_x = c1_np_x;
  c1_op_x = muDoubleScalarCos(c1_op_x);
  c1_pp_x = c1_q4;
  c1_qp_x = c1_pp_x;
  c1_qp_x = muDoubleScalarSin(c1_qp_x);
  c1_mb_A = 165.0 * c1_op_x * c1_qp_x;
  c1_rp_x = c1_mb_A;
  c1_sp_x = c1_rp_x;
  c1_tp_x = c1_sp_x;
  c1_mb_y = c1_tp_x / 2.0;
  c1_up_x = c1_q2 + c1_q3;
  c1_vp_x = c1_up_x;
  c1_vp_x = muDoubleScalarSin(c1_vp_x);
  c1_wp_x = c1_q4;
  c1_xp_x = c1_wp_x;
  c1_xp_x = muDoubleScalarCos(c1_xp_x);
  c1_nb_A = 165.0 * c1_vp_x * c1_xp_x;
  c1_yp_x = c1_nb_A;
  c1_aq_x = c1_yp_x;
  c1_bq_x = c1_aq_x;
  c1_nb_y = c1_bq_x / 2.0;
  c1_cq_x = c1_q2 + c1_q3;
  c1_dq_x = c1_cq_x;
  c1_dq_x = muDoubleScalarCos(c1_dq_x);
  c1_eq_x = c1_q4;
  c1_fq_x = c1_eq_x;
  c1_fq_x = muDoubleScalarCos(c1_fq_x);
  c1_ob_A = 379.0 * c1_dq_x * c1_fq_x;
  c1_gq_x = c1_ob_A;
  c1_hq_x = c1_gq_x;
  c1_iq_x = c1_hq_x;
  c1_ob_y = c1_iq_x / 4.0;
  c1_jq_x = c1_q2 + c1_q3;
  c1_kq_x = c1_jq_x;
  c1_kq_x = muDoubleScalarSin(c1_kq_x);
  c1_lq_x = c1_q4;
  c1_mq_x = c1_lq_x;
  c1_mq_x = muDoubleScalarSin(c1_mq_x);
  c1_pb_A = 379.0 * c1_kq_x * c1_mq_x;
  c1_nq_x = c1_pb_A;
  c1_oq_x = c1_nq_x;
  c1_pq_x = c1_oq_x;
  c1_pb_y = c1_pq_x / 4.0;
  c1_qq_x = c1_q1;
  c1_rq_x = c1_qq_x;
  c1_rq_x = muDoubleScalarSin(c1_rq_x);
  c1_sq_x = (c1_q2 + c1_q3) + c1_q4;
  c1_tq_x = c1_sq_x;
  c1_tq_x = muDoubleScalarCos(c1_tq_x);
  c1_qb_A = 379.0 * c1_tq_x;
  c1_uq_x = c1_qb_A;
  c1_vq_x = c1_uq_x;
  c1_wq_x = c1_vq_x;
  c1_qb_y = c1_wq_x / 4.0;
  c1_xq_x = c1_q2 + c1_q3;
  c1_yq_x = c1_xq_x;
  c1_yq_x = muDoubleScalarSin(c1_yq_x);
  c1_ar_x = (c1_q2 + c1_q3) + c1_q4;
  c1_br_x = c1_ar_x;
  c1_br_x = muDoubleScalarSin(c1_br_x);
  c1_cr_x = c1_q5;
  c1_dr_x = c1_cr_x;
  c1_dr_x = muDoubleScalarSin(c1_dr_x);
  c1_rb_A = 165.0 * c1_br_x * c1_dr_x;
  c1_er_x = c1_rb_A;
  c1_fr_x = c1_er_x;
  c1_gr_x = c1_fr_x;
  c1_rb_y = c1_gr_x / 2.0;
  c1_hr_x = c1_q1;
  c1_ir_x = c1_hr_x;
  c1_ir_x = muDoubleScalarSin(c1_ir_x);
  c1_jr_x = (c1_q2 + c1_q3) + c1_q4;
  c1_kr_x = c1_jr_x;
  c1_kr_x = muDoubleScalarCos(c1_kr_x);
  c1_sb_A = 379.0 * c1_kr_x;
  c1_lr_x = c1_sb_A;
  c1_mr_x = c1_lr_x;
  c1_nr_x = c1_mr_x;
  c1_sb_y = c1_nr_x / 4.0;
  c1_or_x = (c1_q2 + c1_q3) + c1_q4;
  c1_pr_x = c1_or_x;
  c1_pr_x = muDoubleScalarSin(c1_pr_x);
  c1_qr_x = c1_q5;
  c1_rr_x = c1_qr_x;
  c1_rr_x = muDoubleScalarSin(c1_rr_x);
  c1_tb_A = 165.0 * c1_pr_x * c1_rr_x;
  c1_sr_x = c1_tb_A;
  c1_tr_x = c1_sr_x;
  c1_ur_x = c1_tr_x;
  c1_tb_y = c1_ur_x / 2.0;
  c1_vr_x = c1_q2;
  c1_wr_x = c1_vr_x;
  c1_wr_x = muDoubleScalarCos(c1_wr_x);
  c1_xr_x = c1_q5;
  c1_yr_x = c1_xr_x;
  c1_yr_x = muDoubleScalarCos(c1_yr_x);
  c1_as_x = c1_q1;
  c1_bs_x = c1_as_x;
  c1_bs_x = muDoubleScalarSin(c1_bs_x);
  c1_cs_x = c1_q3;
  c1_ds_x = c1_cs_x;
  c1_ds_x = muDoubleScalarSin(c1_ds_x);
  c1_es_x = c1_q4;
  c1_fs_x = c1_es_x;
  c1_fs_x = muDoubleScalarSin(c1_fs_x);
  c1_ub_A = 165.0 * c1_wr_x * c1_yr_x * c1_bs_x * c1_ds_x * c1_fs_x;
  c1_gs_x = c1_ub_A;
  c1_hs_x = c1_gs_x;
  c1_is_x = c1_hs_x;
  c1_ub_y = c1_is_x / 2.0;
  c1_js_x = c1_q2;
  c1_ks_x = c1_js_x;
  c1_ks_x = muDoubleScalarCos(c1_ks_x);
  c1_ls_x = c1_q3;
  c1_ms_x = c1_ls_x;
  c1_ms_x = muDoubleScalarCos(c1_ms_x);
  c1_ns_x = c1_q4;
  c1_os_x = c1_ns_x;
  c1_os_x = muDoubleScalarCos(c1_os_x);
  c1_ps_x = c1_q5;
  c1_qs_x = c1_ps_x;
  c1_qs_x = muDoubleScalarCos(c1_qs_x);
  c1_rs_x = c1_q1;
  c1_ss_x = c1_rs_x;
  c1_ss_x = muDoubleScalarSin(c1_ss_x);
  c1_vb_A = 165.0 * c1_ks_x * c1_ms_x * c1_os_x * c1_qs_x * c1_ss_x;
  c1_ts_x = c1_vb_A;
  c1_us_x = c1_ts_x;
  c1_vs_x = c1_us_x;
  c1_vb_y = c1_vs_x / 2.0;
  c1_ws_x = c1_q1;
  c1_xs_x = c1_ws_x;
  c1_xs_x = muDoubleScalarCos(c1_xs_x);
  c1_ys_x = c1_q5;
  c1_at_x = c1_ys_x;
  c1_at_x = muDoubleScalarSin(c1_at_x);
  c1_wb_A = 165.0 * c1_xs_x * c1_at_x;
  c1_bt_x = c1_wb_A;
  c1_ct_x = c1_bt_x;
  c1_dt_x = c1_ct_x;
  c1_wb_y = c1_dt_x / 2.0;
  c1_et_x = c1_q3;
  c1_ft_x = c1_et_x;
  c1_ft_x = muDoubleScalarCos(c1_ft_x);
  c1_gt_x = c1_q5;
  c1_ht_x = c1_gt_x;
  c1_ht_x = muDoubleScalarCos(c1_ht_x);
  c1_it_x = c1_q1;
  c1_jt_x = c1_it_x;
  c1_jt_x = muDoubleScalarSin(c1_jt_x);
  c1_kt_x = c1_q2;
  c1_lt_x = c1_kt_x;
  c1_lt_x = muDoubleScalarSin(c1_lt_x);
  c1_mt_x = c1_q4;
  c1_nt_x = c1_mt_x;
  c1_nt_x = muDoubleScalarSin(c1_nt_x);
  c1_xb_A = 165.0 * c1_ft_x * c1_ht_x * c1_jt_x * c1_lt_x * c1_nt_x;
  c1_ot_x = c1_xb_A;
  c1_pt_x = c1_ot_x;
  c1_qt_x = c1_pt_x;
  c1_xb_y = c1_qt_x / 2.0;
  c1_rt_x = c1_q4;
  c1_st_x = c1_rt_x;
  c1_st_x = muDoubleScalarCos(c1_st_x);
  c1_tt_x = c1_q5;
  c1_ut_x = c1_tt_x;
  c1_ut_x = muDoubleScalarCos(c1_ut_x);
  c1_vt_x = c1_q1;
  c1_wt_x = c1_vt_x;
  c1_wt_x = muDoubleScalarSin(c1_wt_x);
  c1_xt_x = c1_q2;
  c1_yt_x = c1_xt_x;
  c1_yt_x = muDoubleScalarSin(c1_yt_x);
  c1_au_x = c1_q3;
  c1_bu_x = c1_au_x;
  c1_bu_x = muDoubleScalarSin(c1_bu_x);
  c1_yb_A = 165.0 * c1_st_x * c1_ut_x * c1_wt_x * c1_yt_x * c1_bu_x;
  c1_cu_x = c1_yb_A;
  c1_du_x = c1_cu_x;
  c1_eu_x = c1_du_x;
  c1_yb_y = c1_eu_x / 2.0;
  c1_fu_x = ((c1_q2 + c1_q3) + c1_q4) + c1_q5;
  c1_gu_x = c1_fu_x;
  c1_gu_x = muDoubleScalarSin(c1_gu_x);
  c1_ac_A = 165.0 * c1_gu_x;
  c1_hu_x = c1_ac_A;
  c1_iu_x = c1_hu_x;
  c1_ju_x = c1_iu_x;
  c1_ac_y = c1_ju_x / 4.0;
  c1_ku_x = (c1_q2 + c1_q3) + c1_q4;
  c1_lu_x = c1_ku_x;
  c1_lu_x = muDoubleScalarSin(c1_lu_x);
  c1_bc_A = 379.0 * c1_lu_x;
  c1_mu_x = c1_bc_A;
  c1_nu_x = c1_mu_x;
  c1_ou_x = c1_nu_x;
  c1_bc_y = c1_ou_x / 4.0;
  c1_pu_x = c1_q2 + c1_q3;
  c1_qu_x = c1_pu_x;
  c1_qu_x = muDoubleScalarCos(c1_qu_x);
  c1_ru_x = c1_q2;
  c1_su_x = c1_ru_x;
  c1_su_x = muDoubleScalarCos(c1_su_x);
  c1_tu_x = ((c1_q2 + c1_q3) + c1_q4) - c1_q5;
  c1_uu_x = c1_tu_x;
  c1_uu_x = muDoubleScalarSin(c1_uu_x);
  c1_cc_A = 165.0 * c1_uu_x;
  c1_vu_x = c1_cc_A;
  c1_wu_x = c1_vu_x;
  c1_xu_x = c1_wu_x;
  c1_cc_y = c1_xu_x / 4.0;
  c1_yu_x = ((c1_q2 + c1_q3) + c1_q4) + c1_q5;
  c1_av_x = c1_yu_x;
  c1_av_x = muDoubleScalarSin(c1_av_x);
  c1_dc_A = 165.0 * c1_av_x;
  c1_bv_x = c1_dc_A;
  c1_cv_x = c1_bv_x;
  c1_dv_x = c1_cv_x;
  c1_dc_y = c1_dv_x / 4.0;
  c1_ev_x = (c1_q2 + c1_q3) + c1_q4;
  c1_fv_x = c1_ev_x;
  c1_fv_x = muDoubleScalarSin(c1_fv_x);
  c1_ec_A = 379.0 * c1_fv_x;
  c1_gv_x = c1_ec_A;
  c1_hv_x = c1_gv_x;
  c1_iv_x = c1_hv_x;
  c1_ec_y = c1_iv_x / 4.0;
  c1_jv_x = c1_q2 + c1_q3;
  c1_kv_x = c1_jv_x;
  c1_kv_x = muDoubleScalarCos(c1_kv_x);
  c1_lv_x = ((c1_q2 + c1_q3) + c1_q4) - c1_q5;
  c1_mv_x = c1_lv_x;
  c1_mv_x = muDoubleScalarSin(c1_mv_x);
  c1_fc_A = 165.0 * c1_mv_x;
  c1_nv_x = c1_fc_A;
  c1_ov_x = c1_nv_x;
  c1_pv_x = c1_ov_x;
  c1_fc_y = c1_pv_x / 4.0;
  c1_qv_x = ((c1_q2 + c1_q3) + c1_q4) + c1_q5;
  c1_rv_x = c1_qv_x;
  c1_rv_x = muDoubleScalarSin(c1_rv_x);
  c1_gc_A = 165.0 * c1_rv_x;
  c1_sv_x = c1_gc_A;
  c1_tv_x = c1_sv_x;
  c1_uv_x = c1_tv_x;
  c1_gc_y = c1_uv_x / 4.0;
  c1_vv_x = (c1_q2 + c1_q3) + c1_q4;
  c1_wv_x = c1_vv_x;
  c1_wv_x = muDoubleScalarSin(c1_wv_x);
  c1_hc_A = 379.0 * c1_wv_x;
  c1_xv_x = c1_hc_A;
  c1_yv_x = c1_xv_x;
  c1_aw_x = c1_yv_x;
  c1_hc_y = c1_aw_x / 4.0;
  c1_bw_x = ((c1_q2 + c1_q3) + c1_q4) - c1_q5;
  c1_cw_x = c1_bw_x;
  c1_cw_x = muDoubleScalarSin(c1_cw_x);
  c1_ic_A = 165.0 * c1_cw_x;
  c1_dw_x = c1_ic_A;
  c1_ew_x = c1_dw_x;
  c1_fw_x = c1_ew_x;
  c1_ic_y = c1_fw_x / 4.0;
  c1_gw_x = ((c1_q2 + c1_q3) + c1_q4) + c1_q5;
  c1_hw_x = c1_gw_x;
  c1_hw_x = muDoubleScalarSin(c1_hw_x);
  c1_jc_A = 165.0 * c1_hw_x;
  c1_iw_x = c1_jc_A;
  c1_jw_x = c1_iw_x;
  c1_kw_x = c1_jw_x;
  c1_jc_y = c1_kw_x / 4.0;
  c1_lw_x = ((c1_q2 + c1_q3) + c1_q4) - c1_q5;
  c1_mw_x = c1_lw_x;
  c1_mw_x = muDoubleScalarSin(c1_mw_x);
  c1_kc_A = 165.0 * c1_mw_x;
  c1_nw_x = c1_kc_A;
  c1_ow_x = c1_nw_x;
  c1_pw_x = c1_ow_x;
  c1_kc_y = c1_pw_x / 4.0;
  c1_qw_x = c1_q1;
  c1_rw_x = c1_qw_x;
  c1_rw_x = muDoubleScalarSin(c1_rw_x);
  c1_sw_x = c1_q1;
  c1_tw_x = c1_sw_x;
  c1_tw_x = muDoubleScalarSin(c1_tw_x);
  c1_uw_x = c1_q1;
  c1_vw_x = c1_uw_x;
  c1_vw_x = muDoubleScalarSin(c1_vw_x);
  c1_ww_x = (c1_q2 + c1_q3) + c1_q4;
  c1_xw_x = c1_ww_x;
  c1_xw_x = muDoubleScalarSin(c1_xw_x);
  c1_yw_x = c1_q1;
  c1_ax_x = c1_yw_x;
  c1_ax_x = muDoubleScalarCos(c1_ax_x);
  c1_bx_x = c1_q5;
  c1_cx_x = c1_bx_x;
  c1_cx_x = muDoubleScalarCos(c1_cx_x);
  c1_dx_x = c1_q1;
  c1_ex_x = c1_dx_x;
  c1_ex_x = muDoubleScalarSin(c1_ex_x);
  c1_fx_x = (c1_q2 + c1_q3) + c1_q4;
  c1_gx_x = c1_fx_x;
  c1_gx_x = muDoubleScalarCos(c1_gx_x);
  c1_hx_x = c1_q1;
  c1_ix_x = c1_hx_x;
  c1_ix_x = muDoubleScalarCos(c1_ix_x);
  c1_jx_x = c1_q5;
  c1_kx_x = c1_jx_x;
  c1_kx_x = muDoubleScalarSin(c1_kx_x);
  c1_lx_x = c1_q1;
  c1_mx_x = c1_lx_x;
  c1_mx_x = muDoubleScalarCos(c1_mx_x);
  c1_nx_x = c1_q1;
  c1_ox_x = c1_nx_x;
  c1_ox_x = muDoubleScalarCos(c1_ox_x);
  c1_px_x = c1_q1;
  c1_qx_x = c1_px_x;
  c1_qx_x = muDoubleScalarCos(c1_qx_x);
  c1_rx_x = (c1_q2 + c1_q3) + c1_q4;
  c1_sx_x = c1_rx_x;
  c1_sx_x = muDoubleScalarSin(c1_sx_x);
  c1_tx_x = c1_q1;
  c1_ux_x = c1_tx_x;
  c1_ux_x = muDoubleScalarSin(c1_ux_x);
  c1_vx_x = c1_q1;
  c1_wx_x = c1_vx_x;
  c1_wx_x = muDoubleScalarCos(c1_wx_x);
  c1_xx_x = c1_q5;
  c1_yx_x = c1_xx_x;
  c1_yx_x = muDoubleScalarCos(c1_yx_x);
  c1_ay_x = (c1_q2 + c1_q3) + c1_q4;
  c1_by_x = c1_ay_x;
  c1_by_x = muDoubleScalarCos(c1_by_x);
  c1_cy_x = c1_q1;
  c1_dy_x = c1_cy_x;
  c1_dy_x = muDoubleScalarSin(c1_dy_x);
  c1_ey_x = c1_q5;
  c1_fy_x = c1_ey_x;
  c1_fy_x = muDoubleScalarSin(c1_fy_x);
  c1_gy_x = (c1_q2 + c1_q3) + c1_q4;
  c1_hy_x = c1_gy_x;
  c1_hy_x = muDoubleScalarCos(c1_hy_x);
  c1_iy_x = (c1_q2 + c1_q3) + c1_q4;
  c1_jy_x = c1_iy_x;
  c1_jy_x = muDoubleScalarSin(c1_jy_x);
  c1_ky_x = c1_q5;
  c1_ly_x = c1_ky_x;
  c1_ly_x = muDoubleScalarSin(c1_ly_x);
  c1_J[0] = ((((((425.0 * c1_vf_x * c1_xf_x - c1_o_y) - c1_p_y) - 392.0 *
                c1_mg_x * c1_og_x * c1_qg_x) + c1_q_y) - c1_r_y) - c1_s_y) +
    392.0 * c1_uh_x * c1_wh_x * c1_yh_x;
  c1_J[6] = c1_bi_x * ((((392.0 * c1_di_x + 425.0 * c1_fi_x) + c1_hi_x * (c1_t_y
    + c1_u_y)) + c1_v_y) - c1_w_y);
  c1_J[12] = c1_mj_x * ((c1_x_y + 392.0 * c1_tj_x) + c1_y_y);
  c1_J[18] = c1_dk_x * (c1_ab_y + c1_bb_y);
  c1_J[24] = (((c1_cb_y - c1_db_y) + c1_eb_y) + c1_fb_y) + c1_gb_y;
  c1_J[30] = 0.0;
  c1_J[1] = ((((((c1_hb_y - 425.0 * c1_kn_x * c1_mn_x) - c1_ib_y) - c1_jb_y) -
               c1_kb_y) + c1_lb_y) - 392.0 * c1_to_x * c1_vo_x * c1_xo_x) +
    392.0 * c1_ap_x * c1_cp_x * c1_ep_x;
  c1_J[7] = c1_gp_x * ((((392.0 * c1_ip_x + 425.0 * c1_kp_x) + c1_mp_x *
    (c1_mb_y + c1_nb_y)) + c1_ob_y) - c1_pb_y);
  c1_J[13] = c1_rq_x * ((c1_qb_y + 392.0 * c1_yq_x) + c1_rb_y);
  c1_J[19] = c1_ir_x * (c1_sb_y + c1_tb_y);
  c1_J[25] = (((c1_ub_y - c1_vb_y) - c1_wb_y) + c1_xb_y) + c1_yb_y;
  c1_J[31] = 0.0;
  c1_J[2] = 0.0;
  c1_J[8] = (((c1_ac_y - c1_bc_y) + 392.0 * c1_qu_x) + 425.0 * c1_su_x) -
    c1_cc_y;
  c1_J[14] = ((c1_dc_y - c1_ec_y) + 392.0 * c1_kv_x) - c1_fc_y;
  c1_J[20] = (c1_gc_y - c1_hc_y) - c1_ic_y;
  c1_J[26] = c1_jc_y + c1_kc_y;
  c1_J[32] = 0.0;
  c1_J[3] = 0.0;
  c1_J[9] = -c1_rw_x;
  c1_J[15] = -c1_tw_x;
  c1_J[21] = -c1_vw_x;
  c1_J[27] = c1_xw_x * c1_ax_x;
  c1_J[33] = -c1_cx_x * c1_ex_x - c1_gx_x * c1_ix_x * c1_kx_x;
  c1_J[4] = 0.0;
  c1_J[10] = c1_mx_x;
  c1_J[16] = c1_ox_x;
  c1_J[22] = c1_qx_x;
  c1_J[28] = c1_sx_x * c1_ux_x;
  c1_J[34] = c1_wx_x * c1_yx_x - c1_by_x * c1_dy_x * c1_fy_x;
  c1_J[5] = 1.0;
  c1_J[11] = 0.0;
  c1_J[17] = 0.0;
  c1_J[23] = 0.0;
  c1_J[29] = c1_hy_x;
  c1_J[35] = c1_jy_x * c1_ly_x;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 24);
  for (c1_i9 = 0; c1_i9 < 36; c1_i9++) {
    c1_k[c1_i9] = c1_a[c1_i9];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 25);
  for (c1_i10 = 0; c1_i10 < 6; c1_i10++) {
    c1_b[c1_i10] = c1_x_des[c1_i10] - c1_x[c1_i10];
  }

  c1_c_eml_scalar_eg(chartInstance);
  c1_c_eml_scalar_eg(chartInstance);
  c1_e_threshold(chartInstance);
  for (c1_i11 = 0; c1_i11 < 6; c1_i11++) {
    c1_lc_y[c1_i11] = 0.0;
    c1_i12 = 0;
    for (c1_i13 = 0; c1_i13 < 6; c1_i13++) {
      c1_lc_y[c1_i11] += c1_a[c1_i12 + c1_i11] * c1_b[c1_i13];
      c1_i12 += 6;
    }
  }

  for (c1_i14 = 0; c1_i14 < 36; c1_i14++) {
    c1_b_J[c1_i14] = c1_J[c1_i14];
  }

  c1_pinv(chartInstance, c1_b_J, c1_b_a);
  for (c1_i15 = 0; c1_i15 < 6; c1_i15++) {
    c1_lc_y[c1_i15] += c1_xd_des[c1_i15];
  }

  c1_c_eml_scalar_eg(chartInstance);
  c1_c_eml_scalar_eg(chartInstance);
  for (c1_i16 = 0; c1_i16 < 6; c1_i16++) {
    c1_dq[c1_i16] = 0.0;
  }

  for (c1_i17 = 0; c1_i17 < 6; c1_i17++) {
    c1_dq[c1_i17] = 0.0;
  }

  for (c1_i18 = 0; c1_i18 < 6; c1_i18++) {
    c1_b[c1_i18] = c1_dq[c1_i18];
  }

  for (c1_i19 = 0; c1_i19 < 6; c1_i19++) {
    c1_dq[c1_i19] = c1_b[c1_i19];
  }

  c1_e_threshold(chartInstance);
  for (c1_i20 = 0; c1_i20 < 6; c1_i20++) {
    c1_b[c1_i20] = c1_dq[c1_i20];
  }

  for (c1_i21 = 0; c1_i21 < 6; c1_i21++) {
    c1_dq[c1_i21] = c1_b[c1_i21];
  }

  for (c1_i22 = 0; c1_i22 < 6; c1_i22++) {
    c1_dq[c1_i22] = 0.0;
    c1_i23 = 0;
    for (c1_i24 = 0; c1_i24 < 6; c1_i24++) {
      c1_dq[c1_i22] += c1_b_a[c1_i23 + c1_i22] * c1_lc_y[c1_i24];
      c1_i23 += 6;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -25);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i25 = 0; c1_i25 < 6; c1_i25++) {
    (*c1_b_dq)[c1_i25] = c1_dq[c1_i25];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_UR5Model(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  (void)c1_chartNumber;
  (void)c1_instanceNumber;
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i26;
  real_T c1_b_inData[6];
  int32_T c1_i27;
  real_T c1_u[6];
  const mxArray *c1_y = NULL;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i26 = 0; c1_i26 < 6; c1_i26++) {
    c1_b_inData[c1_i26] = (*(real_T (*)[6])c1_inData)[c1_i26];
  }

  for (c1_i27 = 0; c1_i27 < 6; c1_i27++) {
    c1_u[c1_i27] = c1_b_inData[c1_i27];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_dq, const char_T *c1_identifier, real_T c1_y[6])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_dq), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_dq);
}

static void c1_b_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6])
{
  real_T c1_dv1[6];
  int32_T c1_i28;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv1, 1, 0, 0U, 1, 0U, 1, 6);
  for (c1_i28 = 0; c1_i28 < 6; c1_i28++) {
    c1_y[c1_i28] = c1_dv1[c1_i28];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_dq;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[6];
  int32_T c1_i29;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_dq = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_dq), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_dq);
  for (c1_i29 = 0; c1_i29 < 6; c1_i29++) {
    (*(real_T (*)[6])c1_outData)[c1_i29] = c1_y[c1_i29];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_c_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i30;
  int32_T c1_i31;
  int32_T c1_i32;
  real_T c1_b_inData[36];
  int32_T c1_i33;
  int32_T c1_i34;
  int32_T c1_i35;
  real_T c1_u[36];
  const mxArray *c1_y = NULL;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i30 = 0;
  for (c1_i31 = 0; c1_i31 < 6; c1_i31++) {
    for (c1_i32 = 0; c1_i32 < 6; c1_i32++) {
      c1_b_inData[c1_i32 + c1_i30] = (*(real_T (*)[36])c1_inData)[c1_i32 +
        c1_i30];
    }

    c1_i30 += 6;
  }

  c1_i33 = 0;
  for (c1_i34 = 0; c1_i34 < 6; c1_i34++) {
    for (c1_i35 = 0; c1_i35 < 6; c1_i35++) {
      c1_u[c1_i35 + c1_i33] = c1_b_inData[c1_i35 + c1_i33];
    }

    c1_i33 += 6;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[36])
{
  real_T c1_dv2[36];
  int32_T c1_i36;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv2, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c1_i36 = 0; c1_i36 < 36; c1_i36++) {
    c1_y[c1_i36] = c1_dv2[c1_i36];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_J;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[36];
  int32_T c1_i37;
  int32_T c1_i38;
  int32_T c1_i39;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_J = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_J), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_J);
  c1_i37 = 0;
  for (c1_i38 = 0; c1_i38 < 6; c1_i38++) {
    for (c1_i39 = 0; c1_i39 < 6; c1_i39++) {
      (*(real_T (*)[36])c1_outData)[c1_i39 + c1_i37] = c1_y[c1_i39 + c1_i37];
    }

    c1_i37 += 6;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_UR5Model_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 210, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  c1_b_info_helper(&c1_nameCaptureInfo);
  c1_c_info_helper(&c1_nameCaptureInfo);
  c1_d_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  const mxArray *c1_rhs17 = NULL;
  const mxArray *c1_lhs17 = NULL;
  const mxArray *c1_rhs18 = NULL;
  const mxArray *c1_lhs18 = NULL;
  const mxArray *c1_rhs19 = NULL;
  const mxArray *c1_lhs19 = NULL;
  const mxArray *c1_rhs20 = NULL;
  const mxArray *c1_lhs20 = NULL;
  const mxArray *c1_rhs21 = NULL;
  const mxArray *c1_lhs21 = NULL;
  const mxArray *c1_rhs22 = NULL;
  const mxArray *c1_lhs22 = NULL;
  const mxArray *c1_rhs23 = NULL;
  const mxArray *c1_lhs23 = NULL;
  const mxArray *c1_rhs24 = NULL;
  const mxArray *c1_lhs24 = NULL;
  const mxArray *c1_rhs25 = NULL;
  const mxArray *c1_lhs25 = NULL;
  const mxArray *c1_rhs26 = NULL;
  const mxArray *c1_lhs26 = NULL;
  const mxArray *c1_rhs27 = NULL;
  const mxArray *c1_lhs27 = NULL;
  const mxArray *c1_rhs28 = NULL;
  const mxArray *c1_lhs28 = NULL;
  const mxArray *c1_rhs29 = NULL;
  const mxArray *c1_lhs29 = NULL;
  const mxArray *c1_rhs30 = NULL;
  const mxArray *c1_lhs30 = NULL;
  const mxArray *c1_rhs31 = NULL;
  const mxArray *c1_lhs31 = NULL;
  const mxArray *c1_rhs32 = NULL;
  const mxArray *c1_lhs32 = NULL;
  const mxArray *c1_rhs33 = NULL;
  const mxArray *c1_lhs33 = NULL;
  const mxArray *c1_rhs34 = NULL;
  const mxArray *c1_lhs34 = NULL;
  const mxArray *c1_rhs35 = NULL;
  const mxArray *c1_lhs35 = NULL;
  const mxArray *c1_rhs36 = NULL;
  const mxArray *c1_lhs36 = NULL;
  const mxArray *c1_rhs37 = NULL;
  const mxArray *c1_lhs37 = NULL;
  const mxArray *c1_rhs38 = NULL;
  const mxArray *c1_lhs38 = NULL;
  const mxArray *c1_rhs39 = NULL;
  const mxArray *c1_lhs39 = NULL;
  const mxArray *c1_rhs40 = NULL;
  const mxArray *c1_lhs40 = NULL;
  const mxArray *c1_rhs41 = NULL;
  const mxArray *c1_lhs41 = NULL;
  const mxArray *c1_rhs42 = NULL;
  const mxArray *c1_lhs42 = NULL;
  const mxArray *c1_rhs43 = NULL;
  const mxArray *c1_lhs43 = NULL;
  const mxArray *c1_rhs44 = NULL;
  const mxArray *c1_lhs44 = NULL;
  const mxArray *c1_rhs45 = NULL;
  const mxArray *c1_lhs45 = NULL;
  const mxArray *c1_rhs46 = NULL;
  const mxArray *c1_lhs46 = NULL;
  const mxArray *c1_rhs47 = NULL;
  const mxArray *c1_lhs47 = NULL;
  const mxArray *c1_rhs48 = NULL;
  const mxArray *c1_lhs48 = NULL;
  const mxArray *c1_rhs49 = NULL;
  const mxArray *c1_lhs49 = NULL;
  const mxArray *c1_rhs50 = NULL;
  const mxArray *c1_lhs50 = NULL;
  const mxArray *c1_rhs51 = NULL;
  const mxArray *c1_lhs51 = NULL;
  const mxArray *c1_rhs52 = NULL;
  const mxArray *c1_lhs52 = NULL;
  const mxArray *c1_rhs53 = NULL;
  const mxArray *c1_lhs53 = NULL;
  const mxArray *c1_rhs54 = NULL;
  const mxArray *c1_lhs54 = NULL;
  const mxArray *c1_rhs55 = NULL;
  const mxArray *c1_lhs55 = NULL;
  const mxArray *c1_rhs56 = NULL;
  const mxArray *c1_lhs56 = NULL;
  const mxArray *c1_rhs57 = NULL;
  const mxArray *c1_lhs57 = NULL;
  const mxArray *c1_rhs58 = NULL;
  const mxArray *c1_lhs58 = NULL;
  const mxArray *c1_rhs59 = NULL;
  const mxArray *c1_lhs59 = NULL;
  const mxArray *c1_rhs60 = NULL;
  const mxArray *c1_lhs60 = NULL;
  const mxArray *c1_rhs61 = NULL;
  const mxArray *c1_lhs61 = NULL;
  const mxArray *c1_rhs62 = NULL;
  const mxArray *c1_lhs62 = NULL;
  const mxArray *c1_rhs63 = NULL;
  const mxArray *c1_lhs63 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837572U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825922U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825936U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825996U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eye"), "name", "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857498U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1368190230U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c1_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c1_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c1_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c1_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c1_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c1_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c1_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c1_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c1_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c1_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c1_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c1_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383880894U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c1_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c1_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("pinv"), "name", "name", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826028U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c1_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c1_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c1_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c1_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("svd"), "name", "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826032U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c1_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c1_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c1_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isfinite"), "name", "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c1_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c1_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c1_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c1_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c1_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c1_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/svd.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgesvd"), "name", "name",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgesvd.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826006U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c1_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgesvd.m"),
                  "context", "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_lapack_xgesvd"), "name",
                  "name", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgesvd.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286826010U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c1_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgesvd.m"),
                  "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zsvdc"), "name",
                  "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1295288466U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c1_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c1_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c1_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c1_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c1_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("min"), "name", "name", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c1_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c1_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c1_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c1_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c1_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c1_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c1_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c1_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c1_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("max"), "name", "name", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c1_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "context",
                  "context", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378303184U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c1_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c1_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c1_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
  sf_mex_destroy(&c1_rhs17);
  sf_mex_destroy(&c1_lhs17);
  sf_mex_destroy(&c1_rhs18);
  sf_mex_destroy(&c1_lhs18);
  sf_mex_destroy(&c1_rhs19);
  sf_mex_destroy(&c1_lhs19);
  sf_mex_destroy(&c1_rhs20);
  sf_mex_destroy(&c1_lhs20);
  sf_mex_destroy(&c1_rhs21);
  sf_mex_destroy(&c1_lhs21);
  sf_mex_destroy(&c1_rhs22);
  sf_mex_destroy(&c1_lhs22);
  sf_mex_destroy(&c1_rhs23);
  sf_mex_destroy(&c1_lhs23);
  sf_mex_destroy(&c1_rhs24);
  sf_mex_destroy(&c1_lhs24);
  sf_mex_destroy(&c1_rhs25);
  sf_mex_destroy(&c1_lhs25);
  sf_mex_destroy(&c1_rhs26);
  sf_mex_destroy(&c1_lhs26);
  sf_mex_destroy(&c1_rhs27);
  sf_mex_destroy(&c1_lhs27);
  sf_mex_destroy(&c1_rhs28);
  sf_mex_destroy(&c1_lhs28);
  sf_mex_destroy(&c1_rhs29);
  sf_mex_destroy(&c1_lhs29);
  sf_mex_destroy(&c1_rhs30);
  sf_mex_destroy(&c1_lhs30);
  sf_mex_destroy(&c1_rhs31);
  sf_mex_destroy(&c1_lhs31);
  sf_mex_destroy(&c1_rhs32);
  sf_mex_destroy(&c1_lhs32);
  sf_mex_destroy(&c1_rhs33);
  sf_mex_destroy(&c1_lhs33);
  sf_mex_destroy(&c1_rhs34);
  sf_mex_destroy(&c1_lhs34);
  sf_mex_destroy(&c1_rhs35);
  sf_mex_destroy(&c1_lhs35);
  sf_mex_destroy(&c1_rhs36);
  sf_mex_destroy(&c1_lhs36);
  sf_mex_destroy(&c1_rhs37);
  sf_mex_destroy(&c1_lhs37);
  sf_mex_destroy(&c1_rhs38);
  sf_mex_destroy(&c1_lhs38);
  sf_mex_destroy(&c1_rhs39);
  sf_mex_destroy(&c1_lhs39);
  sf_mex_destroy(&c1_rhs40);
  sf_mex_destroy(&c1_lhs40);
  sf_mex_destroy(&c1_rhs41);
  sf_mex_destroy(&c1_lhs41);
  sf_mex_destroy(&c1_rhs42);
  sf_mex_destroy(&c1_lhs42);
  sf_mex_destroy(&c1_rhs43);
  sf_mex_destroy(&c1_lhs43);
  sf_mex_destroy(&c1_rhs44);
  sf_mex_destroy(&c1_lhs44);
  sf_mex_destroy(&c1_rhs45);
  sf_mex_destroy(&c1_lhs45);
  sf_mex_destroy(&c1_rhs46);
  sf_mex_destroy(&c1_lhs46);
  sf_mex_destroy(&c1_rhs47);
  sf_mex_destroy(&c1_lhs47);
  sf_mex_destroy(&c1_rhs48);
  sf_mex_destroy(&c1_lhs48);
  sf_mex_destroy(&c1_rhs49);
  sf_mex_destroy(&c1_lhs49);
  sf_mex_destroy(&c1_rhs50);
  sf_mex_destroy(&c1_lhs50);
  sf_mex_destroy(&c1_rhs51);
  sf_mex_destroy(&c1_lhs51);
  sf_mex_destroy(&c1_rhs52);
  sf_mex_destroy(&c1_lhs52);
  sf_mex_destroy(&c1_rhs53);
  sf_mex_destroy(&c1_lhs53);
  sf_mex_destroy(&c1_rhs54);
  sf_mex_destroy(&c1_lhs54);
  sf_mex_destroy(&c1_rhs55);
  sf_mex_destroy(&c1_lhs55);
  sf_mex_destroy(&c1_rhs56);
  sf_mex_destroy(&c1_lhs56);
  sf_mex_destroy(&c1_rhs57);
  sf_mex_destroy(&c1_lhs57);
  sf_mex_destroy(&c1_rhs58);
  sf_mex_destroy(&c1_lhs58);
  sf_mex_destroy(&c1_rhs59);
  sf_mex_destroy(&c1_lhs59);
  sf_mex_destroy(&c1_rhs60);
  sf_mex_destroy(&c1_lhs60);
  sf_mex_destroy(&c1_rhs61);
  sf_mex_destroy(&c1_lhs61);
  sf_mex_destroy(&c1_rhs62);
  sf_mex_destroy(&c1_lhs62);
  sf_mex_destroy(&c1_rhs63);
  sf_mex_destroy(&c1_lhs63);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static void c1_b_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs64 = NULL;
  const mxArray *c1_lhs64 = NULL;
  const mxArray *c1_rhs65 = NULL;
  const mxArray *c1_lhs65 = NULL;
  const mxArray *c1_rhs66 = NULL;
  const mxArray *c1_lhs66 = NULL;
  const mxArray *c1_rhs67 = NULL;
  const mxArray *c1_lhs67 = NULL;
  const mxArray *c1_rhs68 = NULL;
  const mxArray *c1_lhs68 = NULL;
  const mxArray *c1_rhs69 = NULL;
  const mxArray *c1_lhs69 = NULL;
  const mxArray *c1_rhs70 = NULL;
  const mxArray *c1_lhs70 = NULL;
  const mxArray *c1_rhs71 = NULL;
  const mxArray *c1_lhs71 = NULL;
  const mxArray *c1_rhs72 = NULL;
  const mxArray *c1_lhs72 = NULL;
  const mxArray *c1_rhs73 = NULL;
  const mxArray *c1_lhs73 = NULL;
  const mxArray *c1_rhs74 = NULL;
  const mxArray *c1_lhs74 = NULL;
  const mxArray *c1_rhs75 = NULL;
  const mxArray *c1_lhs75 = NULL;
  const mxArray *c1_rhs76 = NULL;
  const mxArray *c1_lhs76 = NULL;
  const mxArray *c1_rhs77 = NULL;
  const mxArray *c1_lhs77 = NULL;
  const mxArray *c1_rhs78 = NULL;
  const mxArray *c1_lhs78 = NULL;
  const mxArray *c1_rhs79 = NULL;
  const mxArray *c1_lhs79 = NULL;
  const mxArray *c1_rhs80 = NULL;
  const mxArray *c1_lhs80 = NULL;
  const mxArray *c1_rhs81 = NULL;
  const mxArray *c1_lhs81 = NULL;
  const mxArray *c1_rhs82 = NULL;
  const mxArray *c1_lhs82 = NULL;
  const mxArray *c1_rhs83 = NULL;
  const mxArray *c1_lhs83 = NULL;
  const mxArray *c1_rhs84 = NULL;
  const mxArray *c1_lhs84 = NULL;
  const mxArray *c1_rhs85 = NULL;
  const mxArray *c1_lhs85 = NULL;
  const mxArray *c1_rhs86 = NULL;
  const mxArray *c1_lhs86 = NULL;
  const mxArray *c1_rhs87 = NULL;
  const mxArray *c1_lhs87 = NULL;
  const mxArray *c1_rhs88 = NULL;
  const mxArray *c1_lhs88 = NULL;
  const mxArray *c1_rhs89 = NULL;
  const mxArray *c1_lhs89 = NULL;
  const mxArray *c1_rhs90 = NULL;
  const mxArray *c1_lhs90 = NULL;
  const mxArray *c1_rhs91 = NULL;
  const mxArray *c1_lhs91 = NULL;
  const mxArray *c1_rhs92 = NULL;
  const mxArray *c1_lhs92 = NULL;
  const mxArray *c1_rhs93 = NULL;
  const mxArray *c1_lhs93 = NULL;
  const mxArray *c1_rhs94 = NULL;
  const mxArray *c1_lhs94 = NULL;
  const mxArray *c1_rhs95 = NULL;
  const mxArray *c1_lhs95 = NULL;
  const mxArray *c1_rhs96 = NULL;
  const mxArray *c1_lhs96 = NULL;
  const mxArray *c1_rhs97 = NULL;
  const mxArray *c1_lhs97 = NULL;
  const mxArray *c1_rhs98 = NULL;
  const mxArray *c1_lhs98 = NULL;
  const mxArray *c1_rhs99 = NULL;
  const mxArray *c1_lhs99 = NULL;
  const mxArray *c1_rhs100 = NULL;
  const mxArray *c1_lhs100 = NULL;
  const mxArray *c1_rhs101 = NULL;
  const mxArray *c1_lhs101 = NULL;
  const mxArray *c1_rhs102 = NULL;
  const mxArray *c1_lhs102 = NULL;
  const mxArray *c1_rhs103 = NULL;
  const mxArray *c1_lhs103 = NULL;
  const mxArray *c1_rhs104 = NULL;
  const mxArray *c1_lhs104 = NULL;
  const mxArray *c1_rhs105 = NULL;
  const mxArray *c1_lhs105 = NULL;
  const mxArray *c1_rhs106 = NULL;
  const mxArray *c1_lhs106 = NULL;
  const mxArray *c1_rhs107 = NULL;
  const mxArray *c1_lhs107 = NULL;
  const mxArray *c1_rhs108 = NULL;
  const mxArray *c1_lhs108 = NULL;
  const mxArray *c1_rhs109 = NULL;
  const mxArray *c1_lhs109 = NULL;
  const mxArray *c1_rhs110 = NULL;
  const mxArray *c1_lhs110 = NULL;
  const mxArray *c1_rhs111 = NULL;
  const mxArray *c1_lhs111 = NULL;
  const mxArray *c1_rhs112 = NULL;
  const mxArray *c1_lhs112 = NULL;
  const mxArray *c1_rhs113 = NULL;
  const mxArray *c1_lhs113 = NULL;
  const mxArray *c1_rhs114 = NULL;
  const mxArray *c1_lhs114 = NULL;
  const mxArray *c1_rhs115 = NULL;
  const mxArray *c1_lhs115 = NULL;
  const mxArray *c1_rhs116 = NULL;
  const mxArray *c1_lhs116 = NULL;
  const mxArray *c1_rhs117 = NULL;
  const mxArray *c1_lhs117 = NULL;
  const mxArray *c1_rhs118 = NULL;
  const mxArray *c1_lhs118 = NULL;
  const mxArray *c1_rhs119 = NULL;
  const mxArray *c1_lhs119 = NULL;
  const mxArray *c1_rhs120 = NULL;
  const mxArray *c1_lhs120 = NULL;
  const mxArray *c1_rhs121 = NULL;
  const mxArray *c1_lhs121 = NULL;
  const mxArray *c1_rhs122 = NULL;
  const mxArray *c1_lhs122 = NULL;
  const mxArray *c1_rhs123 = NULL;
  const mxArray *c1_lhs123 = NULL;
  const mxArray *c1_rhs124 = NULL;
  const mxArray *c1_lhs124 = NULL;
  const mxArray *c1_rhs125 = NULL;
  const mxArray *c1_lhs125 = NULL;
  const mxArray *c1_rhs126 = NULL;
  const mxArray *c1_lhs126 = NULL;
  const mxArray *c1_rhs127 = NULL;
  const mxArray *c1_lhs127 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c1_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c1_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_relop"), "name", "name",
                  66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("function_handle"),
                  "dominantType", "dominantType", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1342458382U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c1_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m"), "context",
                  "context", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c1_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c1_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c1_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c1_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c1_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("max"), "name", "name", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m"), "resolved",
                  "resolved", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311262516U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c1_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c1_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c1_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c1_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c1_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c1_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c1_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c1_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c1_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c1_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c1_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c1_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c1_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1381857500U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c1_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c1_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c1_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c1_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c1_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c1_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c1_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c1_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_realmin"), "name", "name",
                  93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658444U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c1_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c1_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c1_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c1_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c1_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c1_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c1_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xscal"), "name", "name",
                  100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c1_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c1_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xscal"),
                  "name", "name", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c1_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c1_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c1_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c1_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c1_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xscal"),
                  "name", "name", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c1_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c1_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c1_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c1_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c1_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xdotc"), "name", "name",
                  112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c1_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c1_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xdotc"),
                  "name", "name", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c1_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "context", "context", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xdot"),
                  "name", "name", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c1_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c1_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c1_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c1_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xdot"),
                  "name", "name", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c1_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "context", "context", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xdotx"),
                  "name", "name", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c1_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c1_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c1_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c1_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xaxpy"), "name", "name",
                  124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xaxpy.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c1_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xaxpy.m"), "context",
                  "context", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c1_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xaxpy.m"), "context",
                  "context", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xaxpy"),
                  "name", "name", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c1_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p"),
                  "context", "context", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c1_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c1_rhs64);
  sf_mex_destroy(&c1_lhs64);
  sf_mex_destroy(&c1_rhs65);
  sf_mex_destroy(&c1_lhs65);
  sf_mex_destroy(&c1_rhs66);
  sf_mex_destroy(&c1_lhs66);
  sf_mex_destroy(&c1_rhs67);
  sf_mex_destroy(&c1_lhs67);
  sf_mex_destroy(&c1_rhs68);
  sf_mex_destroy(&c1_lhs68);
  sf_mex_destroy(&c1_rhs69);
  sf_mex_destroy(&c1_lhs69);
  sf_mex_destroy(&c1_rhs70);
  sf_mex_destroy(&c1_lhs70);
  sf_mex_destroy(&c1_rhs71);
  sf_mex_destroy(&c1_lhs71);
  sf_mex_destroy(&c1_rhs72);
  sf_mex_destroy(&c1_lhs72);
  sf_mex_destroy(&c1_rhs73);
  sf_mex_destroy(&c1_lhs73);
  sf_mex_destroy(&c1_rhs74);
  sf_mex_destroy(&c1_lhs74);
  sf_mex_destroy(&c1_rhs75);
  sf_mex_destroy(&c1_lhs75);
  sf_mex_destroy(&c1_rhs76);
  sf_mex_destroy(&c1_lhs76);
  sf_mex_destroy(&c1_rhs77);
  sf_mex_destroy(&c1_lhs77);
  sf_mex_destroy(&c1_rhs78);
  sf_mex_destroy(&c1_lhs78);
  sf_mex_destroy(&c1_rhs79);
  sf_mex_destroy(&c1_lhs79);
  sf_mex_destroy(&c1_rhs80);
  sf_mex_destroy(&c1_lhs80);
  sf_mex_destroy(&c1_rhs81);
  sf_mex_destroy(&c1_lhs81);
  sf_mex_destroy(&c1_rhs82);
  sf_mex_destroy(&c1_lhs82);
  sf_mex_destroy(&c1_rhs83);
  sf_mex_destroy(&c1_lhs83);
  sf_mex_destroy(&c1_rhs84);
  sf_mex_destroy(&c1_lhs84);
  sf_mex_destroy(&c1_rhs85);
  sf_mex_destroy(&c1_lhs85);
  sf_mex_destroy(&c1_rhs86);
  sf_mex_destroy(&c1_lhs86);
  sf_mex_destroy(&c1_rhs87);
  sf_mex_destroy(&c1_lhs87);
  sf_mex_destroy(&c1_rhs88);
  sf_mex_destroy(&c1_lhs88);
  sf_mex_destroy(&c1_rhs89);
  sf_mex_destroy(&c1_lhs89);
  sf_mex_destroy(&c1_rhs90);
  sf_mex_destroy(&c1_lhs90);
  sf_mex_destroy(&c1_rhs91);
  sf_mex_destroy(&c1_lhs91);
  sf_mex_destroy(&c1_rhs92);
  sf_mex_destroy(&c1_lhs92);
  sf_mex_destroy(&c1_rhs93);
  sf_mex_destroy(&c1_lhs93);
  sf_mex_destroy(&c1_rhs94);
  sf_mex_destroy(&c1_lhs94);
  sf_mex_destroy(&c1_rhs95);
  sf_mex_destroy(&c1_lhs95);
  sf_mex_destroy(&c1_rhs96);
  sf_mex_destroy(&c1_lhs96);
  sf_mex_destroy(&c1_rhs97);
  sf_mex_destroy(&c1_lhs97);
  sf_mex_destroy(&c1_rhs98);
  sf_mex_destroy(&c1_lhs98);
  sf_mex_destroy(&c1_rhs99);
  sf_mex_destroy(&c1_lhs99);
  sf_mex_destroy(&c1_rhs100);
  sf_mex_destroy(&c1_lhs100);
  sf_mex_destroy(&c1_rhs101);
  sf_mex_destroy(&c1_lhs101);
  sf_mex_destroy(&c1_rhs102);
  sf_mex_destroy(&c1_lhs102);
  sf_mex_destroy(&c1_rhs103);
  sf_mex_destroy(&c1_lhs103);
  sf_mex_destroy(&c1_rhs104);
  sf_mex_destroy(&c1_lhs104);
  sf_mex_destroy(&c1_rhs105);
  sf_mex_destroy(&c1_lhs105);
  sf_mex_destroy(&c1_rhs106);
  sf_mex_destroy(&c1_lhs106);
  sf_mex_destroy(&c1_rhs107);
  sf_mex_destroy(&c1_lhs107);
  sf_mex_destroy(&c1_rhs108);
  sf_mex_destroy(&c1_lhs108);
  sf_mex_destroy(&c1_rhs109);
  sf_mex_destroy(&c1_lhs109);
  sf_mex_destroy(&c1_rhs110);
  sf_mex_destroy(&c1_lhs110);
  sf_mex_destroy(&c1_rhs111);
  sf_mex_destroy(&c1_lhs111);
  sf_mex_destroy(&c1_rhs112);
  sf_mex_destroy(&c1_lhs112);
  sf_mex_destroy(&c1_rhs113);
  sf_mex_destroy(&c1_lhs113);
  sf_mex_destroy(&c1_rhs114);
  sf_mex_destroy(&c1_lhs114);
  sf_mex_destroy(&c1_rhs115);
  sf_mex_destroy(&c1_lhs115);
  sf_mex_destroy(&c1_rhs116);
  sf_mex_destroy(&c1_lhs116);
  sf_mex_destroy(&c1_rhs117);
  sf_mex_destroy(&c1_lhs117);
  sf_mex_destroy(&c1_rhs118);
  sf_mex_destroy(&c1_lhs118);
  sf_mex_destroy(&c1_rhs119);
  sf_mex_destroy(&c1_lhs119);
  sf_mex_destroy(&c1_rhs120);
  sf_mex_destroy(&c1_lhs120);
  sf_mex_destroy(&c1_rhs121);
  sf_mex_destroy(&c1_lhs121);
  sf_mex_destroy(&c1_rhs122);
  sf_mex_destroy(&c1_lhs122);
  sf_mex_destroy(&c1_rhs123);
  sf_mex_destroy(&c1_lhs123);
  sf_mex_destroy(&c1_rhs124);
  sf_mex_destroy(&c1_lhs124);
  sf_mex_destroy(&c1_rhs125);
  sf_mex_destroy(&c1_lhs125);
  sf_mex_destroy(&c1_rhs126);
  sf_mex_destroy(&c1_lhs126);
  sf_mex_destroy(&c1_rhs127);
  sf_mex_destroy(&c1_lhs127);
}

static void c1_c_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs128 = NULL;
  const mxArray *c1_lhs128 = NULL;
  const mxArray *c1_rhs129 = NULL;
  const mxArray *c1_lhs129 = NULL;
  const mxArray *c1_rhs130 = NULL;
  const mxArray *c1_lhs130 = NULL;
  const mxArray *c1_rhs131 = NULL;
  const mxArray *c1_lhs131 = NULL;
  const mxArray *c1_rhs132 = NULL;
  const mxArray *c1_lhs132 = NULL;
  const mxArray *c1_rhs133 = NULL;
  const mxArray *c1_lhs133 = NULL;
  const mxArray *c1_rhs134 = NULL;
  const mxArray *c1_lhs134 = NULL;
  const mxArray *c1_rhs135 = NULL;
  const mxArray *c1_lhs135 = NULL;
  const mxArray *c1_rhs136 = NULL;
  const mxArray *c1_lhs136 = NULL;
  const mxArray *c1_rhs137 = NULL;
  const mxArray *c1_lhs137 = NULL;
  const mxArray *c1_rhs138 = NULL;
  const mxArray *c1_lhs138 = NULL;
  const mxArray *c1_rhs139 = NULL;
  const mxArray *c1_lhs139 = NULL;
  const mxArray *c1_rhs140 = NULL;
  const mxArray *c1_lhs140 = NULL;
  const mxArray *c1_rhs141 = NULL;
  const mxArray *c1_lhs141 = NULL;
  const mxArray *c1_rhs142 = NULL;
  const mxArray *c1_lhs142 = NULL;
  const mxArray *c1_rhs143 = NULL;
  const mxArray *c1_lhs143 = NULL;
  const mxArray *c1_rhs144 = NULL;
  const mxArray *c1_lhs144 = NULL;
  const mxArray *c1_rhs145 = NULL;
  const mxArray *c1_lhs145 = NULL;
  const mxArray *c1_rhs146 = NULL;
  const mxArray *c1_lhs146 = NULL;
  const mxArray *c1_rhs147 = NULL;
  const mxArray *c1_lhs147 = NULL;
  const mxArray *c1_rhs148 = NULL;
  const mxArray *c1_lhs148 = NULL;
  const mxArray *c1_rhs149 = NULL;
  const mxArray *c1_lhs149 = NULL;
  const mxArray *c1_rhs150 = NULL;
  const mxArray *c1_lhs150 = NULL;
  const mxArray *c1_rhs151 = NULL;
  const mxArray *c1_lhs151 = NULL;
  const mxArray *c1_rhs152 = NULL;
  const mxArray *c1_lhs152 = NULL;
  const mxArray *c1_rhs153 = NULL;
  const mxArray *c1_lhs153 = NULL;
  const mxArray *c1_rhs154 = NULL;
  const mxArray *c1_lhs154 = NULL;
  const mxArray *c1_rhs155 = NULL;
  const mxArray *c1_lhs155 = NULL;
  const mxArray *c1_rhs156 = NULL;
  const mxArray *c1_lhs156 = NULL;
  const mxArray *c1_rhs157 = NULL;
  const mxArray *c1_lhs157 = NULL;
  const mxArray *c1_rhs158 = NULL;
  const mxArray *c1_lhs158 = NULL;
  const mxArray *c1_rhs159 = NULL;
  const mxArray *c1_lhs159 = NULL;
  const mxArray *c1_rhs160 = NULL;
  const mxArray *c1_lhs160 = NULL;
  const mxArray *c1_rhs161 = NULL;
  const mxArray *c1_lhs161 = NULL;
  const mxArray *c1_rhs162 = NULL;
  const mxArray *c1_lhs162 = NULL;
  const mxArray *c1_rhs163 = NULL;
  const mxArray *c1_lhs163 = NULL;
  const mxArray *c1_rhs164 = NULL;
  const mxArray *c1_lhs164 = NULL;
  const mxArray *c1_rhs165 = NULL;
  const mxArray *c1_lhs165 = NULL;
  const mxArray *c1_rhs166 = NULL;
  const mxArray *c1_lhs166 = NULL;
  const mxArray *c1_rhs167 = NULL;
  const mxArray *c1_lhs167 = NULL;
  const mxArray *c1_rhs168 = NULL;
  const mxArray *c1_lhs168 = NULL;
  const mxArray *c1_rhs169 = NULL;
  const mxArray *c1_lhs169 = NULL;
  const mxArray *c1_rhs170 = NULL;
  const mxArray *c1_lhs170 = NULL;
  const mxArray *c1_rhs171 = NULL;
  const mxArray *c1_lhs171 = NULL;
  const mxArray *c1_rhs172 = NULL;
  const mxArray *c1_lhs172 = NULL;
  const mxArray *c1_rhs173 = NULL;
  const mxArray *c1_lhs173 = NULL;
  const mxArray *c1_rhs174 = NULL;
  const mxArray *c1_lhs174 = NULL;
  const mxArray *c1_rhs175 = NULL;
  const mxArray *c1_lhs175 = NULL;
  const mxArray *c1_rhs176 = NULL;
  const mxArray *c1_lhs176 = NULL;
  const mxArray *c1_rhs177 = NULL;
  const mxArray *c1_lhs177 = NULL;
  const mxArray *c1_rhs178 = NULL;
  const mxArray *c1_lhs178 = NULL;
  const mxArray *c1_rhs179 = NULL;
  const mxArray *c1_lhs179 = NULL;
  const mxArray *c1_rhs180 = NULL;
  const mxArray *c1_lhs180 = NULL;
  const mxArray *c1_rhs181 = NULL;
  const mxArray *c1_lhs181 = NULL;
  const mxArray *c1_rhs182 = NULL;
  const mxArray *c1_lhs182 = NULL;
  const mxArray *c1_rhs183 = NULL;
  const mxArray *c1_lhs183 = NULL;
  const mxArray *c1_rhs184 = NULL;
  const mxArray *c1_lhs184 = NULL;
  const mxArray *c1_rhs185 = NULL;
  const mxArray *c1_lhs185 = NULL;
  const mxArray *c1_rhs186 = NULL;
  const mxArray *c1_lhs186 = NULL;
  const mxArray *c1_rhs187 = NULL;
  const mxArray *c1_lhs187 = NULL;
  const mxArray *c1_rhs188 = NULL;
  const mxArray *c1_lhs188 = NULL;
  const mxArray *c1_rhs189 = NULL;
  const mxArray *c1_lhs189 = NULL;
  const mxArray *c1_rhs190 = NULL;
  const mxArray *c1_lhs190 = NULL;
  const mxArray *c1_rhs191 = NULL;
  const mxArray *c1_lhs191 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p!below_threshold"),
                  "context", "context", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c1_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p!below_threshold"),
                  "context", "context", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303153406U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c1_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p"),
                  "context", "context", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c1_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xaxpy.p"),
                  "context", "context", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xaxpy"),
                  "name", "name", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c1_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "context", "context", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.isaUint"),
                  "name", "name", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/isaUint.p"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c1_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "context", "context", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c1_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "context", "context", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c1_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "context", "context", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c1_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xaxpy.p"),
                  "context", "context", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c1_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c1_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c1_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307658442U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c1_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c1_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825982U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c1_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_eps"), "name", "name", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c1_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c1_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c1_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c1_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_const_nonsingleton_dim"),
                  "name", "name", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c1_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m"),
                  "context", "context", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c1_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c1_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum"),
                  "context", "context", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c1_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c1_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717458U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c1_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c1_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c1_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub"),
                  "context", "context", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_relop"), "name", "name",
                  154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("function_handle"),
                  "dominantType", "dominantType", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m"), "resolved",
                  "resolved", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1342458382U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c1_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c1_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_error"), "name", "name",
                  156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837558U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c1_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825938U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c1_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xrotg"), "name", "name",
                  158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrotg.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c1_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrotg.m"), "context",
                  "context", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c1_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrotg.m"), "context",
                  "context", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xrotg"),
                  "name", "name", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrotg.p"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c1_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrotg.p"),
                  "context", "context", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c1_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrotg.p"),
                  "context", "context", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xrotg"),
                  "name", "name", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrotg.p"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c1_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrotg.p"),
                  "context", "context", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c1_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrotg.p"),
                  "context", "context", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name",
                  164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1388463696U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370017086U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c1_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrotg.p"),
                  "context", "context", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sqrt"), "name", "name", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1343837586U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c1_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrotg.p!eml_ceval_xrotg"),
                  "context", "context", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c1_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xrot"), "name", "name",
                  167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrot.m"), "resolved",
                  "resolved", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c1_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrot.m"), "context",
                  "context", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c1_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xrot.m"), "context",
                  "context", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xrot"),
                  "name", "name", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrot.p"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c1_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrot.p"),
                  "context", "context", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c1_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrot.p!below_threshold"),
                  "context", "context", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c1_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrot.p"),
                  "context", "context", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c1_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xrot.p"),
                  "context", "context", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xrot"),
                  "name", "name", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrot.p"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c1_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrot.p"),
                  "context", "context", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c1_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xrot.p"),
                  "context", "context", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c1_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zsvdc.m"),
                  "context", "context", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xswap"), "name", "name",
                  176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c1_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c1_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xswap"),
                  "name", "name", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c1_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c1_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p!below_threshold"),
                  "context", "context", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c1_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xswap"),
                  "name", "name", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "resolved", "resolved", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c1_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c1_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c1_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286825912U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c1_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c1_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c1_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c1_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c1_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c1_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c1_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xscal"), "name", "name",
                  191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987892U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c1_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c1_rhs128);
  sf_mex_destroy(&c1_lhs128);
  sf_mex_destroy(&c1_rhs129);
  sf_mex_destroy(&c1_lhs129);
  sf_mex_destroy(&c1_rhs130);
  sf_mex_destroy(&c1_lhs130);
  sf_mex_destroy(&c1_rhs131);
  sf_mex_destroy(&c1_lhs131);
  sf_mex_destroy(&c1_rhs132);
  sf_mex_destroy(&c1_lhs132);
  sf_mex_destroy(&c1_rhs133);
  sf_mex_destroy(&c1_lhs133);
  sf_mex_destroy(&c1_rhs134);
  sf_mex_destroy(&c1_lhs134);
  sf_mex_destroy(&c1_rhs135);
  sf_mex_destroy(&c1_lhs135);
  sf_mex_destroy(&c1_rhs136);
  sf_mex_destroy(&c1_lhs136);
  sf_mex_destroy(&c1_rhs137);
  sf_mex_destroy(&c1_lhs137);
  sf_mex_destroy(&c1_rhs138);
  sf_mex_destroy(&c1_lhs138);
  sf_mex_destroy(&c1_rhs139);
  sf_mex_destroy(&c1_lhs139);
  sf_mex_destroy(&c1_rhs140);
  sf_mex_destroy(&c1_lhs140);
  sf_mex_destroy(&c1_rhs141);
  sf_mex_destroy(&c1_lhs141);
  sf_mex_destroy(&c1_rhs142);
  sf_mex_destroy(&c1_lhs142);
  sf_mex_destroy(&c1_rhs143);
  sf_mex_destroy(&c1_lhs143);
  sf_mex_destroy(&c1_rhs144);
  sf_mex_destroy(&c1_lhs144);
  sf_mex_destroy(&c1_rhs145);
  sf_mex_destroy(&c1_lhs145);
  sf_mex_destroy(&c1_rhs146);
  sf_mex_destroy(&c1_lhs146);
  sf_mex_destroy(&c1_rhs147);
  sf_mex_destroy(&c1_lhs147);
  sf_mex_destroy(&c1_rhs148);
  sf_mex_destroy(&c1_lhs148);
  sf_mex_destroy(&c1_rhs149);
  sf_mex_destroy(&c1_lhs149);
  sf_mex_destroy(&c1_rhs150);
  sf_mex_destroy(&c1_lhs150);
  sf_mex_destroy(&c1_rhs151);
  sf_mex_destroy(&c1_lhs151);
  sf_mex_destroy(&c1_rhs152);
  sf_mex_destroy(&c1_lhs152);
  sf_mex_destroy(&c1_rhs153);
  sf_mex_destroy(&c1_lhs153);
  sf_mex_destroy(&c1_rhs154);
  sf_mex_destroy(&c1_lhs154);
  sf_mex_destroy(&c1_rhs155);
  sf_mex_destroy(&c1_lhs155);
  sf_mex_destroy(&c1_rhs156);
  sf_mex_destroy(&c1_lhs156);
  sf_mex_destroy(&c1_rhs157);
  sf_mex_destroy(&c1_lhs157);
  sf_mex_destroy(&c1_rhs158);
  sf_mex_destroy(&c1_lhs158);
  sf_mex_destroy(&c1_rhs159);
  sf_mex_destroy(&c1_lhs159);
  sf_mex_destroy(&c1_rhs160);
  sf_mex_destroy(&c1_lhs160);
  sf_mex_destroy(&c1_rhs161);
  sf_mex_destroy(&c1_lhs161);
  sf_mex_destroy(&c1_rhs162);
  sf_mex_destroy(&c1_lhs162);
  sf_mex_destroy(&c1_rhs163);
  sf_mex_destroy(&c1_lhs163);
  sf_mex_destroy(&c1_rhs164);
  sf_mex_destroy(&c1_lhs164);
  sf_mex_destroy(&c1_rhs165);
  sf_mex_destroy(&c1_lhs165);
  sf_mex_destroy(&c1_rhs166);
  sf_mex_destroy(&c1_lhs166);
  sf_mex_destroy(&c1_rhs167);
  sf_mex_destroy(&c1_lhs167);
  sf_mex_destroy(&c1_rhs168);
  sf_mex_destroy(&c1_lhs168);
  sf_mex_destroy(&c1_rhs169);
  sf_mex_destroy(&c1_lhs169);
  sf_mex_destroy(&c1_rhs170);
  sf_mex_destroy(&c1_lhs170);
  sf_mex_destroy(&c1_rhs171);
  sf_mex_destroy(&c1_lhs171);
  sf_mex_destroy(&c1_rhs172);
  sf_mex_destroy(&c1_lhs172);
  sf_mex_destroy(&c1_rhs173);
  sf_mex_destroy(&c1_lhs173);
  sf_mex_destroy(&c1_rhs174);
  sf_mex_destroy(&c1_lhs174);
  sf_mex_destroy(&c1_rhs175);
  sf_mex_destroy(&c1_lhs175);
  sf_mex_destroy(&c1_rhs176);
  sf_mex_destroy(&c1_lhs176);
  sf_mex_destroy(&c1_rhs177);
  sf_mex_destroy(&c1_lhs177);
  sf_mex_destroy(&c1_rhs178);
  sf_mex_destroy(&c1_lhs178);
  sf_mex_destroy(&c1_rhs179);
  sf_mex_destroy(&c1_lhs179);
  sf_mex_destroy(&c1_rhs180);
  sf_mex_destroy(&c1_lhs180);
  sf_mex_destroy(&c1_rhs181);
  sf_mex_destroy(&c1_lhs181);
  sf_mex_destroy(&c1_rhs182);
  sf_mex_destroy(&c1_lhs182);
  sf_mex_destroy(&c1_rhs183);
  sf_mex_destroy(&c1_lhs183);
  sf_mex_destroy(&c1_rhs184);
  sf_mex_destroy(&c1_lhs184);
  sf_mex_destroy(&c1_rhs185);
  sf_mex_destroy(&c1_lhs185);
  sf_mex_destroy(&c1_rhs186);
  sf_mex_destroy(&c1_lhs186);
  sf_mex_destroy(&c1_rhs187);
  sf_mex_destroy(&c1_lhs187);
  sf_mex_destroy(&c1_rhs188);
  sf_mex_destroy(&c1_lhs188);
  sf_mex_destroy(&c1_rhs189);
  sf_mex_destroy(&c1_lhs189);
  sf_mex_destroy(&c1_rhs190);
  sf_mex_destroy(&c1_lhs190);
  sf_mex_destroy(&c1_rhs191);
  sf_mex_destroy(&c1_lhs191);
}

static void c1_d_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs192 = NULL;
  const mxArray *c1_lhs192 = NULL;
  const mxArray *c1_rhs193 = NULL;
  const mxArray *c1_lhs193 = NULL;
  const mxArray *c1_rhs194 = NULL;
  const mxArray *c1_lhs194 = NULL;
  const mxArray *c1_rhs195 = NULL;
  const mxArray *c1_lhs195 = NULL;
  const mxArray *c1_rhs196 = NULL;
  const mxArray *c1_lhs196 = NULL;
  const mxArray *c1_rhs197 = NULL;
  const mxArray *c1_lhs197 = NULL;
  const mxArray *c1_rhs198 = NULL;
  const mxArray *c1_lhs198 = NULL;
  const mxArray *c1_rhs199 = NULL;
  const mxArray *c1_lhs199 = NULL;
  const mxArray *c1_rhs200 = NULL;
  const mxArray *c1_lhs200 = NULL;
  const mxArray *c1_rhs201 = NULL;
  const mxArray *c1_lhs201 = NULL;
  const mxArray *c1_rhs202 = NULL;
  const mxArray *c1_lhs202 = NULL;
  const mxArray *c1_rhs203 = NULL;
  const mxArray *c1_lhs203 = NULL;
  const mxArray *c1_rhs204 = NULL;
  const mxArray *c1_lhs204 = NULL;
  const mxArray *c1_rhs205 = NULL;
  const mxArray *c1_lhs205 = NULL;
  const mxArray *c1_rhs206 = NULL;
  const mxArray *c1_lhs206 = NULL;
  const mxArray *c1_rhs207 = NULL;
  const mxArray *c1_lhs207 = NULL;
  const mxArray *c1_rhs208 = NULL;
  const mxArray *c1_lhs208 = NULL;
  const mxArray *c1_rhs209 = NULL;
  const mxArray *c1_lhs209 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372589616U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c1_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/pinv.m!eml_pinv"),
                  "context", "context", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c1_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c1_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c1_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c1_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c1_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("min"), "name", "name", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311262518U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c1_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c1_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311522U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c1_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c1_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389311520U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c1_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c1_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c1_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c1_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "context", "context", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372590360U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c1_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c1_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987888U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c1_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375987890U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c1_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs209), "lhs", "lhs",
                  209);
  sf_mex_destroy(&c1_rhs192);
  sf_mex_destroy(&c1_lhs192);
  sf_mex_destroy(&c1_rhs193);
  sf_mex_destroy(&c1_lhs193);
  sf_mex_destroy(&c1_rhs194);
  sf_mex_destroy(&c1_lhs194);
  sf_mex_destroy(&c1_rhs195);
  sf_mex_destroy(&c1_lhs195);
  sf_mex_destroy(&c1_rhs196);
  sf_mex_destroy(&c1_lhs196);
  sf_mex_destroy(&c1_rhs197);
  sf_mex_destroy(&c1_lhs197);
  sf_mex_destroy(&c1_rhs198);
  sf_mex_destroy(&c1_lhs198);
  sf_mex_destroy(&c1_rhs199);
  sf_mex_destroy(&c1_lhs199);
  sf_mex_destroy(&c1_rhs200);
  sf_mex_destroy(&c1_lhs200);
  sf_mex_destroy(&c1_rhs201);
  sf_mex_destroy(&c1_lhs201);
  sf_mex_destroy(&c1_rhs202);
  sf_mex_destroy(&c1_lhs202);
  sf_mex_destroy(&c1_rhs203);
  sf_mex_destroy(&c1_lhs203);
  sf_mex_destroy(&c1_rhs204);
  sf_mex_destroy(&c1_lhs204);
  sf_mex_destroy(&c1_rhs205);
  sf_mex_destroy(&c1_lhs205);
  sf_mex_destroy(&c1_rhs206);
  sf_mex_destroy(&c1_lhs206);
  sf_mex_destroy(&c1_rhs207);
  sf_mex_destroy(&c1_lhs207);
  sf_mex_destroy(&c1_rhs208);
  sf_mex_destroy(&c1_lhs208);
  sf_mex_destroy(&c1_rhs209);
  sf_mex_destroy(&c1_lhs209);
}

static real_T c1_eml_div(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x,
  real_T c1_y)
{
  real_T c1_b_x;
  real_T c1_b_y;
  (void)chartInstance;
  c1_b_x = c1_x;
  c1_b_y = c1_y;
  return c1_b_x / c1_b_y;
}

static void c1_eml_switch_helper(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_pinv(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_A[36],
                    real_T c1_X[36])
{
  int32_T c1_i40;
  int32_T c1_k;
  int32_T c1_b_k;
  real_T c1_x;
  real_T c1_b_x;
  boolean_T c1_b;
  boolean_T c1_b0;
  real_T c1_c_x;
  boolean_T c1_b_b;
  boolean_T c1_b1;
  boolean_T c1_c_b;
  int32_T c1_i41;
  real_T c1_b_A[36];
  real_T c1_V[36];
  real_T c1_s[6];
  real_T c1_U[36];
  int32_T c1_i42;
  real_T c1_S[36];
  int32_T c1_c_k;
  real_T c1_d_k;
  real_T c1_tol;
  int32_T c1_r;
  int32_T c1_e_k;
  int32_T c1_f_k;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_vcol;
  int32_T c1_b_r;
  int32_T c1_d_b;
  int32_T c1_e_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  real_T c1_y;
  real_T c1_b_y;
  real_T c1_z;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_i43;
  real_T c1_b_V[36];
  int32_T c1_i44;
  real_T c1_b_U[36];
  boolean_T exitg1;
  for (c1_i40 = 0; c1_i40 < 36; c1_i40++) {
    c1_X[c1_i40] = 0.0;
  }

  for (c1_k = 1; c1_k < 37; c1_k++) {
    c1_b_k = c1_k;
    c1_x = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_k), 1, 36, 1, 0) - 1];
    c1_b_x = c1_x;
    c1_b = muDoubleScalarIsInf(c1_b_x);
    c1_b0 = !c1_b;
    c1_c_x = c1_x;
    c1_b_b = muDoubleScalarIsNaN(c1_c_x);
    c1_b1 = !c1_b_b;
    c1_c_b = (c1_b0 && c1_b1);
    if (!c1_c_b) {
      c1_eml_error(chartInstance);
    }
  }

  for (c1_i41 = 0; c1_i41 < 36; c1_i41++) {
    c1_b_A[c1_i41] = c1_A[c1_i41];
  }

  c1_eml_xgesvd(chartInstance, c1_b_A, c1_U, c1_s, c1_V);
  for (c1_i42 = 0; c1_i42 < 36; c1_i42++) {
    c1_S[c1_i42] = 0.0;
  }

  for (c1_c_k = 0; c1_c_k < 6; c1_c_k++) {
    c1_d_k = 1.0 + (real_T)c1_c_k;
    c1_S[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c1_d_k),
           1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", c1_d_k), 1, 6, 2, 0) - 1)) - 1] =
      c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c1_d_k), 1, 6, 1, 0) - 1];
  }

  c1_eps(chartInstance);
  c1_tol = 6.0 * c1_S[0] * 2.2204460492503131E-16;
  c1_r = 0;
  c1_e_k = 1;
  exitg1 = false;
  while ((exitg1 == false) && (c1_e_k < 7)) {
    c1_f_k = c1_e_k;
    if (!(c1_S[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_f_k), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_f_k), 1, 6, 2, 0) - 1)) -
          1] > c1_tol)) {
      exitg1 = true;
    } else {
      c1_a = c1_r;
      c1_b_a = c1_a + 1;
      c1_r = c1_b_a;
      c1_e_k++;
    }
  }

  if (c1_r > 0) {
    c1_vcol = 1;
    c1_b_r = c1_r;
    c1_d_b = c1_b_r;
    c1_e_b = c1_d_b;
    if (1 > c1_e_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_e_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_j = 1; c1_j <= c1_b_r; c1_j++) {
      c1_b_j = c1_j;
      c1_y = c1_S[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_b_j), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
                     (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2, 0)
        - 1)) - 1];
      c1_b_y = c1_y;
      c1_z = 1.0 / c1_b_y;
      c1_f_eml_xscal(chartInstance, c1_z, c1_V, c1_vcol);
      c1_c_a = c1_vcol;
      c1_d_a = c1_c_a + 6;
      c1_vcol = c1_d_a;
    }

    for (c1_i43 = 0; c1_i43 < 36; c1_i43++) {
      c1_b_V[c1_i43] = c1_V[c1_i43];
    }

    for (c1_i44 = 0; c1_i44 < 36; c1_i44++) {
      c1_b_U[c1_i44] = c1_U[c1_i44];
    }

    c1_b_eml_xgemm(chartInstance, c1_r, c1_b_V, c1_b_U, c1_X);
  }
}

static void c1_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  int32_T c1_i45;
  static char_T c1_cv0[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 's', 'v', 'd', '_', 'm', 'a', 't', 'r', 'i', 'x', 'W', 'i',
    't', 'h', 'N', 'a', 'N', 'I', 'n', 'f' };

  char_T c1_u[33];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i45 = 0; c1_i45 < 33; c1_i45++) {
    c1_u[c1_i45] = c1_cv0[c1_i45];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static void c1_eml_xgesvd(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_A[36], real_T c1_U[36], real_T c1_S[6], real_T c1_V[36])
{
  int32_T c1_i46;
  real_T c1_b_A[36];
  int32_T c1_i47;
  real_T c1_s[6];
  int32_T c1_i48;
  real_T c1_e[6];
  int32_T c1_i49;
  real_T c1_work[6];
  int32_T c1_i50;
  int32_T c1_i51;
  real_T c1_Vf[36];
  int32_T c1_q;
  int32_T c1_b_q;
  int32_T c1_a;
  int32_T c1_b_a;
  int32_T c1_qp1;
  int32_T c1_c_a;
  int32_T c1_d_a;
  int32_T c1_qm1;
  int32_T c1_b;
  int32_T c1_b_b;
  int32_T c1_c;
  int32_T c1_e_a;
  int32_T c1_c_b;
  int32_T c1_f_a;
  int32_T c1_d_b;
  int32_T c1_qq;
  int32_T c1_e_b;
  int32_T c1_f_b;
  int32_T c1_nmq;
  int32_T c1_g_a;
  int32_T c1_h_a;
  int32_T c1_nmqp1;
  int32_T c1_i52;
  real_T c1_c_A[36];
  real_T c1_nrm;
  real_T c1_absx;
  real_T c1_d;
  real_T c1_y;
  real_T c1_d1;
  int32_T c1_b_qp1;
  int32_T c1_i_a;
  int32_T c1_j_a;
  boolean_T c1_overflow;
  int32_T c1_jj;
  int32_T c1_b_jj;
  int32_T c1_k_a;
  int32_T c1_l_a;
  int32_T c1_b_c;
  int32_T c1_g_b;
  int32_T c1_h_b;
  int32_T c1_c_c;
  int32_T c1_m_a;
  int32_T c1_i_b;
  int32_T c1_n_a;
  int32_T c1_j_b;
  int32_T c1_qjj;
  int32_T c1_i53;
  real_T c1_d_A[36];
  int32_T c1_i54;
  real_T c1_e_A[36];
  real_T c1_t;
  int32_T c1_c_q;
  int32_T c1_o_a;
  int32_T c1_p_a;
  boolean_T c1_b_overflow;
  int32_T c1_ii;
  int32_T c1_b_ii;
  int32_T c1_k_b;
  int32_T c1_l_b;
  int32_T c1_pmq;
  int32_T c1_i55;
  real_T c1_b_e[6];
  real_T c1_b_absx;
  real_T c1_b_d;
  real_T c1_b_y;
  real_T c1_d2;
  int32_T c1_c_qp1;
  int32_T c1_q_a;
  int32_T c1_r_a;
  boolean_T c1_c_overflow;
  int32_T c1_c_ii;
  int32_T c1_d_qp1;
  int32_T c1_s_a;
  int32_T c1_t_a;
  boolean_T c1_d_overflow;
  int32_T c1_c_jj;
  int32_T c1_u_a;
  int32_T c1_v_a;
  int32_T c1_d_c;
  int32_T c1_m_b;
  int32_T c1_n_b;
  int32_T c1_e_c;
  int32_T c1_w_a;
  int32_T c1_o_b;
  int32_T c1_x_a;
  int32_T c1_p_b;
  int32_T c1_qp1jj;
  int32_T c1_i56;
  real_T c1_f_A[36];
  int32_T c1_e_qp1;
  int32_T c1_y_a;
  int32_T c1_ab_a;
  boolean_T c1_e_overflow;
  int32_T c1_d_jj;
  int32_T c1_bb_a;
  int32_T c1_cb_a;
  int32_T c1_f_c;
  int32_T c1_q_b;
  int32_T c1_r_b;
  int32_T c1_g_c;
  int32_T c1_db_a;
  int32_T c1_s_b;
  int32_T c1_eb_a;
  int32_T c1_t_b;
  int32_T c1_i57;
  real_T c1_b_work[6];
  int32_T c1_f_qp1;
  int32_T c1_fb_a;
  int32_T c1_gb_a;
  boolean_T c1_f_overflow;
  int32_T c1_d_ii;
  int32_T c1_m;
  int32_T c1_e_ii;
  int32_T c1_d_q;
  int32_T c1_hb_a;
  int32_T c1_ib_a;
  int32_T c1_u_b;
  int32_T c1_v_b;
  int32_T c1_jb_a;
  int32_T c1_kb_a;
  int32_T c1_lb_a;
  int32_T c1_mb_a;
  int32_T c1_h_c;
  int32_T c1_w_b;
  int32_T c1_x_b;
  int32_T c1_i_c;
  int32_T c1_nb_a;
  int32_T c1_y_b;
  int32_T c1_ob_a;
  int32_T c1_ab_b;
  int32_T c1_g_qp1;
  int32_T c1_pb_a;
  int32_T c1_qb_a;
  boolean_T c1_g_overflow;
  int32_T c1_e_jj;
  int32_T c1_rb_a;
  int32_T c1_sb_a;
  int32_T c1_j_c;
  int32_T c1_bb_b;
  int32_T c1_cb_b;
  int32_T c1_k_c;
  int32_T c1_tb_a;
  int32_T c1_db_b;
  int32_T c1_ub_a;
  int32_T c1_eb_b;
  int32_T c1_i58;
  real_T c1_b_U[36];
  int32_T c1_i59;
  real_T c1_c_U[36];
  int32_T c1_e_q;
  int32_T c1_vb_a;
  int32_T c1_wb_a;
  boolean_T c1_h_overflow;
  int32_T c1_f_ii;
  int32_T c1_xb_a;
  int32_T c1_yb_a;
  int32_T c1_i60;
  int32_T c1_fb_b;
  int32_T c1_gb_b;
  boolean_T c1_i_overflow;
  int32_T c1_g_ii;
  int32_T c1_h_ii;
  int32_T c1_f_q;
  int32_T c1_ac_a;
  int32_T c1_bc_a;
  int32_T c1_hb_b;
  int32_T c1_ib_b;
  int32_T c1_cc_a;
  int32_T c1_dc_a;
  int32_T c1_l_c;
  int32_T c1_jb_b;
  int32_T c1_kb_b;
  int32_T c1_m_c;
  int32_T c1_ec_a;
  int32_T c1_lb_b;
  int32_T c1_fc_a;
  int32_T c1_mb_b;
  int32_T c1_qp1q;
  int32_T c1_h_qp1;
  int32_T c1_gc_a;
  int32_T c1_hc_a;
  boolean_T c1_j_overflow;
  int32_T c1_f_jj;
  int32_T c1_ic_a;
  int32_T c1_jc_a;
  int32_T c1_n_c;
  int32_T c1_nb_b;
  int32_T c1_ob_b;
  int32_T c1_o_c;
  int32_T c1_kc_a;
  int32_T c1_pb_b;
  int32_T c1_lc_a;
  int32_T c1_qb_b;
  int32_T c1_i61;
  real_T c1_b_Vf[36];
  int32_T c1_i62;
  real_T c1_c_Vf[36];
  int32_T c1_i_ii;
  int32_T c1_g_q;
  real_T c1_rt;
  real_T c1_r;
  int32_T c1_mc_a;
  int32_T c1_nc_a;
  int32_T c1_p_c;
  int32_T c1_rb_b;
  int32_T c1_sb_b;
  int32_T c1_q_c;
  int32_T c1_tb_b;
  int32_T c1_ub_b;
  int32_T c1_colq;
  int32_T c1_oc_a;
  int32_T c1_pc_a;
  int32_T c1_r_c;
  int32_T c1_qc_a;
  int32_T c1_rc_a;
  int32_T c1_s_c;
  int32_T c1_vb_b;
  int32_T c1_wb_b;
  int32_T c1_t_c;
  int32_T c1_xb_b;
  int32_T c1_yb_b;
  int32_T c1_colqp1;
  real_T c1_iter;
  real_T c1_tiny;
  real_T c1_snorm;
  int32_T c1_j_ii;
  real_T c1_varargin_1;
  real_T c1_varargin_2;
  real_T c1_b_varargin_2;
  real_T c1_varargin_3;
  real_T c1_x;
  real_T c1_c_y;
  real_T c1_b_x;
  real_T c1_d_y;
  real_T c1_xk;
  real_T c1_yk;
  real_T c1_c_x;
  real_T c1_e_y;
  real_T c1_maxval;
  real_T c1_b_varargin_1;
  real_T c1_c_varargin_2;
  real_T c1_d_varargin_2;
  real_T c1_b_varargin_3;
  real_T c1_d_x;
  real_T c1_f_y;
  real_T c1_e_x;
  real_T c1_g_y;
  real_T c1_b_xk;
  real_T c1_b_yk;
  real_T c1_f_x;
  real_T c1_h_y;
  int32_T c1_sc_a;
  int32_T c1_tc_a;
  int32_T c1_uc_a;
  int32_T c1_vc_a;
  int32_T c1_i63;
  int32_T c1_wc_a;
  int32_T c1_xc_a;
  boolean_T c1_k_overflow;
  int32_T c1_k_ii;
  int32_T c1_yc_a;
  int32_T c1_ad_a;
  int32_T c1_u_c;
  real_T c1_test0;
  real_T c1_ztest0;
  int32_T c1_bd_a;
  int32_T c1_cd_a;
  int32_T c1_v_c;
  real_T c1_kase;
  int32_T c1_qs;
  int32_T c1_b_m;
  int32_T c1_h_q;
  int32_T c1_dd_a;
  int32_T c1_ac_b;
  int32_T c1_ed_a;
  int32_T c1_bc_b;
  boolean_T c1_l_overflow;
  int32_T c1_l_ii;
  real_T c1_test;
  int32_T c1_fd_a;
  int32_T c1_gd_a;
  int32_T c1_w_c;
  int32_T c1_hd_a;
  int32_T c1_id_a;
  int32_T c1_x_c;
  real_T c1_ztest;
  int32_T c1_jd_a;
  int32_T c1_kd_a;
  int32_T c1_ld_a;
  int32_T c1_md_a;
  int32_T c1_y_c;
  real_T c1_f;
  int32_T c1_nd_a;
  int32_T c1_od_a;
  int32_T c1_ab_c;
  int32_T c1_pd_a;
  int32_T c1_qd_a;
  int32_T c1_i64;
  int32_T c1_i_q;
  int32_T c1_rd_a;
  int32_T c1_cc_b;
  int32_T c1_sd_a;
  int32_T c1_dc_b;
  boolean_T c1_m_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  real_T c1_t1;
  real_T c1_b_t1;
  real_T c1_b_f;
  real_T c1_sn;
  real_T c1_cs;
  real_T c1_b_cs;
  real_T c1_b_sn;
  int32_T c1_td_a;
  int32_T c1_ud_a;
  int32_T c1_km1;
  int32_T c1_vd_a;
  int32_T c1_wd_a;
  int32_T c1_bb_c;
  int32_T c1_ec_b;
  int32_T c1_fc_b;
  int32_T c1_cb_c;
  int32_T c1_gc_b;
  int32_T c1_hc_b;
  int32_T c1_colk;
  int32_T c1_xd_a;
  int32_T c1_yd_a;
  int32_T c1_db_c;
  int32_T c1_ic_b;
  int32_T c1_jc_b;
  int32_T c1_eb_c;
  int32_T c1_kc_b;
  int32_T c1_lc_b;
  int32_T c1_colm;
  int32_T c1_ae_a;
  int32_T c1_be_a;
  int32_T c1_j_q;
  int32_T c1_c_m;
  int32_T c1_ce_a;
  int32_T c1_mc_b;
  int32_T c1_de_a;
  int32_T c1_nc_b;
  boolean_T c1_n_overflow;
  int32_T c1_c_k;
  real_T c1_c_t1;
  real_T c1_unusedU0;
  real_T c1_c_sn;
  real_T c1_c_cs;
  int32_T c1_ee_a;
  int32_T c1_fe_a;
  int32_T c1_fb_c;
  int32_T c1_oc_b;
  int32_T c1_pc_b;
  int32_T c1_gb_c;
  int32_T c1_qc_b;
  int32_T c1_rc_b;
  int32_T c1_ge_a;
  int32_T c1_he_a;
  int32_T c1_hb_c;
  int32_T c1_sc_b;
  int32_T c1_tc_b;
  int32_T c1_ib_c;
  int32_T c1_uc_b;
  int32_T c1_vc_b;
  int32_T c1_colqm1;
  int32_T c1_ie_a;
  int32_T c1_je_a;
  int32_T c1_mm1;
  real_T c1_d3;
  real_T c1_d4;
  real_T c1_d5;
  real_T c1_d6;
  real_T c1_d7;
  real_T c1_c_varargin_1[5];
  int32_T c1_ixstart;
  real_T c1_mtmp;
  real_T c1_g_x;
  boolean_T c1_wc_b;
  int32_T c1_ix;
  int32_T c1_b_ix;
  real_T c1_h_x;
  boolean_T c1_xc_b;
  int32_T c1_ke_a;
  int32_T c1_le_a;
  int32_T c1_i65;
  int32_T c1_me_a;
  int32_T c1_ne_a;
  boolean_T c1_o_overflow;
  int32_T c1_c_ix;
  real_T c1_oe_a;
  real_T c1_yc_b;
  boolean_T c1_p;
  real_T c1_b_mtmp;
  real_T c1_scale;
  real_T c1_sm;
  real_T c1_smm1;
  real_T c1_emm1;
  real_T c1_sqds;
  real_T c1_eqds;
  real_T c1_ad_b;
  real_T c1_jb_c;
  real_T c1_shift;
  real_T c1_g;
  int32_T c1_k_q;
  int32_T c1_b_mm1;
  int32_T c1_pe_a;
  int32_T c1_bd_b;
  int32_T c1_qe_a;
  int32_T c1_cd_b;
  boolean_T c1_p_overflow;
  int32_T c1_d_k;
  int32_T c1_re_a;
  int32_T c1_se_a;
  int32_T c1_te_a;
  int32_T c1_ue_a;
  int32_T c1_kp1;
  real_T c1_c_f;
  real_T c1_unusedU1;
  real_T c1_d_sn;
  real_T c1_d_cs;
  int32_T c1_ve_a;
  int32_T c1_we_a;
  int32_T c1_kb_c;
  int32_T c1_dd_b;
  int32_T c1_ed_b;
  int32_T c1_lb_c;
  int32_T c1_fd_b;
  int32_T c1_gd_b;
  int32_T c1_hd_b;
  int32_T c1_id_b;
  int32_T c1_mb_c;
  int32_T c1_jd_b;
  int32_T c1_kd_b;
  int32_T c1_colkp1;
  real_T c1_d_f;
  real_T c1_unusedU2;
  real_T c1_e_sn;
  real_T c1_e_cs;
  int32_T c1_xe_a;
  int32_T c1_ye_a;
  int32_T c1_nb_c;
  int32_T c1_ld_b;
  int32_T c1_md_b;
  int32_T c1_ob_c;
  int32_T c1_nd_b;
  int32_T c1_od_b;
  int32_T c1_pd_b;
  int32_T c1_qd_b;
  int32_T c1_pb_c;
  int32_T c1_rd_b;
  int32_T c1_sd_b;
  int32_T c1_af_a;
  int32_T c1_bf_a;
  int32_T c1_qb_c;
  int32_T c1_e_k;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_i;
  int32_T c1_b_i;
  int32_T c1_rb_c;
  int32_T c1_cf_a;
  int32_T c1_sb_c;
  int32_T c1_td_b;
  int32_T c1_ud_b;
  int32_T c1_df_a;
  int32_T c1_ef_a;
  int32_T c1_tb_c;
  int32_T c1_ff_a;
  int32_T c1_ub_c;
  int32_T c1_vd_b;
  int32_T c1_wd_b;
  int32_T c1_vb_c;
  int32_T c1_xd_b;
  int32_T c1_yd_b;
  int32_T c1_wb_c;
  int32_T c1_gf_a;
  int32_T c1_xb_c;
  int32_T c1_ae_b;
  int32_T c1_be_b;
  int32_T c1_yb_c;
  int32_T c1_ce_b;
  int32_T c1_de_b;
  int32_T c1_hf_a;
  int32_T c1_if_a;
  int32_T c1_ee_b;
  int32_T c1_fe_b;
  int32_T c1_jf_a;
  int32_T c1_kf_a;
  int32_T c1_lf_a;
  int32_T c1_ge_b;
  int32_T c1_he_b;
  int32_T c1_ie_b;
  int32_T c1_je_b;
  int32_T c1_mf_a;
  int32_T c1_ke_b;
  int32_T c1_le_b;
  int32_T c1_me_b;
  int32_T c1_ne_b;
  int32_T c1_nf_a;
  real_T c1_d8;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T exitg5;
  boolean_T guard11 = false;
  for (c1_i46 = 0; c1_i46 < 36; c1_i46++) {
    c1_b_A[c1_i46] = c1_A[c1_i46];
  }

  c1_eml_scalar_eg(chartInstance);
  for (c1_i47 = 0; c1_i47 < 6; c1_i47++) {
    c1_s[c1_i47] = 0.0;
  }

  for (c1_i48 = 0; c1_i48 < 6; c1_i48++) {
    c1_e[c1_i48] = 0.0;
  }

  for (c1_i49 = 0; c1_i49 < 6; c1_i49++) {
    c1_work[c1_i49] = 0.0;
  }

  for (c1_i50 = 0; c1_i50 < 36; c1_i50++) {
    c1_U[c1_i50] = 0.0;
  }

  for (c1_i51 = 0; c1_i51 < 36; c1_i51++) {
    c1_Vf[c1_i51] = 0.0;
  }

  for (c1_q = 1; c1_q < 6; c1_q++) {
    c1_b_q = c1_q;
    c1_a = c1_b_q;
    c1_b_a = c1_a + 1;
    c1_qp1 = c1_b_a;
    c1_c_a = c1_b_q;
    c1_d_a = c1_c_a;
    c1_qm1 = c1_d_a;
    c1_b = c1_qm1 - 1;
    c1_b_b = c1_b;
    c1_c = 6 * c1_b_b;
    c1_e_a = c1_b_q;
    c1_c_b = c1_c;
    c1_f_a = c1_e_a;
    c1_d_b = c1_c_b;
    c1_qq = c1_f_a + c1_d_b;
    c1_e_b = c1_b_q;
    c1_f_b = c1_e_b;
    c1_nmq = 6 - c1_f_b;
    c1_g_a = c1_nmq;
    c1_h_a = c1_g_a + 1;
    c1_nmqp1 = c1_h_a;
    if (c1_b_q <= 5) {
      for (c1_i52 = 0; c1_i52 < 36; c1_i52++) {
        c1_c_A[c1_i52] = c1_b_A[c1_i52];
      }

      c1_nrm = c1_eml_xnrm2(chartInstance, c1_nmqp1, c1_c_A, c1_qq);
      if (c1_nrm > 0.0) {
        c1_absx = c1_nrm;
        c1_d = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_qq), 1, 36, 1, 0) - 1];
        if (c1_d < 0.0) {
          c1_y = -c1_absx;
        } else {
          c1_y = c1_absx;
        }

        c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = c1_y;
        c1_d1 = c1_eml_div(chartInstance, 1.0, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
        c1_d_eml_xscal(chartInstance, c1_nmqp1, c1_d1, c1_b_A, c1_qq);
        c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_qq), 1, 36, 1, 0) - 1] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK
          ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qq), 1, 36, 1, 0) - 1]
          + 1.0;
        c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = -c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1];
      } else {
        c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = 0.0;
      }
    }

    c1_b_qp1 = c1_qp1;
    c1_i_a = c1_b_qp1;
    c1_j_a = c1_i_a;
    if (c1_j_a > 6) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = false;
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_jj = c1_b_qp1; c1_jj < 7; c1_jj++) {
      c1_b_jj = c1_jj;
      c1_k_a = c1_b_jj;
      c1_l_a = c1_k_a;
      c1_b_c = c1_l_a;
      c1_g_b = c1_b_c - 1;
      c1_h_b = c1_g_b;
      c1_c_c = 6 * c1_h_b;
      c1_m_a = c1_b_q;
      c1_i_b = c1_c_c;
      c1_n_a = c1_m_a;
      c1_j_b = c1_i_b;
      c1_qjj = c1_n_a + c1_j_b;
      if (c1_b_q <= 5) {
        if (c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
          for (c1_i53 = 0; c1_i53 < 36; c1_i53++) {
            c1_d_A[c1_i53] = c1_b_A[c1_i53];
          }

          for (c1_i54 = 0; c1_i54 < 36; c1_i54++) {
            c1_e_A[c1_i54] = c1_b_A[c1_i54];
          }

          c1_t = c1_eml_xdotc(chartInstance, c1_nmqp1, c1_d_A, c1_qq, c1_e_A,
                              c1_qjj);
          c1_t = -c1_eml_div(chartInstance, c1_t, c1_b_A
                             [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) + 6 *
                               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0) - 1)) - 1]);
          c1_d_eml_xaxpy(chartInstance, c1_nmqp1, c1_t, c1_qq, c1_b_A, c1_qjj);
        }
      }

      c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_jj), 1, 6, 1, 0) - 1] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK(
        "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qjj), 1, 36, 1, 0) - 1];
    }

    if (c1_b_q <= 5) {
      c1_c_q = c1_b_q;
      c1_o_a = c1_c_q;
      c1_p_a = c1_o_a;
      if (c1_p_a > 6) {
        c1_b_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_b_overflow = false;
      }

      if (c1_b_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
      }

      for (c1_ii = c1_c_q; c1_ii < 7; c1_ii++) {
        c1_b_ii = c1_ii;
        c1_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_ii), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0)
               - 1)) - 1] = c1_b_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1, 6, 1, 0) + 6 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 2, 0) - 1)) - 1];
      }
    }

    if (c1_b_q <= 4) {
      c1_k_b = c1_b_q;
      c1_l_b = c1_k_b;
      c1_pmq = 6 - c1_l_b;
      for (c1_i55 = 0; c1_i55 < 6; c1_i55++) {
        c1_b_e[c1_i55] = c1_e[c1_i55];
      }

      c1_nrm = c1_b_eml_xnrm2(chartInstance, c1_pmq, c1_b_e, c1_qp1);
      if (c1_nrm == 0.0) {
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = 0.0;
      } else {
        c1_b_absx = c1_nrm;
        c1_b_d = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_qp1), 1, 6, 1, 0) - 1];
        if (c1_b_d < 0.0) {
          c1_b_y = -c1_b_absx;
        } else {
          c1_b_y = c1_b_absx;
        }

        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = c1_b_y;
        c1_d2 = c1_eml_div(chartInstance, 1.0, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
        c1_e_eml_xscal(chartInstance, c1_pmq, c1_d2, c1_e, c1_qp1);
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_qp1), 1, 6, 1, 0) - 1] = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qp1), 1, 6, 1, 0) - 1]
          + 1.0;
      }

      c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_q), 1, 6, 1, 0) - 1] = -c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1];
      if (c1_qp1 <= 6) {
        if (c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
          c1_c_qp1 = c1_qp1;
          c1_q_a = c1_c_qp1;
          c1_r_a = c1_q_a;
          if (c1_r_a > 6) {
            c1_c_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_c_overflow = false;
          }

          if (c1_c_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_c_overflow);
          }

          for (c1_c_ii = c1_c_qp1; c1_c_ii < 7; c1_c_ii++) {
            c1_b_ii = c1_c_ii;
            c1_work[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
              "", (real_T)c1_b_ii), 1, 6, 1, 0) - 1] = 0.0;
          }

          c1_d_qp1 = c1_qp1;
          c1_s_a = c1_d_qp1;
          c1_t_a = c1_s_a;
          if (c1_t_a > 6) {
            c1_d_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_d_overflow = false;
          }

          if (c1_d_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_d_overflow);
          }

          for (c1_c_jj = c1_d_qp1; c1_c_jj < 7; c1_c_jj++) {
            c1_b_jj = c1_c_jj;
            c1_u_a = c1_b_jj;
            c1_v_a = c1_u_a;
            c1_d_c = c1_v_a;
            c1_m_b = c1_d_c - 1;
            c1_n_b = c1_m_b;
            c1_e_c = 6 * c1_n_b;
            c1_w_a = c1_qp1;
            c1_o_b = c1_e_c;
            c1_x_a = c1_w_a;
            c1_p_b = c1_o_b;
            c1_qp1jj = c1_x_a + c1_p_b;
            for (c1_i56 = 0; c1_i56 < 36; c1_i56++) {
              c1_f_A[c1_i56] = c1_b_A[c1_i56];
            }

            c1_e_eml_xaxpy(chartInstance, c1_nmq,
                           c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_jj), 1, 6, 1, 0) - 1], c1_f_A,
                           c1_qp1jj, c1_work, c1_qp1);
          }

          c1_e_qp1 = c1_qp1;
          c1_y_a = c1_e_qp1;
          c1_ab_a = c1_y_a;
          if (c1_ab_a > 6) {
            c1_e_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_e_overflow = false;
          }

          if (c1_e_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_e_overflow);
          }

          for (c1_d_jj = c1_e_qp1; c1_d_jj < 7; c1_d_jj++) {
            c1_b_jj = c1_d_jj;
            c1_t = c1_eml_div(chartInstance, -c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK(
              "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_jj), 1, 6, 1, 0)
                              - 1], c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qp1), 1, 6, 1, 0) - 1]);
            c1_bb_a = c1_b_jj;
            c1_cb_a = c1_bb_a;
            c1_f_c = c1_cb_a;
            c1_q_b = c1_f_c - 1;
            c1_r_b = c1_q_b;
            c1_g_c = 6 * c1_r_b;
            c1_db_a = c1_qp1;
            c1_s_b = c1_g_c;
            c1_eb_a = c1_db_a;
            c1_t_b = c1_s_b;
            c1_qp1jj = c1_eb_a + c1_t_b;
            for (c1_i57 = 0; c1_i57 < 6; c1_i57++) {
              c1_b_work[c1_i57] = c1_work[c1_i57];
            }

            c1_f_eml_xaxpy(chartInstance, c1_nmq, c1_t, c1_b_work, c1_qp1,
                           c1_b_A, c1_qp1jj);
          }
        }
      }

      c1_f_qp1 = c1_qp1;
      c1_fb_a = c1_f_qp1;
      c1_gb_a = c1_fb_a;
      if (c1_gb_a > 6) {
        c1_f_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_f_overflow = false;
      }

      if (c1_f_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_f_overflow);
      }

      for (c1_d_ii = c1_f_qp1; c1_d_ii < 7; c1_d_ii++) {
        c1_b_ii = c1_d_ii;
        c1_Vf[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c1_b_ii), 1, 6, 1, 0) + 6 *
               (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_q), 1, 6, 2, 0) - 1)) - 1] =
          c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_ii), 1, 6, 1, 0) - 1];
      }
    }
  }

  c1_m = 6;
  c1_s[5] = c1_b_A[35];
  c1_e[4] = c1_b_A[34];
  c1_e[5] = 0.0;
  for (c1_e_ii = 1; c1_e_ii < 7; c1_e_ii++) {
    c1_b_ii = c1_e_ii;
    c1_U[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_ii), 1, 6, 1, 0) + 29] = 0.0;
  }

  c1_U[35] = 1.0;
  for (c1_d_q = 5; c1_d_q > 0; c1_d_q--) {
    c1_b_q = c1_d_q;
    c1_hb_a = c1_b_q;
    c1_ib_a = c1_hb_a;
    c1_qp1 = c1_ib_a;
    c1_u_b = c1_b_q;
    c1_v_b = c1_u_b;
    c1_nmq = 6 - c1_v_b;
    c1_jb_a = c1_nmq;
    c1_kb_a = c1_jb_a + 1;
    c1_nmqp1 = c1_kb_a;
    c1_lb_a = c1_b_q;
    c1_mb_a = c1_lb_a;
    c1_h_c = c1_mb_a;
    c1_w_b = c1_h_c - 1;
    c1_x_b = c1_w_b;
    c1_i_c = 6 * c1_x_b;
    c1_nb_a = c1_b_q;
    c1_y_b = c1_i_c;
    c1_ob_a = c1_nb_a;
    c1_ab_b = c1_y_b;
    c1_qq = c1_ob_a + c1_ab_b;
    if (c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
      c1_g_qp1 = c1_qp1 + 1;
      c1_pb_a = c1_g_qp1;
      c1_qb_a = c1_pb_a;
      if (c1_qb_a > 6) {
        c1_g_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_g_overflow = false;
      }

      if (c1_g_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_g_overflow);
      }

      for (c1_e_jj = c1_g_qp1; c1_e_jj < 7; c1_e_jj++) {
        c1_b_jj = c1_e_jj;
        c1_rb_a = c1_b_jj;
        c1_sb_a = c1_rb_a;
        c1_j_c = c1_sb_a;
        c1_bb_b = c1_j_c - 1;
        c1_cb_b = c1_bb_b;
        c1_k_c = 6 * c1_cb_b;
        c1_tb_a = c1_b_q;
        c1_db_b = c1_k_c;
        c1_ub_a = c1_tb_a;
        c1_eb_b = c1_db_b;
        c1_qjj = c1_ub_a + c1_eb_b;
        for (c1_i58 = 0; c1_i58 < 36; c1_i58++) {
          c1_b_U[c1_i58] = c1_U[c1_i58];
        }

        for (c1_i59 = 0; c1_i59 < 36; c1_i59++) {
          c1_c_U[c1_i59] = c1_U[c1_i59];
        }

        c1_t = c1_eml_xdotc(chartInstance, c1_nmqp1, c1_b_U, c1_qq, c1_c_U,
                            c1_qjj);
        c1_t = -c1_eml_div(chartInstance, c1_t, c1_U[_SFD_EML_ARRAY_BOUNDS_CHECK
                           ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qq),
                            1, 36, 1, 0) - 1]);
        c1_d_eml_xaxpy(chartInstance, c1_nmqp1, c1_t, c1_qq, c1_U, c1_qjj);
      }

      c1_e_q = c1_b_q;
      c1_vb_a = c1_e_q;
      c1_wb_a = c1_vb_a;
      if (c1_wb_a > 6) {
        c1_h_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_h_overflow = false;
      }

      if (c1_h_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_h_overflow);
      }

      for (c1_f_ii = c1_e_q; c1_f_ii < 7; c1_f_ii++) {
        c1_b_ii = c1_f_ii;
        c1_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_ii), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0)
               - 1)) - 1] = -c1_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1, 6, 1, 0) + 6 *
          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 2, 0) - 1)) - 1];
      }

      c1_U[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_qq), 1, 36, 1, 0) - 1] = c1_U[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_qq), 1, 36, 1, 0) - 1] + 1.0;
      c1_xb_a = c1_b_q;
      c1_yb_a = c1_xb_a - 1;
      c1_i60 = c1_yb_a;
      c1_fb_b = c1_i60;
      c1_gb_b = c1_fb_b;
      if (1 > c1_gb_b) {
        c1_i_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_i_overflow = (c1_gb_b > 2147483646);
      }

      if (c1_i_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_i_overflow);
      }

      for (c1_g_ii = 1; c1_g_ii <= c1_i60; c1_g_ii++) {
        c1_b_ii = c1_g_ii;
        c1_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_ii), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0)
               - 1)) - 1] = 0.0;
      }
    } else {
      for (c1_h_ii = 1; c1_h_ii < 7; c1_h_ii++) {
        c1_b_ii = c1_h_ii;
        c1_U[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_ii), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0)
               - 1)) - 1] = 0.0;
      }

      c1_U[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_qq), 1, 36, 1, 0) - 1] = 1.0;
    }
  }

  for (c1_f_q = 6; c1_f_q > 0; c1_f_q--) {
    c1_b_q = c1_f_q;
    if (c1_b_q <= 4) {
      if (c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
        c1_ac_a = c1_b_q;
        c1_bc_a = c1_ac_a + 1;
        c1_qp1 = c1_bc_a;
        c1_hb_b = c1_b_q;
        c1_ib_b = c1_hb_b;
        c1_pmq = 6 - c1_ib_b;
        c1_cc_a = c1_b_q;
        c1_dc_a = c1_cc_a;
        c1_l_c = c1_dc_a;
        c1_jb_b = c1_l_c - 1;
        c1_kb_b = c1_jb_b;
        c1_m_c = 6 * c1_kb_b;
        c1_ec_a = c1_qp1;
        c1_lb_b = c1_m_c;
        c1_fc_a = c1_ec_a;
        c1_mb_b = c1_lb_b;
        c1_qp1q = c1_fc_a + c1_mb_b;
        c1_h_qp1 = c1_qp1;
        c1_gc_a = c1_h_qp1;
        c1_hc_a = c1_gc_a;
        if (c1_hc_a > 6) {
          c1_j_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_j_overflow = false;
        }

        if (c1_j_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_j_overflow);
        }

        for (c1_f_jj = c1_h_qp1; c1_f_jj < 7; c1_f_jj++) {
          c1_b_jj = c1_f_jj;
          c1_ic_a = c1_b_jj;
          c1_jc_a = c1_ic_a;
          c1_n_c = c1_jc_a;
          c1_nb_b = c1_n_c - 1;
          c1_ob_b = c1_nb_b;
          c1_o_c = 6 * c1_ob_b;
          c1_kc_a = c1_qp1;
          c1_pb_b = c1_o_c;
          c1_lc_a = c1_kc_a;
          c1_qb_b = c1_pb_b;
          c1_qp1jj = c1_lc_a + c1_qb_b;
          for (c1_i61 = 0; c1_i61 < 36; c1_i61++) {
            c1_b_Vf[c1_i61] = c1_Vf[c1_i61];
          }

          for (c1_i62 = 0; c1_i62 < 36; c1_i62++) {
            c1_c_Vf[c1_i62] = c1_Vf[c1_i62];
          }

          c1_t = c1_eml_xdotc(chartInstance, c1_pmq, c1_b_Vf, c1_qp1q, c1_c_Vf,
                              c1_qp1jj);
          c1_t = -c1_eml_div(chartInstance, c1_t,
                             c1_Vf[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_qp1q), 1, 36, 1, 0) - 1]);
          c1_d_eml_xaxpy(chartInstance, c1_pmq, c1_t, c1_qp1q, c1_Vf, c1_qp1jj);
        }
      }
    }

    for (c1_i_ii = 1; c1_i_ii < 7; c1_i_ii++) {
      c1_b_ii = c1_i_ii;
      c1_Vf[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               (real_T)c1_b_ii), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
               "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0)
              - 1)) - 1] = 0.0;
    }

    c1_Vf[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             (real_T)c1_b_q), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 2, 0) - 1))
      - 1] = 1.0;
  }

  for (c1_g_q = 1; c1_g_q < 7; c1_g_q++) {
    c1_b_q = c1_g_q;
    if (c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
      c1_rt = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
      c1_r = c1_eml_div(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
        (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1], c1_rt);
      c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_q), 1, 6, 1, 0) - 1] = c1_rt;
      if (c1_b_q < 6) {
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = c1_eml_div(chartInstance,
          c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1], c1_r);
      }

      if (c1_b_q <= 6) {
        c1_mc_a = c1_b_q;
        c1_nc_a = c1_mc_a;
        c1_p_c = c1_nc_a;
        c1_rb_b = c1_p_c - 1;
        c1_sb_b = c1_rb_b;
        c1_q_c = 6 * c1_sb_b;
        c1_tb_b = c1_q_c;
        c1_ub_b = c1_tb_b;
        c1_colq = c1_ub_b;
        c1_f_eml_xscal(chartInstance, c1_r, c1_U, c1_colq + 1);
      }
    }

    if (c1_b_q < 6) {
      if (c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_q), 1, 6, 1, 0) - 1] != 0.0) {
        c1_rt = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
        c1_r = c1_eml_div(chartInstance, c1_rt, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK
                          ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q),
                           1, 6, 1, 0) - 1]);
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_q), 1, 6, 1, 0) - 1] = c1_rt;
        c1_oc_a = c1_b_q;
        c1_pc_a = c1_oc_a;
        c1_r_c = c1_pc_a;
        c1_qc_a = c1_b_q;
        c1_rc_a = c1_qc_a;
        c1_s_c = c1_rc_a;
        c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c1_r_c + 1)), 1, 6, 1, 0) - 1] =
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c1_s_c + 1)), 1, 6, 1, 0) - 1] * c1_r;
        c1_vb_b = c1_b_q;
        c1_wb_b = c1_vb_b;
        c1_t_c = 6 * c1_wb_b;
        c1_xb_b = c1_t_c;
        c1_yb_b = c1_xb_b;
        c1_colqp1 = c1_yb_b;
        c1_f_eml_xscal(chartInstance, c1_r, c1_Vf, c1_colqp1 + 1);
      }
    }
  }

  c1_iter = 0.0;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  c1_tiny = c1_eml_div(chartInstance, 2.2250738585072014E-308,
                       2.2204460492503131E-16);
  c1_snorm = 0.0;
  for (c1_j_ii = 1; c1_j_ii < 7; c1_j_ii++) {
    c1_b_ii = c1_j_ii;
    c1_varargin_1 = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1, 6, 1, 0) - 1]);
    c1_varargin_2 = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1, 6, 1, 0) - 1]);
    c1_b_varargin_2 = c1_varargin_1;
    c1_varargin_3 = c1_varargin_2;
    c1_x = c1_b_varargin_2;
    c1_c_y = c1_varargin_3;
    c1_b_x = c1_x;
    c1_d_y = c1_c_y;
    c1_b_eml_scalar_eg(chartInstance);
    c1_xk = c1_b_x;
    c1_yk = c1_d_y;
    c1_c_x = c1_xk;
    c1_e_y = c1_yk;
    c1_b_eml_scalar_eg(chartInstance);
    c1_maxval = muDoubleScalarMax(c1_c_x, c1_e_y);
    c1_b_varargin_1 = c1_snorm;
    c1_c_varargin_2 = c1_maxval;
    c1_d_varargin_2 = c1_b_varargin_1;
    c1_b_varargin_3 = c1_c_varargin_2;
    c1_d_x = c1_d_varargin_2;
    c1_f_y = c1_b_varargin_3;
    c1_e_x = c1_d_x;
    c1_g_y = c1_f_y;
    c1_b_eml_scalar_eg(chartInstance);
    c1_b_xk = c1_e_x;
    c1_b_yk = c1_g_y;
    c1_f_x = c1_b_xk;
    c1_h_y = c1_b_yk;
    c1_b_eml_scalar_eg(chartInstance);
    c1_snorm = muDoubleScalarMax(c1_f_x, c1_h_y);
  }

  exitg1 = false;
  while ((exitg1 == false) && (c1_m > 0)) {
    if (c1_iter >= 75.0) {
      c1_b_eml_error(chartInstance);
      exitg1 = true;
    } else {
      c1_sc_a = c1_m;
      c1_tc_a = c1_sc_a - 1;
      c1_b_q = c1_tc_a;
      c1_uc_a = c1_m;
      c1_vc_a = c1_uc_a - 1;
      c1_i63 = c1_vc_a;
      c1_wc_a = c1_i63;
      c1_xc_a = c1_wc_a;
      if (c1_xc_a < 0) {
        c1_k_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_k_overflow = false;
      }

      if (c1_k_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_k_overflow);
      }

      c1_k_ii = c1_i63;
      guard3 = false;
      guard4 = false;
      exitg5 = false;
      while ((exitg5 == false) && (c1_k_ii > -1)) {
        c1_b_ii = c1_k_ii;
        c1_b_q = c1_b_ii;
        if (c1_b_ii == 0) {
          exitg5 = true;
        } else {
          c1_yc_a = c1_b_ii;
          c1_ad_a = c1_yc_a;
          c1_u_c = c1_ad_a;
          c1_test0 = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1,
            6, 1, 0) - 1]) + c1_abs(chartInstance,
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                                      (real_T)(c1_u_c + 1)), 1, 6, 1, 0) - 1]);
          c1_ztest0 = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii),
            1, 6, 1, 0) - 1]);
          c1_eps(chartInstance);
          if (c1_ztest0 <= 2.2204460492503131E-16 * c1_test0) {
            guard4 = true;
            exitg5 = true;
          } else if (c1_ztest0 <= c1_tiny) {
            guard4 = true;
            exitg5 = true;
          } else {
            guard11 = false;
            if (c1_iter > 20.0) {
              c1_eps(chartInstance);
              if (c1_ztest0 <= 2.2204460492503131E-16 * c1_snorm) {
                guard3 = true;
                exitg5 = true;
              } else {
                guard11 = true;
              }
            } else {
              guard11 = true;
            }

            if (guard11 == true) {
              c1_k_ii--;
              guard3 = false;
              guard4 = false;
            }
          }
        }
      }

      if (guard4 == true) {
        guard3 = true;
      }

      if (guard3 == true) {
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_ii), 1, 6, 1, 0) - 1] = 0.0;
      }

      c1_bd_a = c1_m;
      c1_cd_a = c1_bd_a;
      c1_v_c = c1_cd_a;
      if (c1_b_q == c1_v_c - 1) {
        c1_kase = 4.0;
      } else {
        c1_qs = c1_m;
        c1_b_m = c1_m;
        c1_h_q = c1_b_q;
        c1_dd_a = c1_b_m;
        c1_ac_b = c1_h_q;
        c1_ed_a = c1_dd_a;
        c1_bc_b = c1_ac_b;
        if (c1_ed_a < c1_bc_b) {
          c1_l_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_l_overflow = (c1_bc_b < -2147483647);
        }

        if (c1_l_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_l_overflow);
        }

        c1_l_ii = c1_b_m;
        guard2 = false;
        exitg4 = false;
        while ((exitg4 == false) && (c1_l_ii >= c1_h_q)) {
          c1_b_ii = c1_l_ii;
          c1_qs = c1_b_ii;
          if (c1_b_ii == c1_b_q) {
            exitg4 = true;
          } else {
            c1_test = 0.0;
            if (c1_b_ii < c1_m) {
              c1_test = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii), 1, 6, 1, 0)
                               - 1]);
            }

            c1_fd_a = c1_b_q;
            c1_gd_a = c1_fd_a;
            c1_w_c = c1_gd_a;
            if (c1_b_ii > c1_w_c + 1) {
              c1_hd_a = c1_b_ii;
              c1_id_a = c1_hd_a;
              c1_x_c = c1_id_a;
              c1_test += c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)(c1_x_c - 1)), 1, 6,
                1, 0) - 1]);
            }

            c1_ztest = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                               (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_ii),
              1, 6, 1, 0) - 1]);
            c1_eps(chartInstance);
            if (c1_ztest <= 2.2204460492503131E-16 * c1_test) {
              guard2 = true;
              exitg4 = true;
            } else if (c1_ztest <= c1_tiny) {
              guard2 = true;
              exitg4 = true;
            } else {
              c1_l_ii--;
              guard2 = false;
            }
          }
        }

        if (guard2 == true) {
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ii), 1, 6, 1, 0) - 1] = 0.0;
        }

        if (c1_qs == c1_b_q) {
          c1_kase = 3.0;
        } else if (c1_qs == c1_m) {
          c1_kase = 1.0;
        } else {
          c1_kase = 2.0;
          c1_b_q = c1_qs;
        }
      }

      c1_jd_a = c1_b_q;
      c1_kd_a = c1_jd_a + 1;
      c1_b_q = c1_kd_a;
      switch ((int32_T)_SFD_INTEGER_CHECK("", c1_kase)) {
       case 1:
        c1_ld_a = c1_m;
        c1_md_a = c1_ld_a;
        c1_y_c = c1_md_a;
        c1_f = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)(c1_y_c - 1)), 1, 6, 1, 0) - 1];
        c1_nd_a = c1_m;
        c1_od_a = c1_nd_a;
        c1_ab_c = c1_od_a;
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c1_ab_c - 1)), 1, 6, 1, 0) - 1] = 0.0;
        c1_pd_a = c1_m;
        c1_qd_a = c1_pd_a - 1;
        c1_i64 = c1_qd_a;
        c1_i_q = c1_b_q;
        c1_rd_a = c1_i64;
        c1_cc_b = c1_i_q;
        c1_sd_a = c1_rd_a;
        c1_dc_b = c1_cc_b;
        if (c1_sd_a < c1_dc_b) {
          c1_m_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_m_overflow = (c1_dc_b < -2147483647);
        }

        if (c1_m_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_m_overflow);
        }

        for (c1_k = c1_i64; c1_k >= c1_i_q; c1_k--) {
          c1_b_k = c1_k;
          c1_t1 = c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
          c1_b_t1 = c1_t1;
          c1_b_f = c1_f;
          c1_b_eml_xrotg(chartInstance, &c1_b_t1, &c1_b_f, &c1_cs, &c1_sn);
          c1_t1 = c1_b_t1;
          c1_f = c1_b_f;
          c1_b_cs = c1_cs;
          c1_b_sn = c1_sn;
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] = c1_t1;
          if (c1_b_k > c1_b_q) {
            c1_td_a = c1_b_k;
            c1_ud_a = c1_td_a - 1;
            c1_km1 = c1_ud_a;
            c1_f = -c1_b_sn * c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_km1), 1, 6, 1, 0) - 1];
            c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_km1), 1, 6, 1, 0) - 1] =
              c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
              "", (real_T)c1_km1), 1, 6, 1, 0) - 1] * c1_b_cs;
          }

          c1_vd_a = c1_b_k;
          c1_wd_a = c1_vd_a;
          c1_bb_c = c1_wd_a;
          c1_ec_b = c1_bb_c - 1;
          c1_fc_b = c1_ec_b;
          c1_cb_c = 6 * c1_fc_b;
          c1_gc_b = c1_cb_c;
          c1_hc_b = c1_gc_b;
          c1_colk = c1_hc_b;
          c1_xd_a = c1_m;
          c1_yd_a = c1_xd_a;
          c1_db_c = c1_yd_a;
          c1_ic_b = c1_db_c - 1;
          c1_jc_b = c1_ic_b;
          c1_eb_c = 6 * c1_jc_b;
          c1_kc_b = c1_eb_c;
          c1_lc_b = c1_kc_b;
          c1_colm = c1_lc_b;
          c1_b_eml_xrot(chartInstance, c1_Vf, c1_colk + 1, c1_colm + 1, c1_b_cs,
                        c1_b_sn);
        }
        break;

       case 2:
        c1_ae_a = c1_b_q;
        c1_be_a = c1_ae_a - 1;
        c1_qm1 = c1_be_a;
        c1_f = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_qm1), 1, 6, 1, 0) - 1];
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_qm1), 1, 6, 1, 0) - 1] = 0.0;
        c1_j_q = c1_b_q;
        c1_c_m = c1_m;
        c1_ce_a = c1_j_q;
        c1_mc_b = c1_c_m;
        c1_de_a = c1_ce_a;
        c1_nc_b = c1_mc_b;
        if (c1_de_a > c1_nc_b) {
          c1_n_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_n_overflow = (c1_nc_b > 2147483646);
        }

        if (c1_n_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_n_overflow);
        }

        for (c1_c_k = c1_j_q; c1_c_k <= c1_c_m; c1_c_k++) {
          c1_b_k = c1_c_k;
          c1_t1 = c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
          c1_c_t1 = c1_t1;
          c1_unusedU0 = c1_f;
          c1_b_eml_xrotg(chartInstance, &c1_c_t1, &c1_unusedU0, &c1_c_cs,
                         &c1_c_sn);
          c1_t1 = c1_c_t1;
          c1_b_cs = c1_c_cs;
          c1_b_sn = c1_c_sn;
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] = c1_t1;
          c1_f = -c1_b_sn * c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
          c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK
            ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) -
            1] * c1_b_cs;
          c1_ee_a = c1_b_k;
          c1_fe_a = c1_ee_a;
          c1_fb_c = c1_fe_a;
          c1_oc_b = c1_fb_c - 1;
          c1_pc_b = c1_oc_b;
          c1_gb_c = 6 * c1_pc_b;
          c1_qc_b = c1_gb_c;
          c1_rc_b = c1_qc_b;
          c1_colk = c1_rc_b;
          c1_ge_a = c1_qm1;
          c1_he_a = c1_ge_a;
          c1_hb_c = c1_he_a;
          c1_sc_b = c1_hb_c - 1;
          c1_tc_b = c1_sc_b;
          c1_ib_c = 6 * c1_tc_b;
          c1_uc_b = c1_ib_c;
          c1_vc_b = c1_uc_b;
          c1_colqm1 = c1_vc_b;
          c1_b_eml_xrot(chartInstance, c1_U, c1_colk + 1, c1_colqm1 + 1, c1_b_cs,
                        c1_b_sn);
        }
        break;

       case 3:
        c1_ie_a = c1_m;
        c1_je_a = c1_ie_a - 1;
        c1_mm1 = c1_je_a;
        c1_d3 = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_m), 1, 6, 1, 0) - 1]);
        c1_d4 = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_mm1), 1, 6, 1, 0) - 1]);
        c1_d5 = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_mm1), 1, 6, 1, 0) - 1]);
        c1_d6 = c1_abs(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
        c1_d7 = c1_abs(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1]);
        c1_c_varargin_1[0] = c1_d3;
        c1_c_varargin_1[1] = c1_d4;
        c1_c_varargin_1[2] = c1_d5;
        c1_c_varargin_1[3] = c1_d6;
        c1_c_varargin_1[4] = c1_d7;
        c1_ixstart = 1;
        c1_mtmp = c1_c_varargin_1[0];
        c1_g_x = c1_mtmp;
        c1_wc_b = muDoubleScalarIsNaN(c1_g_x);
        if (c1_wc_b) {
          c1_eml_switch_helper(chartInstance);
          c1_ix = 2;
          exitg2 = false;
          while ((exitg2 == false) && (c1_ix < 6)) {
            c1_b_ix = c1_ix;
            c1_ixstart = c1_b_ix;
            c1_h_x = c1_c_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_ix), 1, 5, 1, 0) - 1];
            c1_xc_b = muDoubleScalarIsNaN(c1_h_x);
            if (!c1_xc_b) {
              c1_mtmp = c1_c_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_ix), 1, 5, 1, 0) - 1];
              exitg2 = true;
            } else {
              c1_ix++;
            }
          }
        }

        if (c1_ixstart < 5) {
          c1_ke_a = c1_ixstart;
          c1_le_a = c1_ke_a + 1;
          c1_i65 = c1_le_a;
          c1_me_a = c1_i65;
          c1_ne_a = c1_me_a;
          if (c1_ne_a > 5) {
            c1_o_overflow = false;
          } else {
            c1_eml_switch_helper(chartInstance);
            c1_o_overflow = false;
          }

          if (c1_o_overflow) {
            c1_check_forloop_overflow_error(chartInstance, c1_o_overflow);
          }

          for (c1_c_ix = c1_i65; c1_c_ix < 6; c1_c_ix++) {
            c1_b_ix = c1_c_ix;
            c1_oe_a = c1_c_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_ix), 1, 5, 1, 0) - 1];
            c1_yc_b = c1_mtmp;
            c1_p = (c1_oe_a > c1_yc_b);
            if (c1_p) {
              c1_mtmp = c1_c_varargin_1[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
                _SFD_INTEGER_CHECK("", (real_T)c1_b_ix), 1, 5, 1, 0) - 1];
            }
          }
        }

        c1_b_mtmp = c1_mtmp;
        c1_scale = c1_b_mtmp;
        c1_sm = c1_eml_div(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_m), 1, 6, 1, 0) - 1],
                           c1_scale);
        c1_smm1 = c1_eml_div(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_mm1), 1,
          6, 1, 0) - 1], c1_scale);
        c1_emm1 = c1_eml_div(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_mm1), 1,
          6, 1, 0) - 1], c1_scale);
        c1_sqds = c1_eml_div(chartInstance, c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1,
          6, 1, 0) - 1], c1_scale);
        c1_eqds = c1_eml_div(chartInstance, c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("",
                              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1,
          6, 1, 0) - 1], c1_scale);
        c1_ad_b = c1_eml_div(chartInstance, (c1_smm1 + c1_sm) * (c1_smm1 - c1_sm)
                             + c1_emm1 * c1_emm1, 2.0);
        c1_jb_c = c1_sm * c1_emm1;
        c1_jb_c *= c1_jb_c;
        c1_shift = 0.0;
        guard1 = false;
        if (c1_ad_b != 0.0) {
          guard1 = true;
        } else {
          if (c1_jb_c != 0.0) {
            guard1 = true;
          }
        }

        if (guard1 == true) {
          c1_shift = c1_ad_b * c1_ad_b + c1_jb_c;
          c1_b_sqrt(chartInstance, &c1_shift);
          if (c1_ad_b < 0.0) {
            c1_shift = -c1_shift;
          }

          c1_shift = c1_eml_div(chartInstance, c1_jb_c, c1_ad_b + c1_shift);
        }

        c1_f = (c1_sqds + c1_sm) * (c1_sqds - c1_sm) + c1_shift;
        c1_g = c1_sqds * c1_eqds;
        c1_k_q = c1_b_q;
        c1_b_mm1 = c1_mm1;
        c1_pe_a = c1_k_q;
        c1_bd_b = c1_b_mm1;
        c1_qe_a = c1_pe_a;
        c1_cd_b = c1_bd_b;
        if (c1_qe_a > c1_cd_b) {
          c1_p_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_p_overflow = (c1_cd_b > 2147483646);
        }

        if (c1_p_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_p_overflow);
        }

        for (c1_d_k = c1_k_q; c1_d_k <= c1_b_mm1; c1_d_k++) {
          c1_b_k = c1_d_k;
          c1_re_a = c1_b_k;
          c1_se_a = c1_re_a;
          c1_km1 = c1_se_a;
          c1_te_a = c1_b_k;
          c1_ue_a = c1_te_a + 1;
          c1_kp1 = c1_ue_a;
          c1_c_f = c1_f;
          c1_unusedU1 = c1_g;
          c1_b_eml_xrotg(chartInstance, &c1_c_f, &c1_unusedU1, &c1_d_cs,
                         &c1_d_sn);
          c1_f = c1_c_f;
          c1_b_cs = c1_d_cs;
          c1_b_sn = c1_d_sn;
          if (c1_b_k > c1_b_q) {
            c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)(c1_km1 - 1)), 1, 6, 1, 0) - 1] = c1_f;
          }

          c1_f = c1_b_cs * c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1] + c1_b_sn *
            c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1];
          c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] = c1_b_cs *
            c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] - c1_b_sn *
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1];
          c1_g = c1_b_sn * c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_kp1), 1, 6, 1, 0) - 1];
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_kp1), 1, 6, 1, 0) - 1] = c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK
            ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_kp1), 1, 6, 1, 0) -
            1] * c1_b_cs;
          c1_ve_a = c1_b_k;
          c1_we_a = c1_ve_a;
          c1_kb_c = c1_we_a;
          c1_dd_b = c1_kb_c - 1;
          c1_ed_b = c1_dd_b;
          c1_lb_c = 6 * c1_ed_b;
          c1_fd_b = c1_lb_c;
          c1_gd_b = c1_fd_b;
          c1_colk = c1_gd_b;
          c1_hd_b = c1_b_k;
          c1_id_b = c1_hd_b;
          c1_mb_c = 6 * c1_id_b;
          c1_jd_b = c1_mb_c;
          c1_kd_b = c1_jd_b;
          c1_colkp1 = c1_kd_b;
          c1_b_eml_xrot(chartInstance, c1_Vf, c1_colk + 1, c1_colkp1 + 1,
                        c1_b_cs, c1_b_sn);
          c1_d_f = c1_f;
          c1_unusedU2 = c1_g;
          c1_b_eml_xrotg(chartInstance, &c1_d_f, &c1_unusedU2, &c1_e_cs,
                         &c1_e_sn);
          c1_f = c1_d_f;
          c1_b_cs = c1_e_cs;
          c1_b_sn = c1_e_sn;
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] = c1_f;
          c1_f = c1_b_cs * c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1] + c1_b_sn *
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_kp1), 1, 6, 1, 0) - 1];
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_kp1), 1, 6, 1, 0) - 1] = -c1_b_sn *
            c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) - 1] + c1_b_cs *
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_kp1), 1, 6, 1, 0) - 1];
          c1_g = c1_b_sn * c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_kp1), 1, 6, 1, 0) - 1];
          c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_kp1), 1, 6, 1, 0) - 1] = c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK
            ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_kp1), 1, 6, 1, 0) -
            1] * c1_b_cs;
          if (c1_b_k < 6) {
            c1_xe_a = c1_b_k;
            c1_ye_a = c1_xe_a;
            c1_nb_c = c1_ye_a;
            c1_ld_b = c1_nb_c - 1;
            c1_md_b = c1_ld_b;
            c1_ob_c = 6 * c1_md_b;
            c1_nd_b = c1_ob_c;
            c1_od_b = c1_nd_b;
            c1_colk = c1_od_b;
            c1_pd_b = c1_b_k;
            c1_qd_b = c1_pd_b;
            c1_pb_c = 6 * c1_qd_b;
            c1_rd_b = c1_pb_c;
            c1_sd_b = c1_rd_b;
            c1_colkp1 = c1_sd_b;
            c1_b_eml_xrot(chartInstance, c1_U, c1_colk + 1, c1_colkp1 + 1,
                          c1_b_cs, c1_b_sn);
          }
        }

        c1_af_a = c1_m;
        c1_bf_a = c1_af_a;
        c1_qb_c = c1_bf_a;
        c1_e[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c1_qb_c - 1)), 1, 6, 1, 0) - 1] = c1_f;
        c1_iter++;
        break;

       default:
        if (c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_q), 1, 6, 1, 0) - 1] < 0.0) {
          c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_q), 1, 6, 1, 0) - 1] =
            -c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_q), 1, 6, 1, 0) - 1];
          c1_cf_a = c1_b_q;
          c1_if_a = c1_cf_a;
          c1_rb_c = c1_if_a;
          c1_td_b = c1_rb_c - 1;
          c1_ee_b = c1_td_b;
          c1_sb_c = 6 * c1_ee_b;
          c1_ud_b = c1_sb_c;
          c1_fe_b = c1_ud_b;
          c1_colq = c1_fe_b;
          c1_eml_scalar_eg(chartInstance);
          c1_d8 = -1.0;
          c1_f_eml_xscal(chartInstance, c1_d8, c1_Vf, c1_colq + 1);
        }

        c1_df_a = c1_b_q;
        c1_jf_a = c1_df_a + 1;
        c1_qp1 = c1_jf_a;
        exitg3 = false;
        while ((exitg3 == false) && (c1_b_q < 6)) {
          if (c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_b_q), 1, 6, 1, 0) - 1] <
              c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                "", (real_T)c1_qp1), 1, 6, 1, 0) - 1]) {
            c1_rt = c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
              _SFD_INTEGER_CHECK("", (real_T)c1_b_q), 1, 6, 1, 0) - 1];
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_q), 1, 6, 1, 0) - 1] =
              c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
              "", (real_T)c1_qp1), 1, 6, 1, 0) - 1];
            c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_qp1), 1, 6, 1, 0) - 1] = c1_rt;
            if (c1_b_q < 6) {
              c1_ff_a = c1_b_q;
              c1_lf_a = c1_ff_a;
              c1_tb_c = c1_lf_a;
              c1_vd_b = c1_tb_c - 1;
              c1_ge_b = c1_vd_b;
              c1_ub_c = 6 * c1_ge_b;
              c1_wd_b = c1_ub_c;
              c1_he_b = c1_wd_b;
              c1_colq = c1_he_b;
              c1_xd_b = c1_b_q;
              c1_ie_b = c1_xd_b;
              c1_vb_c = 6 * c1_ie_b;
              c1_yd_b = c1_vb_c;
              c1_je_b = c1_yd_b;
              c1_colqp1 = c1_je_b;
              c1_b_eml_xswap(chartInstance, c1_Vf, c1_colq + 1, c1_colqp1 + 1);
            }

            if (c1_b_q < 6) {
              c1_gf_a = c1_b_q;
              c1_mf_a = c1_gf_a;
              c1_wb_c = c1_mf_a;
              c1_ae_b = c1_wb_c - 1;
              c1_ke_b = c1_ae_b;
              c1_xb_c = 6 * c1_ke_b;
              c1_be_b = c1_xb_c;
              c1_le_b = c1_be_b;
              c1_colq = c1_le_b;
              c1_ce_b = c1_b_q;
              c1_me_b = c1_ce_b;
              c1_yb_c = 6 * c1_me_b;
              c1_de_b = c1_yb_c;
              c1_ne_b = c1_de_b;
              c1_colqp1 = c1_ne_b;
              c1_b_eml_xswap(chartInstance, c1_U, c1_colq + 1, c1_colqp1 + 1);
            }

            c1_b_q = c1_qp1;
            c1_hf_a = c1_b_q;
            c1_nf_a = c1_hf_a + 1;
            c1_qp1 = c1_nf_a;
          } else {
            exitg3 = true;
          }
        }

        c1_iter = 0.0;
        c1_ef_a = c1_m;
        c1_kf_a = c1_ef_a - 1;
        c1_m = c1_kf_a;
        break;
      }
    }
  }

  for (c1_e_k = 1; c1_e_k < 7; c1_e_k++) {
    c1_b_k = c1_e_k;
    c1_S[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_k), 1, 6, 1, 0) - 1] = c1_s[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
  }

  for (c1_j = 1; c1_j < 7; c1_j++) {
    c1_b_j = c1_j;
    for (c1_i = 1; c1_i < 7; c1_i++) {
      c1_b_i = c1_i;
      c1_V[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c1_b_i), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2, 0) - 1))
        - 1] = c1_Vf[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 6, 1, 0) + 6 *
                      (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2, 0) - 1)) - 1];
    }
  }
}

static real_T c1_eml_xnrm2(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[36], int32_T c1_ix0)
{
  real_T c1_y;
  int32_T c1_b_n;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  int32_T c1_c_ix0;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_scale;
  int32_T c1_kstart;
  int32_T c1_a;
  int32_T c1_c;
  int32_T c1_b_a;
  int32_T c1_b_c;
  int32_T c1_c_a;
  int32_T c1_b;
  int32_T c1_kend;
  int32_T c1_b_kstart;
  int32_T c1_b_kend;
  int32_T c1_d_a;
  int32_T c1_b_b;
  int32_T c1_e_a;
  int32_T c1_c_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_absxk;
  real_T c1_t;
  c1_b_n = c1_n;
  c1_b_ix0 = c1_ix0;
  c1_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_ix0 = c1_b_ix0;
  c1_y = 0.0;
  if (c1_c_n < 1) {
  } else if (c1_c_n == 1) {
    c1_b_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_c_ix0), 1, 36, 1, 0) - 1];
    c1_c_x = c1_b_x;
    c1_y = muDoubleScalarAbs(c1_c_x);
  } else {
    c1_realmin(chartInstance);
    c1_scale = 2.2250738585072014E-308;
    c1_kstart = c1_c_ix0;
    c1_a = c1_c_n;
    c1_c = c1_a;
    c1_b_a = c1_c - 1;
    c1_b_c = c1_b_a;
    c1_c_a = c1_kstart;
    c1_b = c1_b_c;
    c1_kend = c1_c_a + c1_b;
    c1_b_kstart = c1_kstart;
    c1_b_kend = c1_kend;
    c1_d_a = c1_b_kstart;
    c1_b_b = c1_b_kend;
    c1_e_a = c1_d_a;
    c1_c_b = c1_b_b;
    if (c1_e_a > c1_c_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_c_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = c1_b_kstart; c1_k <= c1_b_kend; c1_k++) {
      c1_b_k = c1_k;
      c1_d_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_b_k), 1, 36, 1, 0) - 1];
      c1_e_x = c1_d_x;
      c1_absxk = muDoubleScalarAbs(c1_e_x);
      if (c1_absxk > c1_scale) {
        c1_t = c1_scale / c1_absxk;
        c1_y = 1.0 + c1_y * c1_t * c1_t;
        c1_scale = c1_absxk;
      } else {
        c1_t = c1_absxk / c1_scale;
        c1_y += c1_t * c1_t;
      }
    }

    c1_y = c1_scale * muDoubleScalarSqrt(c1_y);
  }

  return c1_y;
}

static void c1_threshold(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_abs(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x)
{
  real_T c1_b_x;
  (void)chartInstance;
  c1_b_x = c1_x;
  return muDoubleScalarAbs(c1_b_x);
}

static void c1_realmin(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_check_forloop_overflow_error(SFc1_UR5ModelInstanceStruct
  *chartInstance, boolean_T c1_overflow)
{
  int32_T c1_i66;
  static char_T c1_cv1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  int32_T c1_i67;
  static char_T c1_cv2[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c1_b_u[23];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  if (!c1_overflow) {
  } else {
    for (c1_i66 = 0; c1_i66 < 34; c1_i66++) {
      c1_u[c1_i66] = c1_cv1[c1_i66];
    }

    c1_y = NULL;
    sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c1_i67 = 0; c1_i67 < 23; c1_i67++) {
      c1_b_u[c1_i67] = c1_cv2[c1_i67];
    }

    c1_b_y = NULL;
    sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c1_y, 14, c1_b_y));
  }
}

static void c1_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_b_x[36])
{
  int32_T c1_i68;
  for (c1_i68 = 0; c1_i68 < 36; c1_i68++) {
    c1_b_x[c1_i68] = c1_x[c1_i68];
  }

  c1_d_eml_xscal(chartInstance, c1_n, c1_a, c1_b_x, c1_ix0);
}

static void c1_b_threshold(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_eml_xdotc(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0)
{
  real_T c1_d;
  int32_T c1_b_n;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_n;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_d_n;
  int32_T c1_d_ix0;
  int32_T c1_d_iy0;
  int32_T c1_e_n;
  int32_T c1_e_ix0;
  int32_T c1_e_iy0;
  int32_T c1_ix;
  int32_T c1_iy;
  int32_T c1_f_n;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_a;
  int32_T c1_b_a;
  c1_b_n = c1_n;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_c_n = c1_b_n;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  c1_d_n = c1_c_n;
  c1_d_ix0 = c1_c_ix0;
  c1_d_iy0 = c1_c_iy0;
  c1_e_n = c1_d_n;
  c1_e_ix0 = c1_d_ix0;
  c1_e_iy0 = c1_d_iy0;
  c1_d = 0.0;
  if (c1_e_n < 1) {
  } else {
    c1_ix = c1_e_ix0;
    c1_iy = c1_e_iy0;
    c1_f_n = c1_e_n;
    c1_b = c1_f_n;
    c1_b_b = c1_b;
    if (1 > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = 1; c1_k <= c1_f_n; c1_k++) {
      c1_d += c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_ix), 1, 36, 1, 0) - 1] * c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_iy), 1, 36, 1, 0) - 1];
      c1_a = c1_ix + 1;
      c1_ix = c1_a;
      c1_b_a = c1_iy + 1;
      c1_iy = c1_b_a;
    }
  }

  return c1_d;
}

static void c1_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0, real_T
  c1_b_y[36])
{
  int32_T c1_i69;
  for (c1_i69 = 0; c1_i69 < 36; c1_i69++) {
    c1_b_y[c1_i69] = c1_y[c1_i69];
  }

  c1_d_eml_xaxpy(chartInstance, c1_n, c1_a, c1_ix0, c1_b_y, c1_iy0);
}

static void c1_c_threshold(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_b_eml_xnrm2(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_x[6], int32_T c1_ix0)
{
  real_T c1_y;
  int32_T c1_b_n;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  int32_T c1_c_ix0;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_scale;
  int32_T c1_kstart;
  int32_T c1_a;
  int32_T c1_c;
  int32_T c1_b_a;
  int32_T c1_b_c;
  int32_T c1_c_a;
  int32_T c1_b;
  int32_T c1_kend;
  int32_T c1_b_kstart;
  int32_T c1_b_kend;
  int32_T c1_d_a;
  int32_T c1_b_b;
  int32_T c1_e_a;
  int32_T c1_c_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_absxk;
  real_T c1_t;
  c1_b_n = c1_n;
  c1_b_ix0 = c1_ix0;
  c1_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_ix0 = c1_b_ix0;
  c1_y = 0.0;
  if (c1_c_n < 1) {
  } else if (c1_c_n == 1) {
    c1_b_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_c_ix0), 1, 6, 1, 0) - 1];
    c1_c_x = c1_b_x;
    c1_y = muDoubleScalarAbs(c1_c_x);
  } else {
    c1_realmin(chartInstance);
    c1_scale = 2.2250738585072014E-308;
    c1_kstart = c1_c_ix0;
    c1_a = c1_c_n;
    c1_c = c1_a;
    c1_b_a = c1_c - 1;
    c1_b_c = c1_b_a;
    c1_c_a = c1_kstart;
    c1_b = c1_b_c;
    c1_kend = c1_c_a + c1_b;
    c1_b_kstart = c1_kstart;
    c1_b_kend = c1_kend;
    c1_d_a = c1_b_kstart;
    c1_b_b = c1_b_kend;
    c1_e_a = c1_d_a;
    c1_c_b = c1_b_b;
    if (c1_e_a > c1_c_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_c_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = c1_b_kstart; c1_k <= c1_b_kend; c1_k++) {
      c1_b_k = c1_k;
      c1_d_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
      c1_e_x = c1_d_x;
      c1_absxk = muDoubleScalarAbs(c1_e_x);
      if (c1_absxk > c1_scale) {
        c1_t = c1_scale / c1_absxk;
        c1_y = 1.0 + c1_y * c1_t * c1_t;
        c1_scale = c1_absxk;
      } else {
        c1_t = c1_absxk / c1_scale;
        c1_y += c1_t * c1_t;
      }
    }

    c1_y = c1_scale * muDoubleScalarSqrt(c1_y);
  }

  return c1_y;
}

static void c1_b_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_b_x[6])
{
  int32_T c1_i70;
  for (c1_i70 = 0; c1_i70 < 6; c1_i70++) {
    c1_b_x[c1_i70] = c1_x[c1_i70];
  }

  c1_e_eml_xscal(chartInstance, c1_n, c1_a, c1_b_x, c1_ix0);
}

static void c1_b_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[6], int32_T
  c1_iy0, real_T c1_b_y[6])
{
  int32_T c1_i71;
  int32_T c1_i72;
  real_T c1_b_x[36];
  for (c1_i71 = 0; c1_i71 < 6; c1_i71++) {
    c1_b_y[c1_i71] = c1_y[c1_i71];
  }

  for (c1_i72 = 0; c1_i72 < 36; c1_i72++) {
    c1_b_x[c1_i72] = c1_x[c1_i72];
  }

  c1_e_eml_xaxpy(chartInstance, c1_n, c1_a, c1_b_x, c1_ix0, c1_b_y, c1_iy0);
}

static void c1_c_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_y[36], int32_T
  c1_iy0, real_T c1_b_y[36])
{
  int32_T c1_i73;
  int32_T c1_i74;
  real_T c1_b_x[6];
  for (c1_i73 = 0; c1_i73 < 36; c1_i73++) {
    c1_b_y[c1_i73] = c1_y[c1_i73];
  }

  for (c1_i74 = 0; c1_i74 < 6; c1_i74++) {
    c1_b_x[c1_i74] = c1_x[c1_i74];
  }

  c1_f_eml_xaxpy(chartInstance, c1_n, c1_a, c1_b_x, c1_ix0, c1_b_y, c1_iy0);
}

static void c1_c_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_b_x[36])
{
  int32_T c1_i75;
  for (c1_i75 = 0; c1_i75 < 36; c1_i75++) {
    c1_b_x[c1_i75] = c1_x[c1_i75];
  }

  c1_f_eml_xscal(chartInstance, c1_a, c1_b_x, c1_ix0);
}

static void c1_eps(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  int32_T c1_i76;
  static char_T c1_cv3[30] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 's', 'v', 'd', '_', 'N', 'o', 'C', 'o', 'n', 'v', 'e', 'r',
    'g', 'e', 'n', 'c', 'e' };

  char_T c1_u[30];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i76 = 0; c1_i76 < 30; c1_i76++) {
    c1_u[c1_i76] = c1_cv3[c1_i76];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static real_T c1_sqrt(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x)
{
  real_T c1_b_x;
  c1_b_x = c1_x;
  c1_b_sqrt(chartInstance, &c1_b_x);
  return c1_b_x;
}

static void c1_c_eml_error(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  int32_T c1_i77;
  static char_T c1_cv4[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c1_u[30];
  const mxArray *c1_y = NULL;
  int32_T c1_i78;
  static char_T c1_cv5[4] = { 's', 'q', 'r', 't' };

  char_T c1_b_u[4];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  for (c1_i77 = 0; c1_i77 < 30; c1_i77++) {
    c1_u[c1_i77] = c1_cv4[c1_i77];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c1_i78 = 0; c1_i78 < 4; c1_i78++) {
    c1_b_u[c1_i78] = c1_cv5[c1_i78];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c1_y, 14, c1_b_y));
}

static void c1_eml_xrotg(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_a,
  real_T c1_b, real_T *c1_b_a, real_T *c1_b_b, real_T *c1_c, real_T *c1_s)
{
  *c1_b_a = c1_a;
  *c1_b_b = c1_b;
  c1_b_eml_xrotg(chartInstance, c1_b_a, c1_b_b, c1_c, c1_s);
}

static void c1_eml_xrot(SFc1_UR5ModelInstanceStruct *chartInstance, real_T c1_x
  [36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_c, real_T c1_s, real_T c1_b_x
  [36])
{
  int32_T c1_i79;
  for (c1_i79 = 0; c1_i79 < 36; c1_i79++) {
    c1_b_x[c1_i79] = c1_x[c1_i79];
  }

  c1_b_eml_xrot(chartInstance, c1_b_x, c1_ix0, c1_iy0, c1_c, c1_s);
}

static void c1_eml_xswap(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_b_x[36])
{
  int32_T c1_i80;
  for (c1_i80 = 0; c1_i80 < 36; c1_i80++) {
    c1_b_x[c1_i80] = c1_x[c1_i80];
  }

  c1_b_eml_xswap(chartInstance, c1_b_x, c1_ix0, c1_iy0);
}

static void c1_d_threshold(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_xgemm(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_k, real_T c1_A[36], real_T c1_B[36], real_T c1_C[36], real_T c1_b_C[36])
{
  int32_T c1_i81;
  int32_T c1_i82;
  real_T c1_b_A[36];
  int32_T c1_i83;
  real_T c1_b_B[36];
  for (c1_i81 = 0; c1_i81 < 36; c1_i81++) {
    c1_b_C[c1_i81] = c1_C[c1_i81];
  }

  for (c1_i82 = 0; c1_i82 < 36; c1_i82++) {
    c1_b_A[c1_i82] = c1_A[c1_i82];
  }

  for (c1_i83 = 0; c1_i83 < 36; c1_i83++) {
    c1_b_B[c1_i83] = c1_B[c1_i83];
  }

  c1_b_eml_xgemm(chartInstance, c1_k, c1_b_A, c1_b_B, c1_b_C);
}

static void c1_e_threshold(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_c_eml_scalar_eg(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_e_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i84;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i84, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i84;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_UR5ModelInstanceStruct *chartInstance;
  chartInstance = (SFc1_UR5ModelInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_f_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_UR5Model, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_UR5Model), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_UR5Model);
  return c1_y;
}

static uint8_T c1_g_emlrt_marshallIn(SFc1_UR5ModelInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0)
{
  int32_T c1_b_n;
  real_T c1_b_a;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_d_ix0;
  int32_T c1_d_a;
  int32_T c1_c;
  int32_T c1_b;
  int32_T c1_b_c;
  int32_T c1_e_a;
  int32_T c1_b_b;
  int32_T c1_i85;
  int32_T c1_f_a;
  int32_T c1_c_b;
  int32_T c1_g_a;
  int32_T c1_d_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  c1_b_n = c1_n;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_d_ix0 = c1_c_ix0;
  c1_d_a = c1_c_n;
  c1_c = c1_d_a;
  c1_b = c1_c - 1;
  c1_b_c = c1_b;
  c1_e_a = c1_c_ix0;
  c1_b_b = c1_b_c;
  c1_i85 = c1_e_a + c1_b_b;
  c1_f_a = c1_d_ix0;
  c1_c_b = c1_i85;
  c1_g_a = c1_f_a;
  c1_d_b = c1_c_b;
  if (c1_g_a > c1_d_b) {
    c1_overflow = false;
  } else {
    c1_eml_switch_helper(chartInstance);
    c1_overflow = (c1_d_b > 2147483646);
  }

  if (c1_overflow) {
    c1_check_forloop_overflow_error(chartInstance, c1_overflow);
  }

  for (c1_k = c1_d_ix0; c1_k <= c1_i85; c1_k++) {
    c1_b_k = c1_k;
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_k), 1, 36, 1, 0) - 1] = c1_c_a * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 36, 1, 0) - 1];
  }
}

static void c1_d_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, int32_T c1_ix0, real_T c1_y[36], int32_T c1_iy0)
{
  int32_T c1_b_n;
  real_T c1_b_a;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_n;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_d_a;
  int32_T c1_ix;
  int32_T c1_e_a;
  int32_T c1_iy;
  int32_T c1_f_a;
  int32_T c1_i86;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_g_a;
  int32_T c1_c;
  int32_T c1_h_a;
  int32_T c1_b_c;
  int32_T c1_i_a;
  int32_T c1_c_c;
  int32_T c1_j_a;
  int32_T c1_k_a;
  c1_b_n = c1_n;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_c_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  if (c1_c_n < 1) {
  } else if (c1_c_a == 0.0) {
  } else {
    c1_d_a = c1_c_ix0 - 1;
    c1_ix = c1_d_a;
    c1_e_a = c1_c_iy0 - 1;
    c1_iy = c1_e_a;
    c1_f_a = c1_c_n - 1;
    c1_i86 = c1_f_a;
    c1_b = c1_i86;
    c1_b_b = c1_b;
    if (0 > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = 0; c1_k <= c1_i86; c1_k++) {
      c1_g_a = c1_iy;
      c1_c = c1_g_a;
      c1_h_a = c1_iy;
      c1_b_c = c1_h_a;
      c1_i_a = c1_ix;
      c1_c_c = c1_i_a;
      c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_c + 1)), 1, 36, 1, 0) - 1] =
        c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_b_c + 1)), 1, 36, 1, 0) - 1] + c1_c_a *
        c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_c_c + 1)), 1, 36, 1, 0) - 1];
      c1_j_a = c1_ix + 1;
      c1_ix = c1_j_a;
      c1_k_a = c1_iy + 1;
      c1_iy = c1_k_a;
    }
  }
}

static void c1_e_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0)
{
  int32_T c1_b_n;
  real_T c1_b_a;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_d_ix0;
  int32_T c1_d_a;
  int32_T c1_c;
  int32_T c1_b;
  int32_T c1_b_c;
  int32_T c1_e_a;
  int32_T c1_b_b;
  int32_T c1_i87;
  int32_T c1_f_a;
  int32_T c1_c_b;
  int32_T c1_g_a;
  int32_T c1_d_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  c1_b_n = c1_n;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_d_ix0 = c1_c_ix0;
  c1_d_a = c1_c_n;
  c1_c = c1_d_a;
  c1_b = c1_c - 1;
  c1_b_c = c1_b;
  c1_e_a = c1_c_ix0;
  c1_b_b = c1_b_c;
  c1_i87 = c1_e_a + c1_b_b;
  c1_f_a = c1_d_ix0;
  c1_c_b = c1_i87;
  c1_g_a = c1_f_a;
  c1_d_b = c1_c_b;
  if (c1_g_a > c1_d_b) {
    c1_overflow = false;
  } else {
    c1_eml_switch_helper(chartInstance);
    c1_overflow = (c1_d_b > 2147483646);
  }

  if (c1_overflow) {
    c1_check_forloop_overflow_error(chartInstance, c1_overflow);
  }

  for (c1_k = c1_d_ix0; c1_k <= c1_i87; c1_k++) {
    c1_b_k = c1_k;
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_k), 1, 6, 1, 0) - 1] = c1_c_a * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 6, 1, 0) - 1];
  }
}

static void c1_e_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[36], int32_T c1_ix0, real_T c1_y[6], int32_T
  c1_iy0)
{
  int32_T c1_b_n;
  real_T c1_b_a;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_n;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_d_a;
  int32_T c1_ix;
  int32_T c1_e_a;
  int32_T c1_iy;
  int32_T c1_f_a;
  int32_T c1_i88;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_g_a;
  int32_T c1_c;
  int32_T c1_h_a;
  int32_T c1_b_c;
  int32_T c1_i_a;
  int32_T c1_c_c;
  int32_T c1_j_a;
  int32_T c1_k_a;
  c1_b_n = c1_n;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_c_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  if (c1_c_n < 1) {
  } else if (c1_c_a == 0.0) {
  } else {
    c1_d_a = c1_c_ix0 - 1;
    c1_ix = c1_d_a;
    c1_e_a = c1_c_iy0 - 1;
    c1_iy = c1_e_a;
    c1_f_a = c1_c_n - 1;
    c1_i88 = c1_f_a;
    c1_b = c1_i88;
    c1_b_b = c1_b;
    if (0 > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = 0; c1_k <= c1_i88; c1_k++) {
      c1_g_a = c1_iy;
      c1_c = c1_g_a;
      c1_h_a = c1_iy;
      c1_b_c = c1_h_a;
      c1_i_a = c1_ix;
      c1_c_c = c1_i_a;
      c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_c + 1)), 1, 6, 1, 0) - 1] = c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)(c1_b_c + 1)), 1, 6, 1, 0)
        - 1] + c1_c_a * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)(c1_c_c + 1)), 1, 36, 1, 0) - 1];
      c1_j_a = c1_ix + 1;
      c1_ix = c1_j_a;
      c1_k_a = c1_iy + 1;
      c1_iy = c1_k_a;
    }
  }
}

static void c1_f_eml_xaxpy(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_n, real_T c1_a, real_T c1_x[6], int32_T c1_ix0, real_T c1_y[36], int32_T
  c1_iy0)
{
  int32_T c1_b_n;
  real_T c1_b_a;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_n;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_d_a;
  int32_T c1_ix;
  int32_T c1_e_a;
  int32_T c1_iy;
  int32_T c1_f_a;
  int32_T c1_i89;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_g_a;
  int32_T c1_c;
  int32_T c1_h_a;
  int32_T c1_b_c;
  int32_T c1_i_a;
  int32_T c1_c_c;
  int32_T c1_j_a;
  int32_T c1_k_a;
  c1_b_n = c1_n;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_c_threshold(chartInstance);
  c1_c_n = c1_b_n;
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  if (c1_c_n < 1) {
  } else if (c1_c_a == 0.0) {
  } else {
    c1_d_a = c1_c_ix0 - 1;
    c1_ix = c1_d_a;
    c1_e_a = c1_c_iy0 - 1;
    c1_iy = c1_e_a;
    c1_f_a = c1_c_n - 1;
    c1_i89 = c1_f_a;
    c1_b = c1_i89;
    c1_b_b = c1_b;
    if (0 > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_k = 0; c1_k <= c1_i89; c1_k++) {
      c1_g_a = c1_iy;
      c1_c = c1_g_a;
      c1_h_a = c1_iy;
      c1_b_c = c1_h_a;
      c1_i_a = c1_ix;
      c1_c_c = c1_i_a;
      c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_c + 1)), 1, 36, 1, 0) - 1] =
        c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_b_c + 1)), 1, 36, 1, 0) - 1] + c1_c_a *
        c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)(c1_c_c + 1)), 1, 6, 1, 0) - 1];
      c1_j_a = c1_ix + 1;
      c1_ix = c1_j_a;
      c1_k_a = c1_iy + 1;
      c1_iy = c1_k_a;
    }
  }
}

static void c1_f_eml_xscal(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_a, real_T c1_x[36], int32_T c1_ix0)
{
  real_T c1_b_a;
  int32_T c1_b_ix0;
  real_T c1_c_a;
  int32_T c1_c_ix0;
  int32_T c1_d_ix0;
  int32_T c1_d_a;
  int32_T c1_i90;
  int32_T c1_e_a;
  int32_T c1_b;
  int32_T c1_f_a;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  c1_b_a = c1_a;
  c1_b_ix0 = c1_ix0;
  c1_b_threshold(chartInstance);
  c1_c_a = c1_b_a;
  c1_c_ix0 = c1_b_ix0;
  c1_d_ix0 = c1_c_ix0;
  c1_d_a = c1_c_ix0 + 5;
  c1_i90 = c1_d_a;
  c1_e_a = c1_d_ix0;
  c1_b = c1_i90;
  c1_f_a = c1_e_a;
  c1_b_b = c1_b;
  if (c1_f_a > c1_b_b) {
    c1_overflow = false;
  } else {
    c1_eml_switch_helper(chartInstance);
    c1_overflow = (c1_b_b > 2147483646);
  }

  if (c1_overflow) {
    c1_check_forloop_overflow_error(chartInstance, c1_overflow);
  }

  for (c1_k = c1_d_ix0; c1_k <= c1_i90; c1_k++) {
    c1_b_k = c1_k;
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_b_k), 1, 36, 1, 0) - 1] = c1_c_a * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 36, 1, 0) - 1];
  }
}

static void c1_b_sqrt(SFc1_UR5ModelInstanceStruct *chartInstance, real_T *c1_x)
{
  if (*c1_x < 0.0) {
    c1_c_eml_error(chartInstance);
  }

  *c1_x = muDoubleScalarSqrt(*c1_x);
}

static void c1_b_eml_xrotg(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  *c1_a, real_T *c1_b, real_T *c1_c, real_T *c1_s)
{
  real_T c1_b_a;
  real_T c1_b_b;
  real_T c1_c_b;
  real_T c1_c_a;
  real_T c1_d_a;
  real_T c1_d_b;
  real_T c1_e_b;
  real_T c1_e_a;
  real_T c1_b_c;
  real_T c1_b_s;
  double * c1_a_t;
  double * c1_b_t;
  double * c1_c_t;
  double * c1_s_t;
  real_T c1_c_c;
  real_T c1_c_s;
  (void)chartInstance;
  c1_b_a = *c1_a;
  c1_b_b = *c1_b;
  c1_c_b = c1_b_b;
  c1_c_a = c1_b_a;
  c1_d_a = c1_c_a;
  c1_d_b = c1_c_b;
  c1_e_b = c1_d_b;
  c1_e_a = c1_d_a;
  c1_b_c = 0.0;
  c1_b_s = 0.0;
  c1_a_t = (double *)(&c1_e_a);
  c1_b_t = (double *)(&c1_e_b);
  c1_c_t = (double *)(&c1_b_c);
  c1_s_t = (double *)(&c1_b_s);
  drotg(c1_a_t, c1_b_t, c1_c_t, c1_s_t);
  c1_c_a = c1_e_a;
  c1_c_b = c1_e_b;
  c1_c_c = c1_b_c;
  c1_c_s = c1_b_s;
  *c1_a = c1_c_a;
  *c1_b = c1_c_b;
  *c1_c = c1_c_c;
  *c1_s = c1_c_s;
}

static void c1_b_eml_xrot(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0, real_T c1_c, real_T c1_s)
{
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  real_T c1_b_c;
  real_T c1_b_s;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  real_T c1_c_c;
  real_T c1_c_s;
  int32_T c1_ix;
  int32_T c1_iy;
  int32_T c1_k;
  real_T c1_temp;
  int32_T c1_a;
  int32_T c1_b_a;
  (void)chartInstance;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_b_c = c1_c;
  c1_b_s = c1_s;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  c1_c_c = c1_b_c;
  c1_c_s = c1_b_s;
  c1_ix = c1_c_ix0;
  c1_iy = c1_c_iy0;
  for (c1_k = 1; c1_k < 7; c1_k++) {
    c1_temp = c1_c_c * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_ix), 1, 36, 1, 0) - 1] + c1_c_s *
      c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_iy), 1, 36, 1, 0) - 1];
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_iy), 1, 36, 1, 0) - 1] = c1_c_c * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("",
      (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_iy), 1, 36, 1, 0) - 1] - c1_c_s
      * c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_ix), 1, 36, 1, 0) - 1];
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_ix), 1, 36, 1, 0) - 1] = c1_temp;
    c1_a = c1_iy + 1;
    c1_iy = c1_a;
    c1_b_a = c1_ix + 1;
    c1_ix = c1_b_a;
  }
}

static void c1_b_eml_xswap(SFc1_UR5ModelInstanceStruct *chartInstance, real_T
  c1_x[36], int32_T c1_ix0, int32_T c1_iy0)
{
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_ix;
  int32_T c1_iy;
  int32_T c1_k;
  real_T c1_temp;
  int32_T c1_a;
  int32_T c1_b_a;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_d_threshold(chartInstance);
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  c1_ix = c1_c_ix0;
  c1_iy = c1_c_iy0;
  for (c1_k = 1; c1_k < 7; c1_k++) {
    c1_temp = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c1_ix), 1, 36, 1, 0) - 1];
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_ix), 1, 36, 1, 0) - 1] = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c1_iy), 1, 36, 1, 0) - 1];
    c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c1_iy), 1, 36, 1, 0) - 1] = c1_temp;
    c1_a = c1_ix + 1;
    c1_ix = c1_a;
    c1_b_a = c1_iy + 1;
    c1_iy = c1_b_a;
  }
}

static void c1_b_eml_xgemm(SFc1_UR5ModelInstanceStruct *chartInstance, int32_T
  c1_k, real_T c1_A[36], real_T c1_B[36], real_T c1_C[36])
{
  int32_T c1_b_k;
  int32_T c1_c_k;
  int32_T c1_a;
  int32_T c1_km1;
  int32_T c1_cr;
  int32_T c1_b_cr;
  int32_T c1_b_a;
  int32_T c1_i91;
  int32_T c1_c_a;
  int32_T c1_i92;
  int32_T c1_d_a;
  int32_T c1_b;
  int32_T c1_e_a;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_ic;
  int32_T c1_b_ic;
  int32_T c1_br;
  int32_T c1_c_cr;
  int32_T c1_ar;
  int32_T c1_f_a;
  int32_T c1_b_br;
  int32_T c1_c_b;
  int32_T c1_c;
  int32_T c1_g_a;
  int32_T c1_d_b;
  int32_T c1_i93;
  int32_T c1_h_a;
  int32_T c1_e_b;
  int32_T c1_i_a;
  int32_T c1_f_b;
  boolean_T c1_b_overflow;
  int32_T c1_ib;
  int32_T c1_b_ib;
  real_T c1_temp;
  int32_T c1_ia;
  int32_T c1_j_a;
  int32_T c1_i94;
  int32_T c1_k_a;
  int32_T c1_i95;
  int32_T c1_l_a;
  int32_T c1_g_b;
  int32_T c1_m_a;
  int32_T c1_h_b;
  boolean_T c1_c_overflow;
  int32_T c1_c_ic;
  int32_T c1_n_a;
  int32_T c1_o_a;
  c1_b_k = c1_k;
  c1_e_threshold(chartInstance);
  c1_c_k = c1_b_k;
  c1_a = c1_c_k;
  c1_km1 = c1_a;
  for (c1_cr = 0; c1_cr < 31; c1_cr += 6) {
    c1_b_cr = c1_cr;
    c1_b_a = c1_b_cr + 1;
    c1_i91 = c1_b_a;
    c1_c_a = c1_b_cr + 6;
    c1_i92 = c1_c_a;
    c1_d_a = c1_i91;
    c1_b = c1_i92;
    c1_e_a = c1_d_a;
    c1_b_b = c1_b;
    if (c1_e_a > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_ic = c1_i91; c1_ic <= c1_i92; c1_ic++) {
      c1_b_ic = c1_ic;
      c1_C[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_b_ic), 1, 36, 1, 0) - 1] = 0.0;
    }
  }

  c1_br = 0;
  for (c1_c_cr = 0; c1_c_cr < 31; c1_c_cr += 6) {
    c1_b_cr = c1_c_cr;
    c1_ar = 0;
    c1_f_a = c1_br + 1;
    c1_br = c1_f_a;
    c1_b_br = c1_br;
    c1_c_b = c1_km1 - 1;
    c1_c = 6 * c1_c_b;
    c1_g_a = c1_br;
    c1_d_b = c1_c;
    c1_i93 = c1_g_a + c1_d_b;
    c1_h_a = c1_b_br;
    c1_e_b = c1_i93;
    c1_i_a = c1_h_a;
    c1_f_b = c1_e_b;
    if (c1_i_a > c1_f_b) {
      c1_b_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_b_overflow = (c1_f_b > 2147483641);
    }

    if (c1_b_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
    }

    for (c1_ib = c1_b_br; c1_ib <= c1_i93; c1_ib += 6) {
      c1_b_ib = c1_ib;
      if (c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ib), 1, 36, 1, 0) - 1] != 0.0) {
        c1_temp = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_ib), 1, 36, 1, 0) - 1];
        c1_ia = c1_ar;
        c1_j_a = c1_b_cr + 1;
        c1_i94 = c1_j_a;
        c1_k_a = c1_b_cr + 6;
        c1_i95 = c1_k_a;
        c1_l_a = c1_i94;
        c1_g_b = c1_i95;
        c1_m_a = c1_l_a;
        c1_h_b = c1_g_b;
        if (c1_m_a > c1_h_b) {
          c1_c_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_c_overflow = (c1_h_b > 2147483646);
        }

        if (c1_c_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_c_overflow);
        }

        for (c1_c_ic = c1_i94; c1_c_ic <= c1_i95; c1_c_ic++) {
          c1_b_ic = c1_c_ic;
          c1_n_a = c1_ia + 1;
          c1_ia = c1_n_a;
          c1_C[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ic), 1, 36, 1, 0) - 1] =
            c1_C[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ic), 1, 36, 1, 0) - 1] + c1_temp *
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ia), 1, 36, 1, 0) - 1];
        }
      }

      c1_o_a = c1_ar + 6;
      c1_ar = c1_o_a;
    }
  }
}

static void init_dsm_address_info(SFc1_UR5ModelInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c1_UR5Model_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(739505933U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1565683081U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(901706100U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3434804950U);
}

mxArray *sf_c1_UR5Model_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("1QfYD2jDd0RrpSY579S9kE");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_UR5Model_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_UR5Model_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c1_UR5Model(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"dq\",},{M[8],M[0],T\"is_active_c1_UR5Model\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_UR5Model_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_UR5ModelInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_UR5ModelInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _UR5ModelMachineNumber_,
           1,
           1,
           1,
           0,
           4,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_UR5ModelMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_UR5ModelMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _UR5ModelMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"x_des");
          _SFD_SET_DATA_PROPS(1,2,0,1,"dq");
          _SFD_SET_DATA_PROPS(2,1,1,0,"xd_des");
          _SFD_SET_DATA_PROPS(3,1,1,0,"q");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,6318);

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c1_x_des)[6];
          real_T (*c1_dq)[6];
          real_T (*c1_xd_des)[6];
          real_T (*c1_q)[6];
          c1_q = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
          c1_xd_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
          c1_dq = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
          c1_x_des = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c1_x_des);
          _SFD_SET_DATA_VALUE_PTR(1U, *c1_dq);
          _SFD_SET_DATA_VALUE_PTR(2U, *c1_xd_des);
          _SFD_SET_DATA_VALUE_PTR(3U, *c1_q);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _UR5ModelMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "rfp0r80RldCO1TerDxtSf";
}

static void sf_opaque_initialize_c1_UR5Model(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_UR5ModelInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
  initialize_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_UR5Model(void *chartInstanceVar)
{
  enable_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_UR5Model(void *chartInstanceVar)
{
  disable_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_UR5Model(void *chartInstanceVar)
{
  sf_gateway_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_UR5Model(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_UR5Model((SFc1_UR5ModelInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_UR5Model();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c1_UR5Model(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c1_UR5Model();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_UR5Model((SFc1_UR5ModelInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_UR5Model(SimStruct* S)
{
  return sf_internal_get_sim_state_c1_UR5Model(S);
}

static void sf_opaque_set_sim_state_c1_UR5Model(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c1_UR5Model(S, st);
}

static void sf_opaque_terminate_c1_UR5Model(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_UR5ModelInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_UR5Model_optimization_info();
    }

    finalize_c1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_UR5Model((SFc1_UR5ModelInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_UR5Model(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c1_UR5Model((SFc1_UR5ModelInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_UR5Model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_UR5Model_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2333141629U));
  ssSetChecksum1(S,(998330753U));
  ssSetChecksum2(S,(3962719223U));
  ssSetChecksum3(S,(580298832U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_UR5Model(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_UR5Model(SimStruct *S)
{
  SFc1_UR5ModelInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_UR5ModelInstanceStruct *)utMalloc(sizeof
    (SFc1_UR5ModelInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_UR5ModelInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_UR5Model;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_UR5Model;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_UR5Model;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_UR5Model;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_UR5Model;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_UR5Model;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_UR5Model;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_UR5Model;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_UR5Model;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_UR5Model;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_UR5Model;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c1_UR5Model_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_UR5Model(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_UR5Model(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_UR5Model(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_UR5Model_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
