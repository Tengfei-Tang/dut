/*******************************************************************************
 * TRANS_NEUT case
 *
 * 
 * Solves equations for 
 *  plasmas density                  Ni
 *  parallel ion velocity            Vi
 *  electron and ion temperatures    Te, Ti
 *  vorticity and potential          U00, phi01
 * and also equations of neutrals:
 *  atom density                     Nn
 *  perpendicular velocity           Vn 
 *  molecule density                 Nm
 *  molecule radial velocity         Vmx
 * 
 * Intended to be run for NZ=1 (i.e. X and Y only) at first,
 * 3D SMBI has also been tested no problem for a short time injection
 *
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>

#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <interpolation.hxx>
#include <bout/globalfield.hxx>
//#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <math.h>

// Vectors
Vector2D b0xcv;    // Curvature term
Field2D  J0;       // Current
Vector3D Vm;
Vector3D grad_perp_Diffi,grad_perp_Diffn,grad_perp_chii,grad_perp_chie;
Field3D grad_para_kappTe,grad_para_kappTi;

//2D Evolving fields from Grid file
Field2D Te_grid,Ti_grid,Ni_grid,U00_grid,phi01_grid;
Field2D Te_exp,Ti_exp,Ne_exp,Ni_exp,U00_exp,phi01_exp;

//2D Psi from Grid file
Field2D psixy;
// 3D evolving fields
Field3D  Te, Ne,Ni, Vi, Ve, Ti, Nn,Tn, Vn, Nm ,Vmx;             // Vi and Vn are paralell velocity
Field3D temp_Ni,temp_Nn,temp_Nm,temp_Ti,temp_Te;                // temporary variables for minivalue protection

// 3D initial profiles Hmode
Field3D Ni0_Hmode,Ti0_Hmode,Te0_Hmode;

// Non-linear coefficientsx
Field3D kappa_Te, kappa_Ti;
Field3D kappa_Te_fl,kappa_Ti_fl;                                // heat flux limit coefficients 

// Thermal Ion and Electron Speed
Field3D V_th_i,V_th_e,V_th_n;

//ionization/CX/rates for Plasmas and Neutrals
Field3D nu_ionz,nu_CX,nu_ionz_n,nu_CX_n,nu_diss,nu_rec;

// Collisional time between electron and ions
Field3D tau_ei;

//ion viscosity
Field3D eta0_i,eta0_n;
Field2D Diff_grid,Chi_grid,Che_grid;
// Density/heat diffusion coefficients of plasmas
Field3D Diff_ni_perp,chi_i_perp,chi_e_perp;
Field3D Diffc_ni_step,chic_i_step,chic_e_step;
Field3D Diffc_ni_Hmode,chic_i_Hmode,chic_e_Hmode;

// average over Y direction for Diffscoefs_Hmode 
Field2D aveY_J,aveY_g11;
Field3D aveY_g11J_ni,aveY_g11J_ti,aveY_g11J_te;
Field3D Gamma_ni_xin_exp,Qi_xin_exp,Qe_xin_exp;

// Classical Diffusion coefficients of neutrals 
Field3D Diffc_nn_perp,Diffc_nn_par,Diffc_nm_perp;
Field3D Diffc_nn_perp_fl;  //flux limited
// Source Terms
Field3D Si_p,Si_CX,S_diss,S_rec;
Field3D S_pi_ext,S_Ee_ext,S_Ei_ext;      // external sources

// 3D total values
Field3D Nit, Tit, Tet, Vit;

// pressures
Field3D  Pe,Pe0,Pi,Pi0,pei,pn,pm;
BoutReal delta;

// parallel pressure gradient
Field3D Grad_par_pei,Grad_par_pn;

// gradient length related quatities
Field3D Grad_par_logNn,Grad_par_Tn;

//pressure gradient length
Field3D dpei_dx,Heavistep;

//radial derivatives for Hmode diffusion coef calculation
Field3D DDX_Ni,DDX_Te,DDX_Ti;
Field3D D2DX2_Ni, D2DX2_Ti, D2DX2_Te;

// Variations for Sheath Boundary condition 
Field3D c_se,Vn_pol,Vn_perp,q_si_kappa,q_se_kappa;
Field3D Gamma_ni,Gamma_nn,Gamma_ni_Diffc,Gamma_nn_Diffc;
Field3D Gamma_ni_wall,Gamma_nn_wall,Gamma_nnw_Diffc,Grad_par_Nn;

// spitzer resistivty
Field3D eta_spitzer;

// Bootstrap current 
Field3D Jpar_BS0;
Field3D L31, L32, L34;
Field3D f31, f32ee, f32ei, f34, ft;
Field3D nu_estar, nu_istar, BSal0, BSal;
BoutReal Aratio;
BoutReal test_f1,test_f2,test_f3;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe,B0;
Field2D I;

// B field vectors
Vector2D B0vec; // B0 field vector

// Max grid number in X Y Z
int NX,NY,NZ;

// width of SOL in Y grid
int jyseps11;

// filter 
BoutReal filter_para;

//Prameter for heat flux limit
Field2D q95;
BoutReal q95_input;
BoutReal q_alpha; 
BoutReal R_major, r_minor,P_input,q_input,Elgt,chi_coef;
BoutReal leak_ni,leak_ti,leak_te,decay_bndry_legth;
bool energy_flux_bndry,leakage_bndry,decay_bndry;
// Initial profile parameters
BoutReal Te_core,Te_edge,Ti_core,Ti_edge,Ni_core,Ni_edge;
BoutReal Initfile_x0,Initfile_w_ped;
BoutReal x0_ped,width_ped,coef_coregrad_ped,coef_corenlarge_ped;
BoutReal density_unit;                   // Number density [m^-3]

// parameters
BoutReal Te_x, Ti_x,Tm_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal Zi;                             // charge number of ion same as ZZ
BoutReal minimum_val;
BoutReal Mm,Mi,Mn;
BoutReal W_ionz,W_diss,W_bind,W_rec;
BoutReal Lbar,tbar,Lp,Lp_crit,Lnn_min;
BoutReal lambda_ei, lambda_ii;
BoutReal N_cubic_Lbar;                   // particle numbers in cubic of Lbar^3
BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
Field3D nu_ei,nu_ii;
BoutReal Nm0,Vm0;

BoutReal Diffc_ni_perp,Difft_ni_perp,chic_i_perp,chit_i_perp;
BoutReal chic_e_perp,chit_e_perp,chic_n_perp;
int t_output;

//Boundaries of 3D-const fueling flux
BoutReal CF_BC_x0,CF_BC_y0,CF_BC_y1,CF_BC_z0,CF_BC_z1;  // SMBI_LFS
BoutReal Sheath_BC_x0;
BoutReal Rate_recycle_wall,Rate_recycle_div,alpha_vn,angle_B_plate,Tn_plate,Vn_th_plate;
BoutReal Lni_wall;                       //Gradient length Ni at wall
int Sheath_width;

//radial global x 
BoutReal x_rela;

// parameters of external sources 
BoutReal x0_extS,width_extS,coef_coregrad_extS;
BoutReal amp_spi_ext,amp_see_ext,amp_sei_ext;

//constant gradient at xin with auto unit
BoutReal dNidx_xin_au,dTidx_xin_au,dTedx_xin_au;
BoutReal psi_xout_y0,psi_axis,psi_bndry;

//step function of diffusion coefficients
BoutReal diffusion_coef_step0, diffusion_coef_step1;
bool  diffusion_coef_step_function;

// diffusion coefficients of Hmode profiles
BoutReal diffusion_coef_Hmode0, diffusion_coef_Hmode1;
bool  diffusion_coef_Hmode;

//Solving Equations options
bool Solving_Eq_Ni, Solving_Eq_Ti, Solving_Eq_Te, Solving_Eq_Vi, Solving_Eq_phi01, Solving_Eq_U00, Solving_Eq_Nn; 

//logical paramter
bool profiles_lasttimestep,load_grid_trans,load_grid_profiles, initial_profile_exp, initial_profile_linear;
bool load_experiment_profiles;
bool initial_profile_Hmode,initial_SOL_edgeval;
bool external_sources,extsrcs_balance_diffusion;
bool noshear, include_curvature,curvature_phi,diamag,energy_flux,Turb_Diff_on;
bool SMBI_LFS, Sheath_BC,Sheath_BC_phi,SBC_particle_recycle,Wall_particle_recycle;
bool Secondary_electron_emission;
bool Ni_Constant_Grad,Ni_Zero_Grad;
bool term_GparkappaTe_GparTe;
bool terms_ionization;
bool terms_recombination,terms_Gradpar_pn,terms_Gradpar_eta0n,terms_cross,terms_exb,include_J0,J1_all,terms_Ve,terms_radial_conv,thermal_force,terms_spitzer_resist;
//***Ion orbit loss**//
bool terms_IOL;
int xout,P_N;
Field3D V_iol,V_iol1;
Field2D M_iol0,F_iol0,E_iol0,M_iol1,F_iol1,E_iol1;
//**IOL**//
bool terms_NnGradpar_Vn,terms_NnGradpar_Vi,terms_Diffcnn_par;
bool terms_Gradperp_diffcoefs;
bool terms_Gradpara_diffcoefs;
bool nlfilter_noisy_data,nlfilter_Gradpar_logNn,nlfilter_Vnpol;
bool BScurrent,spitzer_resist;

const Field3D ret_const_flux_BC(const Field3D &var, const BoutReal value);
void Diag_neg_value (const Field3D &f1,const Field3D &f2,const Field3D &f3, const Field3D &f4);
const Field3D field_larger(const Field3D &f, const BoutReal limit);
const Field3D field_smaller(const Field3D &f, const BoutReal limit);
void SBC_Dirichlet_SWidth1(Field3D &var, const Field3D &value);
void SBC_Dirichlet(Field3D &var, const Field3D &value);
void SBC_Dirichlet_phi(Field3D &var, const Field3D &value);
void SBC_Gradpar(Field3D &var, const Field3D &value);
void SBC_yup_eq(Field3D &var, const Field3D &value);
void SBC_ydown_eq(Field3D &var, const Field3D &value);
void SBC_yup_Grad_par(Field3D &var, const Field3D &value);
void SBC_ydown_Grad_par(Field3D &var, const Field3D &value);
void WallBC_Xout_GradX(Field3D &var, const Field3D &value);
void WallBC_Xout_GradX_len(Field3D &var, BoutReal value);

const Field3D BS_ft(const int index);
const Field3D F31(const Field3D input);
const Field3D F32ee(const Field3D input);
const Field3D F32ei(const Field3D input);

//*******************************************************li2016
BoutReal Bbar;                  // Normalisation constants
Field3D eta;                    // Resistivity profile (1 / S)
BoutReal ixsep;
BoutReal jysep1, jysep2;        //index for x-point on y direction
BoutReal xloc;
int yloc;

Field3D Ti0,N0,Te0;             // number density and temperature
Field3D P0;                     // Pressure
bool iterative_phi0,J_para0;
BoutReal sbc_lambda0, sbc_lambda1, Mu_perp, Mu_para, RF_coef;
Vector3D er00, er0, ersh0,ersh1,Vexb,Vdia,Vdia_e;
Field3D phi0,phi00, phi01, U00, eta_sh, phish0, phish1, dphi01dx, d2phishdx2;
Field3D term_Mu_perp,term_Mu_par,term_J1,term_J0,term_curvature,term_exb,J1_para,J1_phi,J1_Te,J1_Pe,term_gyro;

//******** parallel heat flux *******************//
Field3D heatf_cond_i,heatf_cond_e,heatf_conv_i,heatf_conv_e,heatf_conv_e1,heatf_cond_i_x,heatf_cond_e_x ,heatf_exb_i_x,heatf_exb_e_x,heatf_dia_i_x,heatf_dia_e_x,heatf_eng_i_x,heatf_eng_e_x;
Field3D Ti_grad,Te_grad;
//******* radial particle flux ****************//
Field3D gamma_cond_i_x,gamma_exb_i_x,gamma_dia_i_x;

void Grad2_par2_sh(Field3D &var, const Field3D &value);
Field3D Laplace_par_filterPF (const Field3D &var);
Field3D Laplace_perp_filterPF (const Field3D &var);
Field3D Laplace_perp_PFboundary (const Field3D &var);

bool gyroviscous;
Field3D Dperp2phi0, Dperp2Pi, Gradphi02, bracketphi0P; // Temporary variables for gyroviscous
BoutReal Ugyro;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm_exb, bm_mag; // Bracket method for advection terms
int bracket_method_exb, bracket_method_mag;

//*******************************************************liend

const BoutReal PI = 3.14159265;
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Me = 9.1094e-31;
const BoutReal Mi2 = 2.0*1.6726e-27;  // Ion mass
const BoutReal KB = 1.38065e-23;      // Boltamann constant
const BoutReal eV_K = 11605.0;        // 1eV = 11605K
const BoutReal ee = 1.602e-19;        // ln(Lambda)

BoutReal C_fe_sheat;                  // coefficient of parallel heat transmission in sheat BC 
BoutReal C_fi_sheat;


int physics_init(bool restarting)
{
   
  t_output=1;
  
  output.write("Solving transport equations for Plasmas: Ni, Vi, Ti, Te and Neutrals: Nn Vn \n");

  /////////////// LOAD DATA FROM GRID FILE //////////////

  // Load 2D profiles
  mesh->get(Te_grid, "Te0");         // eV
  mesh->get(Ti_grid, "Ti0");         // eV
  mesh->get(Ni_grid, "Ni0");         
  mesh->get(phi01_grid, "PHI_0");    

  mesh->get(Diff_grid, "diff");      
  mesh->get(Chi_grid, "chi_i");      
  mesh->get(Che_grid, "chi_e");      
  mesh->get(Te_exp, "Teexp");       // eV
  //  Te_exp *=1000;
  mesh->get(Ti_exp, "Tiexp");       // eV
  //  Ti_exp *=1000;
  mesh->get(Ne_exp, "Neexp");       //  10^20/m^3
  //  Ne_exp *=10;
  mesh->get(Ni_exp, "Niexp");       //  10^20/m^3
  //  Ni_exp *=10;

  // Load Psi
  mesh->get(psixy,"psixy");               // unit m^2 T
  mesh->get(psi_axis,"psi_axis");
  mesh->get(psi_bndry,"psi_bndry");
  // Load curvature term
  b0xcv.covariant = false;               // Read contravariant components
  Vm.covariant = false;
  mesh->get(b0xcv, "bxcv");              // mixed units x: T y: m^-2 z: m^-2
  mesh->get(J0, "Jpar0");                // A / m^2

  // Load metrics
  GRID_LOAD(Rxy);                        // Major radius [m]
  GRID_LOAD2(Bpxy, Btxy);                // Poloidal, Toroidal B field [T]
  mesh->get(B0,   "Bxy");                // T
  GRID_LOAD(hthe);                       // Poloidal arc length [m / radian]
  //mesh->get(mesh->dx,   "dpsi");
  //mesh->get(mesh->dx, "dx");           // 1D test only
  mesh->get(I,    "sinty");              // m^-2 T^-1
  mesh->get(NX,"nx");   
  mesh->get(NY,"ny");   

  mesh->get(jyseps11,"jyseps1_1");   

  if(mesh->get(bmag, "bmag")) 
     bmag=1.0;
  if(mesh->get(Lbar, "rmag"))
     Lbar=1.0; 
//***********************************************************li2016
  if(mesh->get(Bbar, "bmag"))            // Typical magnetic field
    Bbar = 1.0;
   mesh->get(ixsep, "ixseps1");
   mesh->get(jysep1, "jyseps1_1");
   mesh->get(jysep2, "jyseps2_2");
   mesh->get(P0, "pressure");            // Pascals
   mesh->get(phi0, "PHI_0");
//*********************************************************liend 

  // Load normalisation values
  //GRID_LOAD(Te_x);
  //GRID_LOAD(Ti_x);
  //GRID_LOAD(Ni_x);
  //GRID_LOAD(bmag);


  /////////////// READ OPTIONS //////////////////////////

  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("trans_neu");
 
  OPTION(options, minimum_val,  1.e-10);    // minimum value limit for densities Ns  

  OPTION(options, NZ,  1);                  // maximum grid number in Z

  OPTION(options, Te_x,  10.0);             // Read in eV 
  OPTION(options, Ti_x,  10.0);             // eV
  OPTION(options, Tm_x,  0.0258);           // eV
  OPTION(options, Ni_x,  1.0);              // in 1.e^20 m^-3
  OPTION(options, density_unit,   1.0e20);  // Number density [m^-3]

  OPTION(options, Lp_crit,  5.e-5);         // in m

  OPTION(options, Solving_Eq_Ni, true);
  OPTION(options, Solving_Eq_Te, true);
  OPTION(options, Solving_Eq_Ti, true);
  OPTION(options, Solving_Eq_Vi, true);
  OPTION(options, Solving_Eq_Nn, true);
  OPTION(options, Solving_Eq_phi01, true);
  OPTION(options, Solving_Eq_U00, true);

  OPTION(options, include_curvature, false);
  OPTION(options, diamag, false);
  OPTION(options, curvature_phi, false);
  OPTION(options, energy_flux, false);
  OPTION(options, noshear,           true);

  OPTION(options, load_grid_trans, false);              // priority I
  OPTION(options, profiles_lasttimestep, false);        // priority I
  OPTION(options, load_grid_profiles, false);           // priority II
  OPTION(options, load_experiment_profiles, false);     // priority III
  OPTION(options, initial_profile_exp, false);          // priority III
  OPTION(options, initial_profile_linear, false);       // priority III
  OPTION(options, initial_profile_Hmode, false);        // priority III
  OPTION(options, initial_SOL_edgeval, false);          // set edge value in SOL region in X-gfile


  OPTION(options, Turb_Diff_on, false);                 // include Turbulent diffusions
  OPTION(options, SMBI_LFS, false); 
  OPTION(options, Sheath_BC, false);
  OPTION(options, Ni_Constant_Grad, false);
  OPTION(options, Ni_Zero_Grad, false);
  OPTION(options, Sheath_BC_phi, false);
  OPTION(options, Sheath_width, 0);
  OPTION(options, SBC_particle_recycle, false);
  OPTION(options, Wall_particle_recycle, false);
  OPTION(options, Rate_recycle_wall, 0.50);
  OPTION(options, Rate_recycle_div, 0.50);
  OPTION(options, Lni_wall, 0.05);                     // in m

  OPTION(options, spitzer_resist,         false);

  OPTION(options, BScurrent,         false);
  OPTION(options, Aratio,            0.35);

  OPTION(options, nlfilter_noisy_data, false);
  OPTION(options, nlfilter_Gradpar_logNn, false);
  OPTION(options, nlfilter_Vnpol, false);
  OPTION(options, filter_para, 0.20);

  OPTION(options, term_GparkappaTe_GparTe, true); 
  OPTION(options, terms_Ve, false);
  OPTION(options, Secondary_electron_emission, false);
  OPTION(options, terms_spitzer_resist,         false);
  OPTION(options, thermal_force, false);
  OPTION(options, terms_radial_conv, false);
  OPTION(options, terms_recombination, false);
  OPTION(options, terms_ionization, true);
  OPTION(options, terms_exb, false);  
  OPTION(options, terms_cross, false);  
  OPTION(options, gyroviscous, false);  
  OPTION(options, include_J0, false);  
  OPTION(options, J1_all, false);  
  OPTION(options, terms_Gradpar_pn, true);  
  OPTION(options, terms_Gradpar_eta0n, true);  
  OPTION(options, terms_Diffcnn_par, false);  
  OPTION(options, terms_NnGradpar_Vn, true);  
  OPTION(options, terms_NnGradpar_Vi, true);  
  OPTION(options, terms_Gradperp_diffcoefs, true);  
  OPTION(options, terms_Gradpara_diffcoefs, true);  

  OPTION(options, delta, 0.0);
  OPTION(options, C_fe_sheat, 2.3);
  OPTION(options, C_fi_sheat, 0.0);
  
  OPTION(options, alpha_vn, 0.);            // a.u. control para. for Vn at SBC
  OPTION(options, angle_B_plate, 3.14/6.);  // Read in rad, angle between B and plate

  OPTION(options, AA, 2.0);
  OPTION(options, ZZ, 1.0);
  Zi=ZZ;                                    // Zi is applied in BS_current from Tianyang

  OPTION(options, Mi, 1.0);                 // Read in Mi
  OPTION(options, Mn, 1.0);                 // Read in Mi
  OPTION(options, Mm, 2.0);                 // Read in Mi
  OPTION(options, W_ionz, 20.0);            // Read in eV
  OPTION(options, W_diss, 4.5);             // in eV
  OPTION(options, W_rec, 4.5);              // in eV

  OPTION(options, W_bind, 0.5);             // in eV

  OPTION(options, dNidx_xin_au, -65.);      // a.u.
  OPTION(options, dTidx_xin_au, -3500.);    // a.u.
  OPTION(options, dTedx_xin_au, -3500.);    // a.u.
  OPTION(options, psi_xout_y0, 0.229543);   //m^2 T
  OPTION(options, Tn_plate, 1. );           //eV
  OPTION(options, Te_core, 1000. );         //eV
  OPTION(options, Te_edge, 10. );           //eV
  OPTION(options, Ti_core, 1000. );         //eV
  OPTION(options, Ti_edge, 10. );           //eV
  OPTION(options, Ni_core, 2. );            //in 1.e^20 m^-3
  OPTION(options, Ni_edge, 0.1 );           //in 1.e^20 m^-3
  OPTION(options, Initfile_x0, 0.4 );       //in a.u.
  OPTION(options, Initfile_w_ped, 0.2 );    //in a.u.
  OPTION(options, x0_ped, 0.90);            //in a.u. relative to normalized psi
  OPTION(options, width_ped, 0.062);        //in a.u.
  OPTION(options, coef_coregrad_ped, 0.01); //in a.u.
  OPTION(options, coef_corenlarge_ped, 18.);//in a.u.


  OPTION(options, q95_input,  5.);          // in a.u. a factor for heat flux limiter
  OPTION(options, q_alpha,  1.0);           //flux-limiting coefficient
  OPTION(options, Lnn_min,  1.e-3);         // in m a factor for neutral density flux limit
  OPTION(options, leakage_bndry, false); //leakage boundary condition
  OPTION(options, leak_ni, 1.e-3); //
  OPTION(options, leak_ti, 1.e-2);  //
  OPTION(options, leak_te,  1.e-4);
  OPTION(options, decay_bndry, false);
  OPTION(options, decay_bndry_legth,0.03);
  OPTION(options,energy_flux_bndry,false);
  OPTION(options, R_major, 1.65); //major radius
  OPTION(options, r_minor, 0.4);  //minor radius
  OPTION(options, Elgt,  1.1);    //elongation
  OPTION(options, P_input, 1e6); // input power in W
  OPTION(options, q_input, 1e5); //input energy flux
  OPTION(options, Vm0,   -500.0);           // Read in m/s
  //Vm0 = 0.0;
  OPTION(options, Nm0,  1.e7);              // Read in in 1.e^20 m^-3

  OPTION(options, diffusion_coef_step_function,  false); 
  OPTION(options, diffusion_coef_step0,  0.1);    // Read in m^2 / s 
  OPTION(options, diffusion_coef_step1,  1.0);    // Read in m^2 / s 
  OPTION(options, diffusion_coef_Hmode,  false); 
  OPTION(options, diffusion_coef_Hmode0,  1.);    // Read in m^2 / s 
  OPTION(options, diffusion_coef_Hmode1,  10.0);  // Read in m^2 / s 

  OPTION(options, Diffc_ni_perp,  0.1);           // Read in m^2 / s 
  OPTION(options, Difft_ni_perp,  1.0);
  OPTION(options, chic_i_perp,    0.4);
  OPTION(options, chit_i_perp,    4.0);
  OPTION(options, chic_e_perp,    0.6);
  OPTION(options, chit_e_perp,    6.0);
  OPTION(options, chic_n_perp, 0.4);
  OPTION(options, chi_coef, 1.0); //modify the coef
  
  //**** Boundaries of Const Flux of fueling *//
  OPTION(options, CF_BC_x0,    1.01);
  OPTION(options, CF_BC_y0,    0.47);
  OPTION(options, CF_BC_y1,    0.53);
  OPTION(options, CF_BC_z0,    0.0);
  OPTION(options, CF_BC_z1,    2.0);              //default smbi in whole Z
  OPTION(options, Sheath_BC_x0, 0.84);
  Diffc_nm_perp=0.01;
   
   // external sources
  OPTION(options, external_sources,  false); 
  OPTION(options, extsrcs_balance_diffusion,  false); 
  OPTION(options, x0_extS, 0.85);                 //in a.u. relative to normalized psi
  OPTION(options, width_extS, 0.062);             //in a.u.
  OPTION(options, coef_coregrad_extS, 0.021);     //in a.u.
  OPTION(options, amp_spi_ext, 1.e-5);            //in N0/s,  N0=1e^20/m^3
  OPTION(options, amp_see_ext, 1.e-4);            //in eV*N0/s
  OPTION(options, amp_sei_ext, 1.e-4);            //in eV*N0/s
 // Ion orbit loss
  OPTION(options, terms_IOL, false);
  OPTION(options, xout,    174);                  // Position in X direction, grid index
  OPTION(options, P_N,    24);                    // grid number of pitch number

  OPTION(options, test_f1,    1.0);               // grid number of pitch number
  OPTION(options, test_f2,    1.0);               // grid number of pitch number
  OPTION(options, test_f3,    1.0);               // grid number of pitch number

  // int bracket_method;
  OPTION(options, bracket_method_exb, 0);
  switch(bracket_method_exb) {
  case 0: {
    bm_exb = BRACKET_STD;
    output << "\tBrackets for ExB: default differencing\n";
    break;
  }
  case 1: {
    bm_exb = BRACKET_SIMPLE;
    output << "\tBrackets for ExB: simplified operator\n";
    break;
  }
  case 2: {
    bm_exb = BRACKET_ARAKAWA;
    output << "\tBrackets for ExB: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm_exb = BRACKET_CTU;
    output << "\tBrackets for ExB: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
    }
  // int bracket_method;
  OPTION(options, bracket_method_mag, 2);
  switch(bracket_method_mag) {
  case 0: {
    bm_mag = BRACKET_STD;
    output << "\tBrackets: default differencing\n";
    break;
  }
  case 1: {
    bm_mag = BRACKET_SIMPLE;
    output << "\tBrackets: simplified operator\n";
    break;
  }
  case 2: {
    bm_mag = BRACKET_ARAKAWA;
    output << "\tBrackets: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm_mag = BRACKET_CTU;
    output << "\tBrackets: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

//*****************************************************li2016
  OPTION(options, xloc,           1.0);               // Position in X direction, normalized Psi
  OPTION(options, yloc,           32);                // Position in Y direction, grid index
  OPTION(options, iterative_phi0,          true);
  OPTION(options, J_para0,          true);
  OPTION(options, sbc_lambda0,           1.0);
  OPTION(options, sbc_lambda1,           1.0);
  OPTION(options, Mu_perp,           1.0);
  OPTION(options, Mu_para,           1.0);
  OPTION(options, RF_coef,           1.0);
//**************************************************liend

//////////// CALCULATE PARAMETERS ///////////////////

 // in order to using formulas of wci and rhos in Gaussian cgs units, 
 // bmag in SI unit Tesla multiply by 1.e4 changes to Gaussian unit gauss 
 // Te_x in unit eV, no need to transfer

 // !!!Be very careful when using quantities calculated below, units transfer needed
  
  Ni_x *= 1.0e14;    // in unit cm^-3 now
  bmag *= 1.0e4;     // in unit gauss now

  output.write("Calculating parameters in Gaussian units and then transfer to SI units \n"); 
  output.write("\tIn Gaussian cgs units:  Ni_x %e cm^-3, Te_x %e eV, bmag %e gauss  \n",Ni_x,Te_x,bmag);
  output.write("\tparameters:  AA=mi/mp %e,  ZZ %e \n",AA,ZZ);

  rho_s = 1.02e2*sqrt(AA*Te_x)/ZZ/bmag;   // unit cm
  wci   = 9.58e3*ZZ*bmag/AA;
  Vi_x = wci * rho_s;                     // in unit cm/s                                     

  fmei  = 1./1836.2/AA;                                                      // no unit
  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);                                      // in unit cm
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));                // unit cm
  nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);                         // unit 1/s
  //nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);   // unit 1/s
  nuiix     = 4.78e-8*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);                // unit 1/s

  // after using formulas of wci and rhos in Gaussian cgs units, 
  // bmag in Gaussian unit gauss divided by 1.e4 changes back to Gaussian SI unit Tesla

  Ni_x /= 1.0e14;                                  // back to unit 1.e^20 m^-3 mow
  bmag /=1.0e4;                                    // back to unit tesla now
  Vi_x /= 1.e2;                                    // back to unit m/s

  tbar = Lbar/Vi_x;                                // in unit s

  N_cubic_Lbar = Ni_x*1.e20*Lbar*Lbar*Lbar;        // unit 1


  output.write("\tIn SI units:  Ni_x %e e^20 m^-3,  Te_x %e eV, bmag %e tesla \n", Ni_x,Te_x,bmag);

  ///////////// PRINT Z INFORMATION /////////////////////
  
  BoutReal hthe0;
  if(GRID_LOAD(hthe0) == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/Lbar);
  }

 if(!include_curvature)
    b0xcv = 0.0;

  if(noshear) {
    if(include_curvature)
      b0xcv.z += I*b0xcv.x;
    mesh->ShiftXderivs = false;
    I = 0.0;
  }
  
  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      if(include_curvature)
	b0xcv.z += I*b0xcv.x;
      I = 0.0;  
      // I disappears from metric
    }
  }

  ///////////// NORMALISE QUANTITIES ////////////////////

  output.write("\tNormalising to Lbar = %e m, tbar %e s, V_Ti %e m/s \n", Lbar, tbar, Vi_x);

   // Normalise geometry
  J0 = MU0*Lbar * J0 / B0; 
  Rxy /= Lbar;
  hthe /= Lbar;
  B0   /= Bbar;
  phi0 /=Lbar/tbar*Lbar*Bbar;
  phi01 = phi0*B0*sbc_lambda1;
  mesh->dx /= Lbar*Lbar*bmag;
  //mesh->dx /= Lbar;
  I   *= Lbar*Lbar*bmag;
  psixy /= Lbar*Lbar*bmag;
  psi_axis /= Lbar*Lbar*bmag;
  psi_bndry /= Lbar*Lbar*bmag;
  psi_xout_y0 /= Lbar*Lbar*bmag;
  // Normalise magnetic field
  b0xcv.x /= bmag;
  b0xcv.y *= Lbar*Lbar;
  b0xcv.z *= Lbar*Lbar;

  Bpxy /= bmag;
  Btxy /= bmag;
  mesh->Bxy  /= bmag;

  // Normalise coefficients
  Lni_wall /= Lbar;
  Lp_crit /= Lbar;
  W_ionz /= Te_x;
  W_diss /= Te_x;
  W_bind /= Te_x;

  Tn_plate /= Te_x;
  Tm_x /= Te_x;
  Te_core /= Te_x;
  Te_edge /= Te_x;
  Ti_core /= Te_x;
  Ti_edge /= Te_x;
  Ni_core /= Ni_x;
  Ni_edge /= Ni_x;
    
  Te_grid /= Te_x;
  Ti_grid /= Te_x;
  Ni_grid =Ni_grid/Ni_x;                // in unit 10^20/m^3 
  phi01_grid /=Lbar/tbar*Lbar*Bbar;  

  Te_exp /= Te_x;
  Ti_exp /= Te_x;
  Ni_exp /= Ni_x;  
  Ne_exp /= Ni_x;  

  Vm0 /= Lbar/tbar;

  Diffc_nm_perp /= Lbar*Lbar/tbar;
  Diffc_ni_perp /= Lbar*Lbar/tbar;
  Difft_ni_perp /= Lbar*Lbar/tbar;
  Diff_grid /= Lbar*Lbar/tbar;
  Chi_grid /= Lbar*Lbar/tbar;
  Che_grid /= Lbar*Lbar/tbar;

 // q_input = P_input/4./PI/PI/r_minor/R_major/sqrt(Elgt)/2.;
 // output.write("dinx: %e \n",q_input);
  q_input /= Lbar/tbar*Ni_x*Te_x*1.6e-19*1.e20; 
  
  chic_i_perp /= Lbar*Lbar/tbar;
  chit_i_perp /= Lbar*Lbar/tbar;
  chic_e_perp /= Lbar*Lbar/tbar;
  chit_e_perp /= Lbar*Lbar/tbar;

  diffusion_coef_step0 /= Lbar*Lbar/tbar;
  diffusion_coef_step1 /= Lbar*Lbar/tbar;
  diffusion_coef_Hmode0 /= Lbar*Lbar/tbar;
  diffusion_coef_Hmode1 /= Lbar*Lbar/tbar;

  amp_spi_ext /= Ni_x/tbar;
  amp_see_ext /= Ni_x*Te_x/tbar;
  amp_sei_ext /= Ni_x*Te_x/tbar;
 
  output.write("\tDiffusion coefficients in unit Lbar^2/tbar : \n");
  output.write("\tDiffc_perp %e, Difft_perp %e \n", Diffc_ni_perp, Difft_ni_perp);
  output.write("\tchic_i_perp %e, chit_i_perp %e \n", chic_i_perp,chit_i_perp);
  output.write("\tchic_e_perp %e, chit_e_perp %e \n", chic_e_perp,chit_e_perp);
  output.write("\tdiff_coef_step0 %e, diff_coef_step1 %e \n", diffusion_coef_step0,diffusion_coef_step1);
  output.write("\tdiff_coef_Hmode0 %e, diff_coef_Hmode1 %e \n", diffusion_coef_Hmode0,diffusion_coef_Hmode1);

//*********************************************li2016
  N0 = Ni_exp;
  Ti0 = Ti_exp;
  Te0 = Te_exp;
//*********************************************liend

  /////////////// CALCULATE METRICS /////////////////

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;  
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry();              // Calculate other metrics
 
  // Set B field vector

  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;
 
  // SET VARIABLE LOCATIONS *******************//

  Ni.setLocation(CELL_CENTRE);
  Ti.setLocation(CELL_CENTRE);
  Te.setLocation(CELL_CENTRE);
  Vi.setLocation(CELL_CENTRE);
  Nn.setLocation(CELL_CENTRE);
  Tn.setLocation(CELL_CENTRE);
  Vn.setLocation(CELL_CENTRE);
  Nm.setLocation(CELL_CENTRE);
  Vmx.setLocation(CELL_CENTRE); 
  U00.setLocation(CELL_CENTRE); 
  phi01.setLocation(CELL_CENTRE); 

  //////////////// BOUNDARIES ///////////////////////
  // Set BOUNDARiES first here, and then apply them every time in physics run/////
 
   Ni.setBoundary("Ni");
   Ti.setBoundary("Ti");
   Te.setBoundary("Te");
   Vi.setBoundary("Vi");
   Nn.setBoundary("Nn");
   Tn.setBoundary("Tn");
   Vn.setBoundary("Vn");
   Nm.setBoundary("Nm");
   //Vmx.setBoundary("Vmx");
   Vm.setBoundary("Vm");
   U00.setBoundary("U00");
   phi01.setBoundary("phi01");
   phi00.setBoundary("phi01");
   phish0.setBoundary("phi01");
   phish1.setBoundary("phi01");
   Vexb.setBoundary("neumann");
   Vdia.setBoundary("neumann"); 
   Vdia_e.setBoundary("neumann"); 
   er0.setBoundary("p");
   er00.setBoundary("p");
   (er0.x).setBoundary("p");
   (er00.x).setBoundary("p");
   ersh0.setBoundary("neumann");
   (ersh0.x).setBoundary("neumann");
   ersh1.setBoundary("neumann");
   (ersh1.x).setBoundary("neumann");

   //Set Boundary for other output variables
   tau_ei.setBoundary("Tn");
   nu_ionz.setBoundary("Tn");
   nu_CX.setBoundary("Tn");
   nu_ionz_n.setBoundary("Tn");
   nu_CX_n.setBoundary("Tn");
   Si_p.setBoundary("Tn");
   Si_CX.setBoundary("Tn");
   S_diss.setBoundary("Tn");
   Vn_perp.setBoundary("Tn");

   grad_perp_Diffi.setBoundary("Tn");
   grad_perp_Diffn.setBoundary("Tn");
   grad_perp_chii.setBoundary("Tn");
   grad_perp_chie.setBoundary("Tn");
   //grad_para_kappTe.setBoundary("Tn");
   //grad_para_kappTi.setBoundary("Tn");

   Grad_par_pei.setBoundary("Tn");
   Grad_par_pn.setBoundary("Tn");
   Grad_par_logNn.setBoundary("Tn");
   Grad_par_Tn.setBoundary("Tn");
   DDX_Ni.setBoundary("Tn");
   DDX_Ti.setBoundary("Tn");
   DDX_Te.setBoundary("Tn");
   D2DX2_Ni.setBoundary("Tn");
   D2DX2_Ti.setBoundary("Tn");
   //q_se_kappa.setBoundary("Tn");
   //q_si_kappa.setBoundary("Tn");
   D2DX2_Te.setBoundary("Tn");

  ///////////// SET EVOLVING VARIABLES //////////////
  //
  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object



  Ni= Te = Ti = 0.0;
  //U00= phi01 = 0.0;
  Vi=0.;
  Nn = Tn= Vn=0.0;       //NB: it is **NECESSARY** to set values for evolving quantities if not read from grid 
  Nm = Vmx = 0.0;
  Ni0_Hmode = Ti0_Hmode = Te0_Hmode=0.0;

  SOLVE_FOR4(Ni, Vi, Te, Ti);
  SOLVE_FOR3(Nn,Tn,Vn);
  //SOLVE_FOR2(Nm,Vmx);
  SOLVE_FOR2(Nm,Vm);
  SOLVE_FOR2(U00,phi01);

  if(profiles_lasttimestep)
    {
     output.write("Initial Profiles are loaded from .txt files of all evolving quatities at last time step \n");

     BoutReal lstimeNi[NX][NY][NZ],lstimeTi[NX][NY][NZ],lstimeTe[NX][NY][NZ],lstimeVi[NX][NY][NZ];
     BoutReal lstimeU00[NX][NY][NZ],lstimephi01[NX][NY][NZ];
     BoutReal lstimeNn[NX][NY][NZ],lstimeNm[NX][NY][NZ],lstimeVmx[NX][NY][NZ];  // Vn Tn are not evloved in recent version

	ifstream pFile1,pFile2,pFile3,pFile4,pFile5,pFile6,pFile7,pFile8,pFile9;
	pFile1.open ("data/lstime_ni.txt", ios::in );
	pFile2.open ("data/lstime_ti.txt", ios::in );
	pFile3.open ("data/lstime_te.txt", ios::in );
	pFile4.open ("data/lstime_vi.txt", ios::in );
	pFile5.open ("data/lstime_nn.txt", ios::in );
	pFile6.open ("data/lstime_nm.txt", ios::in );
	pFile7.open ("data/lstime_vmx.txt", ios::in );
	pFile8.open ("data/lstime_u00.txt", ios::in );
	pFile9.open ("data/lstime_phi01.txt", ios::in );

	output.write("\tCheck point 4 : files are open\n" ); 
	for(int jx=0;jx<NX;jx++) 
	 for(int jy=0;jy<NY;jy++) 
	   for(int jz=0;jz<NZ;jz++) 
             {
	       pFile1 >> lstimeNi[jx][jy][jz]; 
	       pFile2 >> lstimeTi[jx][jy][jz]; 
	       pFile3 >> lstimeTe[jx][jy][jz]; 
	       pFile4 >> lstimeVi[jx][jy][jz]; 
	       pFile5 >> lstimeNn[jx][jy][jz]; 
	       pFile6 >> lstimeNm[jx][jy][jz]; 
	       pFile7 >> lstimeVmx[jx][jy][jz]; 
	       pFile8 >> lstimeU00[jx][jy][jz]; 
	       pFile9 >> lstimephi01[jx][jy][jz]; 
             }
	pFile1.close();
	pFile2.close();
	pFile3.close();
	pFile4.close();
	pFile5.close();
	pFile6.close();
	pFile7.close();
	pFile8.close();
	pFile9.close();

	for(int jx=0;jx<mesh->ngx;jx++) 
	for(int jy=0;jy<mesh->ngy;jy++) 
	for(int jz=0;jz<mesh->ngz;jz++) 
        { 
	  Ni[jx][jy][jz] = lstimeNi[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Ti[jx][jy][jz] = lstimeTi[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Te[jx][jy][jz] = lstimeTe[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Vi[jx][jy][jz] = lstimeVi[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Nn[jx][jy][jz] = lstimeNn[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Nm[jx][jy][jz] = lstimeNm[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  Vmx[jx][jy][jz] = lstimeVmx[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  U00[jx][jy][jz] = lstimeU00[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	  phi01[jx][jy][jz] = lstimephi01[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
	}

    Vm.x=Vmx;
    Vm.y=0.;
    Vm.z=0.;

   if(initial_profile_Hmode) 
     {
      output.write("Restart from one chosen time step, Hmode profiles at early beginning also given\n");
      for (int jx=0;jx<mesh->ngx;jx++)
       {
         for (int jy=0;jy<mesh->ngy;jy++)
	  {
            BoutReal psi_normal = (psixy[jx][jy] - psi_axis)/(psi_bndry-psi_axis);  
            BoutReal x_prime = (x0_ped- psi_normal)/width_ped;
            BoutReal Hmode_tanh=(exp(-x_prime)+(1.+ coef_coregrad_ped*x_prime)*exp(x_prime)*coef_corenlarge_ped)/(exp(-x_prime)+exp(x_prime));
 	   for (int jz=0;jz<mesh->ngz;jz++)
            {            
             Ni0_Hmode[jx][jy][jz] = Ni_edge+Ni_core*(Hmode_tanh - 1.);
             Ti0_Hmode[jx][jy][jz] = Ti_edge+Ti_core*(Hmode_tanh - 1.);
             Te0_Hmode[jx][jy][jz] = Te_edge+Te_core*(Hmode_tanh - 1.);

	    }
	  }
        }
      }   //end initial_Hmode
    }
  else{
       if (load_grid_profiles)
        {
         output.write("\tInitial profiles of Ti Te Ni are loaded from grid file\n" ); 
         Te = Te_grid;
         Ti = Ti_grid;
         Ni = Ni_grid;
         phi01 = phi01_grid*B0*sbc_lambda1;
         phi00=phi01/sbc_lambda1;
         U00 = 0.0;
         Vn=0.0;
	 Vi=0.0;                                             
         Tn = minimum_val;
         Nn= minimum_val;
         Nm= minimum_val;
	 Vmx=0.;
        }  
  
      if (load_experiment_profiles)
       {
         output.write("\tInitial experiment profiles of Ti Te Ni are loaded from grid file\n" ); 
         Te = Te_exp;
         Ti = Ti_exp;
         Ni = Ni_exp;
         Ne = Ne_exp;
         //phi01 = phi01_grid*B0*sbc_lambda1;
         phi01= U00 = 0.0;
         phi00=phi01/sbc_lambda1;
         Vn=0.0;
	 Vi=0.0;                                             
         Tn = minimum_val;
         Nn= minimum_val;
         Nm= minimum_val;
	 Vmx=0.;
        }  

  // Initialization of Profiles 
  // ****** NB: profiles should be initialized after "SOLVE_FOR()" ****//
  S_pi_ext = 0.;
  S_Ei_ext = 0.;
  S_Ee_ext = 0.;

  for (int jx=0;jx<mesh->ngx;jx++)
     {
      x_rela= mesh->GlobalX(jx); 

      //BoutReal x0=0.3,w_ped=0.1;
      BoutReal temp=exp(2.*(x_rela-Initfile_x0)/Initfile_w_ped);
      BoutReal x0_nn=1.02,w_nn=0.05;
      BoutReal temp2=exp(-(x_rela-x0_nn)*(x_rela-x0_nn)/w_nn/w_nn);

      for (int jy=0;jy<mesh->ngy;jy++)
	{
        BoutReal x_psi_l = psixy[jx][jy]-psi_xout_y0;
        BoutReal psi_normal = (psixy[jx][jy] - psi_axis)/(psi_bndry-psi_axis);  
        BoutReal x_prime = (x0_ped- psi_normal)/width_ped;
        BoutReal Hmode_tanh=(exp(-x_prime)+(1.+ coef_coregrad_ped*x_prime)*exp(x_prime)*coef_corenlarge_ped)/(exp(-x_prime)+exp(x_prime));

        BoutReal y_rela=mesh->GlobalY(jy);
        int jy_global = mesh->YGLOBAL(jy);
        BoutReal y0_nn=0.5,wy_nn=0.05;
        BoutReal temp2_y=exp(-(y_rela-y0_nn)*(y_rela-y0_nn)/wy_nn/wy_nn);
 	for (int jz=0;jz<mesh->ngz;jz++)
          {
	    if(!load_grid_profiles)
	      {
                if(initial_profile_exp)
	          {
	            Ni[jx][jy][jz]=Ni_edge+Ni_core/(1.+temp);
                    Ti[jx][jy][jz]=Te_edge+Te_core/(1.+temp);
                    Te[jx][jy][jz]=Ti_edge+Ti_core/(1.+temp);
	           }
                if(initial_profile_linear)
                  {
                   Ni[jx][jy][jz] = Ni_edge+dNidx_xin_au*x_psi_l;
                   Ti[jx][jy][jz] = Ti_edge+dTidx_xin_au*x_psi_l;
                   Te[jx][jy][jz] = Te_edge+dTedx_xin_au*x_psi_l;
                  }
                if(initial_profile_Hmode) 
                  {
                   Ni[jx][jy][jz] = Ni_edge+Ni_core*(Hmode_tanh - 1.);
                   Ti[jx][jy][jz] = Ti_edge+Ti_core*(Hmode_tanh - 1.);
                   Te[jx][jy][jz] = Te_edge+Te_core*(Hmode_tanh - 1.);
                   Ni0_Hmode[jx][jy][jz]=Ni[jx][jy][jz];
                   Ti0_Hmode[jx][jy][jz]=Ti[jx][jy][jz];
                   Te0_Hmode[jx][jy][jz]=Te[jx][jy][jz];
                  }
		if(initial_SOL_edgeval)
                  {
		    if(jy_global<=jyseps11 || jy_global>=NY-1-jyseps11)
		      {
			Ni[jx][jy][jz] = Ni_edge;
		        //Ne[jx][jy][jz] = Ne_edge;
                        Ti[jx][jy][jz] = Ti_edge;
                        Te[jx][jy][jz] = Te_edge;
                       }
                   }
	      }

             Vn[jx][jy][jz]=0.0;
	     Vi[jx][jy][jz]=0.0;                              
	     U00[jx][jy][jz]=0.0;                              
             phi01[jx][jy][jz]=phi0[jx][jy][jz];
             phi00=phi01/sbc_lambda1;
             //phi01*=B0*sbc_lambda1;
             // INITIALIZE
             Tn[jx][jy][jz] = minimum_val;         // NB:can not be 0.0
             Nn[jx][jy][jz] = minimum_val;
             Nm[jx][jy][jz] = minimum_val;
	     Vmx[jx][jy][jz]=0.;

	     if(external_sources)
               {
                S_pi_ext[jx][jy][jz] = amp_spi_ext*(Hmode_tanh - 1.);
                S_Ei_ext[jx][jy][jz] = amp_sei_ext*(Hmode_tanh - 1.);
                S_Ee_ext[jx][jy][jz] = amp_see_ext*(Hmode_tanh - 1.);
               }

           }
        }
     }

  } // end of else if (profiles_lasttimestep)

  //End of Initialization of Profiles


  // Set step functions of diffusion coefficients
  Diffc_ni_step=0.;
  chic_i_step=0.;
  chic_e_step=0.; 
  Diffc_ni_Hmode=0.;
  chic_i_Hmode=0.;
  chic_e_Hmode=0.; 
  
  // calculate radial derivative for Hmode diffusion coefficients
  if (diffusion_coef_Hmode && initial_profile_Hmode)
    {
     output.write("Calculated diffusion coefficients for H-mode profiles");
     DDX_Ni = DDX(Ni0_Hmode);
     mesh->communicate(DDX_Ni);
     DDX_Ni.applyBoundary();

     DDX_Ti = DDX(Ti0_Hmode);
     mesh->communicate(DDX_Ti);
     DDX_Ti.applyBoundary();

     DDX_Te = DDX(Te0_Hmode);
     mesh->communicate(DDX_Te);
     DDX_Te.applyBoundary();

     aveY_J=mesh->averageY(mesh->J);
     aveY_g11=mesh->averageY(mesh->g11);

     aveY_g11J_ni=mesh->averageY((mesh->g11*mesh->J*DDX_Ni).DC());
     aveY_g11J_ti=mesh->averageY((mesh->g11*mesh->J*DDX_Ti).DC());
     aveY_g11J_te=mesh->averageY((mesh->g11*mesh->J*DDX_Te).DC());
    }
  // calculate radial derivative for experimental diffusion coefficients
  if (load_experiment_profiles)
    {
     output.write("Calculated diffusion coefficients for experimental profiles \n");
     DDX_Ni = DDX(Ni);
     mesh->communicate(DDX_Ni);
     DDX_Ni.applyBoundary();

     DDX_Ti = DDX(Ti);
     mesh->communicate(DDX_Ti);
     DDX_Ti.applyBoundary();

     DDX_Te = DDX(Te);
     mesh->communicate(DDX_Te);
     DDX_Te.applyBoundary();

     aveY_J=mesh->J;
     aveY_g11=mesh->g11;

     /*aveY_g11J_ni=(mesh->g11*mesh->J*DDX_Ni).DC();
     aveY_g11J_ti=(mesh->g11*mesh->J*DDX_Ti*Ni).DC();
     aveY_g11J_te=(mesh->g11*mesh->J*DDX_Te*Ni).DC();
     */
     aveY_g11J_ni=(mesh->g11*mesh->J*DDX_Ni).DC();
     aveY_g11J_ti=(mesh->g11*mesh->J*DDX_Ti).DC();
     aveY_g11J_te=(mesh->g11*mesh->J*DDX_Te).DC();
    }
  if (load_experiment_profiles)
    {
     GlobalField2D g_g11(mesh),g_J(mesh);
     GlobalField3D gNi(mesh),gDDX_Ni(mesh),gDDX_Ti(mesh),gDDX_Te(mesh),gGamma_ni_xin(mesh),gQi_xin(mesh),gQe_xin(mesh);
     g_g11.gather(aveY_g11);
     g_J.gather(aveY_J);
     gNi.gather(Ni);
     gDDX_Ni.gather(DDX_Ni);
     gDDX_Ti.gather(DDX_Ti);
     gDDX_Te.gather(DDX_Te);
     if( gNi.dataIsLocal() )
       {
        for(int i=0;i<gNi.xSize();i++)  
          {
           for(int j=0;j<gNi.ySize();j++)
            {
             for(int k=0;k<gNi.zSize();k++)
              {
               gGamma_ni_xin(i,j,k) = -diffusion_coef_Hmode0*g_g11(0,j)*g_J(0,j)*gDDX_Ni(0,j,k);
               gQi_xin(i,j,k) = -diffusion_coef_Hmode0*g_g11(0,j)*g_J(0,j)*gDDX_Ti(0,j,k);
               gQe_xin(i,j,k) = -diffusion_coef_Hmode0*g_g11(0,j)*g_J(0,j)*gDDX_Te(0,j,k);
                }
              }
            }
         }
     Gamma_ni_xin_exp = gGamma_ni_xin.scatter();
     Qi_xin_exp = gQi_xin.scatter();
     Qe_xin_exp = gQe_xin.scatter();
     } //end of if (load_experiment_profiles)

  for (int jx=0;jx<mesh->ngx;jx++)
     {
      int jx_global =  mesh->XGLOBAL(jx);
      // BoutReal dindx = indx/ixsep;
      for (int jy=0;jy<mesh->ngy;jy++)
	{
         int jy_global = mesh->YGLOBAL(jy);
         BoutReal psi_normal = (psixy[jx][jy] - psi_axis)/(psi_bndry-psi_axis);
         for (int jz=0;jz<mesh->ngz;jz++)
          {
            if(psi_normal>1.0) 
              {
               Diffc_ni_step[jx][jy][jz]=diffusion_coef_step1;
               chic_i_step[jx][jy][jz]=diffusion_coef_step1;
               chic_e_step[jx][jy][jz]=diffusion_coef_step1;
              }
            else
              {
               Diffc_ni_step[jx][jy][jz]=diffusion_coef_step0;
               chic_i_step[jx][jy][jz]=diffusion_coef_step0;
               chic_e_step[jx][jy][jz]=diffusion_coef_step0;
              }

	    if(diffusion_coef_Hmode && initial_profile_Hmode)
	      {
	       //BoutReal Gamma_ni_xin = -diffusion_coef_Hmode0*sqrt(mesh->g11[jx][jy])*dNidx_xin_au;
       	       //BoutReal Qi_xin = -diffusion_coef_Hmode0*sqrt(mesh->g11[jx][jy])*dTidx_xin_au;
	       //BoutReal Qe_xin = -diffusion_coef_Hmode0*sqrt(mesh->g11[jx][jy])*dTedx_xin_au;

               BoutReal Gamma_ni_xin = -diffusion_coef_Hmode0*aveY_g11[jx][jy]*dNidx_xin_au;
               BoutReal Qi_xin = -diffusion_coef_Hmode0*aveY_g11[jx][jy]*dTidx_xin_au;
               BoutReal Qe_xin = -diffusion_coef_Hmode0*aveY_g11[jx][jy]*dTedx_xin_au;

	       if(DDX_Ni[jx][jy][jz]>-minimum_val) DDX_Ni[jx][jy][jz]=-minimum_val;
	       if(DDX_Ti[jx][jy][jz]>-minimum_val) DDX_Ti[jx][jy][jz]=-minimum_val;
	       if(DDX_Te[jx][jy][jz]>-minimum_val) DDX_Te[jx][jy][jz]=-minimum_val;

	       // Diffc_ni_Hmode[jx][jy][jz] = -Gamma_ni_xin/sqrt(mesh->g11[jx][jy])/DDX_Ni[jx][jy][jz];
               //chic_i_Hmode[jx][jy][jz] = -Qi_xin/sqrt(mesh->g11[jx][jy])/DDX_Ti[jx][jy][jz];
               //chic_e_Hmode[jx][jy][jz] = -Qe_xin/sqrt(mesh->g11[jx][jy])/DDX_Te[jx][jy][jz];

               Diffc_ni_Hmode[jx][jy][jz] = -Gamma_ni_xin*aveY_J[jx][jy]/aveY_g11J_ni[jx][jy][jz];
               chic_i_Hmode[jx][jy][jz] = -Qi_xin*aveY_J[jx][jy]/aveY_g11J_ti[jx][jy][jz];
               chic_e_Hmode[jx][jy][jz] = -Qe_xin*aveY_J[jx][jy]/aveY_g11J_te[jx][jy][jz];


               if(Diffc_ni_Hmode[jx][jy][jz]>diffusion_coef_Hmode1) Diffc_ni_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;
               if(chic_i_Hmode[jx][jy][jz] > diffusion_coef_Hmode1) chic_i_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;
               if(chic_e_Hmode[jx][jy][jz] > diffusion_coef_Hmode1) chic_e_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;
              }
            if(load_experiment_profiles)
              {
               if(DDX_Ni[jx][jy][jz]>-minimum_val) {DDX_Ni[jx][jy][jz]=-minimum_val;}
               if(DDX_Ti[jx][jy][jz]>-minimum_val) {DDX_Ti[jx][jy][jz]=-minimum_val;}
               if(DDX_Te[jx][jy][jz]>-minimum_val) {DDX_Te[jx][jy][jz]=-minimum_val;}
               if(aveY_g11J_ti[jx][jy][jz]>-minimum_val) {aveY_g11J_ni[jx][jy][jz]=-minimum_val;}
               if(aveY_g11J_ti[jx][jy][jz]>-minimum_val) {aveY_g11J_ti[jx][jy][jz]=-minimum_val;}
               if(aveY_g11J_te[jx][jy][jz]>-minimum_val) {aveY_g11J_te[jx][jy][jz]=-minimum_val;}

               Diffc_ni_Hmode[jx][jy][jz] = abs(-Gamma_ni_xin_exp[jx][jy][jz]/aveY_g11J_ni[jx][jy][jz]);
               chic_i_Hmode[jx][jy][jz] = abs(-Qi_xin_exp[jx][jy][jz]/aveY_g11J_ti[jx][jy][jz]);
               chic_e_Hmode[jx][jy][jz] = abs(-Qe_xin_exp[jx][jy][jz]/aveY_g11J_te[jx][jy][jz]);

   
               //if(jx_global <= ixsep && Diffc_ni_Hmode[jx][jy][jz]>diffusion_coef_Hmode0) {Diffc_ni_Hmode[jx][jy][jz]=diffusion_coef_Hmode0;}
               //if(chic_i_Hmode[jx][jy][jz] > 2.*diffusion_coef_Hmode0) {chic_i_Hmode[jx][jy][jz]=2.0*diffusion_coef_Hmode0;}
               //if(chic_e_Hmode[jx][jy][jz] > 2.*diffusion_coef_Hmode0) {chic_e_Hmode[jx][jy][jz]=2.0*diffusion_coef_Hmode0;}     
               //if(jx_global > ixsep && Diffc_ni_Hmode[jx][jy][jz]>diffusion_coef_Hmode1) {Diffc_ni_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;}
               if(Diffc_ni_Hmode[jx][jy][jz]>diffusion_coef_Hmode1) {Diffc_ni_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;}
               if(chic_i_Hmode[jx][jy][jz] > diffusion_coef_Hmode1) {chic_i_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;}
               if(chic_e_Hmode[jx][jy][jz] > diffusion_coef_Hmode1) {chic_e_Hmode[jx][jy][jz]=diffusion_coef_Hmode1;}
                  
               if(jy_global<=jyseps11 || jy_global>=NY-1-jyseps11)
                 {
                  Diffc_ni_Hmode[jx][jy][jz]=diffusion_coef_Hmode0;
                  chic_i_Hmode[jx][jy][jz]=diffusion_coef_Hmode0;
                  chic_e_Hmode[jx][jy][jz]=diffusion_coef_Hmode0;
                 }
               }
            }
         }
      }
  // end of calculating radial derivative for H-mode or experimental diffusion coefficients
  	  
  if(extsrcs_balance_diffusion)
    {
      D2DX2_Ni = D2DX2(Ni);
      mesh->communicate(D2DX2_Ni);
      D2DX2_Ni.applyBoundary();

      D2DX2_Ti = D2DX2(Ti);
      mesh->communicate(D2DX2_Ti);
      D2DX2_Ti.applyBoundary();

      D2DX2_Te = D2DX2(Te);
      mesh->communicate(D2DX2_Te);
      D2DX2_Te.applyBoundary();

      S_pi_ext= - D2DX2_Ni*Diffc_ni_Hmode*mesh->g11;  
      S_Ei_ext= - D2DX2_Ti*chic_i_Hmode*Ni*mesh->g11;  
      S_Ee_ext= - D2DX2_Te*chic_e_Hmode*Ni*mesh->g11;  
    }

  ///////////// ADD OUTPUT VARIABLES ////////////////
  //
  // Add any other variables to be dumped to file
  
  SAVE_ONCE5(Te_x, Ti_x, Ni_x, Lbar, tbar);    // Normalisation factors
  SAVE_ONCE2(bmag,Vi_x);
  SAVE_ONCE(J0);

  // Set flux limit for kappa
  V_th_e= 4.19e5*sqrt(Te*Te_x);
  V_th_i= 9.79e3*sqrt(Ti*Te_x/AA);
  output.write("\tion thermal velocity: %e -> %e [m/s]\n", min(V_th_i), max(V_th_i));
  output.write("\telectron thermal velocity: %e -> %e [m/s]\n", min(V_th_e), max(V_th_e));
  V_th_e /= Lbar/tbar;
  V_th_i /= Lbar/tbar;
  output.write("\tNormalized ion thermal velocity: %e -> %e [Lbar/tbar]\n", min(V_th_i), max(V_th_i));
  output.write("\tNormalized electron thermal velocity: %e -> %e [Lbar/tbar]\n", min(V_th_e), max(V_th_e));

  //kappa_Te = 3.2*(1./fmei)*(1./tbar/nueix)*(Te^2.5);     // power operator '^' works only for Fields
  //kappa_Ti = 3.9*(1./tbar/nuiix)*(Ti^2.5); 
  nu_ei = nueix*Ni/(Te*sqrt(Te));
  nu_ii = nuiix*Ni/(Ti*sqrt(Ti)); 
  kappa_Te = 3.2*V_th_e*V_th_e/nu_ei/tbar*Ni;
  kappa_Ti = 3.9*V_th_i*V_th_i/nu_ii/tbar*Ni;  
  Field3D kappa_Te_realunits, kappa_Ti_realunits;
  kappa_Te_realunits= kappa_Te*Ni_x*pow(Lbar,2.)/tbar;    // pow function works for BoutReal variables     
  kappa_Ti_realunits= kappa_Ti*Ni_x*pow(Lbar,2.)/tbar;

  output.write("\tion para thermal conductivity: %e -> %e [N0 m^2/s]\n", min(kappa_Ti_realunits), max(kappa_Ti_realunits));
  output.write("\telectron para thermal conductivity: %e -> %e [N0 m^2/s]\n", min(kappa_Te_realunits), max(kappa_Te_realunits));   
  output.write("\tNormalzied ion para thermal conductivity: %e -> %e [N0 Lbar^2/tbar]\n", min(kappa_Ti), max(kappa_Ti));
  output.write("\tNormalized electron para thermal conductivity: %e -> %e [N0 Lbar^2/tbar]\n", min(kappa_Te), max(kappa_Te));  

  q95= q95_input*q_alpha;
  if(q95_input < 0.0){q95 = (abs(hthe*Btxy/Bpxy))*q_alpha;}
  kappa_Te_fl = V_th_e*q95*Lbar*Ni;                                   // Ne=Ni quasineutral
  kappa_Ti_fl = V_th_i*q95*Lbar*Ni;

  kappa_Te *= kappa_Te_fl/(kappa_Te_fl+kappa_Te);
  kappa_Ti *= kappa_Ti_fl/(kappa_Ti_fl+kappa_Ti);

  output.write("\tUsed ion para thermal conductivity: %e -> %e [N0 Lbar^2/tbar]\n", min(kappa_Ti), max(kappa_Ti));
  output.write("\tUsed electron para thermal conductivity: %e -> %e [N0 Lbar^2/tbar]\n", min(kappa_Te), max(kappa_Te));  

  // Ionization rate  depending on Nn for Plasmas 
  //nu_ionz = 3.e4*Nn*Ni_x*Te*Te*Te_x*Te_x/(3.+0.01*Te*Te*Te_x*Te_x);   //Ni_x in unit 1e^20 
  nu_ionz = 0.0;
  if(terms_ionization)
   {
    for (int jx=0;jx<mesh->ngx;jx++)
      {
       for (int jy=0;jy<mesh->ngy;jy++)
         {
          for (int jz=0;jz<mesh->ngz;jz++)
              {
               nu_ionz[jx][jy][jz] = 3.e6*tbar*Nn[jx][jy][jz]*Ni_x*(3.854*exp(1.268e-5-(Te[jx][jy][jz]*Te_x+3.13)/6703.89)*erfc(3.561e-3-(Te[jx][jy][jz]*Te_x+3.13)/47.745)-4.34);   //Ni_x in unit 1e^20
              }
           }
        }
    }
  //nu_CX =  1.e6*Nn*Ni_x*(1.7+0.667*((1.5*Ti*Te_x)^0.333-2.466));    // empirical formula 
  nu_CX =  1.e-14*Nn*Ni_x*density_unit*(1.7+1.9*(((1.5*Ti*Te_x)^0.333)-2.466)/(((150.*Ti*Te_x)^0.333)-2.466));      // empirical formula

  //Dissociation rate of molecules
  nu_diss = 3.e4*Nm*Ni_x*Te*Te*Te_x*Te_x/(3.+0.01*Te*Te*Te_x*Te_x);   //need be corrected

  // recombination rate of ions and electrons
  Field3D lambda_rec=1.5789e5/(abs(Te)*Te_x*1.1604e4);                //Te trasfer unit from eV to Kelvin  
  nu_rec = Ni*Ni_x*5.197*ZZ*sqrt(lambda_rec)*(0.4288+0.5*log(lambda_rec)+0.469/(lambda_rec)^0.333); //Seaton M F 1959 'Radiative recombination of hydrogenic ions'

  output.write("\tionization rate: %e -> %e [1/s]\n", min(nu_ionz), max(nu_ionz));
  output.write("\tcharge exchange rate: %e -> %e [1/s]\n", min(nu_CX), max(nu_CX)); 
  output.write("\tdissociation rate: %e -> %e [1/s]\n", min(nu_diss), max(nu_diss));
  output.write("\trecombination rate: %e -> %e [1/s]\n", min(nu_rec), max(nu_rec));   

  if(spitzer_resist) 
    {
    // Use Spitzer resistivity 
    output.write("");
    output.write("\tSpizter parameters");
    eta_spitzer = 1.03e-4*Zi*lambda_ei*((Te*Te_x)^(-1.5));                 // eta in Ohm-m. 
    output.write("\tSpitzer resistivity: %e -> %e [Ohm m]\n", min(eta_spitzer), max(eta_spitzer));
    eta_spitzer /= MU0 * Lbar * Lbar/tbar;
    output.write("\t -> Lundquist %e -> %e\n", 1.0/max(eta_spitzer), 1.0/min(eta_spitzer));
    dump.add(eta_spitzer, "eta_spitzer", 0);
    eta = eta_spitzer;
    dump.add(eta, "eta", 0);
   }
 
  //Bootstrap current calculated by using Sauter's formula 
  if (BScurrent)
    {
      q95=q95_input;
      Jpar_BS0.setLocation(CELL_YLOW);
      Jpar_BS0.setBoundary("Tn");
      pei= Ni*(Te+Ti);
      Pe = Ni*Te;
      Pi = Ni*Ti;

      nu_estar = 100.*nu_ei * q95*tbar / (V_th_e) / (Aratio*sqrt(Aratio));
      nu_istar = 100.*nu_ii * q95*tbar / (V_th_i) / (Aratio*sqrt(Aratio));
      //nu_estar = 0.012 * N0*Nbar*density/1.e20*Zi*Zi*q95*Lbar/(Te0*Tebar/1000. * Aratio^1.5);
      //nu_istar = 0.012 * N0*Nbar*density/1.e20*Zi*q95*Lbar/(Ti0*Tibar/1000. * Aratio^1.5);

      output.write("Bootstrap current is included: \n");
      output.write("Normalized electron collisionality: nu_e* = %e\n", max(nu_estar));
      output.write("Normalized ion collisionality: nu_i* = %e\n", max(nu_istar));
      ft = BS_ft(100);
      output.write("modified collisional trapped particle fraction: ft = %e\n", max(ft));
      f31 = ft / (1.+(1.-0.1*ft)*sqrt(nu_estar) + 0.5*(1.-ft)*nu_estar/Zi);
      f32ee = ft / (1.+0.26*(1.-ft)*sqrt(nu_estar) + 0.18*(1.-0.37*ft)*nu_estar/sqrt(Zi));
      f32ei = ft / (1.+(1.+0.6*ft)*sqrt(nu_estar) + 0.85*(1.-0.37*ft)*nu_estar*(1.+Zi));
      f34 = ft / (1.+(1.-0.1*ft)*sqrt(nu_estar) + 0.5*(1.-0.5*ft)*nu_estar/Zi);

      L31 = F31(f31) ;
      L32 = F32ee(f32ee)+F32ei(f32ei) ;
      L34 = F31(f34) ;

      BSal0 = - (1.17*(1.-ft))/(1.-0.22*ft-0.19*ft*ft);
      BSal = (BSal0+0.25*(1-ft*ft)*sqrt(nu_istar))/(1.+0.5*sqrt(nu_istar)) + 0.31*nu_istar*nu_istar*ft*ft*ft*ft*ft*ft;
      BSal *= 1./(1.+0.15*nu_istar*nu_istar*ft*ft*ft*ft*ft*ft);

      Jpar_BS0 = L31* DDX(pei)/Pe  + L32*DDX(Te)/Te + L34*DDX(Ti)/(Zi*Te)*BSal;
      Jpar_BS0 *= Field3D( -Rxy*Btxy*Pe*(MU0*KB*Ni_x*density_unit*Te_x*eV_K)/(mesh->Bxy*mesh->Bxy)/(bmag*bmag) ); 
      // NB:J_hat = MU0*Lbar * J / mesh->Bxy;

      mesh->communicate(Jpar_BS0);
      Jpar_BS0.applyBoundary();

      dump.add(Jpar_BS0, "jpar_BS0", 1);
    }

  dump.add(psixy,"psixy",0);

  if (diffusion_coef_Hmode)
    {
     dump.add(aveY_J,"aveY_J",0); 
     dump.add(aveY_g11,"aveY_g11",0); 
     dump.add(aveY_g11J_ni,"aveY_g11J_ni",0); 

     dump.add(DDX_Ni,"DDX_Ni",0); 
     dump.add(DDX_Ti,"DDX_Ti",0); 
     dump.add(DDX_Te,"DDX_Te",0); 
    }
  dump.add(Diffc_ni_step,"Diffc_ni_step",0);      // 0: output only at initial
  dump.add(chic_i_step,"chic_i_step",0);
  dump.add(chic_e_step,"chic_e_step",0);
  dump.add(Diffc_ni_Hmode,"Diffc_ni_Hmode",0); 
  dump.add(chic_i_Hmode,"chic_i_Hmode",0);
  dump.add(chic_e_Hmode,"chic_e_Hmode",0);
  
  if(extsrcs_balance_diffusion)
    {
     dump.add(S_pi_ext,"S_pi_ext",0);
     dump.add(S_Ei_ext,"S_Ei_ext",0);
     dump.add(S_Ee_ext,"S_Ee_ext",0);
    }
  dump.add(Diff_ni_perp,"Diff_ni_perp",1); 
  dump.add(chi_i_perp,"chi_i_perp",1);
  dump.add(chi_e_perp,"chi_e_perp",1);

  dump.add(kappa_Te,"kappa_Te",1);               // 1: output at any output step
  dump.add(kappa_Ti,"kappa_Ti",1);
  dump.add(Diffc_nn_perp,"Diffc_nn_perp",1);

  dump.add(nu_ionz,"nu_ionz",1);
  dump.add(nu_CX,"nu_CX",1);
  dump.add(nu_ionz_n,"nu_ionz_n",1);
  dump.add(nu_CX_n,"nu_CX_n",1);
  dump.add(nu_diss,"nu_diss",1);
  dump.add(nu_rec,"nu_rec",1);

  dump.add(Si_p,"Si_p",1);
  dump.add(Si_CX,"Si_CX",1);
  dump.add(S_diss,"S_diss",1);
  dump.add(S_rec,"S_rec",1);
  dump.add(Grad_par_pei,"Grad_par_pei",1);
  dump.add(Grad_par_pn,"Grad_par_pn",1);
  dump.add(Grad_par_logNn,"Grad_par_logNn",1);

  dump.add(tau_ei,"tau_ei",1);
  dump.add(nu_ei,"nu_ei",1);
  dump.add(nu_ii,"nu_ii",1);

  //dump.add(q95,"q95",0);
  dump.add(Vn_perp,"Vn_perp",1);
  dump.add(q_se_kappa,"q_se_kappa",1);
  dump.add(q_si_kappa,"q_si_kappa",1);
  dump.add(Gamma_ni_Diffc,"Gamma_ni_Diffc",1);
  dump.add(Gamma_nn_Diffc,"Gamma_nn_Diffc",1);

  dump.add(phi00, "phi00",1);
  dump.add(phish0, "phish0",1);
  dump.add(phish1, "phish1",1);
  dump.add(er00.x, "er00x",1);
  dump.add(Vexb, "vexb", 1);
  dump.add(Vdia, "vdia", 1);
  dump.add(Vdia_e, "vdia_e", 1);
  dump.add(J1_para, "J1_para",1);
  dump.add(Ve, "Ve",1);
  dump.add(J1_phi, "J1_phi",1);
  dump.add(J1_Pe, "J1_Pe",1);
  dump.add(J1_Te, "J1_Te",1);
  dump.add(term_J1, "term_J1",1);
  dump.add(term_J0, "term_J0",1);
  dump.add(term_Mu_par, "term_Mu_par",1);
  dump.add(term_Mu_perp, "term_Mu_perp",1);
  dump.add(term_gyro, "term_gyro",1);
  dump.add(term_curvature, "term_curvature",1);
  dump.add(term_exb, "term_exb",1);
  dump.add(er0.x, "er0x",1);
  dump.add(ersh0.x, "ersh0x",1);
  dump.add(ersh1.x, "ersh1x",1);

  //dump.add(Ni,"Ni",1); 
  //dump.add(Ti,"Ti",1); 
  //dump.add(Te,"Te",1); 
  dump.add(c_se,"c_se",1);
 
  dump.add(heatf_cond_i, "heatf_cond_i",1);
  dump.add(heatf_cond_e, "heatf_cond_e",1);
  dump.add(heatf_conv_i, "heatf_conv_i",1);
  dump.add(heatf_conv_e, "heatf_conv_e",1);
  dump.add(heatf_conv_e1, "heatf_conv_e1",1);
  dump.add(heatf_cond_i_x, "heatf_cond_i_x",1);
  dump.add(heatf_cond_e_x, "heatf_cond_e_x",1);
  dump.add(heatf_exb_i_x, "heatf_exb_i_x",1);
  dump.add(heatf_exb_e_x, "heatf_exb_e_x",1);
  dump.add(heatf_dia_i_x, "heatf_dia_i_x",1);
  dump.add(heatf_dia_e_x, "heatf_dia_e_x",1);
  dump.add(heatf_eng_i_x, "heatf_eng_i_x",1);
  dump.add(heatf_eng_e_x, "heatf_eng_e_x",1);
  dump.add(gamma_cond_i_x, "gamma_cond_i_x",1);
  dump.add(gamma_exb_i_x, "gamma_exb_i_x",1);
  dump.add(gamma_dia_i_x, "gamma_dia_i_x",1);
  dump.add(Ti_grad, "Ti_grad",1);
  dump.add(Te_grad, "Te_grad",1);
  dump.add(M_iol0, "M_iol0",1);
  dump.add(F_iol0, "F_iol0",1);
  dump.add(E_iol0, "E_iol0",1);
  dump.add(V_iol, "V_iol",1);
  return(0);
}

int physics_run(BoutReal t)
{
  // Communicate variables
  mesh->communicate(Ni, Vi, Te, Ti);
  //mesh->communicate(Nn,Tn,Vn);
  mesh->communicate(Nn);
  mesh->communicate(Nm,Vmx,Vm);
  // NB: Intermediate variables calculated with Grad operators are all necessary to be communicated
  // after being calculated

//**********************************************************li2016 
  mesh->communicate(phi01,U00);
  phi01.applyBoundary();
  U00.applyBoundary();
//**********************************************************liend

  Ni.applyBoundary();
  Vi.applyBoundary();
  Te.applyBoundary();  
  Ti.applyBoundary();
  Nn.applyBoundary();
  Tn.applyBoundary();
  Vn.applyBoundary();
  Nm.applyBoundary();
  Vm.applyBoundary();
/* if(energy_flux_bndry)
   { 
     if(mesh->firstX())
      {
      for ( int jx=0; jx < 2 ; jx++)
       for ( int jy=0; jy < mesh->ngy ; jy++)
       	 for ( int jz=0; jz < mesh->ngz ; jz++)
	   {
 		Ti[jx][jy][jz]=Ti[jx+1][jy][jz]-q_input*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy])/Ni[jx][jy][jz]/chi_i_perp[jx][jy][jz];
                Te[jx][jy][jz]=Te[jx+1][jy][jz]-q_input*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy])/Ni[jx][jy][jz]/chi_e_perp[jx][jy][jz];
	   }
      }
   }*/  
  //smooth noisies
  if(nlfilter_noisy_data)
    {
      //Ni=nl_filter(Ni,filter_para);
      //Vi=nl_filter(Vi,filter_para);
      Nn=nl_filter(Nn,filter_para);
      Vn=nl_filter_y(Vn,filter_para);
    }
  
  //*****@!!@*****
  // NB: Any value re-assignment should be given HERE ONLY!
  //*****@!!@***** 

  temp_Ni=field_larger(Ni,minimum_val);
  temp_Nn=field_larger(Nn,minimum_val);   // in case divided by zero 
  temp_Nm=field_larger(Nm,minimum_val);   //1/temp_N only used to replay 1/N
  temp_Te=field_larger(Te,minimum_val);   // necessary for sheath BC
  temp_Ti=field_larger(Ti,minimum_val);   // necessary for sheath BC

  Ni=field_larger(Ni,minimum_val);  
  Nn=field_larger(Nn,minimum_val);  
  Nm=field_larger(Nm,minimum_val);
  //Ti=field_larger(Ti,minimum_val); 
  //Te=field_larger(Te,minimum_val);  
 
  Tn=temp_Ti;

  Vmx=Vm.x;
  if(SMBI_LFS)
  {
  Nm=ret_const_flux_BC(Nm, Nm0);
  Vmx=ret_const_flux_BC(Vmx, Vm0); 
  }
  else
  {
  Nm=ret_const_flux_BC(Nm,minimum_val);
  Vmx=ret_const_flux_BC(Vmx,0.0);
  }

  // Update non-linear coefficients on the mesh
  //kappa_Te = 3.2*(1./fmei)*(wci/nueix)*(Tet^2.5);
  //kappa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);
  //kappa_Te = 3.2*(1./fmei)*(1./tbar/nueix)*(temp_Te^2.5);
  //kappa_Ti = 3.9*(1./tbar/nuiix)*(temp_Ti^2.5);  

  // Set flux limit for kappa
  V_th_e= 4.19e5*sqrt(temp_Te*Te_x)*tbar/Lbar;
  V_th_i= 9.79e3*sqrt(temp_Ti*Te_x/AA)*tbar/Lbar;

  nu_ei = nueix*temp_Ni/(temp_Te*sqrt(temp_Te));
  nu_ii = nuiix*temp_Ni/(temp_Ti*sqrt(temp_Ti));
  kappa_Te = 3.2*V_th_e*V_th_e/nu_ei/tbar*temp_Ni;
  kappa_Ti = 3.9*V_th_i*V_th_i/nu_ii/tbar*temp_Ni;  
  kappa_Te_fl = V_th_e*q95*Lbar*temp_Ni;  // Ne=Ni quasineutral
  kappa_Ti_fl = V_th_i*q95*Lbar*temp_Ni;

  kappa_Te *= kappa_Te_fl/(kappa_Te_fl+kappa_Te);
  kappa_Ti *= kappa_Ti_fl/(kappa_Ti_fl+kappa_Ti);
  mesh->communicate(kappa_Te);
  (kappa_Te).applyBoundary("neumann");
  mesh->communicate(kappa_Ti);
  (kappa_Ti).applyBoundary("neumann");

  // Thermal Speed normalized in ion thermal speed V_thi
  c_se = sqrt((temp_Te+temp_Ti)/Mi); 
  mesh->communicate(c_se);
  (c_se).applyBoundary("neumann");

  //parallel heat fluxes to calculate Te,i gradients at SBC
  if (Secondary_electron_emission) 
    {
     q_se_kappa=-c_se*temp_Ni*temp_Te*(2.0/(1.0-delta)+1.0/2.0*log((1.-delta)*(1.-delta)*Mi2/Me/2./PI/(1.+temp_Ti/temp_Te)))/kappa_Te;
     q_si_kappa=-c_se*temp_Ni*2.0*temp_Ti/kappa_Ti;
    }
  else
    {
     q_se_kappa = -C_fe_sheat*temp_Ni*temp_Te*c_se/kappa_Te;         // '-' means out-flowing to Sheath
     q_si_kappa = -C_fi_sheat*temp_Ni*temp_Ti*c_se/kappa_Ti;
    }   
     mesh->communicate(q_se_kappa);
     (q_se_kappa).applyBoundary("neumann");
     mesh->communicate(q_si_kappa);
     (q_si_kappa).applyBoundary("neumann");

   
  if (J1_all)
    { 
     J1_para =  -1.0/(eta*sbc_lambda1*B0)*Grad_par(phi01)+((KB*eV_K*Te_x*tbar)/(eta*Lbar*Lbar*ee*Ni*B0*Bbar))*Grad_par((Te*Ni))+((0.71*KB*eV_K*Te_x*tbar)/(eta*Lbar*Lbar*B0*Bbar*ee))*Grad_par(Te);
     mesh->communicate(J1_para);
     J1_para.applyBoundary("neumann");
     Ve = Vi - J1_para*B0*Bbar*tbar/(MU0*Ni*Ni_x*density_unit*ee*Lbar*Lbar); 
     mesh->communicate(Ve);
     Ve.applyBoundary("neumann");
    }
  if (iterative_phi0)
    {  
     phish0 = sbc_lambda1*0.5*log((Mi2/Me)/(4.0*PI)*(2*temp_Te*Te_x/(temp_Ti*Ti_x+temp_Te*Te_x)))*(temp_Te*Te_x*KB*eV_K)/ee;
     if (Secondary_electron_emission) 
       {
        phish0 = sbc_lambda1*0.5*log((Mi2/Me)*(1-delta)*(1-delta)/(4.0*PI)*(2*temp_Te*Te_x/(temp_Ti*Ti_x+temp_Te*Te_x)))*(temp_Te*Te_x*KB*eV_K)/ee;
       }
     if (J_para0) 
       {
        phish0 -= sbc_lambda1*(temp_Te*Te_x*KB*eV_K/ee)*log(1-J1_para*B0*Bbar/(MU0*Lbar)/(temp_Ni*Ni_x*density_unit*ee*sqrt((temp_Ti*Ti_x+temp_Te*Te_x)*eV_K/Mi2)));
       }
     phish0 /= Lbar/tbar*Lbar*Bbar;
     mesh->communicate(phish0);
     (phish0).applyBoundary("neumann");
      
     phish1 = phish0*RF_coef;
     // SAVE_ONCE(phish0);
     // SAVE_ONCE(phish1);
     (ersh0.x).setLocation(CELL_CENTRE);
     //(ersh0.y).setLocation(CELL_CENTRE);
     ersh0 = -1.0*Grad(phish0);
     mesh->communicate(ersh0.x);
     //mesh->communicate(ersh0.y);
     (ersh0.x).applyBoundary("neumann");
     //(ersh0.y).applyBoundary("neumann");
     ersh0.toCovariant();
     ersh0.x *= Rxy*Bpxy;

     (ersh1.x).setLocation(CELL_CENTRE);
     //(ersh1.y).setLocation(CELL_CENTRE);
     ersh1 = -1.0*Grad(phish1);
     mesh->communicate(ersh1.x);
     //mesh->communicate(ersh1.y);
     (ersh1.x).applyBoundary("neumann");
     //(ersh1.y).applyBoundary("neumann");
     ersh1.toCovariant();
     ersh1.x *= Rxy*Bpxy;
   }

  // Apply Sheath Boundary Condition at outside of Separitrix 
  // NB: SBC applied at ydown (theta=0) and yup (theta=2PI) are different due to fluxes flowing 
  // towards X point and hitting on the plates, thus SBC_ydown(var,-value) while SBC_yup(var,value) 
  // If particle recycling, Ni flux Gamma_ni outflows to divertor plates while Nn flux Gamma_nn inflows 
  // from the plates, Gamma_nn = Gamma_ni * Rate_recycle 
  
  if (Sheath_BC)
    {
     SBC_Dirichlet(Vi, c_se); 
     //if(Ni_Zero_Grad) {SBC_Gradpar(Ni, 0.);}
     SBC_Gradpar(Te, q_se_kappa);
     Te_grad = Grad_par(Te);  //unit W/(m^-2)
     (Te_grad).applyBoundary("neumann");
     SBC_Gradpar(Ti, q_si_kappa);
     Ti_grad = Grad_par(Ti);  //unit W/(m^-2)
     (Ti_grad).applyBoundary("neumann");
//***************************************li2016
     if(Sheath_BC_phi) {SBC_Dirichlet_phi(phi01,phish0);}
     if(Ni_Zero_Grad) {SBC_Gradpar(Ni, 0.);}
     // Grad2_par2_sh(phi01,phish1);
     if(Ni_Constant_Grad)
       {
        for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) 
         {
          for(int jz=0; jz<mesh->ngz; jz++) {
            // Free boundary (constant gradient) density
            Ni(r.ind, mesh->ystart-1,jz) = 2.0*Ni(r.ind, mesh->ystart, jz) - Ni(r.ind, mesh->ystart+1, jz);}
          }
        for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) 
         {
          for(int jz=0; jz<mesh->ngz; jz++) {
            // Free boundary (constant gradient) density
            Ni(r.ind, mesh->yend+1,jz) = 2.0*Ni(r.ind, mesh->yend, jz) - Ni(r.ind, mesh->yend-1, jz);}
          }
        }
//***************************************liend
     }  //end of if(Sheath_BC)


  // no flux limitation for kappa at Sheath BC 
  //kappa_Te *= kappa_Te_fl/(kappa_Te_fl+kappa_Te);
  //kappa_Ti *= kappa_Ti_fl/(kappa_Ti_fl+kappa_Ti);

  // Collisional time ion-electrons
  tau_ei = 1./tbar/(nueix*temp_Ni/(temp_Te^1.5));

  //ion viscosity
  eta0_i=0.96*(1./tbar/nuiix)*temp_Ni*temp_Ti;
  eta0_n=0.96*(1./tbar/nuiix)*temp_Nn*temp_Ti;
  // eta0_i=0.96*(1./tbar/nuiix)*(temp_Ti^2.5);
  // eta0_n=0.96*(1./tbar/nuiix)*(Tn^2.5)*Nn/temp_Nn;
  
  //updata collisional rates 
 
  // Ionization rate  depending on Nn for Plasmas 
  //nu_ionz = 3.e4*tbar*temp_Nn*Ni_x*temp_Te*temp_Te*Te_x*Te_x/(3.+0.01*temp_Te*temp_Te*Te_x*Te_x);   //Ni_x in unit 1e^20 
  // Ionization rate  depending on Nn for Plasmas
  if(terms_ionization)
   {
    for (int jx=0;jx<mesh->ngx;jx++)
     {
      for (int jy=0;jy<mesh->ngy;jy++)
         {
          for (int jz=0;jz<mesh->ngz;jz++)
              {
  nu_ionz[jx][jy][jz] = 3.e6*tbar*temp_Nn[jx][jy][jz]*Ni_x*(3.854*exp(1.268e-5-(temp_Te[jx][jy][jz]*Te_x+3.13)/6703.89)*erfc(3.561e-3-(temp_Te[jx][jy][jz]*Te_x+3.13)/47.745)-4.34);   //Ni_x in unit 1e^20
              }
          }
       }
    }
  nu_CX =  1.e6*tbar*temp_Nn*Ni_x*(1.7+1.9*(((1.5*temp_Ti*Te_x)^0.333)-2.4662)/(((150.0*temp_Ti*Te_x)^0.333)-2.4662));      // empirical formula 
  
  //Diag_neg_value(nu_ionz,nu_CX,Nn,Nm);

  // Ionization rate  depending on Ni for Neutrals 
  nu_ionz_n = 3.e4*tbar*temp_Ni*Ni_x*temp_Te*temp_Te*Te_x*Te_x/(3.+0.01*temp_Te*temp_Te*Te_x*Te_x);    //Ni_x in unit 1e^20 
  //nu_CX_n = 1.e5*tbar*temp_Ni*Ni_x*(1.7+0.667*((1.5*temp_Ti*Te_x)^0.333-2.466));      // empirical formula 
  nu_CX_n =  1.e6*tbar*temp_Ni*Ni_x*(1.7+1.9*(((1.5*temp_Ti*Te_x)^0.333)-2.4662)/(((150.0*temp_Ti*Te_x)^0.333)-2.4662));      // empirical formula 

 //Dissociation rate of molecules
  nu_diss = 3.e4*tbar*temp_Nm*Ni_x*temp_Te*temp_Te*Te_x*Te_x/(3.+0.01*temp_Te*temp_Te*Te_x*Te_x);     //need be corrected
  //nu_diss =0.;

  // recombination rate of ions and electrons
  Field3D lambda_rec=1.5789e5/(temp_Te*Te_x*1.1604e4);  //Te trasfer unit from eV to Kelvin  
  nu_rec = tbar*temp_Ni*Ni_x*5.197*ZZ*sqrt(lambda_rec)*(0.4288+0.5*log(lambda_rec)+0.469/(lambda_rec)^0.333); //Seaton M F 1959 'Radiative recombination of hydrogenic ions'

  // source terms
  Si_p = temp_Ni*nu_ionz;
  Si_CX = temp_Ni*nu_CX;
  S_diss = temp_Ni*nu_diss;  
  S_rec = temp_Ni*nu_rec;

  //Set Boundary and Apply Boundary for output variables
  tau_ei.applyBoundary();
  nu_ionz.applyBoundary();
  nu_CX.applyBoundary();
  nu_ionz_n.applyBoundary();
  nu_CX_n.applyBoundary();
  Si_p.applyBoundary();
  Si_CX.applyBoundary();
  S_diss.applyBoundary();

  if(t==0.)  output.write("\n **** t %e  nu_I %e nu_I_n %e Si_p %e tau_ei %e ****\n",t,max(nu_ionz),max(nu_ionz_n),max(Si_p),max(tau_ei)); 
  
  // Perpensicular 3D Diffusion coefficients of Plasmas, set constant by default
  Diff_ni_perp = Diffc_ni_perp;        
  chi_i_perp = chic_i_perp;
  chi_e_perp = chic_e_perp;
  if(diffusion_coef_step_function)
    {
      Diff_ni_perp = Diffc_ni_step;        
      chi_i_perp = chic_i_step;
      chi_e_perp = chic_e_step;
     }
  if(diffusion_coef_Hmode && initial_profile_Hmode)
    {
      Diff_ni_perp = Diffc_ni_Hmode;        
      chi_i_perp = chic_i_Hmode;
      chi_e_perp = chic_e_Hmode;
     } 
  if(load_experiment_profiles)
    {
      Diff_ni_perp = Diffc_ni_Hmode;        
      chi_i_perp = chic_i_Hmode;
      chi_e_perp = chic_e_Hmode;
     }
  if(load_grid_trans)
    {
      Diff_ni_perp = Diff_grid;        
      chi_i_perp = Chi_grid;
      chi_e_perp = Che_grid;
     }
  (Diff_ni_perp).applyBoundary("neumann");
  (chi_i_perp).applyBoundary("neumann");
  (chi_e_perp).applyBoundary("neumann");

  //parallel particle fluxes to calculate Ni,Nn gradients at SBC
  //Diffc_nn_perp = Tn/Mn/nu_CX_n;   
  Diffc_nn_perp = Tn/Mn/(nu_CX_n+2.*S_diss/temp_Nn);
  if(terms_recombination) Diffc_nn_perp = Tn/Mn/(nu_CX_n+2.*S_diss/temp_Nn+S_rec/temp_Nn);
  V_th_n=sqrt(Tn/Mn);
  Diffc_nn_perp_fl=V_th_n*Lnn_min/Lbar;
  Diffc_nn_perp *= Diffc_nn_perp_fl/(Diffc_nn_perp+Diffc_nn_perp_fl);    // Diffc_nn_perp calculated once for SBC and again later
  
  Diffc_nn_par = Tn/Mn/(nu_CX_n+2.*S_diss/temp_Nn);
  if(terms_recombination) Diffc_nn_par = Tn/Mn/(nu_CX_n+2.*S_diss/temp_Nn+S_rec/temp_Nn);
  Diffc_nn_par *= Diffc_nn_perp_fl/(Diffc_nn_par+Diffc_nn_perp_fl); 

  Gamma_ni_wall = Diff_ni_perp*Ni/Lni_wall;
  Gamma_nn_wall = Rate_recycle_wall*Gamma_ni_wall;
  Gamma_nnw_Diffc = Gamma_nn_wall/Diffc_nn_perp;
  Gamma_ni = Ni*c_se;
  Gamma_nn = Rate_recycle_div*Gamma_ni;

  Gamma_ni_Diffc= -Gamma_ni/(Diff_ni_perp+minimum_val);              // '-' means out-flowing to Sheath
  Gamma_nn_Diffc=  Gamma_nn/Diffc_nn_perp;                           // '+' in-flowing from Sheath

  if(Wall_particle_recycle)
    {
     //WallBC_Xout_GradX(Ni,Gamma_ni_Diffc);                         // '-' means out-flowing to wall
     WallBC_Xout_GradX_len(Ni, -1./Lni_wall);                        //fixed Gradient length (real unit m) at wall
     WallBC_Xout_GradX(Nn,Gamma_nnw_Diffc);                          // '+' means in-flowing from wall
     } 
  if(SBC_particle_recycle)
     {      
      /*
      Vn_perp = Diffc_nn_perp*abs(DDX(Nn))*Rxy*abs(Bpxy)/temp_Nn;
      mesh->communicate(Vn_perp);
      Vn_perp.applyBoundary();
      Vn_th_plate = sqrt(Tn_plate/Mn);
      Vn_perp = field_larger(Vn_perp,Vn_th_plate);
      //Vn_pol = -abs(Vn)*sin(angle_B_plate)+abs(Vn_perp)*cos(angle_B_plate);
      Vn_pol = abs(Vn_perp)*cos(angle_B_plate);   //Assumption Vn_perp >> Vn_par
      Vn_pol = field_larger(Vn_pol,Vn_th_plate);
      if(nlfilter_Vnpol)nl_filter(Vn_pol);
    
      SBC_Dirichlet(Nn, Gamma_nn/Vn_pol);
      */       
      Grad_par_Nn=Gamma_nn/Diffc_nn_par;
      mesh->communicate(Grad_par_Nn);
      (Grad_par_Nn).applyBoundary("neumann");
      SBC_Gradpar(Nn, Grad_par_Nn);

      }

if(leakage_bndry)
     {
      if ( mesh->lastX()){
         for ( int jx=mesh->ngx-2; jx < mesh->ngx ; jx++)
             for ( int jy=0; jy < mesh->ngy ; jy++)
	       for ( int jz=0; jz < mesh->ngz ; jz++){
                 Ni[jx][jy][jz]=(1-leak_ni*sqrt((Te[jx-1][jy][jz]+Ti[jx-1][jy][jz])/Mi)/Diff_ni_perp[jx-1][jy][jz]*mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy]))*Ni[jx-1][jy][jz];
   		 Ti[jx][jy][jz]=(1-leak_ti*sqrt((Te[jx-1][jy][jz]+Ti[jx-1][jy][jz])/Mi)/chi_i_perp[jx-1][jy][jz]*mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy]))*Ti[jx-1][jy][jz];
	         Te[jx][jy][jz]=(1-leak_te*sqrt(Te[jx-1][jy][jz]/Me*Mi2)/chi_e_perp[jx-1][jy][jz]*mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy]))*Te[jx-1][jy][jz];	 
	  } 
	 } 
     }

if(decay_bndry)
	{
	 if(mesh->lastX()){
	    for ( int jx=mesh->ngx-2; jx < mesh->ngx ; jx++)
             for ( int jy=0; jy < mesh->ngy ; jy++)
               for ( int jz=0; jz < mesh->ngz ; jz++){
		Ni[jx][jy][jz]=(1-mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy])/decay_bndry_legth*Lbar)*Ni[jx-1][jy][jz];
		Ti[jx][jy][jz]=(1-mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy])/decay_bndry_legth*Lbar)*Ti[jx-1][jy][jz];
		Te[jx][jy][jz]=(1-mesh->dx[jx-1][jy]/sqrt(mesh->g11[jx-1][jy])/decay_bndry_legth*Lbar)*Te[jx-1][jy][jz];
		}}}


   if(energy_flux_bndry)
   {
     if(mesh->firstX())
      {
      for ( int jx=0; jx < 2 ; jx++)
       for ( int jy=0; jy < mesh->ngy ; jy++)
         for ( int jz=0; jz < mesh->ngz ; jz++)
           {
		 int indy = mesh->YGLOBAL(jy);
          //output.write("dinx: %e   indy: %e  jysep1: %e  jysep2: %e\n", dindx, indy, jysep1, jysep2);
           if ((indy > jysep1) && (indy < jysep2))
	     {
                Ti[jx][jy][jz]=Ti[jx+1][jy][jz]+1.0/chi_coef*q_input*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy])/Ni[jx][jy][jz]/chi_i_perp[jx][jy][jz];
                Te[jx][jy][jz]=Te[jx+1][jy][jz]+1.0/chi_coef*q_input*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy])/Ni[jx][jy][jz]/chi_e_perp[jx][jy][jz];
              }}
      }
   }
//******Calculate minmum velocity for Ion orbit loss********//
  V_iol1 = 0.0;
  F_iol1 = 0.0;
  M_iol1 = 0.0;
  E_iol1 = 0.0;
/*  if (terms_IOL)
    {
     GlobalField3D gphi0(mesh);
     GlobalField2D Bp(mesh),Bt(mesh),gRxy(mesh),gpsixy(mesh),gTi(mesh);
     Bp.gather(Bpxy);
     Bt.gather(Btxy);
     gpsixy.gather(psixy);
     gRxy.gather(Rxy);
     gphi0.gather(phi00);
     gTi.gather(Ti_exp);
     GlobalField2D gM_iol0(mesh),gF_iol0(mesh),gE_iol0(mesh);
     GlobalField3D gV_iol(mesh);
     int ngx =gRxy.xSize();
     int ngy =gRxy.ySize();
     int ngz =gV_iol.ySize();
     BoutReal FS0[xout][ngy],FS1[ngy],fphi0[xout][ngy],fphi1[ngy],fphi[ngx][ngy],B_field[ngx][ngy],B_ratio[xout][ngy][ngy],phi_diff[xout][ngy][ngy],FS[ngx][ngy],R_var[ngx][ngy],phi_var[ngx][ngy][ngz];
     BoutReal consi[P_N],dconsi;
     BoutReal Emin_mid[xout][P_N][ngy],Ti_var[ngx][ngy][ngz];
     BoutReal Emin_final[xout][P_N];
     BoutReal Emin[xout][P_N][ngy][ngy],V0_min[xout][P_N][ngy][ngy],v0p[xout][P_N][ngy][ngy],v0m[xout][P_N][ngy][ngy],ke0m[xout][P_N][ngy][ngy],ke0p[xout][P_N][ngy][ngy];
     BoutReal M_iol[xout][ngy][ngy],F_iol[xout][ngy][ngy],E_iol[xout][ngy][ngy],dF_iol[xout][P_N][ngy][ngy],dM_iol[xout][P_N][ngy][ngy],dE_iol[xout][P_N][ngy][ngy];
     BoutReal Emin_final1[xout][P_N],Gamma_m[xout][P_N],Gamma_f[xout][P_N],Gamma_e[xout][P_N];
   
     if (gRxy.dataIsLocal())
       {
        // Data is available on this processor
        for(int i=0;i<ngx;i++)
          {
           for(int j=0;j<ngy;j++)
            {
             FS[i][j] = gpsixy(i,j);
             fphi[i][j] = 1./sqrt(1.+(Bp(i,j)*Bp(i,j))/(Bt(i,j)*Bt(i,j)));
             B_field[i][j] = sqrt(Bp(i,j)*Bp(i,j)+Bt(i,j)*Bt(i,j));
             R_var[i][j] = gRxy(i,j);
             for(int z=0;z<ngz;z++)
              {
               phi_var[i][j][z] = gphi0(i,j,z);
                }
              }
            }
         for(int m=0;m<xout;m++)
           {
            for(int i=0;i<ngy;i++)
              {
               FS0[m][i]=FS[m][i];       //internal magnetic surface
               fphi0[m][i]=fphi[m][i];
               for(int j=0;j<ngy;j++)
                {
                 FS1[j]=FS[xout][j];       //external magnetic surface
                 fphi1[j]=fphi[xout][j];
                 B_ratio[m][i][j] = B_field[xout][j]/B_field[m][i];
                 for(int z=0;z<ngz;z++)
                  {
                   phi_diff[m][i][j] =  phi_var[m][i][z]-phi_var[xout][j][z];
                   }
                  }
                }
              }
           //output.write("\t consi=  \n");
          for(int k=0;k<P_N;k++)         //Define a variable of pitch angle
            {
             dconsi = 2*0.9999/(P_N-1);
             consi[k] = -0.9999+dconsi*k;
             //output.write("\t %e  \n",consi[k]);
             }
          for(int m=0;m<xout;m++)
            {
             for(int k=0;k<24;k++)
              {
               for(int i=0;i<ngy;i++)
                {
                 for(int j=0;j<ngy;j++)
                  {
                   BoutReal a,b,c;
                   a = (B_ratio[m][i][j]*fphi0[m][i]/fphi1[j]*consi[k])*(B_ratio[m][i][j]*fphi0[m][i]/fphi1[j]*consi[k])-1.0+(1.0-consi[k]*consi[k])*B_ratio[m][i][j];
                   b = (2.0*ee*tbar*bmag*(FS0[m][i]-FS1[j])/(R_var[xout][j]*Mi2*fphi1[j]))*(B_ratio[m][i][j]*fphi0[m][i]/fphi1[j]*consi[k]);
                   c = (ee*tbar*bmag*(FS0[m][i]-FS1[j])/(R_var[xout][j]*Mi2*fphi1[j])*ee*tbar*bmag*(FS0[m][i]-FS1[j])/(R_var[xout][j]*Mi2*fphi1[j]))-2.0*ee*tbar*bmag/Mi2*phi_diff[m][i][j];
                   BoutReal abc = b*b-4.0*a*c;
                   if (abc < 0.0)
                    {
                     v0p[m][k][i][j] = 0.0;
                     v0m[m][k][i][j] = 0.0;
                      }
                   else
                    {
                     v0p[m][k][i][j] = (-b+sqrt(abc))/(2.*a);
                     v0m[m][k][i][j] = (-b-sqrt(abc))/(2.*a);
                      }
*/
//******Minimum reduced energy depends on mass*************//
/*
                   ke0m[m][k][i][j] = 0.5*Mi2*(v0m[m][k][i][j]*v0m[m][k][i][j])*Vi_x*Vi_x/ee;
                   ke0p[m][k][i][j] = 0.5*Mi2*(v0p[m][k][i][j]*v0p[m][k][i][j])*Vi_x*Vi_x/ee;
                   if(v0p[m][k][i][j] >= v0m[m][k][i][j] && v0m[m][k][i][j] >= 0.0)// && Emin[m][i][j][k] >=ke0m[m][i][j][k])
                     {
                      Emin[m][k][i][j] = ke0m[m][k][i][j];
                      V0_min[m][k][i][j] = v0m[m][k][i][j];
                       }
                   else if(v0m[m][k][i][j] >= v0p[m][k][i][j] && v0p[m][k][i][j] >= 0.0)// && Emin[m][i][j][k] >=ke0p[m][i][j][k])
                     {
                      Emin[m][k][i][j] = ke0p[m][k][i][j];
                      V0_min[m][k][i][j] = v0p[m][k][i][j];
                       }
                    else if(v0m[m][k][i][j] >= 0.0 && v0p[m][k][i][j] < 0.0)// && Emin[m][i][j][k] >=ke0m[m][i][j][k])
                      {
                       Emin[m][k][i][j] = ke0m[m][k][i][j];
                       V0_min[m][k][i][j] = v0m[m][k][i][j];
                        }
                     else if(v0p[m][k][i][j] >= 0.0 && v0m[m][k][i][j] < 0.0)// && Emin[m][i][j][k] >=ke0p[m][i][j][k])
                       {
                        Emin[m][k][i][j] = ke0p[m][k][i][j];
                        V0_min[m][k][i][j] = v0p[m][k][i][j];
                        }
                      }
                   }
                }
             }
          for(int m=0;m<xout;m++)
            {
             for(int k=0;k<P_N;k++)
              {
               for(int i=0;i<ngy;i++)
                {
                 for(int j=0;j<ngy;j++)
                  {
                   if(Emin[m][k][i][j] < Emin[m][k][i+1][j])
                    { Emin_mid[m][k][i] = Emin[m][k][i][j];}
                   else {Emin_mid[m][k][i] = Emin[m][k][i+1][j];}
                   }
                 }
               }
             }
           //output.write("\t  Emin_final= \n");
           for(int m=0;m<xout;m++)
            {
             for(int k=0;k<P_N;k++)
              {
               for(int i=0;i<ngy;i++)
                {
                 if(Emin_mid[m][k][i] > Emin_mid[m][k][i+1])
                  {Emin_final[m][k] = Emin_mid[m][k][i];}
                 else {Emin_final[m][k] = Emin_mid[m][k][i+1];}
                 }
                 //output.write("\t m=%d k=%d Emin_final=%e \n",m,k, Emin_final[m][k]);
                }
              }
           for(int m=0;m<ngx;m++)
            {
             for(int j=0;j<ngy;j++)
              {
               gM_iol0(m,j) = 0.0;
               gF_iol0(m,j) = 0.0;
               gE_iol0(m,j) = 0.0;
               }
              }
           for(int m=0;m<xout;m++)
            {
             for(int k=0;k<P_N;k++)
              {
               Emin_final1[m][k] = Emin_final[m][k]/(Ti_x*gTi(m,38));
               Gamma_m[m][k] = gsl_sf_gamma_inc(2.0, Emin_final1[m][k]);
               Gamma_f[m][k] = gsl_sf_gamma_inc(1.5, Emin_final1[m][k]);
               Gamma_e[m][k] = gsl_sf_gamma_inc(2.5, Emin_final1[m][k]);
                }
              }
           for(int j=0;j<ngy;j++)
            {
             for(int m=0;m<xout;m++)
              {
               for(int k=0;k<P_N-1;k++)
                {
                 gM_iol0(m,j) += 1.0/2.0*dconsi*consi[k+1]*(Gamma_m[m][k]+Gamma_m[m][k+1])/2.0;
                 gF_iol0(m,j) += 1.0/2.0*dconsi*(Gamma_f[m][k]+Gamma_f[m][k+1])/2.0/tgamma(1.5);
                 gE_iol0(m,j) += 1.0/2.0*dconsi*(Gamma_e[m][k]+Gamma_e[m][k+1])/2.0/tgamma(2.5);
                  }
                }
              }
           for(int m=0;m<xout;m++)
            {
             for(int i=0;i<ngy;i++)
              {
               for(int j=0;j<ngy;j++)
                {
                 for(int k=0;k<P_N-1;k++)
                  {
                   Emin[m][k][i][j] = Emin[m][k][i][j]/(gTi(m,i)*Ti_x);
                   M_iol[m][i][j] += 1.0/2.0*dconsi*(consi[k+1]*gsl_sf_gamma_inc(2.0, Emin[m][k][i][j])+consi[k+1]*gsl_sf_gamma_inc(2.0, Emin[m][k+1][i][j]))/2.0/tgamma(2.0);
                   F_iol[m][i][j] += 1.0/2.0*dconsi*(gsl_sf_gamma_inc(1.5, Emin[m][k][i][j])+gsl_sf_gamma_inc(1.5, Emin[m][k+1][i][j]))/2.0/tgamma(1.5);
                   E_iol[m][i][j] += 1.0/2.0*dconsi*(gsl_sf_gamma_inc(2.5, Emin[m][k][i][j])+gsl_sf_gamma_inc(2.5, Emin[m][k+1][i][j]))/2.0/tgamma(2.5);
                   dF_iol[m][k][i][j] = dconsi*(gsl_sf_gamma_inc(1.5, Emin[m][k][i][j])+gsl_sf_gamma_inc(1.5, Emin[m][k+1][i][j]))/2.0/tgamma(1.5);
                   dM_iol[m][k][i][j] = dconsi*(consi[k]*gsl_sf_gamma_inc(2.0, Emin[m][k][i][j])+consi[k+1]*gsl_sf_gamma_inc(2.0, Emin[m][k+1][i][j]))/2.0/tgamma(2.0);
                   dE_iol[m][k][i][j] = dconsi*(gsl_sf_gamma_inc(2.5, Emin[m][k][i][j])+gsl_sf_gamma_inc(2.5, Emin[m][k+1][i][j]))/2.0/tgamma(2.5);
                    }
                  }
                }
              }
           for(int j=0;j<gV_iol.ySize();j++)
            {
             for(int k=0;k<gV_iol.zSize();k++)
              {
               for(int m=0;m<ngx;m++)
                {
                 if(m<xout)
                  {
                   gV_iol(m,j,k) = 2.0/sqrt(PI)*gM_iol0(m,j)*sqrt(2.0*ee*gTi(m,j)*Ti_x/Mi2);
                    }
                  else
                   {
                    gV_iol(m,j,k) = 0.0;
                     }
                   }
                 }
               }
           }
   M_iol0 = gM_iol0.scatter();
   F_iol0 = gF_iol0.scatter();
   E_iol0 = gE_iol0.scatter();
   V_iol = gV_iol.scatter();
    
   V_iol = V_iol/Vi_x;
   //F_iol1 = 0.0;
   //M_iol1 = 0.0;
   //E_iol1 = 0.0;
   //if(terms_IOL)
   // {
   F_iol1 = F_iol0;
   M_iol1 = M_iol0;
   E_iol1 = E_iol0;
   //  }
  }
*/
//******Calculate minimum velocity for IOL end*************//

//************************************
  // DENSITY EQUATION
//************************************
  // output.write("Now updata Ni \n"); 
  // ddt(Ni) = 0.; 

  Pe = Te*Ni;
  Pi = Ti*Ni;
  if(!Solving_Eq_Ni)
    {ddt(Ni) = 0.;}
  else
    { 
     ddt(Ni) =  (1.0-F_iol1)*Diff_ni_perp*Delp2(Ni)
     //ddt(Ni) = Diff_grid*Delp2(Ni)
              - Ni*Grad_par(Vi)
              - Vpar_Grad_par(Vi,Ni)
              + Si_p
            ; 
     if(terms_recombination) ddt(Ni) -= S_rec;
     if(terms_exb) ddt(Ni) -= (1.0-F_iol1)*bracket((phi01/sbc_lambda1), Ni, bm_exb);       //ExB drift
     if(curvature_phi) ddt(Ni) -= (1.0-F_iol1)*2.0*Ni*b0xcv*Grad((phi01/sbc_lambda1))/B0;  //ExB drift
     if(diamag)  ddt(Ni) -= (1.0-F_iol1)*2.0*(KB*Ti_x*eV_K/(Zi*ee*Bbar*Lbar*Lbar/tbar))*b0xcv*Grad(Pi)/B0;  //magnetic drift
      
     if(external_sources) ddt(Ni) += S_pi_ext;
     if(terms_Gradperp_diffcoefs) 
       {
        grad_perp_Diffi=Grad_perp(test_f1*Diff_ni_perp);
        grad_perp_Diffi.applyBoundary();
        mesh->communicate(grad_perp_Diffi);  
        ddt(Ni) +=  (1.0-F_iol1)*V_dot_Grad(grad_perp_Diffi,Ni); 
       }
     }

//************************************
  // ELECTRON TEMPERATURE   ---Te---
//************************************

  //output.write("Now updata Te \n");
 
  //ddt(Te)=0.0;
  
  Vector3D V_diff = Diff_ni_perp*Grad_perp(Ni)/Ni;  //diffusitive velocity
  mesh->communicate(V_diff.x); 
  mesh->communicate(V_diff.y); 
  mesh->communicate(V_diff.z); 
  (V_diff.z).applyBoundary("neumann");
  (V_diff.y).applyBoundary("neumann");
  (V_diff.x).applyBoundary("neumann");

  if(!Solving_Eq_Te)
    {ddt(Te)=0.0;} 
  else
    {
     ddt(Te) =   0.6667*kappa_Te*Grad2_par2(Te)/temp_Ni
               //+ 0.6667*Grad_par(kappa_Te)*Grad_par(Te)/temp_Ni
               + 0.6667*chi_e_perp*Delp2(Te)
               - nu_ionz*(Te+0.6667*W_ionz)
               //- 0.6667*nu_diss*(W_diss+W_bind)
               - 2.*fmei*(Te-Ti)/tau_ei
               ;
     if(terms_Ve) ddt(Te) -= (Vpar_Grad_par(Ve,Te) + 0.6667*Te*Grad_par(Ve));
     if(terms_radial_conv) ddt(Te) += 0.6667*Te*Div(V_diff)+Diff_ni_perp*Grad_perp(Ni)*Grad_perp(Te)/Ni;
     if(terms_spitzer_resist) ddt(Te) += 2.0/3.0*(Bbar*Bbar/(MU0*Ni_x*density_unit*Te_x*KB*eV_K))*eta*B0*B0*J1_para*J1_para/Ni;
     if(thermal_force) ddt(Te) += 0.71*2.0/3.0*(tbar*Bbar/(Ni_x*density_unit*ee*Lbar*Lbar*MU0))*Te*B0/Ni*Grad_par(J1_para);
     if(terms_recombination) ddt(Te) += nu_rec*W_rec;
     if(terms_exb) ddt(Te) -= bracket((phi01/sbc_lambda1), Te, bm_exb);           //ExB drift
     if(curvature_phi) ddt(Te) -= 4.0/3.0*Te*b0xcv*Grad((phi01/sbc_lambda1))/B0;  //EXB drift 
     if(diamag)  ddt(Te) += 4.0/3.0*(KB*Te_x*eV_K/(ee*Bbar*Lbar*Lbar/tbar))*Te/temp_Ni*b0xcv*Grad(Pe)/B0;  //magnetic drift
     if(energy_flux)  ddt(Te) += 10.0/3.0*(KB*Te_x*eV_K/(ee*Bbar*Lbar*Lbar/tbar))*Te*b0xcv*Grad(Te)/B0;    //energy term
       
     if(external_sources) ddt(Te) += 0.6667*S_Ee_ext/temp_Ni;
     if(terms_Gradperp_diffcoefs) 
      {
        //grad_perp_chie=Grad_perp(temp_Ni*chi_e_perp);
        grad_perp_chie=Grad_perp(chi_e_perp);
        grad_perp_chie.applyBoundary();
        mesh->communicate(grad_perp_chie);  
        ddt(Te) +=  0.6667*V_dot_Grad(grad_perp_chie,Te); 
       }
     if(terms_Gradpara_diffcoefs) 
      {
        grad_para_kappTe = Grad_par(kappa_Te)*Grad_par(Te);
        grad_para_kappTe.applyBoundary("neumann");
        mesh->communicate(grad_para_kappTe);  
        ddt(Te) += 0.6667*grad_para_kappTe/temp_Ni;
       }
    }

//************************************
  // ION TEMPERATURE   ---Ti---
//************************************
  
  // output.write("Now updata Ti \n");
  if(!Solving_Eq_Ti)
    {ddt(Ti)=0.0;} 
  else
    {
     ddt(Ti) = - Vpar_Grad_par(Vi,Ti)
               - 0.6667*Ti*Grad_par(Vi)
               + (1.0-E_iol1)*0.6667*kappa_Ti*Grad2_par2(Ti)/temp_Ni
               + (1.0-E_iol1)*0.6667*Grad_par(kappa_Ti)*Grad_par(Ti)/temp_Ni
               + (1.0-E_iol1)*0.6667*chi_i_perp*Delp2(Ti)
               - nu_ionz*Ti
               - 0.6667*nu_CX*(Ti-Tn)
               + 2.*fmei*(Te-Ti)/tau_ei
               ;
     if(terms_radial_conv) ddt(Ti) +=0.6667*Ti*Div(V_diff)+Diff_ni_perp*Grad_perp(Ni)*Grad_perp(Ti)/Ni;
     if(terms_recombination)  ddt(Ti) += nu_rec*Ti;
     if(terms_exb) ddt(Ti) -= (1.0-E_iol1)*bracket((phi01/sbc_lambda1), Ti, bm_exb);           //ExB drift
     if(curvature_phi) ddt(Ti) -= (1.0-E_iol1)*4.0/3.0*Ti*b0xcv*Grad((phi01/sbc_lambda1))/B0;  //ExB drift
     if(diamag)  ddt(Ti) -= (1.0-E_iol1)*4.0/3.0*(KB*Ti_x*eV_K/(Zi*ee*Bbar*Lbar*Lbar/tbar))*Ti/temp_Ni*b0xcv*Grad(Pi)/B0;  //magnetic drift
     if(energy_flux)  ddt(Ti) -= (1.0-E_iol1)*10.0/3.0*(KB*Ti_x*eV_K/(Zi*ee*Bbar*Lbar*Lbar/tbar))*Ti*b0xcv*Grad(Ti)/B0;  //energy term
       
     if(external_sources) ddt(Ti) += 0.6667*S_Ei_ext/temp_Ni;
     if(terms_Gradperp_diffcoefs) 
      {
       grad_perp_chii=Grad_perp(chi_i_perp);
       //grad_perp_chii=Grad_perp(chi_i_perp*temp_Ni);
       grad_perp_chii.applyBoundary();
       mesh->communicate(grad_perp_chii);  
       //ddt(Ti) += (1.0-E_iol1)*0.6667*V_dot_Grad(grad_perp_chii,Ti)/temp_Ni;
       ddt(Ti) += (1.0-E_iol1)*0.6667*V_dot_Grad(grad_perp_chii,Ti);
       }
     /*if(terms_Gradpara_diffcoefs) 
       {
        grad_para_kappTi = Grad_par(kappa_Ti)*Grad_par(Ti);
        //grad_para_kappTi.applyBoundary();
        grad_para_kappTi.applyBoundary("neumann");
        mesh->communicate(grad_para_kappTi);  
        ddt(Ti) += 0.6667*grad_para_kappTi/temp_Ni;
       }*/
    }

  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // Turbulent Diffusion Update for Ni Te Ti
  //___________________________________________  

  pei = Te*Ni+Ti*Ni;
  Pe = Te*Ni;
  Pi = Ti*Ni;
  Pe0 = Te0*N0;
  Pi0 = Ti0*N0;
//************************************
  // ION Paralell VELOCITY  ---Vi---
//************************************

  // output.write("Now updata Vi \n");
  Grad_par_pei=Grad_par(pei);
  Grad_par_pei.applyBoundary();
  //mesh->communicate(Grad_par_pei);
  if(!Solving_Eq_Vi)
    {ddt(Vi)=0.0;} 
  else
    {   
     ddt(Vi) = - Vpar_Grad_par(Vi,Vi)
               - Grad_par(pei)/temp_Ni/Mi
               + 2.*0.6667*Grad_par(eta0_i)*Grad_par(Vi)/temp_Ni/Mi
               + 2.*0.6667*eta0_i*Grad2_par2(Vi)/temp_Ni/Mi
               - (nu_ionz+nu_CX)*(Vi-Vn)
                ;             
     if(terms_radial_conv) ddt(Vi) +=Diff_ni_perp*Grad_perp(Ni)*Grad(Vi)/Ni;
     if(terms_exb) ddt(Vi) -= bracket((phi01/sbc_lambda1), Vi, bm_exb);       //EXB drift
     if(terms_IOL) ddt(Vi) +=V_iol1;    
     }

//*************************************************************li2016
//************************************
  // Potential  ---phi01---
//************************************
  if(!Solving_Eq_phi01)
    {ddt(phi01)=0.0;} 
  else
    {   
     ddt(phi01) =(sbc_lambda0*Ni)/(sbc_lambda1*B0)*Laplace_perp(phi01)
                  + (sbc_lambda0*KB*eV_K*Ti_x)/(Lbar/tbar*Lbar*Bbar*Zi*ee*B0)*Laplace_perp(Pi)
                  - sbc_lambda0*U00;
     if (terms_cross) ddt(phi01) += sbc_lambda0/(B0*sbc_lambda1)*Grad(Ni)*Grad(phi01);  
     }
      phi00=phi01/sbc_lambda1;

//************************************
  // Vorticity  ---U00---
//************************************
  if(!Solving_Eq_U00)
    {ddt(U00)=0.0;} 
  else
    {   
     ddt(U00) =  (B0*B0*Bbar*Bbar/(MU0*Mi2*Ni_x*density_unit*Lbar*Lbar/tbar/tbar))*Grad_par(J1_para)
                + (Mu_perp/(Lbar/tbar*Lbar))*Laplace_perp(U00)
                + (Mu_para/(Lbar/tbar*Lbar))*Laplace_par(U00);

     if(include_J0) ddt(U00) +=  (B0*B0*Bbar*Bbar/(MU0*Mi2*Ni_x*density_unit*Lbar*Lbar/tbar/tbar))*Grad_par(J0);
     if(include_curvature)  ddt(U00) += 2.0*KB*Te_x*eV_K/(Mi2*Lbar*Lbar/(tbar*tbar))*Div(pei*b0xcv);  // curvature term
     if(terms_exb) ddt(U00) -=bracket((phi01/sbc_lambda1), U00, bm_exb);      //ExB drift
     if(terms_radial_conv) ddt(U00) +=Diff_ni_perp*Grad_perp(Ni)*Grad(U00)/Ni;
     if(gyroviscous)
      {
       Dperp2phi0 = Field3D(Delp2(phi01/sbc_lambda1));
       Dperp2phi0.applyBoundary("neumann");
       mesh->communicate(Dperp2phi0);

       Dperp2Pi = Field3D(Delp2(Pi));
       Dperp2Pi.applyBoundary("neumann");
       mesh->communicate(Dperp2Pi);

       Gradphi02 = Field3D(Grad_perp(phi01/sbc_lambda1) * Grad_perp(phi01/sbc_lambda1) / (B0 * B0));
       Gradphi02.applyBoundary("neumann");
       mesh->communicate(Gradphi02);

       bracketphi0P = bracket((phi01/sbc_lambda1), Pi, bm_exb);
       bracketphi0P.applyBoundary("neumann");
       mesh->communicate(bracketphi0P);

       Ugyro = KB * Ti_x * eV_K *tbar / (Zi * ee * Bbar * Lbar * Lbar);
       ddt(U00) -= 0.5 * Ugyro * bracket(Pi, Dperp2phi0, bm_exb) / B0;
       ddt(U00) += 0.5 * B0 * bracket(Ni, Gradphi02, bm_exb);
       ddt(U00) += 0.5 * Ugyro * bracket((phi01/sbc_lambda1), Dperp2Pi, bm_exb) / B0;
       ddt(U00) -= 0.5 * Ugyro * Delp2(bracketphi0P) / B0;

       term_gyro = - 0.5 * Ugyro * bracket(Pi, Dperp2phi0, bm_exb) /B0 
                   + 0.5 * B0 * bracket(Ni, Gradphi02, bm_exb) 
                   + 0.5 * Ugyro * bracket((phi01/sbc_lambda1), Dperp2Pi, bm_exb) / B0
                   - 0.5 * Ugyro * Delp2(bracketphi0P) / B0;
       }
    }

//************************************
  // Solving_plasma_Eq  ---End---
//************************************

//////////////////////////////////////////////
///-----calculat some output variables-----///
//////////////////////////////////////////////
  term_J1 =  (B0*B0*Bbar*Bbar/(MU0*Mi2*Ni_x*density_unit*Lbar*Lbar/tbar/tbar))*Grad_par(J1_para);
  term_J0 = (B0*B0*Bbar*Bbar/(MU0*Mi2*Ni_x*density_unit*Lbar*Lbar/tbar/tbar))*Grad_par(J0); 
  term_curvature = 2.0*KB*Te_x*eV_K/(Mi2*Lbar*Lbar/(tbar*tbar))*Div(pei*b0xcv);  // curvature term
  term_exb=-1.0* bracket((phi01/sbc_lambda1), U00, bm_exb);
  term_Mu_perp = (Mu_perp/(Lbar/tbar*Lbar))*Laplace_perp(U00);
  term_Mu_par = (Mu_para/(Lbar/tbar*Lbar))*Laplace_par(U00);

  J1_phi =  -1.0/(eta*sbc_lambda1*B0)*Grad_par(phi01);
  J1_Te = ((0.71*KB*eV_K*Te_x*tbar)/(eta*Lbar*Lbar*B0*Bbar*ee))*Grad_par(Te);
  J1_Pe = ((KB*eV_K*Te_x*tbar)/(eta*Lbar*Lbar*ee*Ni*B0*Bbar))*Grad_par(Pe);

  er0 = KB*Ti_x*eV_K*tbar/(Ni*Bbar*Lbar*Lbar*Zi*ee)*Grad(Pi);
  Vexb = (B0vec^Grad(phi00))/B0/B0;
  Vdia = KB*Ti_x*eV_K*tbar/(Zi*ee*Bbar*Lbar*Lbar)*(B0vec^Grad(Pi))/B0/B0/Ni;
  Vdia_e = -KB*Ti_x*eV_K*tbar/(ee*Bbar*Lbar*Lbar)*(B0vec^Grad(Pe))/B0/B0/Ni;
  mesh->communicate(Vexb);
  Vexb.applyBoundary();
  Vexb.toCovariant();
  mesh->communicate(Vdia);
  Vdia.applyBoundary();
  Vdia.toCovariant();
  mesh->communicate(Vdia_e);
  Vdia_e.applyBoundary();
  Vdia_e.toCovariant();
  mesh->communicate(er0.x);
  // mesh->communicate(er0.y);
  (er0.x).applyBoundary();
  //(er0.y).applyBoundary();
  er0.toCovariant();
  er0.x *= Rxy*Bpxy;

  er00 = -1.0*Grad(phi00);
  mesh->communicate(er00.x);
  // mesh->communicate(er00.y);
  (er00.x).applyBoundary();
  // (er00.y).applyBoundary();
  er00.toCovariant();
  er00.x *= Rxy*Bpxy;

  //Energy flux term
  Vector3D eng_i = B0vec^Grad(Ti);
  Vector3D eng_e = B0vec^Grad(Te);
  eng_i.applyBoundary();
  eng_i.toCovariant();
  mesh->communicate(eng_i);
  eng_e.applyBoundary();
  eng_e.toCovariant();
  mesh->communicate(eng_e);

//******** parallel heat flux *******************//
  heatf_cond_i = -kappa_Ti*Grad_par(Ti)*Lbar/tbar*Te_x*ee*Ni_x*density_unit;  //unit W/(m^-2)
  mesh->communicate(heatf_cond_i);
  (heatf_cond_i).applyBoundary("neumann");
  heatf_cond_e = -kappa_Te*Grad_par(Te)*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_cond_e);
  (heatf_cond_e).applyBoundary("neumann");
  heatf_conv_i = 2.5*Ni*Ti*Vi*Lbar/tbar*Ti_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_conv_i);
  (heatf_conv_i).applyBoundary("neumann");
  heatf_conv_e = 2.5*Ni*Te*c_se*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_conv_e);
  (heatf_conv_e).applyBoundary("neumann");
  heatf_conv_e1 = 2.5*Ni*Te*Ve*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_conv_e1);
  (heatf_conv_e1).applyBoundary("neumann");
//******** perpendicular heat flux *******************//
  heatf_cond_i_x = -chi_i_perp*sqrt(mesh->g11)*(Grad(Ti)).x*Lbar/tbar*Ti_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_cond_i_x);
  (heatf_cond_i_x).applyBoundary("neumann");
  heatf_cond_e_x = -chi_e_perp*sqrt(mesh->g11)*(Grad(Te)).x*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_cond_e_x);
  (heatf_cond_e_x).applyBoundary("neumann");
  heatf_exb_i_x = 2.5*Ni*Ti*Vexb.x*Rxy*Bpxy*Lbar/tbar*Ti_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_exb_i_x);
  (heatf_exb_i_x).applyBoundary("neumann");
  heatf_exb_e_x = 2.5*Ni*Te*Vexb.x*Rxy*Bpxy*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_exb_e_x);
  (heatf_exb_e_x).applyBoundary("neumann");
  heatf_dia_i_x = 2.5*Ni*Ti*Vdia.x*Rxy*Bpxy*Lbar/tbar*Ti_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_dia_i_x);
  (heatf_dia_i_x).applyBoundary("neumann");
  heatf_dia_e_x = 2.5*Ni*Te*(Vdia_e).x*Rxy*Bpxy*Lbar/tbar*Te_x*ee*Ni_x*density_unit;
  mesh->communicate(heatf_dia_e_x);
  (heatf_dia_e_x).applyBoundary("neumann");
  heatf_eng_i_x = 2.5*(Ni*Ni_x*density_unit*Ti*Ti_x/(Zi*B0*B0*Bbar))*eng_i.x*sqrt(mesh->g11)*Ti_x*ee/Lbar;
  (heatf_eng_i_x).applyBoundary("neumann");
  mesh->communicate(heatf_eng_i_x);
  heatf_eng_e_x = 2.5*(Ni*Ni_x*density_unit*Te*Te_x/(B0*Bbar*B0))*eng_e.x*sqrt(mesh->g11)*Te_x*ee/Lbar;
  (heatf_eng_e_x).applyBoundary("neumann");
  mesh->communicate(heatf_eng_e_x);
//******* radial particle flux ****************//
  gamma_cond_i_x = -Diff_ni_perp*sqrt(mesh->g11)*(Grad(Ni)).x*Lbar/tbar*Ni_x*density_unit;
  mesh->communicate(gamma_cond_i_x);
  (gamma_cond_i_x).applyBoundary("neumann");
  gamma_exb_i_x = Ni*Vexb.x*Rxy*Bpxy*Lbar/tbar*Ni_x*density_unit;
  mesh->communicate(gamma_exb_i_x);
  (gamma_exb_i_x).applyBoundary("neumann");
  gamma_dia_i_x = Ni*Vdia.x*Rxy*Bpxy*Lbar/tbar*Ni_x*density_unit;
  mesh->communicate(gamma_dia_i_x);
  (gamma_dia_i_x).applyBoundary("neumann");

//*************************************************************liend   

//************************************
  //Neutral Paralell Velocity   ---Vn---
//************************************
  pn=Nn*Tn;                            
  // Grad_par_pn=Grad_par(pn);
  // Grad_par_pn.applyBoundary();
  // mesh->communicate(Grad_par_pn);
  // Grad_par_Tn=Grad_par(Tn);
  // Grad_par_Tn.applyBoundary();
  // mesh->communicate(Grad_par_Tn);
  temp_Nn=field_larger(Nn,minimum_val);
  Grad_par_logNn=Grad_par(log(temp_Nn));
  Grad_par_logNn.applyBoundary();
  mesh->communicate(Grad_par_logNn);
  if(nlfilter_Gradpar_logNn)nl_filter(Grad_par_logNn,filter_para);
   /*
   ddt(Vn) = - Vpar_Grad_par(Vn,Vn)       // upwind convection
             + 2.*0.6667*eta0_n*Grad2_par2(Vn)/temp_Nn/Mn
             + nu_CX_n*(Vi-Vn)
             - 2.*S_diss*Vn/temp_Nn
             ;    
   if(terms_recombination) ddt(Vn) += S_rec*(Vi-Vn)/temp_Nn;
   if(terms_Gradpar_pn)    ddt(Vn) -= (Grad_par(Tn)/Mn;                 // Part I of Grad_par(Pn) term
                                     + Grad_par(log(temp_Nn))*Tn/Mn);   // Part II of Grad_par(Pn) term
   if(terms_Gradpar_eta0n) ddt(Vn) += 2.*0.6667*Grad_par(log(temp_Nn))*Grad_par(Vn)*eta0_n/temp_Nn/Mn  // Part I of Grad_par(eta0_n) term
    	                             +2.*0.6667*Grad_par(Tn)*Grad_par(Vn)*eta0_n/Tn/temp_Nn/Mn;       // Part II of Grad_par(eta0_n) term
   */
  ddt(Vn)=0.0;          // Vn-->Diffc_nn_par
  //Vn = (Vi*nu_CX_n -Grad_par_Tn/Mn-Grad_par_logNn*Tn/Mn)/(nu_CX_n+2.*S_diss/temp_Nn);
  //Vn = Vi-(Grad_par_Tn/Mn+Grad_par_logNn*Tn/Mn)/nu_CX_n;

//************************************
// Neutral Density       ---Nn---
//************************************

  //output.write("Before Nn \n");
 
  Diffc_nn_perp = Tn/Mn/nu_CX_n;               // Diffc_nn_perp re-calculated   
  V_th_n=sqrt(Tn/Mn);         
  Diffc_nn_perp_fl=V_th_n*Lnn_min/Lbar;
  Diffc_nn_perp *= Diffc_nn_perp_fl/(Diffc_nn_perp+Diffc_nn_perp_fl); 

  if(!Solving_Eq_Nn)
    {ddt(Nn)=0.0;}
  else
    { 
     ddt(Nn)=//- Vpar_Grad_par(Vn,Nn)  
            //- Nn*Grad_par(Vn)  re-written in more terms below
              Diffc_nn_perp*Delp2(Nn)
             - Si_p
             + 2.*S_diss
             ;
     if(terms_recombination) ddt(Nn) += S_rec;
     //if(terms_NnGradpar_Vn) ddt(Nn) += (Nn*Grad2_par2(Tn)/Mn/nu_CX_n
     //                                +Nn*Grad_par(Grad_par_logNn)*Tn/Mn/nu_CX_n
     //				        ); 
     if(terms_NnGradpar_Vi) ddt(Nn) -= Nn*Grad_par(Vi)+Vi*Nn*Grad_par_logNn;    
     if(terms_Diffcnn_par) 
       { 
        Diffc_nn_par.applyBoundary("neumann");
        mesh->communicate(Diffc_nn_par);
        ddt(Nn) += Diffc_nn_par*Nn*Grad_par(Grad_par_logNn)
                  + Nn*Grad_par_logNn*Grad_par(Diffc_nn_par);
        }
     }
  /*   // tested with some zig-zag problems at SOL region
  if(terms_Gradperp_diffcoefs) 
    {
      grad_perp_Diffn=Grad_perp(Diffc_nn_perp);
      grad_perp_Diffn.applyBoundary();
      mesh->communicate(grad_perp_Diffn);  
      ddt(Nn) +=  V_dot_Grad(grad_perp_Diffn,Nn); 
    }
  */ 
  //Nn.applyBoundary();

//************************************
  // Neutral Temperature       ---Tn---
//************************************

  ddt(Tn)=0.0;
   /*
  ddt(Tn)= - Vpar_Grad_par(Vn,Tn)
          - 0.6667*Tn*Grad_par(Vn)
          + 0.6667*chic_n_perp*Delp2(Tn)
       // + (Si_p-2.*S_diss)*Tn/temp_Nn
          + 0.6667*nu_CX_n*(Ti-Tn)
       // + 0.6667*S_diss*W_diss/temp_Nn
          ;
  */
//************************************
  // Molecule Perpendicular Velocity in X ---Vmx---
//************************************
  
  pm=Nm*Tm_x; 
  // ddt(Vmx) =0.;
  Vm.x=Vmx;
  Vm.y=0.;
  Vm.z=0.;
  ddt(Vm) = - V_dot_Grad(Vm,Vm) 
            - Grad(pm)/temp_Nm/Mm
            ;
  /* 
  ddt(Vmx) = - VDDX(Vmx,Vmx)//(mesh->J*sqrt(mesh->g_22))
             - DDX(pm)/temp_Nm/Mm
             ;
  */
//************************************
  // Molecule Density       ---Nm---
//************************************
  //ddt(Nm)=0.0;
  ddt(Nm)=-V_dot_Grad(Vm,Nm)-Nm*Div(Vm) - S_diss;
  /*
  ddt(Nm)=- VDDX(Vmx,Nm)//(mesh->J*sqrt(mesh->g_22))
          - Nm*DDX(Vmx)
        //+ 1.e5*Diffc_nm_perp*Delp2(Nm)
       // - S_diss
          ; 
  */

  //Diagnose whether there are negative value

  //Diag_neg_value(Nn,Nm,Tn,Te);

  //Bootstrap current calculated by using Sauter's formula 
  if (BScurrent)
    {
      q95=q95_input;
      pei= Ni*(Te+Ti);
      Pe = Ni*Te;
      Pi = Ni*Ti;

      nu_estar = 100.*nu_ei * q95*tbar / (V_th_e) / (Aratio*sqrt(Aratio));
      nu_istar = 100.*nu_ii * q95*tbar / (V_th_i) / (Aratio*sqrt(Aratio));

      ft = BS_ft(100);
      f31 = ft / (1.+(1.-0.1*ft)*sqrt(nu_estar) + 0.5*(1.-ft)*nu_estar/Zi);
      f32ee = ft / (1.+0.26*(1.-ft)*sqrt(nu_estar) + 0.18*(1.-0.37*ft)*nu_estar/sqrt(Zi));
      f32ei = ft / (1.+(1.+0.6*ft)*sqrt(nu_estar) + 0.85*(1.-0.37*ft)*nu_estar*(1.+Zi));
      f34 = ft / (1.+(1.-0.1*ft)*sqrt(nu_estar) + 0.5*(1.-0.5*ft)*nu_estar/Zi);

      L31 = F31(f31) ;
      L32 = F32ee(f32ee)+F32ei(f32ei) ;
      L34 = F31(f34) ;

      BSal0 = - (1.17*(1.-ft))/(1.-0.22*ft-0.19*ft*ft);
      BSal = (BSal0+0.25*(1-ft*ft)*sqrt(nu_istar))/(1.+0.5*sqrt(nu_istar)) + 0.31*nu_istar*nu_istar*ft*ft*ft*ft*ft*ft;
      BSal *= 1./(1.+0.15*nu_istar*nu_istar*ft*ft*ft*ft*ft*ft);

      Jpar_BS0 = L31* DDX(pei)/Pe  + L32*DDX(Te)/Te + L34*DDX(Ti)/(Zi*Te)*BSal;
      Jpar_BS0 *= Field3D( -Rxy*Btxy*Pe*(MU0*KB*Ni_x*density_unit*Te_x*eV_K)/(mesh->Bxy*mesh->Bxy)/(bmag*bmag) );
      //NB:   J_hat = MU0*Lbar * J / mesh->Bxy;

      mesh->communicate(Jpar_BS0);
      Jpar_BS0.applyBoundary();

    }

  return(0);
}


void update_turb_effects(Field3D &f, const BoutReal Difft_coef, const BoutReal value_crit)
{
  bool first_time=true;
       Field3D dNi_dx=DDX(Ni);
       BoutReal lowest_Lp=1.;
       int jx_lowest_Lp;
       Heavistep=0.;
       dpei_dx=DDX(pei);
       for (int jx=0;jx<mesh->ngx;jx++)
         {
           for (int jy=0;jy<mesh->ngy;jy++)
	     {
 	       for (int jz=0;jz<mesh->ngz;jz++)
                  {                         
                   if(abs(dpei_dx[jx][jy][jz])<1.e-5) dpei_dx[jx][jy][jz]=1.e-5;
                   if(jx<=1 )pei[jx][jy][jz]=Ni[jx+2][jy][jz]*(Te[jx+2][jy][jz]+Ti[jx+2][jy][jz]);
                   if(jx>=mesh->ngx-2 )pei[jx][jy][jz]=Ni[jx-2][jy][jz]*(Te[jx-2][jy][jz]+Ti[jx-2][jy][jz]);
        	   Lp=pei[jx][jy][jz]/(dpei_dx[jx][jy][jz]);
	           if(abs(Lp)<lowest_Lp) {lowest_Lp =abs(Lp);jx_lowest_Lp=jx;}
                   //if (abs(Lp) <= Lp_crit)
                   if (abs(dNi_dx[jx][jy][jz])>=9.e2)
                     {
                      Heavistep[jx][jy][jz]=1.;
                     }
        	   else Heavistep[jx][jy][jz]=0.;
	           }
               }
       	   }
       ddt(f) += Heavistep*Difft_coef*Delp2(f);
}


//3D Boundaries of constant neutral flux injection
  
const Field3D ret_const_flux_BC(const Field3D &var, const BoutReal value)

{

  Field3D result;

  result.allocate();

    for (int jx=0;jx<mesh->ngx;jx++)
        {

          BoutReal x_glb= mesh->GlobalX(jx); 

           for (int jy=0;jy<mesh->ngy;jy++)
	      {
                 BoutReal y_glb=mesh->GlobalY(jy);

 	        for (int jz=0;jz<mesh->ngz;jz++)
                   {
                     BoutReal z_glb=(BoutReal)jz/(BoutReal)(mesh->ngz-1);    // no parallization in Z   
                      
		     if(x_glb>=CF_BC_x0 && y_glb>=CF_BC_y0 && y_glb<=CF_BC_y1 && z_glb>=CF_BC_z0 && z_glb<=CF_BC_z1) 
		       { 
			 //BoutReal CF_BC_yhalf=0.5*(CF_BC_y0+CF_BC_y1);
                         //BoutReal CF_ywidth=CF_BC_yhalf-CF_BC_y0;
                         //BoutReal CF_exp_decay=exp(-(y_glb-CF_BC_yhalf)*(y_glb-CF_BC_yhalf)/CF_ywidth/CF_ywidth);
			 result[jx][jy][jz]=value;//*CF_exp_decay; 
                                                  
                         }
                     else 
                         result[jx][jy][jz]=var[jx][jy][jz];

	              }
        	}
     	}

  mesh->communicate(result);

  return(result);


}

     

void Diag_neg_value (const Field3D &f1,const Field3D &f2,const Field3D &f3,const Field3D &f4 )
{
for (int jx=0;jx<mesh->ngx;jx++)
   {
     BoutReal x_glb=mesh->GlobalX(jx);
      for( int jy=0;jy<mesh->ngy;jy++)
        {
         BoutReal y_glb=mesh->GlobalY(jy);
         for (int jz=0;jz<mesh->ngz;jz++)
           {
             if(f1[jx][jy][jz]<0.) output.write("Field 1 becomes negative %e at x %e y %e \n",f1[jx][jy][jz],x_glb,y_glb);
             if(f2[jx][jy][jz]<0.) output.write("Field 2 becomes negative %e at x %e y %e \n", f2[jx][jy][jz],x_glb,y_glb);
             if(f3[jx][jy][jz]<0.) output.write("Field 3 becomes negative %e at x %e y %e \n",f3[jx][jy][jz],x_glb,y_glb);
             if(f4[jx][jy][jz]<0.) output.write("Field 4 becomes negative %e at x %e y %e \n", f4[jx][jy][jz],x_glb,y_glb);
            }
         }    
    }


}


const Field3D field_larger(const Field3D &f, const BoutReal limit)

{

  Field3D result;

  result.allocate();

//  #pragma omp parallel for

  for(int jx=0;jx<mesh->ngx;jx++)

    for(int jy=0;jy<mesh->ngy;jy++)

      for(int jz=0;jz<mesh->ngz;jz++)

      {

        if(f[jx][jy][jz] >= limit)

             result[jx][jy][jz] = f[jx][jy][jz];

         else

             result[jx][jy][jz] = limit;

      }

  mesh->communicate(result);

  return(result);

}


const Field3D field_smaller(const Field3D &f, const BoutReal limit)

{

  Field3D result;

  result.allocate();

//  #pragma omp parallel for

  for(int jx=0;jx<mesh->ngx;jx++)

    for(int jy=0;jy<mesh->ngy;jy++)

      for(int jz=0;jz<mesh->ngz;jz++)

      {

        if(f[jx][jy][jz] <= limit)

             result[jx][jy][jz] = f[jx][jy][jz];

         else

             result[jx][jy][jz] = limit;

      }

  mesh->communicate(result);

  return(result);

}



/****************BOUNDARY FUNCTIONS*****************************/

// Sheath Boundary Conditions 

// Linearized

void SBC_Dirichlet_SWidth1(Field3D &var, const Field3D &value) //let the boundary equall to the value next to the boundary

{
  RangeIterator xrup = mesh->iterateBndryUpperY();

  for(xrup.first(); !xrup.isDone(); xrup.next())

    {
       {
          for(int jy=mesh->yend+1-1; jy<mesh->ngy; jy++)
            {
              for(int jz=0; jz<mesh->ngz; jz++) 
                {

                  var[xrup.ind][jy][jz] = value[xrup.ind][mesh->yend][jz];

                }
	    }
       }
    }
   RangeIterator xrdn = mesh->iterateBndryLowerY();

  for(xrdn.first(); !xrdn.isDone(); xrdn.next())
    {
       for(int jy=mesh->ystart-1+1; jy>=0; jy--)

            for(int jz=0; jz<mesh->ngz; jz++) 

                {

                  var[xrdn.ind][jy][jz] = value[xrdn.ind][mesh->ystart][jz];

                 }
       
    }


}
void SBC_Dirichlet(Field3D &var, const Field3D &value) //let the boundary equall to the value next to the boundary

{

  SBC_yup_eq(var, value);    

  SBC_ydown_eq(var, -value);   // Fluxes go towards X point and hit on plates, thus SBC at ydown or theta=0 should be negative 

}

 
void SBC_Dirichlet_phi(Field3D &var, const Field3D &value) //let the boundary equall to the value next to the boundary

{

  SBC_yup_eq(var, value);    

  SBC_ydown_eq(var, value);   // Fluxes go towards X point and hit on plates, thus SBC at ydown or theta=0 should be negative 

}

void SBC_Gradpar(Field3D &var, const Field3D &value)

{

  SBC_yup_Grad_par(var, value);  

  SBC_ydown_Grad_par(var, -value); // Fluxes go towards X point and hit on plates, thus SBC at ydown or theta=0 should be different

}

 

// Boundary to specified Field3D object

void SBC_yup_eq(Field3D &var, const Field3D &value)

{
  RangeIterator xrup = mesh->iterateBndryUpperY();

  for(xrup.first(); !xrup.isDone(); xrup.next())

    {
     // BoutReal x_glb= mesh->GlobalX(xrup->ind); 

      //if(x_glb>=Sheath_BC_x0) 

       {
          for(int jy=mesh->yend+1-Sheath_width; jy<mesh->ngy; jy++)
            {
              for(int jz=0; jz<mesh->ngz; jz++) 
                {

                  var[xrup.ind][jy][jz] = value[xrup.ind][mesh->yend][jz];

                }
	    }
       }
    }

}

 

void SBC_ydown_eq(Field3D &var, const Field3D &value)

{

   RangeIterator xrdn = mesh->iterateBndryLowerY();

  for(xrdn.first(); !xrdn.isDone(); xrdn.next())
    {
      //BoutReal x_glb= mesh->GlobalX(xrdn->ind); 

      //if(x_glb>=Sheath_BC_x0) 

       {
         for(int jy=mesh->ystart-1+Sheath_width; jy>=0; jy--)

            for(int jz=0; jz<mesh->ngz; jz++) 

                {

                  var[xrdn.ind][jy][jz] = value[xrdn.ind][mesh->ystart][jz];

                 }
       }
    }
}      

 

// Boundary gradient to specified Field3D object

void SBC_yup_Grad_par(Field3D &var, const Field3D &value)

{
   RangeIterator xrup = mesh->iterateBndryUpperY();

  for(xrup.first(); !xrup.isDone(); xrup.next())

    {
      //BoutReal x_glb= mesh->GlobalX(xrup->ind); 

      //if(x_glb>=Sheath_BC_x0) 

        {
          for(int jy=mesh->yend+1-Sheath_width; jy<mesh->ngy; jy++)

             for(int jz=0; jz<mesh->ngz; jz++) 

                {

                  var[xrup.ind][jy][jz] = var[xrup.ind][jy-1][jz] + mesh->dy[xrup.ind][jy]*sqrt(mesh->g_22[xrup.ind][jy])*value[xrup.ind][jy][jz];
 
                }
        }
    }

}

 

void SBC_ydown_Grad_par(Field3D &var, const Field3D &value)

{
  RangeIterator xrdn = mesh->iterateBndryLowerY();
  for(xrdn.first(); !xrdn.isDone(); xrdn.next())

    {
      //BoutReal x_glb= mesh->GlobalX(xrdn->ind); 

      //if(x_glb>=Sheath_BC_x0) 

        {
           for(int jy=mesh->ystart-1+Sheath_width; jy>=0; jy--)

              for(int jz=0; jz<mesh->ngz; jz++) 

                {

		   var[xrdn.ind][jy][jz] = var[xrdn.ind][jy+1][jz] - mesh->dy[xrdn.ind][jy]*sqrt(mesh->g_22[xrdn.ind][jy])*value[xrdn.ind][jy][jz];

                }
	}
    }

}      

void WallBC_Xout_GradX(Field3D &var, const Field3D &value)
{
  // NB: input value of Gradient X length in real R space 
  for(int jx=0;jx<mesh->ngx;jx++) 
    {
     if ( mesh->XGLOBAL (jx) > NX - 3 ) 
       {
         for(int jy=0;jy<mesh->ngy;jy++) 
             for(int jz=0;jz<mesh->ngz;jz++) 
                { 
                  var[jx][jy][jz] = var[jx-1][jy][jz] + value[jx][jy][jz]*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy]);  // calculated in BOUT psi coordinate
                 }
         }
     }
}

void WallBC_Xout_GradX_len(Field3D &var, BoutReal value)
{
  // NB: input value of Gradient X length in real R space 
  BoutReal temp;
   for(int jx=0;jx<mesh->ngx;jx++) 
    {
     if ( mesh->XGLOBAL (jx) > NX - 3 ) 
       {
         for(int jy=0;jy<mesh->ngy;jy++) 
             for(int jz=0;jz<mesh->ngz;jz++) 
                { 
                  temp=0.5*value*mesh->dx[jx][jy]/sqrt(mesh->g11[jx][jy]);     // transfer to BOUT psi coordinate
                  var[jx][jy][jz] = var[jx-1][jy][jz]*(1.+temp)/(1.-temp);  
                 }
         }
     }
}


const Field3D BS_ft(const int index)
{
  Field3D result, result1;
  result.allocate();
  result1.allocate();
  result1=0.;
  
  BoutReal xlam, dxlam;
  //dxlam = 1./bmag/index;     // wrong since normalization of bmag is needed
  //dxlam = 1./index;          // right when global max(Bxy)=bmag=Bbar
  dxlam = 1./max(mesh->Bxy)/index;       // It is still not perfect since max(mesh->Bxy) is the maximum value at some core in parallel computation
  //output.write("maxmum normalized Bxy %e \n",max(mesh->Bxy));

  xlam = 0.;

  for(int i=0; i<index; i++)
    {
      result1 += xlam*dxlam/sqrt(1.-xlam*mesh->Bxy);
      xlam += dxlam;
    }
  result = 1.- 0.75*mesh->Bxy*mesh->Bxy * result1;

  return(result);
}

const Field3D F31(const Field3D input)
{
  Field3D result;
  result.allocate();

  result = ( 1 + 1.4/(Zi+1.) ) * input;
  result -= 1.9/(Zi+1.) * input*input;
  result += 0.3/(Zi+1.) * input*input*input;
  result += 0.2/(Zi+1.) * input*input*input*input;

  return(result);
}

const Field3D F32ee(const Field3D input)
{
  Field3D result;
  result.allocate();
  
  result = (0.05+0.62*Zi)/(Zi*(1+0.44*Zi))*(input-input*input*input*input);
  result +=1./(1.+0.22*Zi)*( input*input-input*input*input*input-1.2*(input*input*input-input*input*input*input) );
  result += 1.2/(1.+0.5*Zi)*input*input*input*input;

  return(result);
}

const Field3D F32ei(const Field3D input)
{
  Field3D result;
  result.allocate();
  
  result = -(0.56+1.93*Zi)/(Zi*(1+0.44*Zi))*(input-input*input*input*input);
  result += 4.95/(1.+2.48*Zi)*( input*input-input*input*input*input-0.55*(input*input*input-input*input*input*input) );
  result -= 1.2/(1.+0.5*Zi)*input*input*input*input;

  return(result);
}

//***************************************************li2016
Field3D Laplace_par_filterPF (const Field3D &var)
{
  Field3D result;
  result.allocate();
  result = Laplace_par(var);

  for(int jx=0;jx<mesh->ngx;jx++)
    {
      int indx =  mesh->XGLOBAL(jx);
      BoutReal dindx = indx/ixsep;
      for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
        {
          int indy = mesh->YGLOBAL(jy);
          //output.write("dinx: %e   indy: %e  jysep1: %e  jysep2: %e\n", dindx, indy, jysep1, jysep2);
          if ((dindx < 1.0) && ((indy <= jysep1) || (indy > jysep2)))
            result[jx][jy][jz] = 0.;
        }
    }

  mesh->communicate(result);
  return result;
}

Field3D Laplace_perp_filterPF (const Field3D &var)
{
  Field3D result;
  result.allocate();
  result = Laplace_perp(var);

  for(int jx=0;jx<mesh->ngx;jx++)
    {
      int indx =  mesh->XGLOBAL(jx);
      BoutReal dindx = indx/ixsep;
      for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
        {
          int indy = mesh->YGLOBAL(jy);
          //output.write("dinx: %e   indy: %e  jysep1: %e  jysep2: %e\n", dindx, indy, jysep1, jysep2);
          if ((dindx < 1.0) && ((indy <= jysep1) || (indy > jysep2)))
            result[jx][jy][jz] = 0.;
        }
    }

  mesh->communicate(result);
  return result;
}

Field3D Laplace_perp_PFboundary (const Field3D &var)
{
  Field3D result;
  result.allocate();
  result = Delp2(var);


  for(int jz=0;jz<mesh->ngz;jz++)
  for(int jy=0;jy<mesh->ngy;jy++)
    {
      int indy = mesh->YGLOBAL(jy);
      if ((indy <= jysep1) || (indy > jysep2))
        {
          for(int jx=0;jx<mesh->ngx;jx++)
            {
              int indx1 =  mesh->XGLOBAL(jx);
              BoutReal dindx1 = indx1/ixsep;
              int indx0 =  mesh->XGLOBAL(jx-1);
              BoutReal dindx0 = indx0/ixsep;
              if ((dindx0 < 1.0) && (dindx1 > 1.0))
                {
                  result[jx+1][jy][jz] = 2.0*result[jx+2][jy][jz]-result[jx+3][jy][jz];
                  result[jx][jy][jz] = 2.0*result[jx+1][jy][jz]-result[jx+2][jy][jz];

                  result[jx-2][jy][jz] = 2.0*result[jx-3][jy][jz]-result[jx-4][jy][jz];
                  result[jx-1][jy][jz] = 2.0*result[jx-2][jy][jz]-result[jx-3][jy][jz];
                }
            }
        }
    }

  mesh->communicate(result);
  return result;
}
void Grad2_par2_sh(Field3D &var, const Field3D &value)
{

//  Field3D result;
//  result.allocate();
//  result = Grad2_par2(var)*mesh->dy*mesh->dy*mesh->g_22;
//  result = Grad2_par2(var);

  for(int jx=0; jx<mesh->ngx; jx++)
    {
      RangeIterator xrup = mesh->iterateBndryUpperY();
      BoutReal x_glb= mesh->GlobalX(xrup.ind);
      if (x_glb >= 1.0)
//      if (x[jx][0] >= 1.0)
        {
          for (int jy=mesh->ystart; jy<=mesh->yend; jy++)
            for (int jz=0; jz<mesh->ngz; jz++)
            {
              int mgy = mesh->YGLOBAL(jy);

              if (mgy == yloc-1)
                var[jx][jy][jz] = -var[jx][jy-2][jz]/12.0 + 4.0*var[jx][jy-1][jz]/3.0 - 2.5*var[jx][jy][jz] + 4.0*var[jx][jy+1][jz]/3.0 - value[jx][jx+2][jz]/12.0;

              if (mgy == yloc)
                var[jx][jy][jz] = -var[jx][jy-2][jz]/12.0 + 4.0*var[jx][jy-1][jz]/3.0 - 2.5*var[jx][jy][jz] + 4.0*value[jx][jy+1][jz]/3.0 - value[jx][jx+2][jz]/12.0;

              if (mgy == yloc+1)
                var[jx][jy][jz] = -value[jx][jy-2][jz]/12.0 + 4.0*value[jx][jy-1][jz]/3.0 - 2.5*var[jx][jy][jz] + 4.0*var[jx][jy+1][jz]/3.0 - var[jx][jx+2][jz]/12.0;

              if (mgy == yloc+2)
                var[jx][jy][jz] = -value[jx][jy-2][jz]/12.0 + 4.0*var[jx][jy-1][jz]/3.0 - 2.5*var[jx][jy][jz] + 4.0*var[jx][jy+1][jz]/3.0 - var[jx][jx+2][jz]/12.0;
            }
        }
    }

//  result /= mesh->dy*mesh->dy*mesh->g_22;
  /*
  mesh->communicate(result);
  result.setBoundary("U00");
  result.applyBoundary();
  */
 // return result;
}
//*******************************************************liend



