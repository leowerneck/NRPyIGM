int IllinoisGRMHD_Noble2D_con2prim_onepoint(struct eos_struct &eos, struct output_stats &stats, CCTK_REAL *CONSERVS,CCTK_REAL *ADM_3METRIC, CCTK_REAL *TUPMUNU,CCTK_REAL *TDNMUNU, CCTK_REAL *PRIMS) {

  int gamma_equals2 = 1;
  if (fabs(gamma_th-2.0) > 1.e-10) gamma_equals2 = 0;


  int check=0;
  
  /*
from outputC import *            # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

alpha   = sp.symbols('ADM_3METRIC[LAPSE]',real=True)
betaU   = ixp.zerorank2()
gammaDD = ixp.zerorank2()
for i in range(3):
    betaU[i] = sp.symbols('ADM_3METRIC[SHIFT'+chr(ord('X')+i)+"]",real=True)
    for j in range(i,3):
        gammaDD[i][j] = gammaDD[j][i] = sp.symbols('ADM_3METRIC[G'+chr(ord('X')+i)+chr(ord('X')+j)+"]",real=True)

gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
string = outputC([ gammaUU[0][0],       gammaUU[0][1],       gammaUU[0][2],
                   gammaUU[1][1],       gammaUU[1][2],       gammaUU[2][2], sp.sqrt(gammaDET)],
        ["ADM_3METRIC[GUPXX]","ADM_3METRIC[GUPXY]","ADM_3METRIC[GUPXZ]",
         "ADM_3METRIC[GUPYY]","ADM_3METRIC[GUPYZ]","ADM_3METRIC[GUPZZ]", "ADM_3METRIC[SQRTGAMMA]"],
        filename="returnstring", params="outCverbose=False")
# Replace pow(blah, 2) with (blah)*(blah)
string2 = re.sub('pow\(([^,]+), 2\)', '(\\1)*(\\1)', string); string = string2
with open("NRPy+ADMgammaUU.h", "w") as file:
    file.write(string)
  */
#include "NRPy+ADMgammaUU.h"
  // Compute 4-metric, both g_{\mu \nu} and g^{\mu \nu}.
  // This is for computing T_{\mu \nu} and T^{\mu \nu}. Also the HARM con2prim lowlevel function requires them.
  CCTK_REAL g4dn[4][4];
  /* CONTINUED FROM ABOVE
import re
import BSSN.ADMBSSN_tofrom_4metric as AB4m
AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
for i in range(3):
    for j in range(i,3):
        gammaUU[i][j] = gammaUU[j][i] = sp.symbols('ADM_3METRIC[GUP'+chr(ord('X')+i)+chr(ord('X')+j)+"]",real=True)
AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha,gammaUU)

exprlist = []
namelist = []
for mu in range(4):
    for nu in range(4):
        exprlist.append(AB4m.g4DD[mu][nu])
        namelist.append("g4dn["+str(mu)+"]["+str(nu)+"]")
string = outputC(exprlist,namelist,"returnstring", params="outCverbose=False")
# Replace pow(blah, 2) with (blah)*(blah)
string2 = re.sub('pow\(([^,]+), 2\)', '(\\1)*(\\1)', string); string = string2
with open("NRPy+g4DD.h", "w") as file:
    file.write(string)
  */
#include "NRPy+g4DD.h"

  CCTK_REAL g4up[4][4];
  /* CONTINUED FROM ABOVE g4dn[4][4] definition.
exprlist = []
namelist = []
for mu in range(4):
    for nu in range(4):
        exprlist.append(AB4m.g4UU[mu][nu])
        namelist.append("g4up["+str(mu)+"]["+str(nu)+"]")

string = outputC(exprlist,namelist,"returnstring", params="outCverbose=False")
# Replace pow(blah, 2) with (blah)*(blah)
string2 = re.sub('pow\(([^,]+), 2\)', '(\\1)*(\\1)', string); string = string2
string2 = re.sub('pow\(([^,]+), -2\)', '(1.0 / (\\1)*(\\1))', string); string = string2
with open("NRPy+g4UU.h", "w") as file:
    file.write(string)
  */
#include "NRPy+g4UU.h"

  /* Idea: Apply a rho_star floor. rho_* = alpha*sqrt(gamma)*rho_b*u^0, so minimum rho_star is simply alpha*sqrt(gamma)*rho_b_atm*1.0
     if(CONSERVS[RHOSTAR] < METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6]*rho_b_atm*1.0) {
     stats.failure_checker+=1;
     CONSERVS[RHOSTAR] = METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6]*rho_b_atm*1.0;
     stats.rho_star_fix_applied=1;
     }
  */
  if(CONSERVS[RHOSTAR]>0.0) {
    // Apply the tau floor
    apply_tau_floor(tau_atm,rho_b_atm,Psi6threshold,PRIMS,ADM_3METRIC,stats,eos,  CONSERVS);

    stats.font_fixed=0;
    for(int ii=0;ii<3;ii++) {
      check = harm_primitives_gammalaw_lowlevel(/*index,i,j,k,x,y,z,*/
                                                ADM_3METRIC,CONSERVS,PRIMS,  g4dn,g4up,   stats,eos);
      if(check==0) ii=4;
      else stats.failure_checker+=100000;
    }
  } else {
    stats.failure_checker+=1;
    // Set to atmosphere if rho_star<0.
    //FIXME: FOR GAMMA=2 ONLY:
    PRIMS[RHOB]      =rho_b_atm;
    
    /* Set P = P_cold */
    int polytropic_index = find_polytropic_K_and_Gamma_index(eos, rho_b_atm);
    CCTK_REAL K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
    CCTK_REAL Gamma_ppoly_tab = eos.K_ppoly_tab[polytropic_index];
    PRIMS[PRESSURE]      = K_ppoly_tab*pow(rho_b_atm,Gamma_ppoly_tab);
      
    PRIMS[VX]        =-ADM_3METRIC[SHIFTX];
    PRIMS[VY]        =-ADM_3METRIC[SHIFTY];
    PRIMS[VZ]        =-ADM_3METRIC[SHIFTZ];

    stats.rho_star_fix_applied=1;
  }
  // Enforce limits on primitive variables and recompute conservatives.
  const int already_computed_physical_metric_and_inverse=1;
  IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,PRIMS,stats,eos,ADM_3METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS);

  return check;
}
