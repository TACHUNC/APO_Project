* Course: ADVANCED PROCESS OPTIMIZATION (CENG 70003)1
* Author: TACHUN,CHEN SHUNHUA,JIAO

SETS
 i          functional group
 j          group contributional factor
 IC         integer Cut /1*5/
 dym(IC)    dynamic integer Cut
;

Parameters
 Property(i<,j<)  Group contributional table
 integercut(IC,i) Intger cut storage
* The rest are parameters used for result displacement
 number_of_groups(i,IC)
 objective_value(IC)
 melting_temperature(IC)
 boiling_temperature(IC)
 relative_energy_difference(IC)
 liquid_density(IC)
 heat_capacity(IC)
;

* Reading data from excel
$onEmbeddedCode Connect:
- ExcelReader:
    file: /Users/chen/Desktop/APO/contribution_table_v6.xlsx
    symbols:
        - name: Property
          range: Sheet1!B1:U24
          
- GAMSWriter:
    symbols: all 
$offEmbeddedCode


SCALARS
 Tavg   average operating temperature[K]
 T_abs  average temperature of the absorption column[K]                 /313/
 T_des  average temperature of the desorption column[K]                 /393/
 T_b0   universal constant for calculating Tb[~]                        /244.5165/
 T_m0   universal constant for calculating Tm[~]                        /143.5706/
 V_m0   universal constant for calculating Vm[~]                        /0.016/
 dD_co2 Hansen solution parameter of dispersion for co2[MPa^0.5]        /15.6/
 dP_co2 Hansen solution parameter of polarity for co2[MPa^0.5]          /5.2/
 dH_co2 Hansen solution parameter of hydrogen bonds for co2[MPa^0.5]    /5.8/
 R0_co2 interation radius for co2 [MPa^0.5]                             /4.0/
 a      wighting constant for RED in the objective function[~]          /0/
 b      weighting constant for Cp in the objective function[~]          /0/
* The reason that c is negative is that we want to maximize it, not minimize it.
 c      weighting constant for rho in the objective function[~]         /-1/
;

* Average operating temperature is the mean of the T_abs and T_dess
Tavg = (T_abs+T_des)/2; 


Variables
 obj    Objective Varaible
* Variables needed for the calculation of RED
 dD     Hansen solution parameter of dispersion
 dP     Hansen solution parameter of polarity
 dH     Hansen solution parameter of hydrogen bonds
* Variables needed for the calculation of Cp 
 alpha  nominator of the acentric factor
 beta_v denominator of the acentric factor
 omega  acentric factor
;

INTEGER VARIABLES
 n(i)   number of functional group i 
;

POSITIVE VARIABLES
 RED    relative energy difference[~]
 Tm     melting temperature[K]
 Tb     boiling temperature[K]
* Variables needed for the calculation of rho
 Vm     liquid molar volume[cc kmol^-1]
 Mw     molecular weight[g mol^-1]
 rho    liquid density[g cm^-3]
* Variables needed for the calculation of Cp
 Pc     critical pressure[bar]
 Tc     critical temperature[K]
 Tbr    reduced boiling temperature[~]
 Tavgr  reduced average operating temperature[~]
 Cp_a   ideal gas heat capacity[J mol^-1 K^-1]
 Cp_l   liquid heat capacity[J mol^-1 K^-1]
;

* Solvent's melting temperature must be less than T_abs
Tm.UP = T_abs;
* Solvent's boiling temperature must be higher than T_des
Tb.LO = T_des;

EQUATIONS
 MOBJ           Multi-Objective function                (Weighted sum method)
 OCT            Octet rule                              (O. Odele & S. Macchietto)
 BOND(i)        Bonding rule
 int_cut(IC)    integer cut
 Tm_eq          Melting temperature GC equation         (HeHukkerikar et al)
 Tb_eq          Boiling temperature GC equation         (HeHukkerikar et al)
* Equations needed for the calculation of rho
 Vm_eq          Liquid molar volume GC equation         (HeHukkerikar et al)
 Mw_eq          Molecular weight euqation
 rho_eq         Liquid density equation
* Equations needed for the calculation of RED
 dD_eq          Hansen solubility of dispersion         (HeHukkerikar et al)
 dP_eq          Hansen solubility of polarity           (HeHukkerikar et al)
 dH_eq          Hansen solubility of hydrogen bonds     (HeHukkerikar et al)
 RED_eq         Relative energy difference              (C. Hansen)
 RED_ub         Boundary conditions for affinity        (C.M. Hansen & K. Skaarup)
* Equations needed for the calculation of Cp
 Pc_eq          Critical pressure GC equation           (K.G. Joback & R.C. Reid)
 Tc_eq          Critical temperature GC equation        (K.G. Joback & R.C. Reid)
 Tbr_eq         Reduced boiling temperature             (N.V. Sahinidis et al)
 Tavgr_eq       Reduced average operating temperature   (N.V. Sahinidis et al)
 alpha_eq       Nominator of acentric factor            (N.V. Sahinidis et al)
 beta_eq        Denominator of acentric factor          (N.V. Sahinidis et al)
 omega_eq       Acentic factor                          (N.V. Sahinidis et al)
 Cpa_eq         Ideal gas heat capacity at T_avgr       (K.G. Joback & R.C. Reid)
 Cpl_eq         Liquid heat capacity at T_avgr          (J.S. Rowlinson)
* structural constraints for molecule      
 ester_ub_eq
 ether_ub_eq
 amine_ub_eq
 carbon_ub_eq
 carbon_lb_eq
;

* Final objective to be minimized
MOBJ..      obj     =E= a*RED + b*Cp_l + c*rho;
* Octet rule to make sure that there are no free electrons in the molecule
OCT..       0       =E= sum(i,n(i)*(2-Property(i,'Val'))) - 2;
* bonding rule to make sure all the bonds are feasible
ALIAS(i,ii);
BOND(i)..   0       =G= n(i)*(Property(i,'Val') - 1)+ 2 - sum(ii, n(ii));
* Equation needed for the calculation of Tm (melting temperature)
Tm_eq..     Tm      =E= T_m0 * log(sum(i,n(i)*Property(i,'Tm')));
* Equation needed for the calculation of Tb (boiling temperature)
Tb_eq..     Tb      =E= T_b0 * log(sum(i,n(i)*Property(i,'Tb')));
* Equations needed for the calculation of RED (relative energy difference)
dD_eq..     dD      =E= sum(i,n(i)*Property(i,'dD'));
dP_eq..     dP      =E= sum(i,n(i)*Property(i,'dP'));
dH_eq..     dH      =E= sum(i,n(i)*Property(i,'dH'));
RED_eq..    RED     =E= sqrt(4*power(dD-dD_co2,2) + power(dP-dP_co2,2) + power(dH - dH_co2,2)) / R0_co2;
RED_ub..    RED     =L= 1;
* Equations needed for the calculation of rho (liquid density)
Vm_eq..     Vm      =E= V_m0 + sum(i,n(i)*Property(i,'Vm'));
Mw_eq..     Mw      =E= sum(i,n(i)*Property(i,'Mr'));
rho_eq..    rho     =E= Mw * (10**(-3)) / Vm;
* Equations needed for the calculation of Cp (liquid heat capacity)
Pc_eq..     Pc      =E= 1/(0.113 + 0.0032*sum(i,n(i)*Property(i,'a')) - sum(i,n(i)*Property(i,'Pc')))**2;
Tc_eq..     Tc      =E= Tb/(0.584 + 0.965*sum(i,n(i)*Property(i,'Tc')) - (sum(i,n(i)*Property(i,'Tc')))**2);
Tbr_eq..    Tbr     =E= Tb/Tc;
Tavgr_eq..  Tavgr   =E= Tavg/Tc;
alpha_eq..  alpha   =E= -5.97214 - log(Pc/1.013) + 6.09648/Tbr + 1.28862*log(Tbr) - 0.169347*(Tbr**6);
beta_eq..   beta_v  =E= 15.2518 - 15.6875/Tbr - 13.4721*log(Tbr) + 0.43577*(Tbr**6);
omega_eq..  omega   =E= alpha/beta_v;
Cpa_eq..    Cp_a    =E= sum(i,n(i)*Property(i,'CP_a')) - 37.93 + (sum(i,n(i)*Property(i,'CP_b'))+0.21)*Tavg + (sum(i,n(i)*Property(i,'CP_c'))-3.91*(10**(-4)))*(Tavg**2) + (sum(i,n(i)*Property(i,'CP_d'))+2.06*(10**(-7)))*(Tavg**3);
Cpl_eq..    Cp_l    =E= 1/4.1868 * (Cp_a + 8.314*(1.45 + 0.45/(1-Tavgr) + 0.25*omega*(17.11 + 25.2*(((1-Tavgr)**(1/3))/Tavgr) + 1.742/(1-Tavgr))));
* Structural constraints for molecule design
ester_ub_eq..   1 =G= sum(i,n(i)*Property(i,'N_ester'));
ether_ub_eq..   2 =G= sum(i,n(i)*Property(i,'N_ether'));
amine_ub_eq..   2 =G= sum(i,n(i)*Property(i,'N_amine'));
carbon_ub_eq..  6 =G= sum(i,n(i)*Property(i,'N_carbon'));
carbon_lb_eq..  2 =L= sum(i,n(i)*Property(i,'N_carbon'));
* integer cut to filter repeated choices
int_cut(IC)$ (dym(IC)).. sum(i,abs(integercut(IC,i)-n(i)))=g=1;


Model APO_Project /all/;

OPTION SYSOUT=ON
option MIP=cplex;
option MINLP=BARON;
option optcr=0;
option optca=1E-5;

APO_Project.OPTFILE=1;

* initial guess (methanol)
integercut(IC,i)    = 0;
integercut(IC,'CH2') = 1;
integercut(IC,'CH2NH2(amine)')= 1;
integercut(IC,'OH') = 1;
n.l('CH2NH2(amine)') = 1;
n.l('CH2')= 1;
n.l('OH') = 1;

dym(IC)   = No;

Tm.l      = T_m0 * log(sum(i,n.l(i)*Property(i,'Tm')));
Tb.l      = T_b0 * log(sum(i,n.l(i)*Property(i,'Tb')));
Tc.l      = Tb.l/(0.584 + 0.965*sum(i,n.l(i)*Property(i,'Tc')) - (sum(i,n.l(i)*Property(i,'Tc')))**2);
Mw.l      = sum(i,n.l(i)*Property(i,'Mr'));
Vm.l      = V_m0 + sum(i,n.l(i)*Property(i,'Vm'));
Pc.l      = 1/(0.113 + 0.0032*sum(i,n.l(i)*Property(i,'a')) - sum(i,n.l(i)*Property(i,'Pc')))**2;
rho.l     = Mw.l * (10**6) / Vm.l;
Tbr.l     = Tb.l/Tc.l;
Tavgr.l   = Tavg/Tc.l;
alpha.l   = -5.97214 - log(Pc.l/1.013) + 6.09648/Tbr.l + 1.28862*log(Tbr.l) - 0.169347*(Tbr.l**6);
beta_v.l  = 15.2518 - 15.6875/Tbr.l - 13.4721*log(Tbr.l) + 0.43577*(Tbr.l**6);

ALIAS(IC,cc);
LOOP(cc,       
        SOLVE APO_Project MINIMIZE obj USING MINLP;
        integercut(cc,i) = n.l(i);
* store the properties for later displacement
        number_of_groups(i,cc)          = n.l(i);
        heat_capacity(cc)               = Cp_l.l;
        liquid_density(cc)              = rho.l;
        objective_value(cc)             = obj.l;
        melting_temperature(cc)         = Tm.l;
        boiling_temperature(cc)         = Tb.l;
        relative_energy_difference(cc)  = RED.l;        
        dym(cc)=YES;
);

display number_of_groups, objective_value, melting_temperature, boiling_temperature, relative_energy_difference, liquid_density, heat_capacity;












