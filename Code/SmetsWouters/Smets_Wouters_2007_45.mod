@#ifndef SW_ME
@#define SW_ME = 0
@#endif
@#ifndef MH_REPLIC
@#define MH_REPLIC = 5000
@#endif
@#ifndef N_DRAWS
@#define N_DRAWS = 500
@#endif
/*
 * Smets and Wouters (2007), "Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach",
 * AER 97(3), 586-606, in the Dynare 4.5+ implementation of Johannes Pfeifer
 * (github.com/JohannesPfeifer/DSGE_mod, Smets_Wouters_2007/Smets_Wouters_2007_45.mod, commit 26c21f7),
 * adapted for the reliability toolbox. Two variants, chosen by the macro SW_ME (run_sw_variant.m
 * passes it on the command line, so this file is never edited):
 *
 *   SW_ME = 0  BASELINE: the published model, seven structural shocks, seven noise-free observables.
 *   SW_ME = 1  CALIBRATED MEASUREMENT ERROR on output growth and its components (dy, dc, dinve), on
 *              real wage growth dw and on hours labobs; none on inflation and the interest rate.
 *              The five errors are extra shocks me_* in the observation equations (the toolbox does
 *              not use Dynare's M_.H); their standard deviations are FIXED (not estimated) and set
 *              in MATLAB code below from the data:
 *                  SW_ME_SHARE   share of the observable's sample variance (default 0.10, i.e. an
 *                                s.d. of sqrt(0.10) = 32% of the series s.d.), or
 *                  SW_ME_SD      one absolute s.d. for all five (overrides SW_ME_SHARE)
 *              read with getenv when the driver runs, so the size can be changed without editing.
 *
 * In both variants the 36 parameters of the replication file (7 shock s.d. + 29 structural) are
 * ESTIMATED: posterior mode by csminwel (mode_compute = 4) started at the published mode
 * (usmodel_mode.mat), then a SHORT random-walk Metropolis-Hastings chain (MH_REPLIC draws per chain,
 * 2 chains, 20% burn-in). The chain is for illustration only, not a converged posterior. Fixed:
 * ctou = 0.025, clandaw = 1.5, cg = 0.18, curvp = curvw = 10 (as in the replication file).
 * Data: usmodel_data.mat (230 quarters, 1947Q3-2004Q4; get_sw_data.m downloads it).
 *
 * Without measurement error the model is square and noise-free: the shocks are recovered almost
 * exactly (population R^2 of rrf and of the gap close to 1), and the four differenced observables
 * put four transfer-function zeros at z = 1, so the steady-state Riccati equation has no
 * stabilising solution (rho(A - KC) -> 1); the toolbox then falls back to a slow fixed-point
 * iteration and reports the Riccati flag (cell 19) as not converged.
 *
 * States of interest: the natural (flex-price) real rate rrf and the output gap defined as the
 * linear combination y - yf (column 4 of M_.CellList_States), with fixed-lag smoothers m = 1, 4, 20
 * and, for rrf, 4- and 20-quarter changes, as for the UC3 trends.
 *
 * Original header notes of the replication file (kept):
 *  - The consumption Euler equation in the paper, equation (2), premultiplies the risk premium process
 *    by c_3; in the code this prefactor is 1, so b here is c_3*eps^b.
 *  - Equation (8) misses a (1+cbetabar*cgamma) in i_2 in front of q_t (Del Negro/Schorfheide 2012).
 *  - Chib/Ramamurthy (2010): the mode reported in the paper is not the mode.
 *  - At the prior mean [cmap,crhopinf] and [cmaw,crhow] are pairwise collinear.
 *  - Table 1B: rho_ga follows N(0.5,0.25^2) in the code, not B(0.5,0.2^2).
 *
 * Copyright (C) 2007-2013 Frank Smets and Raf Wouters, (C) 2013-15 Johannes Pfeifer (GPL v3+).
 */

var labobs      ${lHOURS}$      (long_name='log hours worked')
    robs        ${FEDFUNDS}$    (long_name='Federal funds rate')
    pinfobs     ${dlP}$         (long_name='Inflation')
    dy          ${dlGDP}$       (long_name='Output growth rate')
    dc          ${dlCONS}$      (long_name='Consumption growth rate')
    dinve       ${dlINV}$       (long_name='Investment growth rate')
    dw          ${dlWAG}$       (long_name='Wage growth rate')
    ewma        ${\eta^{w,aux}}$ (long_name='Auxiliary wage markup moving average variable')
    epinfma     ${\eta^{p,aux}}$ (long_name='Auxiliary price markup moving average variable')
    zcapf       ${z^{flex}}$    (long_name='Capital utilization rate flex price economy')
    rkf         ${r^{k,flex}}$  (long_name='rental rate of capital flex price economy')
    kf          ${k^{s,flex}}$  (long_name='Capital services flex price economy')
    pkf         ${q^{flex}}$    (long_name='real value of existing capital stock flex price economy')
    cf          ${c^{flex}}$    (long_name='Consumption flex price economy')
    invef       ${i^{flex}}$    (long_name='Investment flex price economy')
    yf          ${y^{flex}}$    (long_name='Output flex price economy')
    labf        ${l^{flex}}$    (long_name='hours worked flex price economy')
    wf          ${w^{flex}}$    (long_name='real wage flex price economy')
    rrf         ${r^{flex}}$    (long_name='real interest rate flex price economy')
    mc          ${\mu_p}$       (long_name='gross price markup')
    zcap        ${z}$           (long_name='Capital utilization rate')
    rk          ${r^{k}}$       (long_name='rental rate of capital')
    k           ${k^{s}}$       (long_name='Capital services')
    pk          ${q}$           (long_name='real value of existing capital stock')
    c           ${c}$           (long_name='Consumption')
    inve        ${i}$           (long_name='Investment')
    y           ${y}$           (long_name='Output')
    lab         ${l}$           (long_name='hours worked')
    pinf        ${\pi}$         (long_name='Inflation')
    w           ${w}$           (long_name='real wage')
    r           ${r}$           (long_name='nominal interest rate')
    a           ${\varepsilon_a}$       (long_name='productivity process')
    b           ${c_2*\varepsilon_t^b}$ (long_name='Scaled risk premium shock')
    g           ${\varepsilon^g}$       (long_name='Exogenous spending')
    qs          ${\varepsilon^i}$       (long_name='Investment-specific technology')
    ms          ${\varepsilon^r}$       (long_name='Monetary policy shock process')
    spinf       ${\varepsilon^p}$       (long_name='Price markup shock process')
    sw          ${\varepsilon^w}$       (long_name='Wage markup shock process')
    kpf         ${k^{flex}}$            (long_name='Capital stock flex price economy')
    kp          ${k}$           (long_name='Capital stock')
    ;

varexo ea       ${\eta^a}$      (long_name='productivity shock')
    eb          ${\eta^b}$      (long_name='risk premium shock')
    eg          ${\eta^g}$      (long_name='Spending shock')
    eqs         ${\eta^i}$      (long_name='Investment-specific technology shock')
    em          ${\eta^m}$      (long_name='Monetary policy shock')
    epinf       ${\eta^{p}}$    (long_name='Price markup shock')
    ew          ${\eta^{w}}$    (long_name='Wage markup shock')
@#if SW_ME == 1
    me_dy       ${\eta^{me,y}}$ (long_name='measurement error output growth')
    me_dc       ${\eta^{me,c}}$ (long_name='measurement error consumption growth')
    me_dinve    ${\eta^{me,i}}$ (long_name='measurement error investment growth')
    me_dw       ${\eta^{me,w}}$ (long_name='measurement error wage growth')
    me_lab      ${\eta^{me,l}}$ (long_name='measurement error hours')
@#endif
        ;

parameters curvw ${\varepsilon_w}$  (long_name='Curvature Kimball aggregator wages')
    cgy         ${\rho_{ga}}$       (long_name='Feedback technology on exogenous spending')
    curvp       ${\varepsilon_p}$   (long_name='Curvature Kimball aggregator prices')
    constelab   ${\bar l}$          (long_name='steady state hours')
    constepinf  ${\bar \pi}$        (long_name='steady state inflation rate')
    constebeta  ${100(\beta^{-1}-1)}$ (long_name='time preference rate in percent')
    cmaw        ${\mu_w}$           (long_name='coefficient on MA term wage markup')
    cmap        ${\mu_p}$           (long_name='coefficient on MA term price markup')
    calfa       ${\alpha}$          (long_name='capital share')
    czcap       ${\psi}$            (long_name='capacity utilization cost')
    csadjcost   ${\varphi}$         (long_name='investment adjustment cost')
    ctou        ${\delta}$          (long_name='depreciation rate')
    csigma      ${\sigma_c}$        (long_name='risk aversion')
    chabb       ${\lambda}$         (long_name='external habit degree')
    ccs         ${d_4}$             (long_name='Unused parameter')
    cinvs       ${d_3}$             (long_name='Unused parameter')
    cfc         ${\phi_p}$          (long_name='fixed cost share')
    cindw       ${\iota_w}$         (long_name='Indexation to past wages')
    cprobw      ${\xi_w}$           (long_name='Calvo parameter wages')
    cindp       ${\iota_p}$         (long_name='Indexation to past prices')
    cprobp      ${\xi_p}$           (long_name='Calvo parameter prices')
    csigl       ${\sigma_l}$        (long_name='Frisch elasticity')
    clandaw     ${\phi_w}$          (long_name='Gross markup wages')
    crdpi       ${r_{\Delta \pi}}$  (long_name='Unused parameter')
    crpi        ${r_{\pi}}$         (long_name='Taylor rule inflation feedback')
    crdy        ${r_{\Delta y}}$    (long_name='Taylor rule output growth feedback')
    cry         ${r_{y}}$           (long_name='Taylor rule output level feedback')
    crr         ${\rho}$            (long_name='interest rate persistence')
    crhoa       ${\rho_a}$          (long_name='persistence productivity shock')
    crhoas      ${d_2}$             (long_name='Unused parameter')
    crhob       ${\rho_b}$          (long_name='persistence risk premium shock')
    crhog       ${\rho_g}$          (long_name='persistence spending shock')
    crhols      ${d_1}$             (long_name='Unused parameter')
    crhoqs      ${\rho_i}$          (long_name='persistence risk premium shock')
    crhoms      ${\rho_r}$          (long_name='persistence monetary policy shock')
    crhopinf    ${\rho_p}$          (long_name='persistence price markup shock')
    crhow       ${\rho_w}$          (long_name='persistence wage markup shock')
    ctrend      ${\bar \gamma}$     (long_name='net growth rate in percent')
    cg          ${\frac{\bar g}{\bar y}}$     (long_name='steady state exogenous spending share')
    ;

// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

// estimated parameters initialisation (replaced by the mode file in the estimation command)
calfa=.24;
cbeta=.9995;
csigma=1.5;
cfc=1.5;
cgy=0.51;

csadjcost= 6.0144;
chabb=    0.6361;
cprobw=   0.8087;
csigl=    1.9423;
cprobp=   0.6;
cindw=    0.3243;
cindp=    0.47;
czcap=    0.2696;
crpi=     1.488;
crr=      0.8762;
cry=      0.0593;
crdy=     0.2347;

crhoa=    0.9977;
crhob=    0.5799;
crhog=    0.9957;
crhols=   0.9928;
crhoqs=   0.7165;
crhoas=1;
crhoms=0;
crhopinf=0;
crhow=0;
cmap = 0;
cmaw  = 0;

constelab=0;

constepinf=0.7;
constebeta=0.7420;
ctrend=0.3982;

model(linear);
//deal with parameter dependencies; taken from usmodel_stst.mod
#cpie=1+constepinf/100;         %gross inflation rate
#cgamma=1+ctrend/100 ;          %gross growth rate
#cbeta=1/(1+constebeta/100);    %discount factor

#clandap=cfc;                   %fixed cost share/gross price markup
#cbetabar=cbeta*cgamma^(-csigma);   %growth-adjusted discount factor in Euler equation
#cr=cpie/(cbeta*cgamma^(-csigma));  %steady state gross real interest rate
#crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou); %R^k_{*}: steady state rental rate
#cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));      %steady state real wage
#cikbar=(1-(1-ctou)/cgamma);        %(1-k_1) in equation LOM capital, equation (8)
#cik=(1-(1-ctou)/cgamma)*cgamma;    %i_k: investment-capital ratio
#clk=((1-calfa)/calfa)*(crk/cw);    %labor to capital ratio
#cky=cfc*(clk)^(calfa-1);           %k_y: steady state output ratio
#ciy=cik*cky;                       %investment-output ratio
#ccy=1-cg-cik*cky;                  %consumption-output ratio
#crkky=crk*cky;                     %z_y=R_{*}^k*k_y
#cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy; %W^{h}_{*}*L_{*}/C_{*} used in c_2 in equation (2)
#cwly=1-crk*cky;                    %unused parameter

#conster=(cr-1)*100;                %steady state federal funds rate ($\bar r$)

// flexible economy

          [name='FOC labor with mpl expressed as function of rk and w, flex price economy']
	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      [name='FOC capacity utilization, flex price economy']
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
          [name='Firm FOC capital, flex price economy']
	      rkf =  (wf)+labf-kf ;
          [name='Definition capital services, flex price economy']
	      kf =  kpf(-1)+zcapf ;
          [name='Investment Euler Equation, flex price economy']
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          [name='Arbitrage equation value of capital, flex price economy']
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      [name='Consumption Euler Equation, flex price economy']
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      [name='Aggregate Resource Constraint, flex price economy']
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      [name='Aggregate Production Function, flex price economy']
          yf = cfc*( calfa*kf+(1-calfa)*labf +a );
          [name='Wage equation, flex price economy']
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
          [name='Law of motion for capital, flex price economy (see header notes)']
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy
          [name='FOC labor with mpl expressed as function of rk and w, SW Equation (9)']
	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      [name='FOC capacity utilization, SW Equation (7)']
          zcap =  (1/(czcap/(1-czcap)))* rk ;
          [name='Firm FOC capital, SW Equation (11)']
	      rk =  w+lab-k ;
          [name='Definition capital services, SW Equation (6)']
	      k =  kp(-1)+zcap ;
          [name='Investment Euler Equation, SW Equation (3)']
	      inve = (1/(1+cbetabar*cgamma))* (inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          [name='Arbitrage equation value of capital, SW Equation (4)']
          pk = -r+pinf(1)-0*b
                + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b
                + (crk/(crk+(1-ctou)))*rk(1)
                + ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      [name='Consumption Euler Equation, SW Equation (2)']
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1)
                + (1/(1+chabb/cgamma))*c(+1)
                +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1))
                - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
	      [name='Aggregate Resource Constraint, SW Equation (1)']
          y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
          [name='Aggregate Production Function, SW Equation (5)']
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
          [name='New Keynesian Phillips Curve, SW Equation (10)']
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1)
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ;
	      [name='Wage Phillips Curve, SW Equation (13), with (12) plugged for mu_w']
          w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w)
               + 1*sw ;
          [name='Taylor rule, SW Equation (14)']
	      r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(y-yf)
               +crdy*(y-yf-y(-1)+yf(-1))
               +crr*r(-1)
               +ms  ;
          [name='Law of motion for productivity']
	      a = crhoa*a(-1)  + ea;
          [name='Law of motion for risk premium']
	      b = crhob*b(-1) + eb;
          [name='Law of motion for spending process']
	      g = crhog*(g(-1)) + eg + cgy*ea;
	      [name='Law of motion for investment specific technology shock process']
          qs = crhoqs*qs(-1) + eqs;
          [name='Law of motion for monetary policy shock process']
	      ms = crhoms*ms(-1) + em;
          [name='Law of motion for price markup shock process']
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
          [name='Law of motion for wage markup shock process']
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew;
          [name='Law of motion for capital, SW Equation (8) (see header notes)']
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

// measurement equations (SW_ME = 1: measurement error on dy, dc, dinve, dw, labobs)
@#if SW_ME == 1
[name='Observation equation output']
dy=y-y(-1)+ctrend + me_dy;
[name='Observation equation consumption']
dc=c-c(-1)+ctrend + me_dc;
[name='Observation equation investment']
dinve=inve-inve(-1)+ctrend + me_dinve;
[name='Observation equation real wage']
dw=w-w(-1)+ctrend + me_dw;
[name='Observation equation hours worked']
labobs = lab + constelab + me_lab;
@#else
[name='Observation equation output']
dy=y-y(-1)+ctrend;
[name='Observation equation consumption']
dc=c-c(-1)+ctrend;
[name='Observation equation investment']
dinve=inve-inve(-1)+ctrend;
[name='Observation equation real wage']
dw=w-w(-1)+ctrend;
[name='Observation equation hours worked']
labobs = lab + constelab;
@#endif
[name='Observation equation inflation']
pinfobs = 1*(pinf) + constepinf;
[name='Observation equation interest rate']
robs =    1*(r) + conster;
end;

steady_state_model;
dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs = constepinf;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
labobs = constelab;
end;

shocks;
var ea;
stderr 0.4618;
var eb;
stderr 1.8513;
var eg;
stderr 0.6090;
var eqs;
stderr 0.6017;
var em;
stderr 0.2397;
var epinf;
stderr 0.1455;
var ew;
stderr 0.2089;
@#if SW_ME == 1
// measurement errors: placeholders, set from the data below (SW_ME_SHARE / SW_ME_SD)
var me_dy;    stderr 0.1;
var me_dc;    stderr 0.1;
var me_dinve; stderr 0.1;
var me_dw;    stderr 0.1;
var me_lab;   stderr 0.1;
@#endif
end;

@#if SW_ME == 1
// CALIBRATED measurement-error standard deviations from the data (native MATLAB statements).
// Default: variance share SW_ME_SHARE = 0.10 of each observable's sample variance (full sample of
// usmodel_data.mat); SW_ME_SD sets one absolute s.d. for all five instead. They are not in
// estimated_params, so estimation keeps M_.Sigma_e of the me_* shocks as set here; only the seven
// structural shock s.d. are estimated.
me_share = str2double(getenv('SW_ME_SHARE'));  if isnan(me_share), me_share = 0.10; end
me_sd_abs = str2double(getenv('SW_ME_SD'));
me_tmp = load('usmodel_data.mat');
me_names = {'me_dy', 'me_dc', 'me_dinve', 'me_dw', 'me_lab'};
me_series = {me_tmp.dy, me_tmp.dc, me_tmp.dinve, me_tmp.dw, me_tmp.labobs};
for me_i = 1:5
    if isnan(me_sd_abs), me_sd = sqrt(me_share) * std(me_series{me_i}); else, me_sd = me_sd_abs; end
    M_.Sigma_e(strcmp(M_.exo_names, me_names{me_i}), strcmp(M_.exo_names, me_names{me_i})) = me_sd^2;
    fprintf('measurement error %-9s s.d. = %.4f (series s.d. %.3f)\n', me_names{me_i}, me_sd, std(me_series{me_i}));
end
M_.sw_me = struct('share', me_share, 'sd_abs', me_sd_abs, 'names', {me_names}, 'sd', sqrt(diag(M_.Sigma_e(8:12, 8:12)))');
clear me_tmp me_series me_i me_sd
@#else
M_.sw_me = struct('share', 0, 'sd_abs', NaN, 'names', {{}}, 'sd', []);
@#endif

estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
stderr ea,0.4618,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eb,0.1818513,0.025,5,INV_GAMMA_PDF,0.1,2;
stderr eg,0.6090,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eqs,0.46017,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr em,0.2397,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr epinf,0.1455,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr ew,0.2089,0.01,3,INV_GAMMA_PDF,0.1,2;
crhoa,.9676 ,.01,.9999,BETA_PDF,0.5,0.20;
crhob,.2703,.01,.9999,BETA_PDF,0.5,0.20;
crhog,.9930,.01,.9999,BETA_PDF,0.5,0.20;
crhoqs,.5724,.01,.9999,BETA_PDF,0.5,0.20;
crhoms,.3,.01,.9999,BETA_PDF,0.5,0.20;
crhopinf,.8692,.01,.9999,BETA_PDF,0.5,0.20;
crhow,.9546,.001,.9999,BETA_PDF,0.5,0.20;
cmap,.7652,0.01,.9999,BETA_PDF,0.5,0.2;
cmaw,.8936,0.01,.9999,BETA_PDF,0.5,0.2;
csadjcost,6.3325,2,15,NORMAL_PDF,4,1.5;
csigma,1.2312,0.25,3,NORMAL_PDF,1.50,0.375;
chabb,0.7205,0.001,0.99,BETA_PDF,0.7,0.1;
cprobw,0.7937,0.3,0.95,BETA_PDF,0.5,0.1;
csigl,2.8401,0.25,10,NORMAL_PDF,2,0.75;
cprobp,0.7813,0.5,0.95,BETA_PDF,0.5,0.10;
cindw,0.4425,0.01,0.99,BETA_PDF,0.5,0.15;
cindp,0.3291,0.01,0.99,BETA_PDF,0.5,0.15;
czcap,0.2648,0.01,1,BETA_PDF,0.5,0.15;
cfc,1.4672,1.0,3,NORMAL_PDF,1.25,0.125;
crpi,1.7985,1.0,3,NORMAL_PDF,1.5,0.25;
crr,0.8258,0.5,0.975,BETA_PDF,0.75,0.10;
cry,0.0893,0.001,0.5,NORMAL_PDF,0.125,0.05;
crdy,0.2239,0.001,0.5,NORMAL_PDF,0.125,0.05;
constepinf,0.7,0.1,2.0,GAMMA_PDF,0.625,0.1;//20;
constebeta,0.7420,0.01,2.0,GAMMA_PDF,0.25,0.1;//0.20;
constelab,1.2918,-10.0,10.0,NORMAL_PDF,0.0,2.0;
ctrend,0.3982,0.1,0.8,NORMAL_PDF,0.4,0.10;
cgy,0.05,0.01,2.0,NORMAL_PDF,0.5,0.25;
calfa,0.24,0.01,1.0,NORMAL_PDF,0.3,0.05;
end;

varobs dy dc dinve labobs pinfobs dw robs;

// ReliabilityToolbox on the path (run_sw_variant.m adds it; this covers a manual run from
// Code/SmetsWouters with the repo layout Code/SmetsWouters -> Code/ReliabilityToolbox)
if ~exist('compute_state_reliability_diagnostics', 'file') && exist('../ReliabilityToolbox', 'dir')
    addpath('../ReliabilityToolbox');   // not genpath: octave_compat/ is Octave-only
end

set_dynare_seed(1);

// Estimation: posterior mode by csminwel started at the published mode, then a SHORT MH chain
// (MH_REPLIC draws per chain, 2 chains, 20% burn-in; illustrative, not converged). Sample
// (first_obs = 1: all 230 quarters, 1947Q3-2004Q4), presample, initialisation and proposal scale as
// in the replication file; Smets and Wouters (2007) report estimates for 1966Q1-2004Q4
// (first_obs = 71 with presample = 4).
// KNOWN GAP: the mode found here is not exactly a mode of the likelihood the toolbox evaluates in
// model_assessment (its "mode condition" is about 3 posterior s.d. in the smoke tests). Matching
// the initialisation (lik_init = 1, presample = 0) did not remove the gap (2 s.d.) and left
// csminwel with a non-positive-definite Hessian; newrat (mode_compute = 5) took an hour.
estimation(optim=('MaxIter',1000),datafile=usmodel_data,mode_file=usmodel_mode,mode_compute=4,first_obs=1, presample=4,lik_init=2,prefilter=0,mh_replic=@{MH_REPLIC},mh_nblocks=2,mh_jscale=0.20,mh_drop=0.2, nograph);

// Smoother and filter at the posterior mean (also sets the bookkeeping the toolbox needs:
// bayestopt_.mf, smoother_var_list, smoother_restrict_columns); filtered_vars gives oo_.UpdatedVariables.
calib_smoother(datafile=usmodel_data,first_obs=1,kalman_algo=2,prefilter=0,filtered_vars);

original_estim_params_=estim_params_;

stoch_simul(nograph);

save tmp.mat

// toolbox bookkeeping (as the UC3 / HLW .mod files)
M_.observable_names = dataset_.name ;
M_.data = dataset_info.rawdata(~any(isnan(dataset_info.rawdata')),: );
// States to diagnose: {name, smoother lags m (0 = filter), change horizons h, definition}. The output
// gap is not a model variable: column 4 defines it as y - yf. Lags and changes as for UC3.
M_.CellList_States    = {'ygap', [0 1 4 20], [],     'y - yf';
                         'rrf',  [0 1 4 20], [4 20], ''};
M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
M_.var_list_ = var_list_;
M_.verbose = false;
nstates = length(oo_.dr.state_var);
M_.X_InitialDistribution.X0 = zeros(nstates,1);
M_.X_InitialDistribution.B0 = eye(nstates)*100;
oo_.bayestopt_smoother_var_list=bayestopt_.smoother_var_list;
oo_.bayestopt_smoother_restrict_columns=bayestopt_.smoother_restrict_columns;
oo_.bayestopt_mf=bayestopt_.mf;
oo_.dr.qz_crit = 1+1e-6;

save tmp.mat

// Staged saves as in the UC3 .mod: (1) after Dynare, (2) after the prior stage, (3) final;
// save_results_split keeps the draw sets in SW_results_prior.mat / SW_results_posterior.mat.
save_results_split('SW_results.mat');

tic;
[oo_.prior_function_results,     prior_draws,     prior_info]     = run_reliability_draws('prior',     @{N_DRAWS}, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;
save_results_split('SW_results.mat');    % staged save (2)
tic;
[oo_.posterior_function_results, posterior_draws, posterior_info] = run_reliability_draws('posterior', @{N_DRAWS}, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;

// Direct (single) run of the diagnostics engine at the posterior mean
tic;
[LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, TV_state_estimate, TV_surprises, TV_log_llhood] = compute_state_reliability_diagnostics(M_.observable_names, M_.CellList_States, oo_, M_, M_.data);
ric_str = 'NOT converged'; if ModelMatrices.KF_ric_flag.converged, ric_str = 'converged'; end
fprintf('engine pass at the posterior mean: %.1f s; Riccati %s (method %s, rho(A-KC) = %.6f)\n', toc, ric_str, ModelMatrices.KF_ric_flag.method, ModelMatrices.KF_ric_flag.rho_cl);

save_results_split('SW_results.mat');    % staged save (3): final
