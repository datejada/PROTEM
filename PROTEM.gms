$Title PROTEM   PROactive Transmission Expansion Model

$onText
For more details please refer to Chapter 3, of the following book:
AAA

Contributed by
Dra. Sonja Wogrin          , IEEE Senior Member, email: 
Dr.  Salvador Pineda       , IEEE Senior Member, email: 
Dr.  Diego A. Tejada-Arango, IEEE        Member, email: 

We do request that publications derived from the use of the developed GAMS code
explicitly acknowledge that fact by citing
AAA
DOI: BBB
$offText

* general options
$EOLCOM //
$OnEmpty OnMulti OffListing
OPTION optcr    =    0 ; // tolerance to solve MIP until IntGap < OptcR
OPTION threads  =   -1 ; // number of cores
;
* Indices, sets, parameters, and variables
SETS
g       "Generating unit index                   "
t       "Time period index                       "
n       "Node index                              "
GN(g,n) "Set of generating units at a given node "
LN(n,n) "Set of possible transmission lines      "
;

ALIAS(n,nn);

PARAMETERS
SB               "Base power (MW)                                                       "
C        (g    ) "Linear cost parameter of generating unitg(€/MWh)                      "
I        (g    ) "Annualized investment cost of generating unitg(€/MW)                  "
IT       (n,n  ) "Annualized transmission investment cost of line between n and nn(€/MW)"
X        (n,n  ) "Reactance of line between n and nn(p.u.)                              "
RHO      (g,t  ) "Capacity factor of generating unit g and time t (p.u.)                "
S0       (t,n  ) "Demand slope at nodenand timet(MW^2/€)                                "
D0       (t,n  ) "Demand intercept at nodenand timet(MW)                                "
PHI      (g    ) "Conjectured-price response of generating unitg(€/MW2)                 "
MLB_p    (g,t  ) "Auxiliary large constants used for linearization                      "
MUB_p    (g,t  ) "Auxiliary large constants used for linearization                      "
MLB_beta (g,t  ) "Auxiliary large constants used for linearization                      "
MUB_beta (g,t  ) "Auxiliary large constants used for linearization                      "
MLB_f    (t,n,n) "Auxiliary large constants used for linearization                      "
MUB_f    (t,n,n) "Auxiliary large constants used for linearization                      "
MLB_gamma(t,n,n) "Auxiliary large constants used for linearization                      "
MUB_gamma(t,n,n) "Auxiliary large constants used for linearization                      "
MLB_d    (t,n  ) "Auxiliary large constants used for linearization                      "
MLB_alpha(t,n  ) "Auxiliary large constants used for linearization                      "
;
FREE VARIABLES
ul_of          "Upper level objective function     "
dummy          "Dummy       objective function     "
f      (t,n,n) "Power flow at time t from node n to node nn(MW)"
theta  (t,n  ) "Voltage angle at timetand noden(rad)"
lambda (t,n  ) "Electricity price at time t and node n (€/MWh)"
kappa  (t,n,n) "Dual of definition of power flow (MW)"
;
POSITIVE VARIABLES
p      (g,t  ) "Output of generating unit g in time t             (MW)"
p_bar  (g    ) "Capacity investment of generating unitg           (MW)"
l      (n,n  ) "Capacity investment in line from node n to node nn(MW)"
d      (t,n  ) "Satisfied demand in time t                        (MW)"
betaLB (g,t  ) "Dual of lower bound on production                 (MW)"
betaUB (g,t  ) "Dual of upper bound on production                 (MW)"
alphaLB(t,n  ) "Dual of lower bound on demand                     (MW)"
gammaLB(t,n,n) "Dual of lower bound on power flow                 (MW)"
gammaUB(t,n,n) "Dual of upper bound on power flow                 (MW)"
;
BINARY VARIABLES
bLB_p(g,t  ) "Auxiliary binary variables used for linearization ({0,1})"
bUB_p(g,t  ) "Auxiliary binary variables used for linearization ({0,1})"
bLB_f(t,n,n) "Auxiliary binary variables used for linearization ({0,1})"
bUB_f(t,n,n) "Auxiliary binary variables used for linearization ({0,1})"
bLB_d(t,n  ) "Auxiliary binary variables used for linearization ({0,1})"
;
* Constraints and Model definition
EQUATIONS
eDummy      "Dummy objective function                            "
eObjFunUL   "Upper level objective function                 (10a)"
eDemBal     "Demand balance                                 (8 ) "
eUBProd     "Upper bound generating unit production         (3b) "
eLBProd     "Lower bound generating unit production         (3c) "
eLBFlow     "Lower bound      power flow                    (7c) "
eUBFlow     "Upper bound      power flow                    (7d) "
eDFFlow     "Definition bound power flow                    (7e) "
eLBDemd     "Lower bound demand                             (5b) "
edL_dd      "Derivative of Lagrangian with respect to d     (5a) "
edL_dp      "Derivative of Lagrangian with respect to p     (3a) "
edL_dp_bar  "Derivative of Lagrangian with respect to p_bar (3b) "
edL_dtheta  "Derivative of Lagrangian with respect to theta (7a) "
edL_df      "Derivative of Lagrangian with respect to f     (7b) "
eComLBd     "complementarity                                (5d) "
eComLBp     "complementarity                                (3f) "
eComUBp     "complementarity                                (3g) "
eComLBgamma "complementarity                                (7g) "
eComUBgamma "complementarity                                (7h) "
;

eDummy.. dummy =e= 0 ;

eObjFunUL.. ul_of =e=
             + SUM[(n,nn)$[LN(n,nn)], IT(n,nn)*l    (n,nn)]
             + SUM[   g             , I (g   )*p_bar(g   )]
             + SUM[(g,t)            , C (g   )*p    (g,t )]
;

eDemBal(t,n)..
 + SUM[g $[GN(g ,n )],p(g,t    )]
 + SUM[nn$[LN(nn,n )],f(t,nn,n )]
 - SUM[nn$[LN(n ,nn)],f(t,n ,nn)]
=e=
 + d(t,n)
;

eUBProd(g,t).. RHO(g,t)*p_bar(g)=g= p(g,t) ;
eLBProd(g,t)..p(g,t)=g= 0                 ;
eLBDemd(t,n)..d(t,n)=g= 0                 ;

eLBFlow(t,n,nn)$[LN(n,nn)]..0 =g= -f(t,n,nn)-l(n,nn)                                  ;
eUBFlow(t,n,nn)$[LN(n,nn)]..0 =g= +f(t,n,nn)-l(n,nn)                                  ;
eDFFlow(t,n,nn)$[LN(n,nn)].. f(t,n,nn)         =e=[theta(t,n)-theta(t,nn)]*SB/X(n,nn) ;

edL_dd    (t,n).. lambda(t,n)-[D0(t,n)-d(t,n)]/S0(t,n)-alphaLB(t,n)                       =e= 0 ;
edL_dp    (g,t).. C(g)-SUM[n$[GN(g,n)],lambda(t,n)]+PHI(g)*p(g,t)-betaLB(g,t)+betaUB(g,t) =e= 0 ;
edL_dp_bar(  g).. I(g)-SUM[t          ,betaUB(g,t)]                                       =e= 0 ;
edL_dtheta(t,n)..+SUM[nn$[LN(nn,n )],SB*kappa(t,nn,n )/X(nn,n )]
                 +SUM[nn$[LN(n ,nn)],SB*kappa(t,n ,nn)/X(n ,nn)]                          =e= 0 ;
edL_df(t,n,nn)$[LN(n,nn)].. + lambda (t,n   ) - lambda (t,  nn)
                            - gammaLB(t,n,nn) - gammaUB(t,n,nn) + kappa(t,n,nn)           =e= 0 ;


eComLBd    (t,n   )           ..           d    (t,n)        *alphaLB(t,n) =e= 0 ;
eComLBp    (g,t   )           ..           p    (g,t)        *betaLB (g,t) =e= 0 ;
eComUBp    (g,t   )           .. [RHO(g,t)*p_bar(g  )-p(g,t)]*betaUB (g,t) =e= 0 ;
eComLBgamma(t,n,nn)$[LN(n,nn)]..[l(n,nn)+f(t,n,nn)]*gammaLB(t,n,nn)        =e= 0 ;
eComUBgamma(t,n,nn)$[LN(n,nn)]..[l(n,nn)-f(t,n,nn)]*gammaUB(t,n,nn)        =e= 0 ;

MODEL  mPROTEM_LL_NLP
/
* objective function
 eDummy ,
* constraints
 eDemBal, eLBDemd, eLBProd, eUBProd, eLBFlow, eUBFlow, eDFFlow,
* derivatives of Lagrangian
 edL_dd , edL_dp , edL_dp_bar, edL_dtheta , edL_df,
* complementarity conditions
 eComLBd, eComLBp, eComUBp   , eComLBgamma, eComUBgamma
/
;

MODEL  mPROTEM_LL_MCP
/
* constraints
 eLBDemd.alphaLB, eDemBal.lambda, eLBProd.betaLB, eUBProd.betaUB,
 eLBFlow.gammaLB, eUBFlow.gammaUB, eDFFlow.kappa,
* derivatives of Lagrangian
 edL_dd.d, edL_dp.p , edL_dp_bar.p_bar, edL_dtheta.theta , edL_df.f
/
;

MODEL  mPROTEM_UL_NLP  / mPROTEM_LL_NLP -eDummy +eObjFunUL / ;
MODEL  mPROTEM_UL_MPEC / mPROTEM_LL_MCP +eObjFunUL / ;

* Input data
SETS
g       /g01*g02/
t       /t01    /
n       /n01*n02/
GN(g,n) /g01.n01, g02.n02/
LN(n,n) /n01.n02         /

;
TABLE tGDATA(g,*) 'generator data'
        LinCost InvCost
*       [€/MWh]  [€/MW]
    g01   35       20
    g02   30       50
;
TABLE tLDATA(n,nn,*) 'line data'
           X      InvCost
*        [p.u.]  [€/MW]
 n01.n02  0.1      4
;
TABLE tRHODATA(t,g) 'capacity factor [p.u.]'
        g01     g02
t01    1.00     1.00
;
TABLE tDEMDATA(n,t,*) 'demand [MW]'
        Intercept  Slope
*        [MW]    [MW^2/€]
n01.t01   10      0.01
n02.t01   10      0.01
;

* load data to parameters
SB       = 100 ;
C  (g  ) = tGDATA  (g,'LinCost'  ) ;
I  (g  ) = tGDATA  (g,'InvCost'  ) ;
IT (n,nn)= tLDATA  (n,nn,'InvCost') ;
X  (n,nn)= tLDATA  (n,nn,'X') ;
RHO(g,t) = tRHODATA(t,g          ) ;
D0 (t,n) = tDEMDATA(n,t,'Intercept') ;
S0 (t,n) = tDEMDATA(n,t,'Slope') ;

PHI(g) = 1/ SUM[n$[GN(g,n)], S0('t01',n)];
PHI(g) = 0;

* Big-M values
MLB_p    (g,t  )=1000;
MUB_p    (g,t  )=1000;
MLB_beta (g,t  )=1000;
MUB_beta (g,t  )=1000;
MLB_f    (t,n,n)=1000;
MUB_f    (t,n,n)=1000;
MLB_gamma(t,n,n)=1000;
MUB_gamma(t,n,n)=1000;
MLB_d    (t,n  )=1000;
MLB_alpha(t,n  )=1000;

* solve model
OPTIONS NLP=KNITRO;
OPTIONS MCP=PATH;
OPTIONS MIP=GUROBI;


solve mPROTEM_UL_NLP using NLP minimizing ul_of ;
execute_unload 'mPROTEM_UL_NLP.gdx';

solve mPROTEM_UL_MPEC using MPEC minimizing ul_of ;
execute_unload 'mPROTEM_UL_MPEC.gdx';

l.fx(n,nn) $[LN(n,nn)] = l.l(n,nn) ;

solve mPROTEM_LL_NLP using NLP minimizing dummy ;
* save all in gdx format
execute_unload 'mPROTEM_LL_NLP.gdx';

solve mPROTEM_LL_MCP using MCP;
* save all in gdx format
execute_unload 'mPROTEM_LL_MCP.gdx';



$onlisting
