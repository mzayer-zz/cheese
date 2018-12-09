* Active Distribution Planning
* These should be multiplied by 8760 to account for a whole year (if they are hourly)
* Distribution Network's Income from selling energy to consumers is not modeled yet.
* LossCost has problems. it's variables and parameters are not determined yet.
Sets
t each period /t1*t5/
N  all nodes /n1*n24/
N_ext(N) load-points /n1*n11/
N_can(N)   future load-points /n12*n24/

L all lines /l1*l38/
L_ex(L)  existing-Lines /l1*l9/
L_can(L)   candidate-Lines /l10*l38/  ;

Alias(nn,n) ;

Parameters
Loss
Network_data
LossCost
LOLEcost
SubOC
Demand(n,t)
 ;
table alter_param(Psi,attributes) ;



Variables
obj_func
Inv_Cost
Opr_Cost
LossCost

Revenue
  ;
BINARY Variables
XI(L_can)
XR(L_ex)
Z1(n)
Z2(n)
 ;

equations
TotCost Total Cost
eq_InvCost
eq_OprCost
eq_LossCost
KCL1
KCL2
KVL1
KVL21
KVL22
eq_sub11
eq_sub12
eq_sub21
eq_sub22
eq_OnlySub
eq_powerflow1
eq_powerflow2
eq_powerflow3
eq_powerflow4
eq_replacement
eq_installation
eq_OnlySubIns1
eq_OnlySubIns2
eq_OnlyRep1
eq_OnlyRep2
eq_Vmax
eq_Vmin
Radiality1
Radiality2

;

Variables
obj_func
InvCost
OprCost
LossCost
B(n,nn)

;
TotCost.. obj_func =e= InvCost + OprCost + LossCost  ;

eq_InvCost..  InvCost =e= sum((n,t), Z1(n)*node_data(n,'sub1') + Z2(n)*node_data(n,'sub2') ) + sum((L,t),XR(L_ex,t)*Line_data(L_ex,"RC")
+ XI(L_can,t)*Line_data(L_can,"IC")) ;

eq_OprCost.. OprCost =e= sum((n,t),Z1*node_data(n,'sub1')*SubOC + Z2*node_data(n,'sub2')*SubOC + r(n,t)*LOLEcost) ;

eq_LossCost.. LossCost =e= sum(L_ext, X*Lf*Ploss*LossCost) ;

*KCL:    r(n,t) is load not served.
KCL1(n).. sum(L$(Line_data(L,'tbus')=ord(n)),PF(L,t)) + Pprod(n,t) =e= sum(L$(Line_data(L,'fbus')=ord(n)),PF(L,t)) + Demand(n,t) - r(n,t) ;
KCL2(n).. sum(L$(Line_data(L,'tbus')=ord(n)),QF(L,t)) + Qprod(n,t) =e= sum(L$(Line_data(L,'fbus')=ord(n)),QF(L,t)) + Demand(n,t) - r(n,t) ;

*KVL:   constraints are designed in a way to be relaxed when line is not in use.
KVL1.. sum(n,A(n,L_ex)*V(n,t)) =e= PF(L_ex,t)*Line_data(L_ex,"R") + QF(L_ex,t)*Line_data(L_ex,"X") ;
KVL21.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Line_data(L_can,"R") - QF(L_can,t)*Line_data(L_can,"X") =g= (X(L_can,t) - 1)*M  ;
KVL22.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Line_data(L_can,"R") - QF(L_can,t)*Line_data(L_can,"X") =l= (1 - X(L_can,t))*M  ;

*power flow limitations
eq_powerflow1(L_ext,t).. PF(L_ex,t) =l= PFLmax *(1 - XR(L_ex,t)) + PFRmax(XR(L_ex,t))  ;
eq_powerflow2(L_ext,t).. PF(L_ex,t) =g= -PFLmax *(1 - XR(L_ex,t)) - PFRmax(XR(L_ex,t)) ;
eq_powerflow3(L_can,t).. PF(L_can,t) =l= PFLmax *XI(L_can,t) ;
eq_powerflow4(L_can,t).. PF(L_can,t) =g= -PFLmax *XI(L_can,t) ;

*to model Substation placing. Ppmax is zero for non-substations.
eq_sub11(n,t).. Pprod(n,t) =l= Z1*Psmax ;
eq_sub12(n,t).. Qprod(n,t) =l= Z1*Qsmax ;
eq_sub21(n,t).. Pprod(n,t) =l= Z2*Psmax ;
eq_sub22(n,t).. Qprod(n,t) =l= Z2*Qsmax ;
eq_OnlySub(n,t).. Z1(n,t) + Z2(n,t) =l= 1 ;



*to ensure only one replacement or installation
eq_replacement(L_ex,t).. XR(L_ex,t) =g= XR(L_ex,t-1) ;
eq_installation(L_can,t).. XI(L_can,t) =g= XI(L_can,t-1) ;
eq_OnlySubIns1(n).. sum(t,Z1(n,t)) =l= 1 ;
eq_OnlySubIns2(n).. sum(t,Z2(n,t)) =l= 1 ;
eq_OnlyRep1(n).. sum(t,XR(L_ex,t)) =l= 1 ;
eq_OnlyRep2(n).. sum(t,XI(L_can,t)) =l= 1 ;

*Node Voltage limitations
eq_Vmax(n,t).. V(n,t) =l= Node_data("Vmax") ;
eq_Vmin(n,t).. V(n,t) =g= Node_data("Vmin") ;

*mesh not allowed (N(t) is the bumber of nodes at each period):
Radiality1..exLineNumber + sum(L, XI(L_can,t)) =e= N(t) - 1 ;
Radiality2(L_can)$(Line_data(L,'tbus')=ord(n)).. B(n,nn) + B(nn,n)