* Active Distribution Planning
* These should be multiplied by 8760 to account for a whole year (if they are hourly)
* LossCost has problems. it's variables and parameters are not determined yet.
* I can simply neglect Q. just a Q<=0.8P would suffice (to lessen the run-time)
* everything should be per unit, yet nothing is.
* Voltage of reference bus
* Radiality3 is wrong.
* Revenue is price of energy times (Loads - (loads not served))
* DN should buy its power (Pprod and Qprod) from a wholesale. **it is not modeled yet**
Sets
y each period /y1*y10/
N  all nodes /n1*n24/
N_ex(N) load-points /n5*n12/
N_can(N)   future load-points /n3,n4,n13*n24/
S(N) all substation nodes /n1*n4/
S_ex(S)  existing substation nodes /n1, n2 /
S_can(S) future substation nodes /n3, n4 /
np(N) nodes of only load /n1,n2,n5*n24/
L all lines /l1*l38/
L_ex(L)  existing-Lines /l1*l9/
L_can(L)   candidate-Lines /l10*l38/
nodePars node parameters /Psmax1,Qsmax1,NIC1,Psmax2,Qsmax2,NIC2,Vmax,Vmin,Psmax,Qsmax,PsmaxR,QsmaxR,NRC,state/
LinePars Line parameters /fbus,tbus,R,X,PFmax,QFmax,Rr,Xr,PFRmax,QFRmax,LRC,R1,X1,PFI1max,QFI1max,LIC1,R2,X2,PFI2max,QFI2max,LIC2,R3,X3,PFI3max,QFI3max,LIC3,InitState,Length/

;
Parameters
number_of_nodes(y) /
y1       11
y2       12
y3       13
y4       15
y5       18
y6       20
y7       22
y8       23
y9       23
y10      23 /
numn(n) /
n1  1
n2  2
n3  3
n4  4
n5  5
n6  6
n7  7
n8  8
n9  9
n10 10
n11 11
n12 12
n13 13
n14 14
n15 15
n16 16
n17 17
n18 18
n19 19
n20 20
n21 21
n22 22
n23 23
n24 24
/
YY(y) numerical value of each year
        / y1      1
          y2      2
          y3      3
          y4      4
          y5      5
          y6      6
          y7      7
          y8      8
          y9      9
          y10     10 /
node_data(n,nodePars)
Line_data(L,LinePars)
LoadP_data(n,y)
LoadQ_data(n,y)
A(L,N)
 ;
$ call gdxxrw node_data.xlsx par node_data rng=sheet1!A1:N25 rdim=1 cdim=1
$GDXIN node_data.gdx
$load node_data
$GDXIN
$ call gdxxrw Line_data.xlsx par Line_data rng=sheet1!A1:AC39 rdim=1 cdim=1
$GDXIN Line_data.gdx
$load Line_data
$GDXIN
$ call gdxxrw LoadP.xlsx par LoadP_data rng=sheet1!A1:K25 rdim=1 cdim=1
$GDXIN LoadP.gdx
$load LoadP_data
$GDXIN
$ call gdxxrw LoadQ.xlsx par LoadQ_data rng=sheet1!A1:K25 rdim=1 cdim=1
$GDXIN LoadQ.gdx
$load LoadQ_data
$GDXIN
$ call gdxxrw Inc_matrix.xlsx par A rng=sheet1!A1:Y39 rdim=1 cdim=1
$GDXIN Inc_matrix.gdx
$load A
$GDXIN

Scalars M big number /1000/
Vref KV /20/
Lf Loss factor /0.35/
LossCost_C dolar per Mwh /60/
LOLEcost value of loss load dolar per Mwh /30000/
SubOC Subsationg operation cost dolar per MVA /2000/
r_inf inflation rate /0.02/
roi rate of interest /.1/
LoadFactor /0.8/
;


Variables
DNObjFunc
InvCost(y)
OprCost(y)
LossCost(y)
PF(L,y)
QF(L,y)
Pprod(n,y)
Qprod(n,y)
rp(n,y)
rq(n,y)
V(n,y) ;

positive variables
Pprod(n,y)
Qprod(n,y)
rp(n,y)
rq(n,y)
V(n,y)
Revenue(y)
uu(l,y)
vv(l,y)
PFlin(l,y)
 ;
BINARY Variables
XI1(L_can,y)
XI2(L_can,y)
XI3(L_can,y)
XR(L_ex,y)
Z1(S_can,y)
Z2(S_can,y)
ZR(np,y)
b(l,y)
 ;

equations
DNTotCost Total Cost
eq_InvCost
eq_OprCost
eq_LossCost
*eq_Revenue
eq_lin1
eq_lin2
eq_lin3
eq_lin4

KCL

KVL1
KVL21
KVL22
KVL31
KVL32
KVL41
KVL42

eq_powerflow1
eq_powerflow2
eq_Qflow

eq_sub_ins
*eq_sub21
*eq_sub22
eq_OnlySub
eq_exSub1
eq_exSub2

eq_Line_rep
eq_only_ins
eq_Line_ins1
eq_Line_ins2
eq_Line_ins3
*eq_OnlyRep1
*eq_OnlyRep2
eq_ZOnce1
eq_ZOnce2
eq_OnlySubIns
eq_ZR
eq_ZROnce

eq_referenceV
*eq_Vmax
*eq_Vmin

Radiality1
Radiality2


;


DNTotCost.. DNObjFunc =e= sum(y, ( InvCost(y) + OprCost(y) + LossCost(y))/((1+roi)**YY(y)))  ;

eq_InvCost(y)..  InvCost(y) =e= sum(S_can, (Z1(S_can,y)-Z1(S_can,y-1))*node_data(S_can,'NIC1') + (Z2(S_can,y)-Z2(S_can,y-1))*node_data(S_can,'NIC2'))
+sum(np,(ZR(np,y)-ZR(np,y-1))*node_data(np,'NRC'))+ sum(L_ex,(XR(L_ex,y)-XR(L_ex,y-1))*Line_data(L_ex,"LRC"))
+ sum(L_can,(XI1(L_can,y)-XI1(L_can,y-1))*Line_data(L_can,"LIC1"))
+ sum(L_can,(XI2(L_can,y)-XI2(L_can,y-1))*Line_data(L_can,"LIC2"))
+ sum(L_can,(XI3(L_can,y)-XI3(L_can,y-1))*Line_data(L_can,"LIC3")) ;

eq_OprCost(y).. OprCost(y) =e= sum(S_can,Z1(S_can,y)*node_data(S_can,'Psmax1')*SubOC + Z2(S_can,y)*node_data(S_can,'PSmax2')*SubOC)
+ sum(S_ex,node_data(S_ex,'Psmax')*SubOC) +  sum(n, rp(n,y)*LOLEcost) ;

eq_LossCost(y).. LossCost(y) =e= sum(L_ex, Lf*PFlin(L_ex,y)*Line_data(L_ex,'R')*LossCost_C)
+ sum(L_can, Lf*PFlin(L_can,y)*Line_data(L_can,'R')*LossCost_C) ;
*eq_Salvage.. SalvageCost =e= sum(S_ex,XR(S_ex,'y7')) + sum(S_can,XI(S_can,'y7'))
*eq_Revenue(y).. Revenue(y) =e= price_en(y)*sum(n,LoadP_data(n,y) - rp(n,y)) ;
*KCL:    r(n,t) is load not served.

eq_lin1(l,y).. PFlin(l,y) =e= uu(l,y) + vv(l,y);
eq_lin2(l,y).. PF(l,y) =e= uu(l,y) - vv(l,y) ;
eq_lin3(l,y).. uu(l,y) =l= M*b(l,y) ;
eq_lin4(l,y).. vv(l,y) =l= M*(1-b(l,y)) ;

KCL(n,y).. Pprod(n,y) + rp(n,y) +  sum(L$(Line_data(L,'tbus')=numn(n)),PF(L,y)) =e=  sum(L$(Line_data(L,'fbus')=numn(n)),PF(L,y)) + LoadP_data(n,y) ;

*KVL:   constraints are designed in a way to be relaxed when line is not in use.
* A /v1 is missing in KVL
KVL1(L_ex,y).. sum(n,A(L_ex,n)*V(n,y)) =e= -(PF(L_ex,y)*Line_data(L_ex,"R") + QF(L_ex,y)*Line_data(L_ex,"X"))/Vref ;
KVL21(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R1") + QF(L_can,y)*Line_data(L_can,"X1"))/Vref =g= (XI1(L_can,y) - 1)*M  ;
KVL22(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R1") + QF(L_can,y)*Line_data(L_can,"X1"))/Vref =l= (1 - XI1(L_can,y))*M  ;
KVL31(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R2") + QF(L_can,y)*Line_data(L_can,"X2"))/Vref =g= (XI2(L_can,y) - 1)*M  ;
KVL32(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R2") + QF(L_can,y)*Line_data(L_can,"X2"))/Vref =l= (1 - XI2(L_can,y))*M  ;
KVL41(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R3") + QF(L_can,y)*Line_data(L_can,"X3"))/Vref =g= (XI3(L_can,y) - 1)*M  ;
KVL42(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R3") + QF(L_can,y)*Line_data(L_can,"X3"))/Vref =l= (1 - XI3(L_can,y))*M  ;

*power flow limitations
eq_powerflow1(L_ex,y).. PFlin(L_ex,y) =l= Line_data(L_ex,'PFmax') *(1 - XR(L_ex,y)) + Line_data(L_ex,'PFRmax')*XR(L_ex,y)  ;
eq_powerflow2(L_can,y).. PFlin(L_can,y) =l= Line_data(L_can,'PFI1max')*XI1(L_can,y) + Line_data(L_can,'PFI2max')*XI2(L_can,y)
+ Line_data(L_can,'PFI3max')*XI3(L_can,y) ;
eq_Qflow(L,y).. QF(L,y) =e= 0.75*PFlin(L,y) ;

*to model Substation placing. Ppmax is zero for non-substations.

eq_sub_ins(S_can,y).. Pprod(S_can,y) =l=  Z1(S_can,y)*node_data(S_can,'Psmax1') + Z2(S_can,y)*node_data(S_can,'Psmax2') ;
eq_OnlySub(S_can,y).. Z1(S_can,y) + Z2(S_can,y) =l= 1 ;
*same model for existing substations replacement.
eq_exSub1(np,y).. Pprod(np,y) =l= (1 - ZR(np,y))*node_data(np,'Psmax') + ZR(np,y)*node_data(np,'PsmaxR')  ;
*eq_exSub2(S_ex,y).. Qprod(S_ex,y) =l= (1 - ZR(S_ex,y))*node_data(S_ex,'Qsmax') + ZR(S_ex,y)*node_data(S_ex,'QsmaxR')  ;
eq_exSub2(n,y).. Qprod(n,y) =e= 0.75*Pprod(n,y)  ;

*to ensure only one replacement or installation
eq_Line_rep(L_ex,y).. XR(L_ex,y) =g= XR(L_ex,y-1) ;
eq_Line_ins1(L_can,y).. XI1(L_can,y) =g= XI1(L_can,y-1) ;
eq_Line_ins2(L_can,y).. XI2(L_can,y) =g= XI2(L_can,y-1) ;
eq_Line_ins3(L_can,y).. XI3(L_can,y) =g= XI3(L_can,y-1) ;
eq_only_ins(L_can,y).. XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y) =l= 1 ;

*eq_OnlyRep2(L_can).. sum(t,XI(L_can,y)) =l= 1 ;

eq_ZOnce1(S_can,y).. Z1(S_can,y) =g= Z1(S_can,y-1) ;
eq_ZOnce2(S_can,y).. Z2(S_can,y) =g= Z2(S_can,y-1) ;
eq_OnlySubIns(S_can,y).. Z1(S_can,y) + Z2(S_can,y) =l= 1 ;

eq_ZR(np,y).. ZR(np,y) =l= node_data(np,'state') ;
eq_ZROnce(np,y).. ZR(np,y) =g= ZR(np,y-1) ;


*Node Voltage limitations
eq_referenceV(S,y).. V(S,y) =e= Vref ;
*eq_Vmax(np,y).. V(np,y) =l= Node_data(np,"Vmax") ;
*eq_Vmin(np,y).. V(np,y) =g= Node_data(np,"Vmin") ;


Radiality1(y).. card(L_ex) + sum(L_can, XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y)) =e= number_of_nodes(y) -card(S_ex) ;
Radiality2(n,y).. sum(L_can$(Line_data(L_can,"tbus")=ord(n)), XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y))
+ sum(L_ex$(Line_data(L_ex,"tbus")=ord(n)), Line_data(L_ex,'InitState')) =l= 1 ;

display LoadP_data
model DNPlanning /all/ ;
option optcr = 0.001 ;
solve DNPlanning minimizing DNObjFunc using mip;


