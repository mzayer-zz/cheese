* Active Distribution Planning
* These should be multiplied by 8760 to account for a whole year (if they are hourly)

Sets
t each period /t1*t5/
N  all nodes /n1*n24/
N_ext(N) load-points /n1*n11/
N_can(N)   future load-points /n12*n24/

L all lines /l1*l38/
L_ex(L)  existing-Lines /l1*l9/
L_can(L)   candidate-Lines /l10*l38/  ;

Parameters
Loss
Network_data
LossCost
LOLEcost
SubOC
Demand(n,t)
 ;
table alter_param(Psi,attributes) ;


Binary Variables
x(Psi)
y(L);

Variables
Inv_Cost
Opr_Cost
LossCost
Revenue
  ;
BINARY Variables
XI
XR
Z1
Z2
 ;

equations
TotCost Total Cost
Inv_Cos
KVL1
;

TotCost.. obj_func =e= Inv_Cost + Opr_Cost + Loss_Cost

Inv_Cost..  InvCost =e= sum((n,t), Z1(n)*node_data(n,'sub1') + Z2(n)*node_data(n,'sub2') ) + sum((L,t),XR(L_ex,t)*Network_data(L_ex,"RC")
+ XI(L_can,t)*Network_data(L_can,"IC"))

Opr_Cost =e= sum((n,t),Z1*node_data(n,'sub1')*SubOC + Z2*node_data(n,'sub2')*SubOC + r(n,t)*LOLEcost)

Loss_Cost =e= sum(L_ext, X*Lf*Ploss*LossCost
*KCL:    r(n,t) is load not served.
sum(L,PF(L,t))$(tbus=n) + Pprod(n,t) =e= sum(L,PF(L,t))$(fbus=n) + Demand(n,t) - r(n,t)
sum(L,QF(L,t))$(tbus=n) + Qprod(n,t) =e= sum(L,QF(L,t))$(fbus=n) + Demand(n,t) - r(n,t)

*to model Substation placing. Ppmax is zero for non-substations.
Sprod(n,t) =l= Z1*Ssmax
Sprod(n,t) =l= Z2*Ssmax
Z1(n,t) + Z2(n,t) =l= 1
*KVL:   constraints are designed in a way to be relaxed when line is not in use.
0 <= PF(L_ex,t) <= PFLmax *(1 - XR(L_ex,t)) + PFRmax(XR(L_ex,t))
0 <= PF(L_can,t) <= PFLmax *XI(L_can,t)
XR(L_ex,t) =g= XR(L_ex,t-1)
XI(L_can,t) =g= XI(L_can,t-1)
sum(t,Z1(n,t)) =l= 1
sum(t,Z2(n,t)) =l= 1
sum(t,XR(L_ex,t)) =l= 1
sum(t,XI(L_can,t)) =l= 1

KVL1.. sum(n,A(n,L_ex)*V(n,t)) =e= PF(L_ex,t)*Network_data(L_ex,"R") + QF(L_ex,t)*Network_data(L_ex,"X")
KVL21.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Network_data(L_can,"R") - QF(L_can,t)*Network_data(L_can,"X") =g= (X(L_can,t) - 1)*M
KVL22.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Network_data(L_can,"R") - QF(L_can,t)*Network_data(L_can,"X") =l= (1 - X(L_can,t))*M

V(n,t) =l= Node_data("Vmax")
V(n,t) =g= Node_data("Vmin")

*mesh not allowed (N(t) is the bumber of nodes at each period):
sum(L, y(L,t)) =l= N(t)
