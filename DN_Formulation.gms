* Active Distribution Planning

Sets
t each period /t1*t5/
N_exitsting load-points /n1*n11/

N_future   future load-points /n12*n24/
L all lines /l1*l38/
L_ex(L)  existing-Lines /l1*l9/

L_can(L)   candidate-Lines /l10*l38/  ;

Parameters
IC(Psi)
Network_data
Demand(n,t)
 ;
table alter_param(Psi,attributes) ;


Binary Variables
x(Psi)
y(L);

Variables
InvCost
OpCost
LossCost
Revenue
  ;
BINARY Variables
XI
XR
 ;

equations
TotCost Total Cost
KVL1
;

obj_func =e= Inv_Cost + Opr_Cost + Loss_Cost + LOLE_Cost

Inv_Cost =e= sum((n,t), XS(sub)*IC(sub) ) sum((L,t),XR(L_ex,t)*Network_data(L_ex,"RC")
+ XI(L_can,t)*Network_data(L_can,"IC"))

Opr_Cost =e= Ssmax(sub)*OC(sub) + r(n,t)*LOLE_Cost

Loss_Cost =e= X(line)*Lf*Ploss*Cost(loss)
*KCL:    r(n,t) is load not served.
sum(L,PF(L,t))$(tbus=n) + Pprod(n,t) =e= sum(L,PF(L,t))$(fbus=n) + Demand(n,t) + r(n,t)
sum(L,QF(L,t))$(tbus=n) + Qprod(n,t) =e= sum(L,QF(L,t))$(fbus=n) + Demand(n,t) + r(n,t)

*KVL:   constraints are designed in a way to be relaxed when line is not in use.
0 <= PF(L_ex,t) <= PFLmax *(1 - XR(L_ex,t)) + PFRmax(XR(L_ex,t))
0 <= PF(L_can,t) <= PFLmax *XI(L_can,t)
XR(L_ex,t) =g= XR(L_ex,t-1)
XI(L_can,t) =g= XI(L_can,t-1)

0 <= QL <= y(L,t)*QLmax
KVL1.. sum(n,A(n,L_ex)*V(n,t)) =e= PF(L_ex,t)*Network_data(L_ex,"R") + QF(L_ex,t)*Network_data(L_ex,"X")
KVL21.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Network_data(L_can,"R") + QF(L_can,t)*Network_data(L_can,"X") =g= (X(L_can,t) - 1)*M
KVL22.. sum(n,A(n,L_can)*V(n,t)) - PF(L_can,t)*Network_data(L_can,"R") + QF(L_can,t)*Network_data(L_can,"X") =l= (1 - X(L_can,t))*M

V(n,t) =l= Node_data("Vmax")
V(n,t) =g= Node_data("Vmin")

*mesh not allowed (N(t) is the bumber of nodes at each period):
sum(L, y(L,t)) =l= N(t)
