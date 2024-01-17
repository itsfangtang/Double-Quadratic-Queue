# Double-Quadratic-Queue

Appendix A: Mathematical DQQ model
This section introduces the proposed Double Quadratic Queue (DQQ) model. By adopting the fluid queue derivation process, the model considers the two distinct stages and differentiates itself from a typical single-stage fluid queue model. The fluid queue modeling approach is a method used to analyze and model the behavior of queues in a dynamic setting, taking into account varying arrival and discharge rates. Newell's approach combines Taylor expansion, calculus, and geometric representation to derive an algebraic expression for the queue length at any given time. The DQQ model incorporates both geometric representation and algebraic expression. This fluid queue process helps visualize the relationships between arrival and discharge rates, queue length, and the passage of time. The ultimate goal of the subsequent derivation process is to obtain an algebraic expression for the queue length or cumulative change, which can be employed to study the queue's or system's behavior analytically.
Table 1 summarizes the notations used in the DQQ model. 
Table 1 Notations and definitions used in this paper
Notations	Definitions
π(t) 	net flow rate at time t
Q(t) 	ridership rate at time t, dQ(t)/dt=π(t)
t_0 	start time of disruption, t_0=0
t_1 	time instance with maximum negative net flow rate
t_2 	time instance with maximum cumulative ridership change
t_3 	time instance with maximum positive net flow rate
t_4 	end time of recovery process, where reaching a new equilibrium condition
D 	disruption duration, D=t_2-t_0
R 	recovery duration, R=t_4-t_2
P 	entire duration of two normal/equilibriums, P=t_4-t_0
m 	Relative disruption duration ratio, m=(t_2-t_0)/(t_4-t_0 )=D/P=D/(D+R)
r 	recovery to disruption (ROD) ratio, r=R/D,  which is the proportion of recovery in relation to the disruption, indicating how well the system recovers from the disruption
Q_max 	maximum ridership changes at time t_2
∆_e 	permanent loss at end time t_4
ρ 	parameter in the single quadratic-form net flow rate function
γ 	parameter in the cubic-form net flow rate function
α,β 	parameters in the double quadratic-form net flow rate function

To maintain generality, we normalize both time and ridership together such that their cumulative changes from the baseline condition at the beginning of the disruption are set to 0, i.e., t_0=0. By doing so, the entire duration of the disruption and recovery can be consistently represented as:
P=t_4 	(1)
D=t_2=mP 	(2)
r=R/D=(1-m)/m 	(3)

Typically, parameters t_1, t_2, t_3, t_4, Q_max, ∆_e are observable directly from time series datasets. In the DQQ model, we assume the net flow rate, π(t), can be represented as double quadratic forms, as shown in Figure 1, which has three intersection points, t_0, t_2, and t_4, with the horizontal axis. 
 
Figure 1 Graphical illustration of DQQ model for a single wave in transit cumulative ridership changes

The first-stage net flow rate, π_D (t), could be approximated by Taylor’s Theorem with k=2 at t=t_1:
π_D (t)=π_D (t_1 )+π_D'(t_1 )∙(t-t_1 )+(〖π_D〗^'' (t_1 ))/2 〖(t-t_1)〗^2 	(4) 

Considering the parabola opens downwards, we have:
α=(〖π_D〗^'' (t_1 ))/2 	(5) 

At the first stage of DQQ model, π_D (t_0 )=π_D (t_2 )=0.
π_D (t_0 )=π_D (t_1 )-α(t_0-t_1 )^2 	(6)
π_D (t_2 )=π_D (t_1 )-α(t_2-t_1 )^2 	(7)

We have 
π_D (t_1 )=α(t_0-t_1 )^2 	(8)
t_0+t_2=2t_1 	(9)

By substituting Eq. (8) and Eq. (9) to Eq. (4), the factored form of π_D (t) can be represented by:
π_D (t)=α(t-t_0 )(t-t_2 ) 	(10) 

Similarly, considering Taylor’s Theorem with k=2 at t=t_3, we let β=(〖π_R〗^'' (t_1 ))/2, then have π_R (t)=β(t-t_2 )(t-t_4 ).
Thus, the DQQ model can be represented as Eq. (4). This paper aims to calibrate the unknown second-order gradient parameters, α and β.
π(t)={█(π_D (t)=α(t-t_0 )(t-t_2 ),   t_0<t<t_2@π_R (t)=β(t-t_2 )(t-t_4 ),t_2<t<t_4 )┤ 	(4) 

As depicted in Figure 1(a), the disruption occurs, and ridership starts to decline at time t_0. At time t_1, the flow loss rate reaches its peak. Consequently, the transit system experiences the greatest passenger loss, Q_max, at time t_2, as illustrated in Figure 2(b). Following t_2, transit ridership begins to recover. At t_3, the net flow recovery rate attains its maximum value.  Eventually, the transit system reaches a new normal/equilibrium at t_4, where the ∆_e represents the gap between this new equilibrium and the initial equilibrium status. The duration of disruption and recovery processes determines the horizontal temporal interval and the vertical steepness of the adaptivity evolution curve, i.e. the relative location of roots in the DQQ model. Generally, the durations and the permanent loss can be used to assess the system’s recovery capability. 
As shown in Figure 1, the boundary conditions can be represented as Eq. (7)-(14). 
(dπ_D (t_1 ))/dt=0 	(5)
(dπ_R (t_3 ))/dt=0 	(6)
Q(t_0 )=0 	(7)
Q(t_2 )=Q_max 	(8)
Q(t_4 )=∆_e 	(9)

(1) When t_0<t<t_2, the disruption process is represented by π_D (t). 
The first order derivative of net flow rate in disruption process, π_D (t), is obtained:
(dπ_D (t))/dt=α(2t-t_2 ) 	(10)

Substituting Eq. (4) and α≠0 leads to the following symmetrical relationships between t_0 and t_2:
t_2=2t_1 	(11)

Virtual queue length or cumulative flow change at time t, Q(t) can be obtained:
Q(t)=∫_(t_0)^t▒〖π_D (t) 〗 dτ 	(12)

By substituting π_D (t) function in Eq. (4), Q(t)  can be represented in terms of t_2 and α:
Q(t)=α/6 t^2 (2t-3t_2 ) 	(13)
As stated above, the maximum queue length is achieved at time t2:
Q(t_2 )=α/6 [2〖t_2〗^3-3〖t_2〗^3 ]=-α/6 〖t_2〗^3 	(14)

The Eq. (14) dedicates that Q_max can be represented by the time instance t_2 and the curvature parameter α. Since Q_max and t_2 can be observed from real-life data, the curvature rate α can be calibrated.

(2) When t_2<t<t_4,  the recovery process is represented by π_R (t).
The first order derivative of π_R (t) is derived as:
(dπ_R (t))/dt=β(2t-t_2-t_4 ) 	(15)

By substituting Eq. (6) and β≠0, t_2 and t_4 are symmetric about t_3:
t_2+t_4=2t_3 	(16)

The queue length Q(t) during the recovery process is expressed as:
Q(t)=Q(t_2 )+∫_(t_2)^t▒〖π_R (t) 〗 dτ 	(17)

Substituting Eq. (14) and π_R (t) function in Eq. (4) leads to Q(t)  expression:
Q(t)=-α/6 〖t_2〗^3-β〖〖(t-t〗_2)〗^2 [(t_4-t_2)/2-(t-t_2)/3] 	(18)

The final permanent loss is determined at time t_4, i.e. ∆_e=Q(t_4):
Q(t_4 )=-α/6 〖t_2〗^3-β/6 〖〖(t_4-t〗_2)〗^3 	(19)

Since ∆_e can be observed from real-life data, and the α can be obtained using Eq. (14), the β can be calibrated.
Especially, when the ridership is fully recovered at the end of the adaptative evolution, i.e. ∆_e=0. Substituting Eq. (19) leads to the internal relationship between the α and β in terms of parameter m or r.
β=α〖(m/(m-1))〗^3=-αr^(-3) 	(20)

3.2 Comparison among single-stage quadratic model, cubic model, and DQQ model
It is intriguing to compare the proposed DQQ model with the single-stage quadratic queue model and the cubic queue model, as all three models fall within the overarching modeling framework of fluid queues. As demonstrated in Figure 2, these three models share similar notations for various critical time instances during the recovery-to-disruption process. Table 2 presents a comparison of the boundary conditions and derivation results across these three distinct formulations. The single-stage quadratic queue model is incapable of describing the new equilibrium point with the net flow rates at t_4, and the single-stage cubic queue model can only depict the fast recovery process with m>1/2, accompanied by possible permanent loss. In contrast, the DQQ model can characterize a severe, prolonged disruption with a feasible range of m between 0 and 1 by permitting a first-order gradient jump at t_2. Meanwhile, the single-stage models are continuous and smooth throughout their entire domain. 
                    
(1)                                                                                  (2)
Figure 2 Graphical illustration of net flow rate functions in (1) singe quadratic form and (2) cubic form
                 
(1)                                                                                  (2)
Figure 2 Graphical illustration of net flow rate functions in (1) singe quadratic form and (2) cubic form with different values of m

Table 2 Comparison in different forms of net flow rate
Net flow rate functions	Boundary condition	Derivation results	Remark
Single quadratic forms
π(t)=ρ(t-t_0 )(t-t_2 )  
	(dπ_D (t_1 ))/dt=0, 
Q(t_0 )=0, 
Q(t_2 )=Q_max, 
Q(t_4 )=∆_e 
	Q_max=-ρ/6 〖t_2〗^3,
∆_e=ρ/6 [2〖t_4〗^3-3t_2 〖t_4〗^2 ] 	m∈[2/3,1] 
(1) m=2/3, full recovery, i.e. ∆_e=0;
(2) m=1, disruption process only.
Cubic forms
π(t)=γ(t-t_0 )(t-t_2 )(t-t_4 ) 	(dπ_D (t_1 ))/dt=0,
(dπ_R (t_3 ))/dt=0, 
Q(t_0 )=0, 
Q(t_2 )=Q_max, 
Q(t_4 )=∆_e 
	Q_max=P^4/12 γm^3 (2-m), 
∆_e=P^4/12 γ(2m-1),
∆_e=Q_max  (2m-1)/((2-m)m^3 ) 	m∈[1/2,1] 
(1) m=1/2, full recovery, i.e. ∆_e=0;
(2) m=1, disruption process only.
Double quadratic forms
π(t)={█(π_D (t),   t_0<t<t_2@π_R (t),   t_2<t<t_4 )┤
		Q_max=-α/6 〖t_2〗^3,
∆_e=-α/6 〖t_2〗^3-β/6 〖〖(t_4-t〗_2)〗^3 	m∈[0,1] 
(1) m=0, similar to Newell’s PAQ model;
(2) m=1/2, if π(t_1 )=π(t_3 ), full recovery, i.e. ∆_e=0;
(3) m=1, disruption process only.
![image](https://github.com/FangTang999/Double-Quadratic-Queue/assets/38580581/18dba412-6956-496a-9ef7-6696d178db00)
