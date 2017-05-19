Program analytic_Shock_tube

!
!    Purpose:
!    To Solve the shock tube problem analyticaly.
!
!    Record of revisions:
!        Date               Programer          Description of change
!        ====               =========          =====================
!      06/30/2014           M. S. Samara            Original code
!
!

Implicit none

Real :: x0 = 0.0
Real :: t = 2.0

Real :: rho_4 = 1.0
Real :: P_4 = 1.0*101325
Real :: u_4 = 0.0
Real :: a4
Real :: T4

Real :: rho_1 = 0.125
Real :: P_1 = 0.1*101325
Real :: u_1 = 0.0
Real :: a1
Real :: T1

Real :: rho_3
Real :: P_3
Real :: u_3
Real :: T3
Real :: a3

Real :: rho_2
Real :: P_2
Real :: u_2
Real :: a2
Real :: T2

Real :: gamma = 1.4
Real :: R = 273.16


Real :: x1
Real :: x2
Real :: x3
Real :: x4

Real :: A
Real :: B
Real :: C
Real :: Z


Real :: A2
Real :: B2

Real :: P21
Real :: T21

Integer,PARAMETER :: nx = 20
Integer :: iit
Integer :: nit=10000

Real, Dimension ( nx ):: x
Real, Dimension ( nx ):: u
Real, Dimension ( nx ):: p
Real, Dimension ( nx ):: rho

a1 = ((gamma*P_1)/rho_1)**0.5
a4 = ((gamma*P_4)/rho_4)**0.5

nloop: Do iit = 1, nit
  P21=0.01*iit
  A = (gamma-1)*(a1/a4)*(P21-1)
  B = ((gamma*2)+((gamma+1)*(P21-1)))**0.5
  C = (1-(A/(((gamma*2)**0.5)*B)))**((-2*gamma)/(gamma-1))
  Z = abs((P21*C)-(P_4/P_1))
  Write (*,*) '#try number', iit, '          P2/P1 = ',P21
if(Z < 0.01) exit nloop

end do nloop

P_2 = P21*P_1
P_3 = P_2

u_3 = ((2*a2)/(gamma-1))*(1-(P_3/P_4)**((gamma-1)/(2*gamma)))
u_2 = u_3

rho_3 = ((P_3/P_4)**(1/gamma))*rho_4

A2  = ((gamma-1)/(gamma+1))*P21
B2  = ((gamma-1)/(gamma+1))*(1/P21)
T21 = (1+A2)/(1+B2)

T1  = P_1/(R*rho_1)
T3  = P_3/(R*rho_3)
T4  = P_4/(R*rho_4)

T2  = T21*T1

rho_2 = P_2/(R*T2)
a3 = ((gamma*P_3)/rho_3)**0.5

x1 = x0 - a1*t
x3 = x0 + ((u_3-a3)*t)
x2 = x0 + u_2*t
x4 = x0 + a4*t

Write (*,*) 'Dimensional'
Write (*,*) 'P_1',P_1,'rho_1',rho_1,'u_1',u_1
Write (*,*) 'P_2',P_2,'rho_2',rho_2,'u_2',u_2
Write (*,*) 'P_3',P_3,'rho_3',rho_3,'u_3',u_3
Write (*,*) 'P_4',P_4,'rho_4',rho_4,'u_4',u_4

!Write (*,*) 'non-Dimensional'
!Write (*,*) 'P_R',P_r/(rho_l*al*al),'rho_R',rho_r/rho_l,'u_R',u_r/al
!Write (*,*) 'P_1',P_1/(rho_l*al*al),'rho_1',rho_1/rho_l,'u_1',u_1/al
!Write (*,*) 'P_2',P_2/(rho_l*al*al),'rho_2',rho_2/rho_l,'u_2',u_2/al
!Write (*,*) 'P_l',P_l/(rho_l*al*al),'rho_l',rho_l/rho_l,'u_l',u_l/al

End program analytic_Shock_tube
