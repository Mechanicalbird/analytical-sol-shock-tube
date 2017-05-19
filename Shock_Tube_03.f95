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

Real :: A1
Real :: B1
Real :: C1
Real :: Z1


Real :: A2
Real :: B2

Real :: P21
Real :: T21

Integer,PARAMETER :: nx = 20
Integer :: iit
Integer :: nit=10000000

Real, Dimension ( nx ):: x
Real, Dimension ( nx ):: u
Real, Dimension ( nx ):: p
Real, Dimension ( nx ):: rho

a1 = ((gamma*P_1)/rho_1)**0.5
a4 = ((gamma*P_4)/rho_4)**0.5

nloop: Do iit = 1, nit
  P21=0.0001*iit
  A1 = (gamma-1)*(a1/a4)*(P21-1)
  B1 = ((gamma*2)+((gamma+1)*(P21-1)))**0.5
  C1 = (1-(A1/(((gamma*2)**0.5)*B1)))**((-2*gamma)/(gamma-1))
  Z1 = abs((P21*C1)-(P_4/P_1))
  Write (*,*) '#try number', iit, '          P2/P1 = ',P21
if(Z1 < 0.01) exit nloop

end do nloop

End program analytic_Shock_tube
