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
Real :: rho_l = 1.0
Real :: P_l = 1.0*101325
Real :: u_l = 0.0
Real :: al 


Real :: rho_r = 0.125
Real :: P_r = 0.1*101325
Real :: u_r = 0.0


Real :: rho_1
Real :: P_1
Real :: u_1

Real :: rho_2
Real :: P_2
Real :: u_2
Real :: a2

Real :: Ms 

Real :: gamma = 1.4
Real :: R = 273.16
Real :: ar

Real :: x1
Real :: x2
Real :: x3
Real :: x4

Real :: A
Real :: B
Real :: C
Real :: Z

Integer,PARAMETER :: nx = 20
Integer :: iit
Integer :: nit=100000

Real, Dimension ( nx ):: x
Real, Dimension ( nx ):: u
Real, Dimension ( nx ):: p
Real, Dimension ( nx ):: rho

al = ((gamma*P_l)/rho_l)**0.5
ar = ((gamma*P_r)/rho_r)**0.5

nloop: Do iit = 1, nit
Ms=0.5+0.001*iit
  A = ((2*gamma)/(gamma+1))*(Ms**2)
  B = (A - ((gamma-1)/(gamma+1)))**((gamma-1)/(2*gamma))
  C = (1 - ((P_r/P_l)*B))
  Z = abs((1/Ms)+(al*((gamma+1)/(gamma-1))*C)-Ms)

	Write (*,*) '#try number', iit, '          Ms = ',Ms 
  
  if(Z < 0.01) exit nloop

end do nloop

P_1 = ((((2*gamma)/(gamma+1))*(Ms**2))-((gamma-1)/(gamma+1)))*P_r
rho_1 = rho_r/(((2*(gamma+1))/(Ms**2))+((gamma-1)/(gamma+1)))
u_1 = (2/(gamma+1))*(Ms-(1/Ms))
u_2 = u_1
P_2 = P_1
rho_2 = ((P_2/P_l)**(1/gamma))*rho_l

a2 = al-((u_2*(gamma-1))/2)

x1 = x0 - al*t
x2 = x0 + ((u_2-a2)*t)
x3 = x0 + u_2*t
x4 = x0 + Ms*t

Write (*,*) 'Dimensional'
Write (*,*) 'P_R',P_r,'rho_R',rho_r,'u_R',u_r
Write (*,*) 'P_1',P_1,'rho_1',rho_1,'u_1',u_1
Write (*,*) 'P_2',P_2,'rho_2',rho_2,'u_2',u_2
Write (*,*) 'P_l',P_l,'rho_l',rho_l,'u_l',u_l

Write (*,*) 'non-Dimensional'
Write (*,*) 'P_R',P_r/(rho_l*al*al),'rho_R',rho_r/rho_l,'u_R',u_r/al
Write (*,*) 'P_1',P_1/(rho_l*al*al),'rho_1',rho_1/rho_l,'u_1',u_1/al
Write (*,*) 'P_2',P_2/(rho_l*al*al),'rho_2',rho_2/rho_l,'u_2',u_2/al
Write (*,*) 'P_l',P_l/(rho_l*al*al),'rho_l',rho_l/rho_l,'u_l',u_l/al

End program analytic_Shock_tube
