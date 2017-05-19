      Program analytic_Shock_tube

C
C    Purpose:
c    To Solve the shock tube problem analyticaly.
c
c    Record of revisions:
c        Date               Programer          Description of change
c        ====               =========          =====================
c      06/30/2014           M. S. Samara            Original code
c
c

      Implicit none

      Real ::x0 = 0.0
      Real ::t = 1.00109*1.0e-03
      Real ::rho_4 = 1.0
      Real ::P_4 = 1.0*101325
      Real ::u_4 = 0.0
      Real ::a4
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
      
      
      Real :: AA
      Real :: BB
      
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
      
      rho_3 = ((P_3/P_4)**(1/gamma))*rho_4
      
      AA  = ((gamma-1)/(gamma+1))*P21
      BB  = ((gamma-1)/(gamma+1))*(1/P21)
      T21 = (1+AA)/(1+BB)
      
      T1  = P_1/(R*rho_1)
      T3  = P_3/(R*rho_3)
      T4  = P_4/(R*rho_4)
      
      T2  = T21*T1
      
      rho_2 = P_2/(R*T2)
      a3 = ((gamma*P_3)/rho_3)**0.5

      a2 = ((gamma*P_2)/rho_2)**0.5

      u_3 = ((2*a4)/(gamma-1))*(1-(P_3/P_4)**((gamma-1)/(2*gamma)))
      u_2 = u_3

      x1 = x0 - ((((((gamma+1)/(2*gamma))*(P21-1))+1)**0.5)*a1*t)
      x3 = x0 + ((a4+(((gamma+1)/2)*u_3))*t)
      x2 = x0 - u_2*t
      x4 = x0 + a4*t
      
      Write (*,*) 'Dimensional'
      Write (*,*) 'P_1',P_1,'rho_1',rho_1,'u_1',u_1
      Write (*,*) 'P_2',P_2,'rho_2',rho_2,'u_2',u_2
      Write (*,*) 'P_3',P_3,'rho_3',rho_3,'u_3',u_3
      Write (*,*) 'P_4',P_4,'rho_4',rho_4,'u_4',u_4

      Write (*,*) 'non-Dimensional'
      Write (*,*) 'P1',P_1/P_4,'rho1',rho_1/rho_4,'u1',u_1/a4
      Write (*,*) 'P2',P_2/P_4,'rho2',rho_2/rho_4,'u2',u_2/a4
      Write (*,*) 'P3',P_3/P_4,'rho3',rho_3/rho_4,'u3',u_3/a4
      Write (*,*) 'P4',P_4/P_4,'rho4',rho_4/rho_4,'u4',u_4/a4

      Write (*,*) 'X-Values'
      Write (*,*) 'x1',x1
      Write (*,*) 'x2',x2
      Write (*,*) 'x3',x3
      Write (*,*) 'x4',x4

      End program analytic_Shock_tube
