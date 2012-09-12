!
! Definition of the interpolation grid 
      !
  
     


!!!! 
Module golf

      !Double Precision, Parameter :: kmin=0.002 ! valeur de ki mimum
      Integer, Parameter :: Ntot=813   ! Total number of points
      Integer, Parameter :: NetaI=4878 ! 6*4884 length of EtaI
      Integer, Dimension(32) :: Nk2 
 

CONTAINS

!Compute the effective strain ratefor a rank-2 array
     Subroutine epseff(eps,eps_e,l,m)  
       Double Precision, Dimension(:,:,:), Intent(in) :: eps
       Double Precision, Dimension(:), Intent(out) :: eps_e
       Integer, Intent(in) :: l,m
       Integer :: i,j
       
       eps_e(:)=0
       
       Do i=1,3
         Do j=1,3
           eps_e(:)=eps_e(:)+eps(:,i,j)**2
         End Do
       End Do
       eps_e = 0.5*eps_e
     End Subroutine epseff


!Get the nonrelative viscosity for a rank-3 array of deviatoric stress
!First index spatial. 
     Subroutine  nrvisc(W,eps0,nu,m)

              
       Implicit None
       Integer, Intent(in) :: m
       Double Precision, Intent(out), Dimension(m) :: nu
       Double Precision, Intent(in), Dimension(m)  :: eps0
       Double Precision, Intent(in), Dimension(m,6) :: W        ! Glen law parameters
       

       ! W(1)=B0 Glen's fluidity parameter
       ! W(2)=n Glen's law exponent
       ! W(3)=Q1 Energy activation for Tc<Tl
       ! W(4)=Q2 Energy activation for Tc>Tl
       ! W(5)=T0 temperature reference for B0
       ! W(6)=Tl temperature for Q1->Q2
       !R is the gas constant J mol**-1 K**-1
       Double Precision, Dimension(m) :: Q
       Double Precision, Parameter :: Tzero=273.15,R=8.313
       Double Precision :: DT
       Integer :: i
!f2py intent(in) Tc,W
!f2py intent(out) Bg
      
       
       Do i=1,m
          If (W(i,6) .GE. W(i,5)) Then 
             Q(i)=W(i,3)
             Else 
                Q(i)=W(i,4)
          End If
       End Do

       nu=0.5*W(:,1)*exp(-Q/R*(1/W(:,5)-1/W(:,6)))*eps0**((1-W(:,2))/W(:,2))

     End Subroutine nrvisc
 

    Subroutine cmat(evs,angles,m,eta36)
    Implicit None
       
       Integer, Intent(in) :: m
       
       !f2py2 depend(m) evs, angles, eta36
       Double Precision, Dimension(4878) :: etaI
       Double Precision, Intent(in), Dimension(m,3) :: evs, angles
       Double Precision, Intent(out), Dimension(m,6,6) :: eta36
       Integer :: i
!       Double Precision, Intent(in), Dimension(:,:) :: evs, angles
!       Double Precision, Intent(out), Dimension(:,:,:) :: eta36
!       Integer :: i,m
!       Double Precision, Dimension(6,6) :: e36 
       Double Precision, Dimension(3) :: ai,Angle

!Read in from file. There is a bug in f2py that prevents etaI from
!being read in beforehand.
       Open(1,file='040010010.Va')
       Do i=1,813
          READ( 1, '(6(e14.8))' ) etaI(6*(i-1)+1:6*(i-1)+6)
       End Do
 

       !m=Size(evs,1)
       Do i=1,m
         Angle=angles(i,1:3)
         ai=evs(i,1:3)
         call OPILGGE_ai_nl(evs(i,1:3),angles(i,1:3),eta36(i,1:6,1:6),etaI)
       End Do
    End Subroutine cmat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Non-relative viscosities and Temperature dependancy      !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
     Subroutine  BGlenT(Tc,W,Bg)

              
       Implicit None

       Double Precision,Intent(out) :: Bg
       Double Precision, Intent(in) :: Tc                     ! Temperature en d Celsius
       Double Precision, Intent(in), Dimension(7) :: W        ! Glen law parameters
       ! W(1)=B0 Glen's fluidity parameter
       ! W(2)=n Glen's law exponent
       ! W(3)=Q1 Energy activation for Tc<Tl
       ! W(4)=Q2 Energy activation for Tc>Tl
       ! W(5)=T0 temperature reference for B0
       ! W(6)=Tl temperature for Q1->Q2
       Double Precision, parameter :: Tzero=273.15
       Double Precision :: Q, DT
!f2py intent(in) Tc,W
!f2py intent(out) Bg
       
       Q= W(3)
       If (Tc.GT.W(6)) Q=W(4)
       Q=Q/W(7)
       DT=-1./(Tc+Tzero)+1./(W(5)+Tzero)

      Bg=W(1)*Exp(Q*DT)

     End Subroutine BGlenT
    
       Subroutine ViscGene(eta6,Angle,eta36)
       
        
       
       Implicit None
       Double Precision, Intent(in), Dimension(3) :: Angle    ! Euler Angles        
       Double Precision, Intent(in), Dimension(6) :: eta6     ! 6 Viscosities of the
                                                  ! matrice law

       Double Precision, Intent(out), Dimension(6,6) :: eta36  ! Viscosity matrix in 
                                                   ! the reference frame
        
       Double Precision :: p,t,o,st,ct
       Double Precision, Dimension(3,3) :: Q
       Double Precision :: dQk 
       Double Precision, Dimension(6) :: coef
       Integer, Dimension(6) :: ik
       Integer, Dimension(6) :: jk
       Integer :: k,m,n
       Integer :: i,j 
      coef=(/1.,1.,1.,.0,.0,.0/)
      ik=(/1,2,3,1,2,3/)
      jk=(/1,2,3,2,3,1/)


!  Angle = Phi, Theta, Omega
       
       p=Angle(1)
       t=Angle(2)
       o=Angle(3)

! terms of the rotation matrix from RG to RO

        ct = cos(t)
        Q(3,3) = ct
        Q(1,3) = sin(o)
        Q(2,3) = cos(o)
        Q(3,1) = sin(p)
        Q(3,2) = cos(p)

        Q(1,1) = Q(3,2)*Q(2,3) - Q(3,1)*Q(1,3)*ct 
        Q(1,2) = Q(3,1)*Q(2,3) + Q(3,2)*Q(1,3)*ct 
        Q(2,1) = - Q(3,2)*Q(1,3) - Q(3,1)*Q(2,3)*ct 
        Q(2,2) = - Q(3,1)*Q(1,3) + Q(3,2)*Q(2,3)*ct 

        st = sin(t)
        Q(1,3) = Q(1,3)*st
        Q(2,3) = Q(2,3)*st
        Q(3,1) = Q(3,1)*st
        Q(3,2) = -Q(3,2)*st

! 36  terms of the Voigt matrix in the reference frame 
        Do m=1,6
          Do n=1,6
            eta36(m,n) = 0.
          End Do
        End Do
        Do k=1,3
          dQk=Q(k,1)*Q(k,1)+Q(k,2)*Q(k,2)+Q(k,3)*Q(k,3)
          Do m=1,6
            Do n=1,6
              eta36(m,n)=eta36(m,n)+eta6(k)*Q(k,ik(n))*Q(k,jk(n))*&
               (Q(k,ik(m))*Q(k,jk(m))-1./3.*coef(m)*dQk)    
              eta36(m,n) = eta36(m,n) - 2./3.*eta6(k+3)*coef(m)* &
               Q(k,ik(n))*Q(k,jk(n))
            End Do
          End Do

          eta36(1,1)=eta36(1,1)+eta6(k+3)*Q(k,1)*Q(k,1)*2.
          eta36(1,4)=eta36(1,4)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(1,6)=eta36(1,6)+eta6(k+3)*Q(k,1)*Q(k,3)

          eta36(2,2)=eta36(2,2)+eta6(k+3)*Q(k,2)*Q(k,2)*2.
          eta36(2,4)=eta36(2,4)+eta6(k+3)*Q(k,2)*Q(k,1)
          eta36(2,5)=eta36(2,5)+eta6(k+3)*Q(k,2)*Q(k,3)
          
          eta36(3,3)=eta36(3,3)+eta6(k+3)*Q(k,3)*Q(k,3)*2.
          eta36(3,5)=eta36(3,5)+eta6(k+3)*Q(k,3)*Q(k,2)
          eta36(3,6)=eta36(3,6)+eta6(k+3)*Q(k,3)*Q(k,1)

          eta36(4,1)=eta36(4,1)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,2)=eta36(4,2)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,4)=eta36(4,4)+eta6(k+3)*(dQk-Q(k,3)*Q(k,3))*0.5
          eta36(4,5)=eta36(4,5)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(4,6)=eta36(4,6)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5

          eta36(5,2)=eta36(5,2)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,3)=eta36(5,3)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,4)=eta36(5,4)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(5,5)=eta36(5,5)+eta6(k+3)*(dQk-Q(k,1)*Q(k,1))*0.5
          eta36(5,6)=eta36(5,6)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5

          eta36(6,1)=eta36(6,1)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,3)=eta36(6,3)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,4)=eta36(6,4)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5
          eta36(6,5)=eta36(6,5)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5
          eta36(6,6)=eta36(6,6)+eta6(k+3)*(dQk-Q(k,2)*Q(k,2))*0.5
        End Do
         
         End Subroutine ViscGene


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   sens = 1 : tri des ki pour avoir k1 < k2 < k3
!   sens =-1 : tri des visc apres calcul avec k1 < k2 < k3 
!
!
       Subroutine triki(ki0,ki,visc,ordre,sens)

       
       
       Implicit None
       
       Double Precision, Dimension(3) ::  ki0,ki
       Double Precision, Dimension(6) :: visc,b
       Double Precision :: a
       Integer :: sens,i,j
       Integer, Dimension(3) :: ordre 
       
       

!
!    Passage pour trier les ki
!
       If (sens.EQ.1) Then
         Do i=1,3
           ki(i)=ki0(i)
           ordre(i)=i
         End Do
         Do j=2,3
          a=ki(j)
          Do i=j-1,1,-1
          If (ki(i).LE.a) Goto 20
          ki(i+1)=ki(i)
          ordre(i+1)=ordre(i)
          End Do
  20      Continue
         ki(i+1)=a
         ordre(i+1)=j
         End Do
!
!   Passage pour remettre les viscosite dans le bon ordre
!
       ElseIf (sens.EQ.-1) Then

         Do i=1,6
         b(i)=visc(i)
         End Do
         Do i=1,3
          visc(ordre(i))=b(i)
          visc(ordre(i)+3)=b(i+3)
        End Do

       Else
         Write(*,*)'triki.f : sens <> 1 ou -1' 
       Stop
       End If
       End Subroutine triki

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Interolation of Q
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     Subroutine InterP(t,x,Q,Ip)

      
!
      Implicit None
      Double Precision, Dimension(3) :: x,Q
      Double Precision :: t,d12,d23
      Double Precision, Intent(out) :: Ip

!
      d12=x(2)-x(1)
      d23=x(3)-x(2)
      Ip=Q(1)*(x(2)-t)*(x(3)-t)/((d12+d23)*d12)
      Ip=Ip+Q(2)*(t-x(1))*(x(3)-t)/(d12*d23)
      Ip=Ip+Q(3)*(t-x(1))*(t-x(2))/((d12+d23)*d23)


      
     End Subroutine InterP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Quadratic interpolation of a quantity Q defined on nine points 
!
!         y3 ---     Q7 ---- Q8 ------ Q9
!                    |        |         | 
!                    |        |         |
!         y2 ---     Q4 ----- Q5 ----- Q6
!                    |        |         | 
!                    |        |         |
!         y1 ---     Q1 ----- Q2 ----- Q3
!
!                    |        |         |
!                   x1        x2       x3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      Subroutine InterQ9(x,y,xi,yi,Q,Ip)

      
!
      Implicit None
      Double Precision, Intent(in), Dimension(3) ::  xi,yi
      Double Precision, Intent(in), Dimension(9) ::  Q
      Double Precision, Intent(out)  :: Ip
      !Double Precision ::  InterP
      Double Precision ::  x,y
      Double Precision, Dimension(3) :: a
      Integer  i

!
         Call InterP(x,xi,Q(1),a(1))
         Call InterP(x,xi,Q(4),a(2))
         Call InterP(x,xi,Q(7),a(3))
         Call InterP(y,yi,a, Ip)

      End Subroutine InterQ9

   Subroutine ReadVa(NetaI)
      Integer, Intent(in) :: NetaI
      !f2py depends(NetaI) etaI 
      !Double Precision, Dimension(:), Intent(out) :: etaI
      Double Precision, Dimension(4878) :: aetaI 

       Open(1,file='040010010.Va')
       Do i=1,813
          READ( 1, '(6(e14.8))' ) aetaI(6*(i-1)+1:6*(i-1)+6)
       End Do
       Open(2,file='ViscGrid.v')      
       Do i=1,4878   
        Write(2,'(e14.8)') aetaI(i)
       End Do
   End Subroutine ReadVa
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Orthotropic Polycrystalline Ice Law from LGGE            !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   give the Voigt viscosity matrix eta36 expressed 
!                                    in the general reference frame
!     Si = Aij dj 
!     Aij = eta36(i,j) i=1,6; j=1,6 non-symetric matrix 
!      where S=(S11,S22,S33,S12,S23,S31)      
!      and   d=(d11,d22,d33,2d12,2d23,2d31)      
!         as in Elmer
!
       
       Subroutine OPILGGE_ai_nl(ai,Angle,eta36, etaI)
      

!      
       
       Implicit None
       Double Precision, Intent(in), Dimension(3) :: ai       ! Fabric eigenvalues 
       Double Precision, Dimension(3) :: ki       ! Texture parameters 
       Double Precision, Intent(in), Dimension(3) :: Angle    ! Euler Angles               
       Double Precision, Dimension(4878) :: etaI ! Grid Viscosities 
       Double Precision, Intent(out), Dimension(6,6) :: eta36
       Double Precision, Dimension(6) :: eta6
       Double Precision :: aplusa
       Integer :: i
      



            Do i=1,3
              ki(i)=ai(i)
            End do
            Do i=1,3
             ki(i)=min(Max(ki(i),0.0),1.0)
            End do
            aplusa=ki(1)+ki(2)+ki(3)
            If (aplusa.GT.1.0) then
               !  write(*,*) 'depasse 1 dans R2R0',aplusa,ai(1),ai(2)
                    do i=1,3
                     ki(i)=ki(i)/aplusa
                    end do
            endif
       ! Value of the 6 relative viscosities in the orthotropic frame

       Call ViscMat_ai(ki,eta6,etaI)

       ! Viscosities in the reference frame

       Call ViscGene(eta6,Angle,eta36)

       End Subroutine OPILGGE_ai_nl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc									cccc
!ccccc       subroutine de calcul des viscosites Analytiques            cccc
!ccccc       Les viscosites ont ete calculees sur une grille            cccc
!ccccc       avec le sous-prog makeEtaI.f                               cccc
!ccccc         Elles passent dans etaI(6*814)=EtaI(4884)                cccc
!ccccc 									cccc      
!ccccc      !!! En entree a1,a2,a3 les valeurs propres du tenseur       cccc
!ccccc                                          d'orientation           cccc
!ccccc 									cccc      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
       Subroutine ViscMat_ai(ai0,eta6,etaI)
       
       
       
       Implicit None
       Integer, Parameter :: NetaI=4878 ! 6*4884 length of EtaI
       Double Precision, Intent(in), Dimension(3) :: ai0
       Double Precision, Intent(out), Dimension(6) :: eta6
       Double Precision, Intent(in), Dimension(NetaI) :: etaI
       Double Precision, Dimension(3) :: a1i,a2i
       !Double Precision, Intent :: IQ9
       Double Precision, Dimension(3) :: ai 
       Double Precision, Dimension(9) :: etaN
       Double Precision :: Delta
       Double Precision ::  a1,a2
       Double Precision, parameter ::  UnTier = 0.3333333333333333333333333333333333333333
       Integer, Dimension(3) :: ordre
       Integer :: i,j,n
       Integer :: ik1,ik2
       Integer :: N4,N5,N6
       Integer, Parameter :: Ndiv=30    ! Ndiv+2 Number of points along ik1      
       Nk2 =  (/ -1,  46,  93, 139, 183, 226, 267, 307, 345, 382,&
              417, 451, 483, 514, 543, 571, 597, 622, 645, 667,&
              687, 706, 723, 739, 753, 766, 777, 787, 795, 802,&
              807, 811/) 
  
       Delta = UnTier / Ndiv
!
! tri des ki 
! pour avoir k1 < k2 < k3 
!
!
        Call triki(ai0,ai,eta6,ordre,1)
!
        a1 = ai(1)
        a2 = ai(2)
!
! Position de a1,a2 dans EtaI
! calcul des indices ik1,ik2      
! ik1,ik2 indice du noeud en haut a droite du pt (a1,a2) 
!
         ik1 = Int((a1 + Delta)/Delta) + 1
         ik2 = Int((a2 + Delta)/Delta) + 1

! Si ik1 + 2ik2 -3(Ndiv+1) >0 => on est sur la frontiere a2=a3
! ( a1+2a2=1)
!  => ik1=ik1-1 sauf si ik1=2 , ik2=ik2-1         
!         
         If ((ik1+2*ik2-3*(Ndiv+1)).GE.0) Then
          If ((ik1.NE.2).And.((ik1+2*ik2-3*(Ndiv+1)).NE.0).And.  &
          (abs((a1-Delta*(ik1-1))/a1).GT.1.0E-5))  ik1=ik1-1
          ik2=ik2-1 
         End If
         If (ik1.EQ.1) ik1=2
!
! Indice N4,N5 et N6 des points 4,5 et 6 dans EtaI
!  
         
         N4 = Nk2(ik1-1) + ik2 - ik1 + 3         
         N5 = Nk2(ik1) + ik2 - ik1 + 2         
         N6 = Nk2(ik1+1) + ik2 - ik1 + 1         

!
! Remplissage etaN(1 a 9)
!  7 - 8 - 9  
!  4 - 5 - 6 
!  1 - 2 - 3
!
       Do i=1,3
       a1i(i)=Delta*(ik1-3+i)
       a2i(i)=Delta*(ik2-3+i)
       End Do
!
       Do n=1,6
         Do i=1,3              
           etaN(1+(i-1)*3) = etaI(6*(N4-3+i)+n)
           etaN(2+(i-1)*3) = etaI(6*(N5-3+i)+n)
           etaN(3+(i-1)*3) = etaI(6*(N6-3+i)+n)
         End Do
!
! interpolation sur le Q9  
!
!        If ( (a1 < a1i(1)) .OR. &
!             (a1 > a1i(3)) .OR. &
!             (a2 < a2i(1)) .OR. &
!             (a2 > a2i(3)) ) Then
!           write(*,*)a1,a1i
!           write(*,*)a2,a2i
!           Stop
!         End If
            
         Call InterQ9(a1,a2,a1i,a2i,etaN,eta6(n))

       End Do
!
! tri des eta
!
        Call triki(ai0,ai,eta6,ordre,-1)

       Return 
       End Subroutine Viscmat_ai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine IBOF(a2,a4)

       
       Implicit None
       Double Precision,dimension(6),intent(in):: a2  
       Double Precision,dimension(9),intent(out):: a4  
       Double Precision:: a_11,a_22,a_33,a_12,a_13,a_23
       Double Precision:: b_11,b_22,b_12,b_13,b_23
       Double Precision:: aPlusa

       Double Precision,dimension(21) :: vec
       Double Precision,dimension(3,21) :: Mat
       Double Precision,dimension(6) :: beta
       Double Precision :: Inv2,Inv3
       Integer :: i,j      

       
       a_11=a2(1)
       a_22=a2(2)
       a_33=a2(3)
       a_12=a2(4)
       a_23=a2(5)
       a_13=a2(6)

!       a_11=MIN(MAX(a2(1),0.0),1.0)
!       a_22=MIN(MAX(a2(2),0.0),1.0)
!       aPlusa=a_11+a_22
!       If(aPlusa.GT.1.0) then
!               a_11=a_11/aPlusa
!               a_22=a_22/aPlusa
!       EndIf
!       a_33=1.0-a_11-a_22
!       If(a2(4).GT.0.0) then
!               a_12=Min(a2(4),0.50)
!       Else
!               a_12=Max(a2(4),-0.50)
!       Endif
!       If(a2(5).GT.0.0) then
!               a_23=Min(a2(5),0.50)
!       Else
!               a_23=Max(a2(5),-0.50)
!       Endif
!       If(a2(6).GT.0.0) then
!               a_13=Min(a2(6),0.50)
!       Else
!               a_13=Max(a2(6),-0.50)
!       Endif
       
       

!       write(*,*) 'a2'
!       write(*,*) a2

       !Coefficients 

      Mat(1,1)=0.217774509809788e+02
      Mat(1,2)=-.297570854171128e+03
      Mat(1,3)=0.188686077307885e+04
      Mat(1,4)=-.272941724578513e+03
      Mat(1,5)=0.417148493642195e+03
      Mat(1,6)=0.152038182241196e+04
      Mat(1,7)=-.137643852992708e+04
      Mat(1,8)=-.628895857556395e+03
      Mat(1,9)=-.526081007711996e+04
      Mat(1,10)=-.266096234984017e+03
      Mat(1,11)=-.196278098216953e+04
      Mat(1,12)=-.505266963449819e+03
      Mat(1,13)=-.110483041928547e+03
      Mat(1,14)=0.430488193758786e+04
      Mat(1,15)=-.139197970442470e+02
      Mat(1,16)=-.144351781922013e+04
      Mat(1,17)=-.265701301773249e+03
      Mat(1,18)=-.428821699139210e+02
      Mat(1,19)=-.443236656693991e+01
      Mat(1,20)=0.309742340203200e+04
      Mat(1,21)=0.386473912295113e+00
      Mat(2,1)=-.514850598717222e+00
      Mat(2,2)=0.213316362570669e+02
      Mat(2,3)=-.302865564916568e+03
      Mat(2,4)=-.198569416607029e+02
      Mat(2,5)=-.460306750911640e+02
      Mat(2,6)=0.270825710321281e+01
      Mat(2,7)=0.184510695601404e+03
      Mat(2,8)=0.156537424620061e+03
      Mat(2,9)=0.190613131168980e+04
      Mat(2,10)=0.277006550460850e+03
      Mat(2,11)=-.568117055198608e+02
      Mat(2,12)=0.428921546783467e+03
      Mat(2,13)=0.142494945404341e+03
      Mat(2,14)=-.541945228489881e+04
      Mat(2,15)=0.233351898912768e+02
      Mat(2,16)=0.104183218654671e+04
      Mat(2,17)=0.331489412844667e+03
      Mat(2,18)=0.660002154209991e+02
      Mat(2,19)=0.997500770521877e+01
      Mat(2,20)=0.560508628472486e+04
      Mat(2,21)=0.209909225990756e+01
      Mat(3,1)=0.203814051719994e+02
      Mat(3,2)=-.283958093739548e+03
      Mat(3,3)=0.173908241235198e+04
      Mat(3,4)=-.195566197110461e+03
      Mat(3,5)=-.138012943339611e+03
      Mat(3,6)=0.523629892715050e+03
      Mat(3,7)=0.859266451736379e+03
      Mat(3,8)=-.805606471979730e+02
      Mat(3,9)=-.468711180560599e+04
      Mat(3,10)=0.889580760829066e+01
      Mat(3,11)=-.782994158054881e+02
      Mat(3,12)=-.437214580089117e+02
      Mat(3,13)=0.112996386047623e+01
      Mat(3,14)=0.401746416262936e+04
      Mat(3,15)=0.104927789918320e+01
      Mat(3,16)=-.139340154288711e+03
      Mat(3,17)=-.170995948015951e+02
      Mat(3,18)=0.545784716783902e+00
      Mat(3,19)=0.971126767581517e+00
      Mat(3,20)=0.141909512967882e+04
      Mat(3,21)=0.994142892628410e+00

       
       ! calcul des invariants
       Inv2=0.50*(1.0-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2.0*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
       Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
       
     ! polynome complet de degre 5 des 2 invariants.
         vec(1)=1.0
         vec(2)=Inv2
         vec(3)=vec(2)*vec(2)
         vec(4)=Inv3
         vec(5)=vec(4)*vec(4)
         vec(6)=vec(2)*vec(4)
         vec(7)=vec(3)*vec(4)
         vec(8)=vec(2)*vec(5)
         vec(9)=vec(2)*vec(3)
         vec(10)=vec(5)*vec(4)
         vec(11)=vec(9)*vec(4)
         vec(12)=vec(3)*vec(5)
         vec(13)=vec(2)*vec(10)
         vec(14)=vec(3)*vec(3)
         vec(15)=vec(5)*vec(5)
         vec(16)=vec(14)*vec(4)
         vec(17)=vec(12)*vec(2)
         vec(18)=vec(12)*vec(4)
         vec(19)=vec(2)*vec(15)
         vec(20)=vec(14)*vec(2)
         vec(21)=vec(15)*vec(4)

       ! calcul des beta_bar (cf annexe C Chung)
       ! attention beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
       !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

       ! calcul des trois beta en fonction du polynome
         beta(:)=0.0
         Do i=1,3
          Do j=1,21
            beta(i)=beta(i)+Mat(i,j)*vec(j)
          End do
         End do
          
       ! calcul des 3 autres pour avoir la normalisation
         beta(4)=3.0*(-1.0/7.0+beta(1)*(1.0/7.0+4.0*Inv2/7.0+8.0*Inv3/3.0)/5.0- &
                  beta(2)*(0.20-8.0*Inv2/15.0-14.0*Inv3/15.0)- &
                  beta(3)*(1.0/35.0-24.0*Inv3/105.0-4.0*Inv2/35.0+ &
                  16.0*Inv2*Inv3/15.0+8.0*Inv2*Inv2/35.0))/5.0

         beta(5)=6.0*(1.0-0.20*beta(1)*(1.0+4.0*Inv2)+ &
                  7.0*beta(2)*(1.0/6.0-Inv2)/5.0- &
                  beta(3)*(-0.20+2.0*Inv3/3.0+4.0*Inv2/5.0- &
                  8.0*Inv2*Inv2/5.0))/7.0

         beta(6)=-4.0*beta(1)/5.0-7.0*beta(2)/5.0- &
                   6.0*beta(3)*(1.0-4.0*Inv2/3.0)/5.0

        ! pour avoir les beta_bar
        Do i=1,6
         beta(i)=beta(i)/3.0
        End do
         beta(2)=beta(2)/2.0
         beta(5)=beta(5)/2.0
         beta(6)=beta(6)/2.0

        !! calcul des 5 b=a.a
        b_11=a_11*a_11+a_12*a_12+a_13*a_13
        b_22=a_22*a_22+a_12*a_12+a_23*a_23
        b_12=a_11*a_12+a_12*a_22+a_13*a_23
        b_13=a_11*a_13+a_12*a_23+a_13*a_33
        b_23=a_12*a_13+a_22*a_23+a_23*a_33

        !Calcul des 9 termes de a4

        a4(1)=3.0*beta(4)+6.0*beta(5)*a_11+3.0*beta(1)*a_11*a_11+&
         6.0*beta(2)*b_11+6.0*beta(6)*a_11*b_11+3.0*beta(3)*b_11*b_11
        a4(2)=3.0*beta(4)+6.0*beta(5)*a_22+3.0*beta(1)*a_22*a_22+&
         6.0*beta(2)*b_22+6.0*beta(6)*a_22*b_22+3.0*beta(3)*b_22*b_22

        a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2.0*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4.0*a_12*b_12)+&
         beta(3)*(b_11*b_22+2.0*b_12*b_12)


         a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2.0*a_12*a_13)+beta(2)*b_23+&
          beta(6)*(a_11*b_23+a_23*b_11+2.0*(a_12*b_13+a_13*b_12))+beta(3)*&
          (b_11*b_23+2.0*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2.0*a_12*a_23)+beta(2)*b_13+&
          beta(6)*(a_22*b_13+a_13*b_22+2.0*(a_12*b_23+a_23*b_12))+beta(3)*&
          (b_22*b_13+2.0*b_12*b_23)


         a4(6)=3.0*beta(5)*a_13+3.0*beta(1)*a_11*a_13+3.0*beta(2)*b_13+&
          3.0*beta(6)*(a_11*b_13+a_13*b_11)+3.0*beta(3)*b_11*b_13
         a4(7)=3.0*beta(5)*a_12+3.0*beta(1)*a_11*a_12+3.0*beta(2)*b_12+&
          3.0*beta(6)*(a_11*b_12+a_12*b_11)+3.0*beta(3)*b_11*b_12
         a4(8)=3.0*beta(5)*a_23+3.0*beta(1)*a_22*a_23+3.0*beta(2)*b_23+&
          3.0*beta(6)*(a_22*b_23+a_23*b_22)+3.0*beta(3)*b_22*b_23
         a4(9)=3.0*beta(5)*a_12+3.0*beta(1)*a_22*a_12+3.0*beta(2)*b_12+&
          3.0*beta(6)*(a_22*b_12+a_12*b_22)+3.0*beta(3)*b_22*b_12


          !modif fab 01/09/04
          ! si a11=0 => on annule les termes de a4
          !IF(a_11.EQ.00) then
                 ! write(*,*) 'a11 = 0'
           !       a4(1)=0.0
           !       a4(3)=0.0
           !       a4(4)=0.0
           !       a4(5)=0.0
           !       a4(6)=0.0
           !       a4(7)=0.0
           !       a4(9)=0.0
           !Endif
                  
          

         End 
!**************************************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine de calcul des 3 : - valeurs propres de a2
!                              - angles d'Euler
!
!--------------------------------------------------------
      subroutine R2Ro(a2,dim,ai,Euler)


      Implicit None

      Double Precision,dimension(6),intent(in) :: a2
      Integer :: dim               ! dimension  (2D-3D)
      Double Precision,dimension(3),intent(out) :: ai,Euler
      Double Precision :: disc
      Double Precision, parameter :: Prec=1.e-08
      Real, parameter :: Pi = 3.141592653589793
      Double Precision :: a,b,c,aplusa
      Complex :: Sol(3)
      Integer :: typeSol
      Integer :: i

      SELECT CASE(dim)

! ******** case 2D **********
      CASE(2)
           disc = (a2(1)+a2(2))*(a2(1)+a2(2)) - &
              4.*(a2(1)*a2(2)-a2(4)*a2(4))

           IF (abs(disc).LT.Prec) disc=0.0
           
           IF (disc.LT.0.0) THEN
              print *, 'Pb to calculate eigenvalues of a2',disc
              stop
           END IF

            disc = sqrt(disc)
        
            ai(2) = 0.5 * (a2(1)+a2(2)+disc)
            ai(1) = 0.5 * (a2(1)+a2(2)-disc)
            ai(3) = 1.0-ai(1)-ai(2)

            Euler=0.0
            Euler(1) = atan2(a2(4),(a2(1)-ai(2)))


            Do i=1,3
             ai(i)=Max(Min(ai(i),1.0),0.0)
            End do

! *********** case 3D **************
      CASE(3)
        ! a2 in the orthotropic frame.
        IF ((abs(a2(4)).LT.prec).and.(abs(a2(5)).LT.prec) &
           .and.(abs(a2(6)).LT.prec) ) THEN
              ai(1:3)=a2(1:3)
              Euler=0.
        ! coord 3 = eigenvector      
        ELSE IF ((abs(a2(5)).LT.prec).and.(abs(a2(6)).LT.prec)) THEN

           disc = (a2(1)+a2(2))*(a2(1)+a2(2)) - &
              4.*(a2(1)+a2(2)+a2(4)*a2(4))
              

           IF (disc.LT.0.0) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)

            
            ! 2 eigenvalues
            ai(2) = 0.5 * (a2(1)+a2(2)+disc)
            ai(1) = 0.5 * (a2(1)+a2(2)-disc)

            ! eigenvalues beetwen 0 and 1

            Do i=1,2
             ai(i)=Max(ai(i),0.0)
            End do
            aplusa=ai(1)+ai(2)
            If (aplusa.GT.1.0) then
                    do i=1,2
                     ai(i)=ai(i)/aplusa
                    end do
            endif

            ai(3) = 1.0-ai(1)-ai(2)

            Euler=0.
            Euler(1) = atan2(a2(4),a2(1)-ai(2))
                
         ! coord 1 = eigenvector      
         ELSE IF ((abs(a2(4)).LT.prec).and.(abs(a2(6)).LT.prec)) THEN

           disc = (a2(2)+a2(3))*(a2(2)+a2(3)) - &
              4.*(a2(2)+a2(3)+a2(5)*a2(5))

           IF (disc.LT.0.0) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)
            
            ai(2) = 0.5 * (a2(1)+a2(3)+disc)
            ai(3) = 0.5 * (a2(1)+a2(3)-disc)
        
            Do i=2,3
             ai(i)=Max(ai(i),0.0)
            End do
            aplusa=ai(2)+ai(3)
            If (aplusa.GT.1.0) then
                    do i=2,3
                     ai(i)=ai(i)/aplusa
                    end do
            endif

            ai(1) = 1.0-ai(2)-ai(3)

            Euler=0.
            Euler(2) = atan2(a2(5),a2(2)-ai(3))

         ! coord 2 = eigenvector      
         ELSE IF ((abs(a2(4)).LT.prec).and.(abs(a2(5)).LT.prec)) THEN
         
           disc = (a2(1)+a2(3))*(a2(1)+a2(3)) - &
              4.*(a2(1)+a2(3)+a2(6)*a2(6))

           IF (disc.LT.0.0) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

           disc = sqrt(disc)

            ai(3) = 0.5 * (a2(1)+a2(3)+disc)
            ai(1) = 0.5 * (a2(1)+a2(3)-disc)
        
            Do i=1,3,2
             ai(i)=Max(ai(i),0.0)
            End do
           aplusa=ai(1)+ai(3)
            If (aplusa.GT.1.0) then
                    do i=1,3,2
                     ai(i)=ai(i)/aplusa
                   end do
            endif

            ai(2) = 1.0-ai(1)-ai(3)

           
            Euler(1) = -0.5*Pi
            Euler(2) = atan2(a2(6),a2(1)-ai(3))
            Euler(3)= 0.5*Pi

         ! general case
        ELSE

          print *, 'Eigenvalues of a2, general 3D, not yet implemented, sorry'
          stop 
                
         
        END IF

      CASE DEFAULT
        write(*,*) 'dimens',dim
        print *, 'Eigenvalues of a2, dimension not available'
      
      END SELECT
      

      RETURN
      END
        


End Module golf
