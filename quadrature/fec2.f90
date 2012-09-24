Module fec2

Subroutine LebedevQuad(funvals,leb_data,integral)
   Double Precision, Intent(in), Dimension(:) :: funvals,leb_data
   Double Precision, Intent(out) :: integral
   Integer :: i

      integral=Dot_Product(funvals,leb_data)
      
End Subroutine LebedevQuad

Subroutine GetA(density,leb_data,m,n,A)
   Double Precision, Intent(in), Dimension(:,:) :: density
   Double Precision, Intent(in), Dimension(:) :: leb_data
   Double Precision, Intent(out), Dimension(:,6) :: A
   Integer, Intent(in) :: m,n
   
   Double Precision, Dimension(n) :: x,y,z


   x=sin(leb_data(:,1))*cos(leb_data(:,2))
   y=sin(leb_data(:,2))*sin(leb_data(:,1))
   z=cos(leb_data(:,1))



   ! $OMP PARALLEL DO
   Do i=1:m
      funvals(i,:)=x(:)*x(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,1))
      
      funvals(i,:)=x(:)*y(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,2))
 
      funvals(i,:)=x(:)*z(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,3))
 
      funvals(i,:)=y(:)*y(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,4))
 
       funvals(i,:)=y(:)*z(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,5))
       
      funvals(i,:)=z(:)*z(:)*density(i,:)
      Call LebedevQuad(funvals,leb_data,A(i,6))
   End Do

   

Subroutine dsyar(A,Q,W,m)
   Double Precision, Intent(in), Dimension(:,:) :: A
   Double Precision, Intent(out), Dimension(:,:,:) :: Q
   Double Precision, Intent(out), Dimension(:,:) :: W
   Integer, Intent(in) :: m
   Integer :: j
   Double Precision, Dimension(3,3) :: Aj

   Do j=1:m
      Aj(1,1)=A(j,1)
      Aj(1,2)=A(j,2)
      Aj(1,3)=A(j,3)
      Aj(2,2)=A(j,4)
      Aj(2,3)=A(j,5)
      Aj(3,3)=A(j,6)
      Call DSYEVJ3(Aj,Q,W)
   End DO
End Subroutine dsyar


!This subroutine (C) 2006 Joachim Kopp (LGPL)
!www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
!Joachim Kopp
!Efficient numerical diagonalization of hermition 3x3 matrix
!Int. J. Mod. Phys. C 19 (2008) 523-548
!arXiv.org: physics/0610206
!This is more efficient that LAPACK for this purpose.
* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
      END SUBROUTINE
* End of subroutine DSYEVJ3

   Double Precision Subroutine RotSym(W,Q,A)
!     This rotates a Hermitian diagonalized matrix (3x1 array of eigenvals) back 
!     into the undiagonalized reference frame. Q has to be unitary. 
!     W=Q^T A Q ; A=Q A Q^T
!     Returns 6x1 array A.
      
      Double Precision, Intent(in), Dimension(3,3) :: Q
      Double Precision, Intent(in), Dimension(3) :: W
      Double Precision, Intent(out), Dimension(6) :: A
      
      Qt=transpose(Q)
      !Saves a touch of time by using the fact that W is diagonal.
      !6 dot products instead of 6**3

      A(1)=DOT_PRODUCT(Q(1,:),W(1)*Qt(:,1))
      A(2)=DOT_PRODUCT(Q(1,:),W(2)*Qt(:,2))
      A(3)=DOT_PRODUCT(Q(1,:),W(3)*Qt(:,3))
      A(4)=DOT_PRODUCT(Q(2,:),W(2)*Qt(:,2))
      A(5)=DOT_PRODUCT(Q(5,:),W(2)*Qt(:,3))
      A(6)=DOT_PRODUCT(Q(6,:),W(3)*Qt(:,3))

   End Subroutine RotSym




C          ********************************************************
C
      DOUBLE PRECISION FUNCTION RD(X,Y,Z,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES AN INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE SECOND KIND,
C          RD(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                -1/2     -1/2     -3/2
C                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,
C
C          WHERE X AND Y ARE NONNEGATIVE, X + Y IS POSITIVE, AND Z IS
C          POSITIVE.  IF X OR Y IS ZERO, THE INTEGRAL IS COMPLETE.
C          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
C          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
C          SERIES TO FIFTH ORDER.  REFERENCE: B. C. CARLSON, COMPUTING
C          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
C          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
C          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
C          MARCH 1, 1980..
C
C          CHECK: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
C          = 3 / DSQRT(X * Y * Z), WHERE X, Y, AND Z ARE POSITIVE.
C
      INTEGER IERR,PRINTR
      DOUBLE PRECISION C1,C2,C3,C4,EA,EB,EC,ED,EF,EPSLON,ERRTOL,LAMDA
      DOUBLE PRECISION LOLIM,MU,POWER4,SIGMA,S1,S2,UPLIM,X,XN,XNDEV
      DOUBLE PRECISION XNROOT,Y,YN,YNDEV,YNROOT,Z,ZN,ZNDEV,ZNROOT
C          INTRINSIC FUNCTIONS USED: DABS,DMAX1,DMIN1,DSQRT
C
C          PRINTR IS THE UNIT NUMBER OF THE PRINTER.
      DATA PRINTR/6/
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ** (2/3).
C          UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
C          MINIMUM) ** (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
C          IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
C          NEVER BE CHOSEN SMALLER THAN 1.D-5.
C
C          ACCEPTABLE VALUES FOR:   LOLIM      UPLIM
C          IBM 360/370 SERIES   :   6.D-51     1.D+48
C          CDC 6000/7000 SERIES :   5.D-215    2.D+191
C          UNIVAC 1100 SERIES   :   1.D-205    2.D+201
C
C          WARNING: IF THIS PROGRAM IS CONVERTED TO SINGLE PRECISION,
C          THE VALUES FOR THE UNIVAC 1100 SERIES SHOULD BE CHANGED TO
C          LOLIM = 1.E-25 AND UPLIM = 2.E+21 BECAUSE THE MACHINE
C          EXTREMA CHANGE WITH THE PRECISION.
C
      DATA LOLIM/6.D-51/, UPLIM/1.D+48/
C
C          ON INPUT:
C
C          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RD(X,Y,Z).
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    4.D-18
C                           3.D-3    3.D-15
C                           1.D-2    4.D-12
C                           3.D-2    3.D-9
C                           1.D-1    4.D-6
C
C          ON OUTPUT:
C
C          X, Y, Z, AND ERRTOL ARE UNALTERED.
C
C          IERR IS THE RETURN ERROR CODE:
C               IERR = 0 FOR NORMAL COMPLETION OF THE SUBROUTINE,
C               IERR = 1 FOR ABNORMAL TERMINATION.
C
C          ********************************************************
C          WARNING: CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
C          EXPENSE OF ROBUSTNESS.
C
      IF (DMIN1(X,Y) .LT. 0.D0) GO TO 100
      IF (DMIN1(X+Y,Z) .LT. LOLIM) GO TO 100
      IF (DMAX1(X,Y,Z) .LE. UPLIM) GO TO 112
  100 WRITE(PRINTR,104)
  104 FORMAT(1H0,42H*** ERROR - INVALID ARGUMENTS PASSED TO RD)
      WRITE(PRINTR,108) X,Y,Z
  108 FORMAT(1H ,4HX = ,D23.16,4X,4HY = ,D23.16,4X,4HZ = ,D23.16)
      IERR = 1
      GO TO 124
C
  112 IERR = 0
      XN = X
      YN = Y
      ZN = Z
      SIGMA = 0.D0
      POWER4 = 1.D0
C
  116 MU = (XN + YN + 3.D0 * ZN) * 0.2D0
      XNDEV = (MU - XN) / MU
      YNDEV = (MU - YN) / MU
      ZNDEV = (MU - ZN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 120
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      SIGMA = SIGMA + POWER4 / (ZNROOT * (ZN + LAMDA))
      POWER4 = POWER4 * 0.25D0
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      GO TO 116
C
  120 C1 = 3.D0 / 14.D0
      C2 = 1.D0 / 6.D0
      C3 = 9.D0 / 22.D0
      C4 = 3.D0 / 26.D0
      EA = XNDEV * YNDEV
      EB = ZNDEV * ZNDEV
      EC = EA - EB
      ED = EA - 6.D0 * EB
      EF = ED + EC + EC
      S1 = ED * (- C1 + 0.25D0 * C3 * ED - 1.5D0 * C4 * ZNDEV * EF)
      S2 = ZNDEV * (C2 * EF + ZNDEV * (- C3 * EC + ZNDEV * C4 * EA))
      RD = 3.D0 * SIGMA + POWER4 * (1.D0 + S1 + S2) / (MU * DSQRT(MU))
C
  124 RETURN
      END














End Module fec2



