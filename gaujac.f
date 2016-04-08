      SUBROUTINE GAUJAC(A,B,N,X,W,WK)
C=======================================================================
C  GAUSS-JACOBI ABSCISSAS X(i) AND WEIGHTS W(i); i=1,N
C  WK(i) - WORKING ARRAY OF DIMENSION.GE.(N+1)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (EPS=1.D-14,ITRMAX=150)
      DIMENSION X(*),W(*),WK(*)
C
      AN=DBLE(N)
      AB=A+B
      ANAB=2.D0*AN+AB
      P1=(A-B)*AN/ANAB
      P2=-AN
      TMP=DSQRT(ANAB+1.D0)
      TMP1=2.D0*DSQRT(AN*(AN+A)*(AN+B)*(AN+AB)/(ANAB-1.D0))/ANAB
      P3=TMP*TMP1
      WW=TMP/TMP1
      DO 1 i=N,1,-1
        IF(i.EQ.N) THEN
          AAN=A/AN
          BAN=B/AN
          R1=(1.D0+A)*(2.78D0/(4.D0+AN*AN)+0.768D0*AAN/AN)
          R2=1.D0+1.48D0*AAN+0.96D0*BAN+0.452D0*AAN*AAN+0.83D0*AAN*BAN
c  The following prescription for R1 and R2 seems to work well even
c  when the previous one fails - ???
c         R1=2.D0*(A+1.D0)
c         R2=(AN+1.D0)*(AN+AB+2.D0)
          Z=1.D0-R1/R2
        ELSE IF(i.EQ.N-1) THEN
          R1=(4.1D0+A)/((1.D0+A)*(1.D0+0.156D0*A))
          R2=1.D0+0.06D0*(AN-8.D0)*(1.D0+0.12D0*A)/AN
          R3=1.D0+0.012D0*B*(1.D0+0.25D0*DABS(A))/AN
          Z=Z-(1.D0-Z)*R1*R2*R3
        ELSE IF(i.EQ.N-2) THEN
          R1=(1.67D0+0.28D0*A)/(1.D0+0.37D0*A)
          R2=1.D0+0.22D0*(AN-8.D0)/AN
          R3=1.D0+8.D0*B/((6.28D0+B)*AN*AN)
          Z=Z-(X(N)-Z)*R1*R2*R3
        ELSE IF(i.EQ.2) THEN
          R1=(1.D0+0.235D0*B)/(0.766D0+0.119D0*B)
          R2=1.D0/(1.D0+0.639D0*(AN-4.D0)/(1.D0+0.71D0*(AN-4.D0)))
          R3=1.D0/(1.D0+20.D0*A/((7.5D0+A)*AN*AN))
          Z=Z+(Z-X(4))*R1*R2*R3
        ELSE IF(i.EQ.1) THEN
          R1=(1.D0+0.37D0*B)/(1.67D0+0.28D0*B)
          R2=1.D0/(1.D0+0.22D0*(AN-8.D0)/AN)
          R3=1.D0/(1.D0+8.D0*A/((6.28D0+A)*AN*AN))
          Z=Z+(Z-X(3))*R1*R2*R3
        ELSE
          Z=3.D0*(X(i+1)-X(i+2))+X(i+3)
        ENDIF
      itr=0
      DZ=0.D0
 99   itr=itr+1
      IF(itr.GT.ITRMAX) STOP ' *** GAUJAC ITERATIONS ERROR'
      Z1=Z
      Z=Z1-DZ
      CALL JACPOL(Z,WK(1),N)
      PP=((P1+P2*Z)*WK(N+1)+P3*WK(N))/(1.D0-Z*Z)
      DZ=WK(N+1)/PP
      IF(DABS(DZ).GT.EPS) GOTO 99
      X(i)=Z
      W(i)=WW/(WK(N)*PP)
 1    CONTINUE
C
      RETURN 
      END
C=======================================================================