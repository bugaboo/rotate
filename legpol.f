      SUBROUTINE LEGPOM(NMAX)
C=======================================================================
C  NORMALIZED LEGENDRE POLYNOMIALS 
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NDIM=2999)
      DIMENSION P2(2:NDIM),P3(2:NDIM),PL(0:*)
      COMPLEX*16 CX,CPL(0:*)
      SAVE C0,C2,P2,P3
C
      IF(NMAX.GT.NDIM) STOP ' *** LEGPOM DIMENSION ERROR'
      C0=DSQRT(0.5D0)
      C2=DSQRT(1.5D0)
      DO 1 i=2,NMAX
        AI=DBLE(i)
        P2(i)=2.D0*DSQRT(AI*AI-0.25D0)/AI
        P3(i)=(1.D0-1.D0/AI)*DSQRT((AI+0.5D0)/(AI-1.5D0))
 1    CONTINUE
      RETURN 
C     --------------------
      ENTRY LEGPOL(X,PL,N)
C     --------------------
      PL(0)=C0
      PL(1)=C2*X
      DO i=2,N
        PL(i)=P2(i)*X*PL(i-1)-P3(i)*PL(i-2)
      ENDDO
      RETURN 
C     -----------------------
      ENTRY CLEGPOL(CX,CPL,N)
C     -----------------------
      CPL(0)=C0
      CPL(1)=C2*CX
      DO i=2,N
        CPL(i)=P2(i)*CX*CPL(i-1)-P3(i)*CPL(i-2)
      ENDDO
C
      RETURN 
      END
C=======================================================================