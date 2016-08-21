      SUBROUTINE DVRLEG(NDVR,X,W,T,KHAM,KSYM,H,NDIM)
C=======================================================================
C
C  LEGENDRE DISCRETE VARIABLE REPRESENTATION
C
C  NDVR      - NUMBER OF DVR QUADRATURE POINTS; IF IABS(KSYM)=1, THEN
C              NDVR MUST BE EVEN AND ON OUTPUT IT IS REPLACED BY NDVR/2  
C  X(i),W(i) - GAUSS-LEGENDRE ABSCISSAS AND WEIGHTS
C  T(n,i)    - FBR-DVR TRANSFORMATION MATRIX; n=POLYNOMIAL, i=POINT
C  KHAM      - 0 (1) - DON'T (DO) CALCULATE H(ij)
C  KSYM      - SYMMETRY KEY:
C              0 - IF BOTH EVEN AND ODD SOLUTIONS ARE NEEDED
C              +1 (-1) - IF ONLY EVEN (ODD) SOLUTIONS ARE NEEDED
C  H(ij)     - KINETIC PART OF THE DVR HAMILTONIAN IN SSM
C  NDIM      - DIMENSION PARAMETER, NDIM.GE.NDVR+1
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (EPS=1.D-14,ITRMAX=30,PI=3.1415926535897932385D0)
      DIMENSION X(*),W(*),T(NDIM,*),H(*)
C
C  Calculation of X(i), W(i), and T(n,i) 
C
      NDVR1=NDVR+1
      AN=DBLE(NDVR)
      P2=-AN
      TMP=DSQRT(AN+0.5D0)
      TMP1=AN/DSQRT(AN-0.5D0)
      P3=TMP*TMP1
      WW=2.D0*TMP/TMP1
      M=NDVR1/2
      DO 1 i=1,M
        Z=DCOS(PI*(DBLE(i)-0.25D0)/(AN+0.5D0))
        itr=0
        DZ=0.D0
 99     itr=itr+1
        IF(itr.GT.ITRMAX) STOP ' *** DVRLEG ITERATIONS ERROR'
        Z1=Z
        Z=Z1-DZ
        CALL LEGPOL(Z,T(1,i),NDVR)
        PP=(P2*Z*T(NDVR1,i)+P3*T(NDVR,i))/(1.D0-Z*Z)
        DZ=T(NDVR1,i)/PP
        IF(DABS(DZ).GT.EPS) GOTO 99
        X(i)=-Z
        W(i)=WW/(T(NDVR,i)*PP)
        iN=NDVR-i+1
        X(iN)=Z
        W(iN)=W(i)
        TMP=DSQRT(W(i))
        S=1.D0
        DO n=1,NDVR
          T(n,iN)=TMP*T(n,i)
          T(n,i)=S*T(n,iN)
          S=-S
        ENDDO
 1    CONTINUE
      IF(KHAM.EQ.0) RETURN
C
C  Calculation of H(ij)
C
      ij=0
      DO 2 j=1,NDVR
        DO 3 i=1,j
          ij=ij+1
          H(ij)=0.D0
          DO n=1,NDVR
            En=DBLE((n-1)*n)
            H(ij)=H(ij)+T(n,i)*En*T(n,j)
          ENDDO
 3      CONTINUE
 2    CONTINUE
      IF(KSYM.EQ.0) RETURN
      CALL SYMDVR(NDVR,H,KSYM)
C
      RETURN 
      END
C=======================================================================