      FUNCTION POLJ(A,B,N,X)
C=======================================================================
C  NORMALIZED JACOBI POLYNOMIAL
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NDIM=199)
      DIMENSION P1(2:NDIM),P2(2:NDIM),P3(2:NDIM) ,PL(0:NDIM)
      AB=A+B
C
C --- If standard function DGAMMA(X) is available then
c      TMP=DGAMMA(AB+2.D0)/DGAMMA(A+1.D0)/DGAMMA(B+1.D0)
C --- else
c      TMP=DEXP(GAMMLN(AB+2.D0)-GAMMLN(A+1.D0)-GAMMLN(B+1.D0))
C --- or
      TMP=DEXP(DLGAMMA(AB+2.D0)-DLGAMMA(A+1.D0)-DLGAMMA(B+1.D0))
C
      C0=DSQRT(TMP/2.D0**(AB+1.D0))
      TMP=0.5D0*DSQRT((AB+3.D0)/(A+1.D0)/(B+1.D0))*C0
      C1=(A-B)*TMP
      C2=(AB+2.D0)*TMP
      PL(0) = C0
      PL(1) = C1 + C2 * X
      DO i=2,N
        AI=DBLE(i)
        AIAB=2.D0*AI+AB
        TMP=AI*(AI+A)*(AI+B)*(AI+AB)
        TMP1=0.5D0*DSQRT((AIAB*AIAB-1.D0)/TMP)        
        P1(i)=(A*A-B*B)/(AIAB-2.D0)*TMP1
        P2(i)=AIAB*TMP1
        AI1=AI-1.D0
        TMP1=AI1*(AI1+A)*(AI1+B)*(AI1+AB)
        P3(i)=AIAB/(AIAB-2.D0)*DSQRT(TMP1/TMP*(AIAB+1.D0)/(AIAB-3.D0))
        PL(i) = (P1(i)+P2(i)*X)*PL(i-1)-P3(i)*PL(i-2)
      ENDDO
      POLJ = PL(N)
      RETURN
C     --------------------
      END
