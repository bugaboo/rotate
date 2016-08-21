      SUBROUTINE DVRSL(NDVR,X,W,T,DVR0,DVR1,DVR2,RHO,PIL,PIR,WK,NDVRD)
C=======================================================================
C  SLOW VARIABLE DISCRETIZATION FOR SCATTERING IN FURTHER SECTORS: 
C  LEGENDRE POLYNOMIALS
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*),T(NDVRD,*),DVR0(*),DVR1(*),DVR2(*),RHO(*),
     &          PIL(*),PIR(*),WK(NDVRD,*)
C
      CALL LEGPOM(NDVR)
      CALL DVRLEG(NDVR,X,W,T,1,0,DVR2,NDVRD)
      DELTA=NDVR*NDVR/(4*NDVR*NDVR-1.D0)
C
C  Left (PIL) and right (PIR) end DVR basis amplitudes
C
      DO n=1,NDVR
        RHO(n)=DSQRT(n-0.5D0)
      ENDDO
      DO i=1,NDVR
        PIL(i)=0.D0
        PIR(i)=0.D0
        SS=1.D0
        DO n=1,NDVR
          TMP=T(n,i)*RHO(n)
          PIL(i)=PIL(i)+TMP*SS
          PIR(i)=PIR(i)+TMP
          SS=-SS
        ENDDO
      ENDDO
C
C  Kinetic matrices DVR0, DVR1, and DVR2
C
      DO k=1,NDVR
        RHO(k)=DSQRT(2*k-1.D0)
        RHO(NDVR+k)=k/DSQRT(4*k*k-1.D0)
        DO i=1,NDVR
          WK(i,k)=0.D0
        ENDDO
      ENDDO
      DO n=1,NDVR
        DO k=n-1,1,-2
          TMP=RHO(n)*RHO(k)
          DO i=1,NDVR
            WK(i,k)=WK(i,k)+T(n,i)*TMP
          ENDDO
        ENDDO
      ENDDO
      ij=0
      DO j=1,NDVR
        DO i=1,j
          ij=ij+1
          SUM=0.D0
          DO k=1,NDVR
            SUM=SUM+WK(i,k)*WK(j,k)
          ENDDO
          DVR0(ij)=SUM
          SUM=0.D0
          DO k=1,NDVR-1
            SUM=SUM+RHO(NDVR+k)*(WK(i,k)*WK(j,k+1)+WK(i,k+1)*WK(j,k))
          ENDDO
          DVR1(ij)=SUM
          DVR2(ij)=DVR0(ij)-DVR2(ij)
        ENDDO
      ENDDO
C
C  Weight matrix RHO
C
      ij=0
      DO j=1,NDVR
        DO i=1,j
          ij=ij+1
          RHO(ij)=DELTA*T(NDVR,i)*T(NDVR,j)
        ENDDO
        RHO(ij)=X(j)**2+RHO(ij)
      ENDDO
C
      RETURN
      END
C=======================================================================
