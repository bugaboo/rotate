      SUBROUTINE RPROP(ETA1,ETA2,CESVD,CFSVD1,CFSVD2,CRMAT,CDMAT)
      USE STARK
C=======================================================================
C  Propagating R(ETA) from ETA1 to ETA2 and constructing D(ETA2,ETA1)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CESVD(*),CFSVD1(NSVD,*),CFSVD2(NSVD,*),
     &          CRMAT(NBAS,*),CDMAT(NBAS,*)
      ALLOCATABLE CR11(:,:),CR22(:,:),CR12(:,:)
      COMMON /DIM_STARK/NDVRX,NDVRE,NDVRSE,NBAS,NSVD /KEY_STARK/KRUN
C  Used variables:
C  NBAS, NSVD, EIGSVD,      
C
C  Matrices R11, R22, and R12
C
      ALLOCATE(CR11(NBAS,NBAS),CR22(NBAS,NBAS),CR12(NBAS,NBAS))
      DO mu=1,NBAS
        DO nu=1,mu
          CTMP1=0.D0
          CTMP2=0.D0
          DO n=1,NSVD
            CTMP1=CTMP1+CESVD(n)*CFSVD1(n,nu)*CFSVD1(n,mu)
            CTMP2=CTMP2+CESVD(n)*CFSVD2(n,nu)*CFSVD2(n,mu)
          ENDDO
          CR11(nu,mu)=CTMP1
          CR22(nu,mu)=CTMP2
        ENDDO
        DO nu=1,NBAS
          CTMP=0.D0
          DO n=1,NSVD
            CTMP=CTMP+CESVD(n)*CFSVD1(n,nu)*CFSVD2(n,mu)
          ENDDO
          CR12(nu,mu)=CTMP
        ENDDO
      ENDDO
      DO mu=1,NBAS
        DO nu=mu+1,NBAS
          CR11(nu,mu)=CR11(mu,nu)
          CR22(nu,mu)=CR22(mu,nu)
        ENDDO
      ENDDO
C
C  Propagation R(ETA1)->R(ETA2)
C
      IF(ETA1.EQ.0.D0) THEN
        DO mu=1,NBAS
          DO nu=1,NBAS
            CRMAT(nu,mu)=CR22(nu,mu)
          ENDDO
        ENDDO
        GOTO 10
      ENDIF
C --- Matrix A=R(ETA1)+SDIR*R11 to be inverted
      SDIR=DSIGN(1.D0,ETA2-ETA1)
      DO mu=1,NBAS
        DO nu=1,NBAS
          CRMAT(nu,mu)=CRMAT(nu,mu)+SDIR*CR11(nu,mu)
        ENDDO
      ENDDO
C --- Inversion of A and calculation of B=A*R12
      CALL LAPEQNCG(NBAS,CRMAT,NBAS,CR12,CR11)
C --- Calculation of R(ETA2)=SDIR*R22-R21*B
      DO mu=1,NBAS
        DO nu=1,mu
          CTMP=0.D0
          DO i=1,NBAS
            CTMP=CTMP+CR12(i,nu)*CR11(i,mu)
          ENDDO
          CRMAT(nu,mu)=SDIR*CR22(nu,mu)-CTMP
        ENDDO
      ENDDO
      DO mu=1,NBAS
        DO nu=mu+1,NBAS
          CRMAT(nu,mu)=CRMAT(mu,nu)
        ENDDO
      ENDDO
      IF(KRUN.EQ.1) GOTO 10
C
C  Matrix D(ETA2,ETA1)
C
      TMP=-SDIR*DSQRT(ETA1/ETA2)
      DO mu=1,NBAS
        DO nu=1,NBAS
          CR11(nu,mu)=CRMAT(nu,mu)-SDIR*CR22(nu,mu)
          CR22(nu,mu)=TMP*CR12(mu,nu)
        ENDDO
      ENDDO
      CALL LAPEQNCG(NBAS,CR11,NBAS,CR22,CDMAT)
C
 10   DEALLOCATE(CR11,CR22,CR12)
      RETURN
      END
C=======================================================================
