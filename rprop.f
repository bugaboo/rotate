      SUBROUTINE RPROP(SDIR, ISEC, CE, CRMAT, CDMAT)
      USE STARK
C=======================================================================
C  Propagating R(ETA) from RAD1 to RAD2 and constructing D(RAD2,RAD1)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION CRMAT(NBAS,*),CDMAT(NBAS,*)
      ALLOCATABLE CR11(:,:),CR22(:,:),CR12(:,:), CER(:)
C  Used variables:
C  NBAS, NSVD, EIGSVD,      
C
C  Matrices R11, R22, and R12
C
      ALLOCATE(CR11(NBAS,NBAS),CR22(NBAS,NBAS),CR12(NBAS,NBAS))
      ALLOCATE(CER(NSVD))
      DO i = 1, NSVD
        CER(i) = 0.5D0 / (DCMPLX(ESVD(i, ISEC), 0.D0) - CE)
      ENDDO
      DO mu=1,NBAS
        DO nu=1,mu
          CTMP1=0.D0
          CTMP2=0.D0
          DO n=1,NSVD
            CTMP1=CTMP1+CER(n)*CFSVD1(n,nu,ISEC)*CFSVD1(n,mu,ISEC)
            CTMP2=CTMP2+CER(n)*CFSVD2(n,nu,ISEC)*CFSVD2(n,mu,ISEC)
          ENDDO
          CR11(nu,mu)=CTMP1
          CR22(nu,mu)=CTMP2
        ENDDO
        DO nu=1,NBAS
          CTMP=0.D0
          DO n=1,NSVD
            CTMP=CTMP+CER(n)*CFSVD1(n,nu,ISEC)*CFSVD2(n,mu,ISEC)
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
C  Propagation R(RAD1)->R(RAD2)
C
      IF(ISEC.EQ.1) THEN
        DO mu=1,NBAS
          DO nu=1,NBAS
            CRMAT(nu,mu)=CR22(nu,mu)
          ENDDO
        ENDDO
        GOTO 10
      ENDIF
C --- Matrix A=R(RAD1)+SDIR*R11 to be inverted
      SDIR = DSIGN(1.D0, SDIR)
      DO mu=1,NBAS
        DO nu=1,NBAS
          CRMAT(nu,mu)=CRMAT(nu,mu)+SDIR*CR11(nu,mu)
        ENDDO
      ENDDO
C --- Inversion of A and calculation of B=A*R12
      CALL LAPEQNCG(NBAS,CRMAT,NBAS,CR12,CR11)
C --- Calculation of R(RAD2)=SDIR*R22-R21*B
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
C  Matrix D(RAD2,RAD1)
C
      TMP=-SDIR*DSQRT(RAD1/RAD2)
      DO mu=1,NBAS
        DO nu=1,NBAS
          CR11(nu,mu)=CRMAT(nu,mu)-SDIR*CR22(nu,mu)
          CR22(nu,mu)=TMP*CR12(mu,nu)
        ENDDO
      ENDDO
      CALL LAPEQNCG(NBAS,CR11,NBAS,CR22,CDMAT)
C
 10   DEALLOCATE(CR11,CR22,CR12,CER)
      RETURN
      END
C=======================================================================
