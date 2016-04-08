      SUBROUTINE CPOTMAT(NTET,NPHI,RAD,CV,NANG,L,M,LMAX,X0)
      USE OMP_LIB
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      INTEGER, DIMENSION(NANG) :: L,M
      DIMENSION CV(NANG,NANG)
      ALLOCATABLE X(:,:), W(:,:), T(:), XT(:), WT(:)
C
      ALLOCATE(X(NTET,0:2*LMAX),W(NTET,0:2*LMAX), T(2*NTET))
      ALLOCATE(XT(NTET), WT(NTET))
      PI=2*DACOS(0.D0)
      DO i=0,2*LMAX
        CALL JACPOM(DBLE(i)/2.D0,DBLE(i)/2.D0,NTET)
        CALL GAUJAC(DBLE(i)/2.D0,DBLE(i)/2.D0,NTET,XT,WT,T)
        DO j=1,NTET
          X(j,i) = XT(j)
          W(j,i) = WT(j)
        ENDDO
      ENDDO
    
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NANG,NPHI,NTET,L,M,CV,X,W,PI,
!$OMP& X0, RAD)
!$OMP DO
      DO i = 1,NANG
        DO j = 1,NANG
        CSUMP = (0.D0, 0.D0)
        DO k = 1,NPHI
          PHI = DBLE(2*k-1)*PI/DBLE(NPHI)
          M1 = ABS(M(i))
          M2 = ABS(M(j))
C          PRINT *, M1, M2, M(i), M(j)
          CTMP = CDEXP((0.D0, 1.D0) * DBLE(M(j) - M(i)) * PHI)
          IF (M(i).LT.0) THEN
            CTMP = (-1.D0)**M1 * CTMP
          ENDIF
          IF (M(j).LT.0) THEN
            CTMP = (-1.D0)**M2 * CTMP
          ENDIF
          CSUMT = (0.D0, 0.D0)
          DO n = 1,NTET
      CSUMT=CSUMT+W(n,M1+M2)*POLJ(DBLE(M1),DBLE(M1),L(i)-M1,X(n,M1+M2))
     &             *POLJ(DBLE(M2),DBLE(M2),L(j)-M2,X(n,M1+M2)) 
     &            *POT3D(RAD,DACOS(X(n,M1+M2)),PHI,X0)
C     & *DSIN(DACOS(X(n,M1+M2)))*CDEXP((0,-1.D0)*PHI)*dsqrt(0.375D0/PI)
          ENDDO
          CSUMP = CSUMP + CTMP*CSUMT
        ENDDO
        CV(i,j) = CSUMP/DBLE(NPHI)
      ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      DEALLOCATE(X,W,XT,WT,T)
      RETURN
      END
