      SUBROUTINE SVDEVP(RADL, RADR, RAD, H, EIG, ANGEIG)
      USE STARK
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION :: H(NSVD,NSVD), RAD(NDVR), ANGEIG(NBAS,NDVR), EIG(NSVD)
      ALLOCATABLE :: OVLP(:, :), DVRK(:), WK(:,:), DVR0(:), DVR1(:),
     &                DVR2(:), RHO(:), VK(:)
      INTEGER INFO

C --- Kinetic part
      NDVRD = NDVR + 1
      NDVRT = NDVR * NDVRD / 2 
      ALLOCATE(DVRK(NDVRT), DVR0(NDVRT), DVR1(NDVRT), DVR2(NDVRT))
      ALLOCATE(RHO(NDVRT), WK(NDVRD, NDVR))
      CALL DVRSL(NDVR,XR,WR,TR,DVR0,DVR1,DVR2,RHO,PIL,PIR,WK,NDVRD)
      A0 = (RADL + RADR) ** 2 / (RADR - RADL) ** 2 / 2.D0
      A1 = (RADL + RADR) / (RADR - RADL)
      A2 = 0.5D0
      ij=0
      DO j=1,NDVR
        DO i=1,j
          ij=ij+1
          DVRK(ij) = A0 * DVR0(ij) + A1 * DVR1(ij) + A2 * DVR2(ij)
          DVRK(ij) = DVRK(ij) / RAD(i) / RAD(j)
          ii = NBAS * (i - 1)
          jj = NBAS * (j - 1)
          DO nu = 1, NBAS
            DO mu = 1, NBAS
              H(ii + nu, jj + mu) = H(ii + nu, jj + mu) * DVRK(ij)
C Symmetric H              
C              if (i .NE. j) then
C                H(jj + nu, ii + mu) = H(jj + nu, ii + mu) * DVRK(ij)
C              endif
            ENDDO
C --- Angular eigenvalues            
            IF (i .EQ. j) THEN
              H(ii + nu, jj + nu) = H(ii + nu, jj + nu) + ANGEIG(nu, i)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C --- Diagonalization       
      SR = (RADR - RADL) / 2.D0
      NVK = 20 * NSVD
      ALLOCATE(VK(NVK))
C      CALL LAPEIGRS (0, NBAS, V, E, WK)
      CALL DSYEV('V','U',NSVD,H,NSVD,EIG,VK,NVK,INFO)
      IF(INFO.NE.0)  WRITE(*,*) "ANGEVP DSYEV ERROR",INFO
C --- Normalization (unnecessary?)
      DO i = 1, NSVD
        TMP = 0.D0
        DO j = 1, NSVD
          TMP = TMP + H(i, j) ** 2
        ENDDO
        TMP = DSQRT(SR * TMP)
        DO j = 1, NSVD
          H(i, j) = H(i, j) / TMP
        ENDDO
      ENDDO
C --- Conversion
      DO n = 1, NSVD
        DO i = 1, NSVD
          H(i, n) = H(i, n) / RAD(MOD(i, NBAS))
        ENDDO
      ENDDO
      DEALLOCATE(DVR0, DVR1, DVR2)
      END SUBROUTINE
