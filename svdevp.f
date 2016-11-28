      SUBROUTINE SVDEVP(RADL, RADR, RAD, H, EIG, ANGEIG)
*   Purpose
C =================
C 
C   Calculates eigenvalues and eigenvectors of SVD eigenvalue 
C   problem in the sector between RADL and RADR.
C
* Arguments
C ================
C   
C RADL    (input) left end of the sector
C RADR    (input) right end of the sector. RADR >= RADL
C RAD     (input) array, dimension (NDVR)
C H       (input/output) array, dimension (NSVD, NSVD). On entry has
C          to be the correct overlap matrix, on exit contains converted 
C          normalized eigenvectors of SVD problem.
C EIG     (output) array, dimension(NSVD). On exit contains eigenvalues
C          of SVD problem.
C ANGEIG  (input) array, dimension(NBAS). Eigenvalues of angular part.
C======================================================================
      USE STARK
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION :: H(NSVD,NSVD), RAD(NDVR), ANGEIG(NBAS,NDVR), EIG(NSVD)

      ALLOCATABLE :: OVLP(:, :), DVRK(:), VK(:)
      INTEGER :: INFO

C --- Kinetic part
      ALLOCATE(DVRK(NDVRT))
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
      DEALLOCATE(DVRK)

C --- Diagonalization       
      SR = (RADR - RADL) / 2.D0
      NVK = 20 * NSVD
      ALLOCATE(VK(NVK))
C      CALL LAPEIGRS (0, NBAS, V, E, WK)
      CALL DSYEV('V','U',NSVD,H,NSVD,EIG,VK,NVK,INFO)
      IF(INFO.NE.0)  WRITE(*,*) "SVDEVP DSYEV ERROR",INFO
C --- Normalization (unnecessary?)
C      DO i = 1, NSVD
C        TMP = 0.D0
C        DO j = 1, NSVD
C          TMP = TMP + H(i, j) ** 2
C        ENDDO
C        TMP = DSQRT(SR * TMP)
C        DO j = 1, NSVD
C          H(i, j) = H(i, j) / TMP
C        ENDDO
C      ENDDO
C --- Conversion
      DO n = 1, NSVD
        DO i = 1, NSVD
          H(i, n) = H(i, n) / RAD(MOD(i, NBAS))
        ENDDO
      ENDDO
      END SUBROUTINE
