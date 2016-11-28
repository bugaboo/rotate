      PROGRAM TEST_SVD_SECTOR
      USE STARK

      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      ALLOCATABLE OVLP(:,:), EIG(:), EANG(:, :), RAD(:)
      ALLOCATABLE :: EIGANG(:,:), VECANG(:,:,:), VECANGL(:,:)
      
      OMEGA = 1.D0
      CALL STARKI(10)
      
      X0 = 0.D0
      RADL = 0.D0
      RADR = RMAX
      ALLOCATE(RAD(NDVR))
      DO i = 1, NDVR
        RAD(i) = SR * XR(i) + (RADR + RADL) / 2.D0
      ENDDO
C --- Angular EVPs      
      ALLOCATE(EANG(NBAS,NDVR), VECANG(NBAS,NBAS,NDVR))
      DO i = 1, NDVR
        CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RAD(i),
     &            EANG(:,i),VECANG(:,:,i))
      ENDDO
C --- Overlap matrix
      ALLOCATE (OVLP(NSVD,NSVD))
      DO i = 1, NDVR
        DO j = 1, NDVR
          DO nu = 1, NBAS
            DO mu = 1, NBAS
              tnu = NBAS * (i - 1) + nu
              tmu = NBAS * (j - 1) + mu
              OVLP(tnu, tmu) = 0.D0
              DO k = 1, NBAS
                OVLP(tnu, tmu) = OVLP(tnu, tmu) + VECANG(k, nu, i) *
     &                VECANG(k, mu, j)                
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C --- SVDEVP      
      ALLOCATE(EIG(NSVD))
      CALL SVDEVP(RADL, RADR, RAD, OVLP, EIG, EANG)
 
      OPEN(1, file = 'out')
      DO i = 1, NSVD
        write(1, *) EIG(i)
      ENDDO
      CLOSE(1)
      END
