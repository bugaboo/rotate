      PROGRAM TEST_SVD

      USE STARK

      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      ALLOCATABLE OVLP(:,:), EIG(:), EANG(:, :), RAD(:)
      
      OMEGA = 1.D0
      CALL STARKI(10)
      ALLOCATE (RAD(NDVR))
      DO i = 1, NDVR
        RAD(i) = (XR(i) + 1.D0) * RMAX / 2.D0 
      ENDDO
      ALLOCATE(OVLP(NSVD, NSVD), EANG(NBAS,NDVR), EIG(NSVD))
      OVLP = 0.D0
      DO i = 1, NDVR
        DO j = 1, NDVR
          DO nu = 1, NBAS
            OVLP(NBAS * (i - 1) + nu, NBAS * (j - 1) + nu) = 1.D0 
          ENDDO
        ENDDO
      ENDDO
      PRINT *, MODEL, NDVR
      DO i = 1, NDVR
        DO nu = 1, NBAS
          EANG(nu, i) = L(nu) * (L(nu) + 1) / 2.D0 / RAD(i)**2 - OMEGA 
     &          * M(nu) + POT3D(RAD(i), 0.D0, 0.D0, 0.D0)
        ENDDO
      ENDDO
      CALL SVDEVP(0.D0, RMAX, RAD, OVLP, EIG, EANG)
      OPEN(1, file = 'out')
      DO i = 1, NSVD
        write(1, *) EIG(i)
      ENDDO
      CLOSE(1)
      END
