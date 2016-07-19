      SUBROUTINE SECTOR(ISEC)
     
      USE STARK
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)

      ALLOCATABLE :: RAD(:), OVLP(:,:), ESVD(:)
  
      RADL = WSEC * (ISEC - 1)
      RADR = WSEC * ISEC

      ALLOCATE(RAD(NDVR))
      DO i = 1, NDVR
        RAD(i) = SR * XR(i) + (RADR + RADL) / 2.D0
      ENDDO
      ALLOCATE(ANGEIGVAL(NBAS,NDVR), ANGEIGVEC(NBAS,NBAS,NDVR))
      DO i = 1, NDVR
        CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RAD(i),
     &            ANGEIGVAL(:,i),ANGEIGVEC(:,:,i))
      ENDDO
      ALLOCATE (OVLP(NSVD,NSVD))
      DO i = 1, NDVR
        DO j = 1, NDVR
          DO nu = 1, NBAS
            DO mu = 1, NBAS
              tnu = NBAS * (i - 1) + nu
              tmu = NBAS * (j - 1) + mu
              OVLP(tnu, tmu) = 0.D0
              DO k = 1, NBAS
                OVLP(tnu, tmu) = OVLP(tnu, tmu) + ANGEIGVEC(k, nu, i) *
     &                ANGEIGVEC(k, mu, j)                
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(ESVD(NSVD))
      CALL SVDEVP(RADL, RADR, RAD, OVLP, ESVD, ANGEIGVAL)
      DEALLOCATE (OVLP, RAD, ANGEIGVAL, ANGEIGVEC, ESVD)
      END SUBROUTINE
