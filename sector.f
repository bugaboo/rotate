      SUBROUTINE SECTOR(ISEC)
     
      USE STARK
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)

      ALLOCATABLE :: RAD(:), OVLP(:,:), EIGANGL(:), EIGANGR(:)
      ALLOCATABLE :: EIGANG(:,:), VECANG(:,:,:), VECANGL(:,:)
      ALLOCATABLE :: VECANGR(:,:)
  
      RADL = RSEC(ISEC - 1)
      RADR = RSEC(ISEC)

      ALLOCATE(RAD(NDVR))
      DO i = 1, NDVR
        RAD(i) = SR * XR(i) + (RADR + RADL) / 2.D0
      ENDDO
C --- Angular EVPs      
      ALLOCATE(EIGANG(NBAS,NDVR), VECANG(NBAS,NBAS,NDVR))
      DO i = 1, NDVR
        CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RAD(i),
     &            EIGANG(:,i),VECANG(:,:,i))
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
      CALL SVDEVP(RADL, RADR, RAD, OVLP, ESVD(:, ISEC), EIGANG)

      ALLOCATE (EIGANGL(NBAS),VECANGL(NBAS,NBAS))
      CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RADL,EIGANGL,VECANGL)
      ALLOCATE (EIGANGR(NBAS),VECANGR(NBAS,NBAS))
      CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RADR,EIGANGR,VECANGR)
C --- Surface amplitudes of R-matrix eigenfunctions. Boundary overlap      
      DO n = 1, NSVD
        DO nu = 1, NBAS
          TMP1 = 0.D0
          TMP2 = 0.D0
          DO i = 1, NDVR
            DO mu = 1, NBAS
              DO k = 1, NBAS
                TMP1 = TMP1 + PIL(i) * VECANG(k, mu, i) * VECANGL(k,nu) 
     &                * OVLP(NBAS * (i-1) + mu, n)
                TMP2 = TMP2 + PIR(i) * VECANG(k, mu, i) * VECANGR(k,nu) 
     &                * OVLP(NBAS * (i-1) + mu, n)
              ENDDO
            ENDDO
          ENDDO
          TMP1 = TMP1 * RADL
          TMP2 = TMP2 * RADR
          IF (ISEC .LE. MSEC) THEN 
            CFSVD1(n, nu, ISEC) = DCMPLX(TMP1)
            CFSVD2(n, nu, ISEC) = DCMPLX(TMP2)
          ELSE
            CFSVD1(n, nu, ISEC) = DCMPLX(TMP2)
            CFSVD2(n, nu, ISEC) = DCMPLX(TMP1)
          ENDIF
        ENDDO
      ENDDO

          
      DEALLOCATE (OVLP, RAD, EIGANG, VECANG)
      END SUBROUTINE
