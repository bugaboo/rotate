      MODULE STARK
      
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)

      NAMELIST /INF_STARK/ITRMAX,EPS,NSEC,MSEC,KSYM,MODEL,RMAX,NDVR,LMAX
     &          ,NTET,NPHI 
      ALLOCATABLE :: L(:), M(:), PIL(:), PIR(:), VK1(:), VK2(:)
      ALLOCATABLE :: RSEC(:), XR(:), WR(:), TR(:,:), DVRR(:)
      ALLOCATABLE :: ANGEIGVAL(:,:), ANGEIGVEC(:,:), FSVD1(:,:,:),
     &               FSVD2(:,:,:)
      INTEGER :: NSVD, NDVR, NBAS
      COMMON /POT_C/MODEL
C      SAVE L, M, RSEC, XR, WR, TR, DVRR, ANGEIGVAL, ANGEIGVEC, PIL, PIR
C      SAVE RMAX, SR, WSEC, FSVD1, FSVD2
C      SAVE NDVR, NBAS, NSVD
      SAVE

      CONTAINS
      SUBROUTINE STARKI(IWF0)
C=======================================================================
C  
C  --------------------------------------------------------------------
C  STARKI initializes DVR arrays used by STARKEF and STARKA:
C    IWF0  - output channel for the eigenfunction.
C  STARKEF calculates the SS eigenvalue E(F) or its inverse F(E):
C    KEF     - KEF=1 - E(F), KEF=2 - F(E);
C    CE0,CF0 - initial guess (INPUT) and final result (OUTPUT) 
C              for E(F) (KEF=1) or F(E) (KEF=2);
C    IFAIL   - IFAIL=0 - normal return, IFAIL=1 - Newton method failed,
C              IFAIL=2 - too many iterations;
C    IITR    - output channel for intermediate iterations.
C  STARKA calculates the SS eigenfunction and TMD amplitude A(k_t):
C    CE0     - SS eigenvalue E(F);
C    CF0     - electric field F;
C    NKT     - number of points in transverse momentum k_t;
C    AKT     - values of k_t;
C    CTMD    - TMD amplitudes A(k_t).
C  ADDITIONAL OUTPUT FOR IWF0.NE.0: file 'stark_wf' containing the SS
C    eigenfunction.
C-----------------------------------------------------------------------
C  Externally defined common block /PI_C/PI is needed.
C  Temporary channels 11 and 12 are reserved.
C-----------------------------------------------------------------------
      USE ANGSYM
C
C  Input information
C
      OPEN(1,FILE='inf')
      READ(1,INF_STARK)
      CLOSE(1)


C
C  DVR arrays
C
C --- R
      NT = NDVR + 1
      ALLOCATE(XR(NDVR),WR(NDVR),TR(NT,NT),DVRR(NDVR*(NT)/2))
      CALL LEGPOM(NDVR)
      CALL DVRLEG(NDVR,XR,WR,TR,1,0,DVRR,NT)
C --- Left (PIL) and right (PIR) end DVR basis amplitudes
      ALLOCATE (PIL(NDVR), PIR(NDVR), VK1(NDVR), VK2(NDVR))
      DO n=1,NDVR
        TMP=DSQRT(n-0.5D0)
        VK1(n)=TMP
        VK2(n)=TMP
      ENDDO
      DO n=2,NDVR,2
        VK1(n)=-VK1(n)
      ENDDO
      DO i=1,NDVR
        PIL(i)=0.D0
        PIR(i)=0.D0
        DO n=1,NDVR
          PIL(i)=PIL(i)+TR(n,i)*VK1(n)
          PIR(i)=PIR(i)+TR(n,i)*VK2(n)
        ENDDO
      ENDDO
      DEALLOCATE(VK1, VK2)
C
C  All angular EVP
C
C --- Sectors info
      WSEC = RMAX / DBLE(NSEC)
      SR = WSEC / 2.D0
      ALLOCATE(RSEC(0:NSEC))
      DO i = 0, NSEC
        RSEC(i) = i * WSEC
      ENDDO
C --- Angular basis
      CALL ANGBAS(KSYM, L, M, LMAX, NBAS) 
      NSVD = NBAS * NDVR
      ALLOCATE (FSVD1(NBAS, NSVD, NSEC), FSVD2(NBAS, NSVD, NSEC))
      
      RETURN
      END SUBROUTINE

      SUBROUTINE STARKEF()

      DO i = 1, NSEC
        CALL SECTOR(i)
      ENDDO
      END SUBROUTINE
      END MODULE
