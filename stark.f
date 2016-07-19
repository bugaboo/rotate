      MODULE STARK
      
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)

      DIMENSION AKT(*),CTMD(*)
      ALLOCATABLE XX(:),WX(:),TX(:,:),DVRX(:),XSI(:),
     &            CEIGX(:),CVECX(:,:),
     &            XE(:,:),WE(:,:),PILE(:,:),PIRE(:,:),DVRE(:,:)
      ALLOCATABLE CRMATL(:,:),CRMATR(:,:),CRAS(:),CVK(:),CWK(:,:)
      ALLOCATABLE CXSIM(:,:),CDMAT(:,:),CNMAT(:,:),
     &            CMAT1(:,:),CMAT2(:,:),CDPSI(:,:),CPSIAS(:),
     &            CFAS(:),CVN(:),CWN(:,:),CVK1(:),CVK2(:)
      NAMELIST /INF_STARK_OLD/KPOT,MZ,NCH,NDVRX,ES,NDVRE,ETAS,NSEC,MSEC,
     &                    KROT,FAC1ST,FACMAX,ITRMAX,EPS
      COMMON /PI_C/PI
      COMMON /POT_STARK/KPOT /DIM_STARK/NDVRX,NDVRE,NDVRSE,NCH,NSVD
     &       /ES_STARK/ES,SCALE /EFM_STARK/CE,CF,AMM,CALF
     &       /COF_STARK/COF(4) /KEY_STARK/KRUN,KWF
      SAVE MZ,AMZ,ETAS,ETAC,NSEC,MSEC,KROT,FAC1ST,FACMAX,ITRMAX,EPS,IWF
      SAVE XX,WX,TX,DVRX,XSI,CEIGX,CVECX,XE,WE,PILE,PIRE,DVRE

      NAMELIST /INF_STARK/ITRMAX,EPS,NSEC,RMAX,NDVR,LMAX,NTET,NPHI
      ALLOCATABLE :: L(:), M(:), PIL(:), PIR(:), VK(:)
      ALLOCATABLE :: RSEC(:), XR(:), WR(:), TR(:,:), DVRR(:)
      ALLOCATABLE :: ANGEIGVAL(:,:), ANGEIGVEC(:,:)
      SAVE L, M, RSEC, XR, WR, TR, DVRR, ANGEIGVAL, ANGEIGVEC, PIL, PIR
      SAVE RMAX, SR, WSEC

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

      AMZ=DBLE(MZ)
      AMM=0.5D0*AMZ*AMZ
      COF(1)=0.25D0*(1.D0-AMZ*AMZ)
      ETAC=DBLE(NSEC)*ETAS
      NDVRSX=(NDVRX*(NDVRX+1))/2
      NDVRSE=(NDVRE*(NDVRE+1))/2
      NSVD=NCH*NDVRE


C
C  DVR arrays
C
C --- R
      NT = NDVR + 1
      ALLOCATE(XR(NDVR),WR(NDVR),TR(NT,NT),DVRR(NDVR*(NT)/2))
      CALL LEGPOM(NDVR)
      CALL DVRLEG(NDVR,XR,WR,TR,1,0,DVRR,NT)
C --- Left (PIL) and right (PIR) end DVR basis amplitudes
      ALLOCATE (PIL(NDVR), PIR(NDVR), VK(NDVR))
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
          PIL(i)=PIL(i)+T(n,i)*VK1(n)
          PIR(i)=PIR(i)+T(n,i)*VK2(n)
        ENDDO
      ENDDO
      DEALLOCATE(VK)
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
C --- ANGEVPs
      ALLOCATE(ANGEIGVAL(NSEC,NBAS), ANGEIGVEC(NSEC,NBAS))
      DO i = 1, NSEC
        DO j = 1, NDVR
          RAD = RSEC(i) - (NDVR - j) * WSEC / DBLE(NDVR) 
          CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RAD,
     &              ANGEIGVAL(i),ANGEIGVEC(i))
        ENDDO
      ENDDO
C --- XSI
      ALLOCATE(XX(NDVRX),WX(NDVRX),TX(NDVRX+1,NDVRX+1),DVRX(NDVRSX),
     &         XSI(NDVRX),CEIGX(NCH),CVECX(NDVRX,NCH))
      CALL LAGPOM(AMZ,NDVRX)
      CALL DVRLAG(AMZ,NDVRX,XX,WX,TX,1,DVRX,NDVRX+1)
      SCALE=1.D0/DSQRT(-2.D0*ES)
      DO i=1,NDVRX
        XSI(i)=SCALE*XX(i)
      ENDDO
