C***********************************************************************
      PROGRAM MAIN
C
C  SPS 3D expansion
C
C***********************************************************************
      USE ANGSYM
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(PI=3.141592653589793238462643D0)
      PARAMETER(LANMAX=28)
      ALLOCATABLE CK(:),CE(:),CPHI(:,:),CVEC(:,:),CS(:,:), CEL(:) 
      ALLOCATABLE CZER(:,:)
      ALLOCATABLE L(:),M(:),PIR(:)
      NAMELIST /INF3D/MODEL,RADA,NDVR,LMAX,LAN,MAN,KSYM,NTET,NPHI,
     & 		KEYA,AMIN,AMAX,NUMA
      COMMON /POT_C/MODEL
      CHARACTER(LEN=30) :: SPSNAME, FILENAME, TMPNAME, SYMNAME
C
C  Input parameters
C
      OPEN(1,file='inf')
      READ(1,INF3D)
      CLOSE(1)
      IF(LAN.GT.LANMAX) STOP ' *** LANMAX ERROR'
      
C
C  Calculating SPS
C
      CALL ANGBAS(KSYM,L,M,LMAX,NANG,LAN,MAN)
      NSPS = 0
      DO i= 1, NANG
        NSPS = NSPS + 2 * NDVR + L(i)
      ENDDO
      ALLOCATE(CK(NSPS),CE(NSPS),CVEC(NDVR*NANG,NSPS))
      ALLOCATE(PIR(NDVR),CPHI(NANG,NSPS),CS(NANG,NANG),CEL(NANG))
      CALL SPS3D(LMAX,NANG,L,M,RADA,NDVR,NTET,NPHI,PIR,CK,
     &	    CE,CVEC,CPHI,NB,NA,NOI)       
      
      IF(NB.GT.0) WRITE(*,70) DREAL(CE(NB))
      
C  Saving eigenvalues and eigenvectors

      IF (KSYM .NE. 2) THEN
        FILENAME = SYMNAME('L',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN)
        PRINT *, FILENAME
        OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &          ACCESS = 'direct', RECL = NSPS * 16)
        WRITE(1, REC = 1) CK
        CLOSE(1)
        FILENAME = SYMNAME('C',MODEL,KSYM,NDVR,LMAX,INT(RADA),MAN)
        OPEN(UNIT = 1, FILE=FILENAME, FORM = 'unformatted', 
     &          ACCESS = 'direct', RECL = NSPS * NANG * 16)
        WRITE(1, REC = 1) CPHI
        CLOSE(1)
      ENDIF

C  Scattering matrix and length      
      
      CAK=0.2D0
      CALL SSUM3D(L,RADA,NSPS,NANG,CK,CPHI,CAK,CS)
      CALL CDETS(CS,CEL,CDET,NANG)
      PRINT *, "Sum det", CDET, AIMAG(CDLOG(CDET))/2.D0
      CALL SPRO(RADA,NSPS,CK,CAK,CDETP)
      CDETP = CDETP*CDEXP(-(0.D0,2.D0)*CAK*RADA*(NANG-1))
      PRINT *, "Prod det", CDETP, AIMAG(CDLOG(CDETP))/2.D0
      PRINT *, "Scattering length", SLENGTH(RADA,NSPS,NANG,CK)
      LTMP = UCHECK(CS, NANG)
      
C  Printing eigenvalues

      WRITE(SPSNAME,"(A7,I0,A1,I0,A1,I2.2,A1,I0,'_'I0)")
     & "spseig_",KSYM,'_',NDVR,'_',LMAX,'_',INT(RADA),MAN
      OPEN(1,FILE=SPSNAME)
      DO n=1,NSPS
        write(1,77) CK(n),CE(n)
      ENDDO
      CLOSE(1)

      DEALLOCATE(PIR,CK,CE,CVEC,CPHI,CS,CEL)
      
 70   FORMAT(' ground state energy = ',E19.12)
 77   FORMAT(4(E19.12,1X))
 
      CONTAINS 
        FUNCTION UCHECK(CS, NANG)
          LOGICAL UCHECK
          DIMENSION CS(NANG, NANG)

          UCHECK = .FALSE.
          DIFMAX = 0.D0
          DO i = 1, NANG
            DO j = 1, NANG
              CTMP = 0.D0
              DO k = 1, NANG
                CTMP = CTMP + CS(k,i) * DCONJG(CS(k,j))
              ENDDO
              IF (i.EQ.j) CTMP = CTMP - 1.D0
              DIFMAX = MAX(DIFMAX, CDABS(CTMP))
              IF (CDABS(CTMP) .GT. 1.D-6) THEN
                PRINT *, "NONUNITARY S MATRIX", i, j, CTMP
              ENDIF
            ENDDO
          ENDDO
          print *, "max diff: ", DIFMAX
          UCHECK = .TRUE.
          RETURN 
        END FUNCTION
      END PROGRAM
