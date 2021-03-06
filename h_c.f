      PROGRAM H_ANG

      USE ANGSYM
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(PI=3.141592653589793238462643D0)
      COMMON /POT_C/MODEL

      CHARACTER(LEN = 30) :: EIGNAME
      LOGICAL :: CONS_INPUT
      ALLOCATABLE :: L(:), M(:), CV(:,:), E(:),CZ(:),CW(:),WK(:),IK(:)
      ALLOCATABLE :: IFAIL(:)
      X0 = 2.D0
      NEIG = 100
      OMEGA = 1.D0
C      EIGNAME = 'eigenvalues'
      
      IF (.NOT.CONSOLE_INPUT(RMAX, NR, LMAX, KSYM, NTET, NPHI, MODEL, X0
     &  ,OMEGA)) STOP 'Arguments stop'
      CALL ANGBAS(KSYM, L, M, LMAX, NANG) 
      NEIG = MIN(NEIG, NANG)
C      DO i = 1, NANG
C        PRINT *, L(i), M(i), KSYM
C      ENDDO
      NW = NANG * NANG + 2
      ALLOCATE (CV(NANG, NANG), E(NANG), CZ(NEIG), CW(NW), WK(NANG*7))
      ALLOCATE (IK(5 * NANG), IFAIL(NANG))
      write(*,*) EIGNAME(MODEL, KSYM, X0, LMAX, OMEGA), X0
      open(1, FILE=EIGNAME(MODEL, KSYM, X0, LMAX, OMEGA))
      TNR = RMAX / DBLE(NR)
      DO i = 1, NR
        IF (MOD(i, NR / 1) .EQ. 0) PRINT *, i
        CV = 0.D0
        R = 0.5D0 + DBLE(i) * TNR
        CALL CPOTMAT(NTET,NPHI,R,CV,NANG,L,M,LMAX,X0)
        if (i .eq. 1) then
          do j = 1, nang
            do k = 1, nang
              write(*, '(g17.10,1x)', ADVANCE = 'NO') 
     &            DREAL(cv(j, k))              
              write(*,50)L(j), M(j), L(k), M(k)
            enddo
          enddo
        endif
        DO j = 1, NANG
          CV(j, j) = CV(j,j) + DBLE(L(j)*(L(j)+1))/2.D0/R/R - OMEGA*M(j)
        ENDDO
        CALL ZHEEVX('N', 'I', 'U', NANG, CV, NANG, 0, 1, 1, NEIG, 0.D0,
     &      NEIG, E, CZ, 1, CW, NW, WK, IK, IFAIL, INFO)
        IF (INFO .GT. 0) STOP 'ZHEEVX INFO>0'
        WRITE(1, 77, ADVANCE = 'NO') R
        DO j = 1, NEIG - 1 
          WRITE(1,77, ADVANCE = 'NO') E(j)
        ENDDO
        WRITE(1, 77) E(NEIG)
      ENDDO
        
 77   FORMAT((E19.12,1X))
 50   FORMAT(I0,' ',I0,' ',I0,' ',I0)

      CONTAINS 
        FUNCTION CONSOLE_INPUT(RMAX, NR, LMAX, KSYM, NTET, NPHI, MODEL,
     &    X0, OMEGA)
        REAL*8, INTENT(OUT) :: RMAX, X0, OMEGA
        INTEGER, INTENT(OUT) :: NR, LMAX, KSYM
        LOGICAL CONSOLE_INPUT
        CHARACTER(LEN=20) :: STR
        
        CONSOLE_INPUT = .FALSE.
        IF (IARGC() .LT. 9) RETURN
        CALL GETARG(1, STR)
        READ(STR,*) RMAX
        CALL GETARG(2, STR)
        READ(STR,*) NR
        CALL GETARG(3, STR)
        READ(STR,*) LMAX
        CALL GETARG(4, STR)
        READ(STR,*) KSYM
        CALL GETARG(5, STR)
        READ(STR,*) NTET
        CALL GETARG(6, STR)
        READ(STR,*) NPHI
        CALL GETARG(7, STR)
        READ(STR,*) MODEL
        CALL GETARG(8, STR)
        READ(STR,*) X0
        CALL GETARG(9, STR)
        READ(STR,*) OMEGA
        CONSOLE_INPUT = .TRUE.
        RETURN
        END FUNCTION
      END PROGRAM

