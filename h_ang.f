      PROGRAM H_ANG

      USE ANGSYM
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      PARAMETER(PI=3.141592653589793238462643D0)
      COMMON /POT_C/MODEL

      CHARACTER(LEN = 30) :: EIGNAME
      LOGICAL :: CONS_INPUT
      ALLOCATABLE :: L(:), M(:), V(:,:), E(:),Z(:),WI(:),WK(:)
      ALLOCATABLE :: IFAIL(:)
      X0 = 2.D0
      NEIG = 3
      OMEGA = 1.D0
      
      IF (.NOT.CONSOLE_INPUT(RMAX, NR, LMAX, KSYM, NTET, NPHI, MODEL, X0
     &  ,OMEGA)) STOP 'Arguments stop'

      CALL ANGBAS(KSYM, L, M, LMAX, NBAS) 
      NEIG = MIN(NEIG, NBAS)

      IW = (NBAS + 3) * NBAS
      ALLOCATE (V(NBAS, NBAS), E(NBAS), WK(IW))
      open(1, FILE=EIGNAME(MODEL, KSYM, X0, LMAX, OMEGA))
      TNR = RMAX / DBLE(NR)
      DO i = 1, NR
C        IF (MOD(i, NR / 10 + 1) .EQ. 0) PRINT *, 'Progress:', DBLE(i)/NR
        V = 0.D0
        R = DBLE(i) * TNR
        CALL POTMAT(NTET,NPHI,R,V,NBAS,L,M,LMAX,X0)
        DO j = 1, NBAS
          V(j, j) = V(j,j) + DBLE(L(j)*(L(j)+1))/2.D0/R/R - OMEGA*M(j)
        ENDDO

        CALL LAPEIGRS (0, NBAS, V, E, WK) 
c        CALL DSYEV ('N', 'U', NBAS, V, NBAS, E, WK, IW, INFO)    

        WRITE(1, 77, ADVANCE = 'NO') R
        DO j = 1, NEIG - 1 
          WRITE(1,77, ADVANCE = 'NO') E(j)
        ENDDO
        WRITE(1, 77) E(NEIG)
        WRITE(*, '(I0,1X,I0,1X)', ADVANCE = 'NO') NTET, NPHI
        DO j = 1, NEIG - 1 
          WRITE(*,77, ADVANCE = 'NO') E(j)
        ENDDO
        WRITE(*, 77) E(NEIG)
      ENDDO

C      CALL ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RMAX,E,V)
      DO i = 1, NEIG
C        PRINT *, E(i)
      ENDDO
      
      DEALLOCATE(L, M, V, E, WK) 

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

