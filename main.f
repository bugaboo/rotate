      program main
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      NAMELIST /INF_MAIN/ CE0, F0, F1, NX0, OMEGA
      COMMON /PI_C/PI

C   Input information      
      PI=3.1415926535897932384626433832795D0
      OPEN(1,FILE='inf')
      READ(1,INF_MAIN)
      CLOSE(1)
       
C   Initialization of STARK       
      CONTAINS 
        SUBROUTINE CONSOLE_INPUT(CE0, F0, F1, NX0, OMEGA)
        REAL*8, INTENT(OUT) :: RMAX, X0, OMEGA, F0, F1
        INTEGER, INTENT(OUT) :: MODEL, KSYM, NX0
        COMPLEX*16, INTENT(OUT) :: CE0
        CHARACTER(LEN=20) :: STR
        
        CONSOLE_INPUT = .FALSE.
        IF (IARGC() .LT. 5) RETURN
        CALL GETARG(1, STR)
        READ(STR,*) CE0
        CALL GETARG(2, STR)
        READ(STR,*) F0
        CALL GETARG(3, STR)
        READ(STR,*) F1
        CALL GETARG(4, STR)
        READ(STR,*) NX0
        CALL GETARG(5, STR)
        READ(STR,*) OMEGA
        END SUBROUTINE
      end program
