      program main
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      NAMELIST /INF_MAIN/ KSYM, MODEL, X0, OMEGA
        
      CONTAINS 
        SUBROUTINE CONSOLE_INPUT(RMAX, KSYM, MODEL, X0, OMEGA)
        REAL*8, INTENT(OUT) :: RMAX, X0, OMEGA
        INTEGER, INTENT(OUT) :: MODEL, KSYM
        CHARACTER(LEN=20) :: STR
        
        CONSOLE_INPUT = .FALSE.
        IF (IARGC() .LT. 5) RETURN
        CALL GETARG(1, STR)
        READ(STR,*) RMAX
        CALL GETARG(2, STR)
        READ(STR,*) KSYM
        CALL GETARG(3, STR)
        READ(STR,*) MODEL
        CALL GETARG(4, STR)
        READ(STR,*) X0
        CALL GETARG(5, STR)
        READ(STR,*) OMEGA
        END SUBROUTINE
      end program
