      FUNCTION EIGNAME(MODE, MODEL, KSYM, X0, LMAX, OMEGA) RESULT(RES)
C ---------------------------------------------------------------      
C  Filename for storing eigenvectors and eigenvalues
C  Modes: V - eigenvectors, E - eigenvalues
C  Filename format: eigen(vec,val)LMAX_KSYM_MODEL_10*X0_10*OMEGA
C ---------------------------------------------------------------
      CHARACTER(30) :: RES, TMP
      INTEGER, INTENT(IN) :: MODEL, KSYM, LMAX
      REAL*8, INTENT(IN) :: X0, OMEGA
      CHARACTER(2), INTENT(IN) :: MODE

      IF (MODE .EQ. 'V') THEN
        TMP = 'vec'
      ELSE IF (MODE .EQ. 'E') THEN
        TMP = 'val'
      ELSE
        STOP ' *** WRONG MODE FOR EIGNAME'
      END IF
      WRITE(RES, 78) TMP,LMAX,KSYM,MODEL,INT(10.D0*X0), INT(10.D0*OMEGA)
      RETURN
 78   FORMAT('eigen',A3,I0,'_',I0,'_',I0,'_',I0,'_',I0)
      END
