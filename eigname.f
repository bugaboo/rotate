      FUNCTION EIGNAME(MODEL, KSYM, X0, LMAX, OMEGA) RESULT(RES)
C  Filename for storing eigenvectors and eigenvalues
C  Modes: C - eigenvectors, L - eigenvalues
      CHARACTER(30) :: RES, TMP
      INTEGER, INTENT(IN) :: MODEL, KSYM, LMAX
      REAL*8, INTENT(IN) :: X0, OMEGA
      WRITE(RES, 78) LMAX, KSYM, MODEL, INT(1.d1*X0), INT(1.d1*OMEGA)
      RETURN
 78   FORMAT('eigenval',I0,'_',I0,'_',I0,'_',I0,'_',I0)
      END
