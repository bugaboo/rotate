      SUBROUTINE SYMDVR(NDVR,H,KSYM)
C=======================================================================
C  SYMMETRIZATION OF A DVR HAMILTONIAN
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(*)
C
      IF(MOD(NDVR,2).NE.0) STOP ' *** SYMDVR NDVR ERROR'
      IF(IABS(KSYM).NE.1)  STOP ' *** SYMDVR KSYM ERROR'
      S=DBLE(KSYM)
      NDVRS=NDVR/2
      ij=0
      DO j=1,NDVRS
        ijN=((NDVR-j)*(NDVR-j+1))/2
        DO i=1,j
          ij=ij+1
          ijN=ijN+1
          H(ij)=H(ij)+S*H(ijN)
        ENDDO
      ENDDO
      NDVR=NDVRS
C
      RETURN 
      END
C=======================================================================