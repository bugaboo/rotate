      SUBROUTINE SVDEVP(ISEC,NSVD,CEIG,CVEC)
      USE STARK
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION :: CEIG(*), CVEC(NSVD,*)
      ALLOCATABLE :: OVLP(:, :)

C --- Overlap matrix
      ALLOCATE (OVLP(NSVD, NSVD))
      CALL GET_OVLP(ISEC, NSVD, OVLP)
C --- Kinetic part
      ALLOCATE (CA(NSVD, NSVD))
      DO i = 0, NDVR - 1
        DO j = 0, NDVR - 1
          DO n = 1, NANG
            DO m = 1, NANG
              CA(i*NANG+n, j*NANG+m) = DCMPLX(TK(i,j) * 
     &                                OVLP(i*NANG+n, j*NANG+m))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL LAPEIGCG(1,N,CA,CEIG,CVEC)
      END SUBROUTINE
