      SUBROUTINE ANGEVP(L,M,NBAS,LMAX,NTET,NPHI,X0,OMEGA,RAD,EIG,VEC)
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      DIMENSION :: L(NBAS), M(NBAS), EIG(NBAS), VEC(NBAS,NBAS)
      ALLOCATABLE :: WK(:)
      
      CALL POTMAT(NTET, NPHI, RAD, VEC, NBAS, L, M, LMAX, X0)
      DO j = 1, NBAS
        VEC(j, j) = VEC(j,j)+DBLE(L(j)*(L(j)+1))/2.D0/RAD/RAD-OMEGA*M(j)
      ENDDO
      
      NWK = 20 * NBAS
      ALLOCATE(WK(NWK))
C      CALL LAPEIGRS (0, NBAS, V, E, WK)
      CALL DSYEV('V','U',NBAS,VEC,NBAS,EIG,WK,NWK,info)
      IF(info.NE.0)  WRITE(*,*) "ANGEVP DSYEV ERROR",info

      END SUBROUTINE