      MODULE ANGSYM
      CONTAINS
      SUBROUTINE ANGBAS(KSYM, L, M, LMAX, NANG)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL :: SX, SY, SZ
      INTEGER, INTENT(OUT) :: NANG
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: L, M

      n = 1

C  KSYM=0 No symmetry          
      SELECT CASE(KSYM)
      CASE(0)
        NANG = (LMAX + 1)**2
        ALLOCATE(L(NANG), M(NANG))
        DO i = 0, LMAX
          DO k = -i, i
            L(n) = i
            M(n) = k
            n = n + 1
          ENDDO
        ENDDO
C  KSYM= 1 Inversion oz symmetry
      CASE(-1, 1)
        SZ = (KSYM .GT. 0)
        NANG = 0
        DO i = 0, LMAX
          DO j = -i, i
            IF (SYMZ(i,j,SZ)) THEN
              NANG = NANG + 1
            ENDIF
          ENDDO
        ENDDO
        ALLOCATE(L(NANG), M(NANG))
        DO i = 0, LMAX
          DO j = -i, i
            IF (SYMZ(i,j,SZ)) THEN
              L(n) = i
              M(n) = j
              n = n + 1
            ENDIF
          ENDDO
        ENDDO
C  KSYM=2 Spherical symmetry        
      CASE(2)
        NANG = 1
        ALLOCATE(L(1),M(1))
        L(1)=LAN
        M(1)=0
      END SELECT
      RETURN
      END SUBROUTINE
      
      FUNCTION SYMX(L, M, B)
        LOGICAL :: SYMX, B
        INTEGER :: L, M                    
        SYMX = MOD(M, 2) .EQ. 0
        IF (M.GE.0 .NEQV. B) THEN
          SYMX = (.NOT.SYMX)
        ENDIF
        RETURN
      END FUNCTION
      
      FUNCTION SYMY(L, M, B)
        LOGICAL :: SYMY, B
        INTEGER :: L, M
        SYMY = (M.GE.0)
        IF (.NOT.B) SYMY = (.NOT.SYMY)
        RETURN
      END FUNCTION
C      
      FUNCTION SYMZ(L, M, B)
        LOGICAL :: SYMZ, B
        INTEGER :: L, M
        SYMZ = (MOD(L-M,2).EQ.0)
        IF (.NOT.B) SYMZ = (.NOT.SYMZ)
        RETURN
      END FUNCTION
      END MODULE ANGSYM
