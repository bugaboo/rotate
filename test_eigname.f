      PROGRAM TEST_EIGNAME
      character(30) :: EIGNAME
      character(2) :: a
      a = 'V'
      IF (eigname(a, 1, 2, 3.D0, 4, 5.D2) .NE. 'eigenvec4_2_1_30_5000') 
     &    STOP '*** EIGNAME TEST FAILED' 
      end 
