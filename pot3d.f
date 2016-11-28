      FUNCTION POT3D(RAD,TET,PHI,X0)
      USE STARK

C=======================================================================
C  Potential energy V(r) for one-channel SPS EVP in the R problem
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C      COMMON /POT_C/MODEL
C
      R = DSQRT(RAD*RAD + 2*RAD*X0*DSIN(TET)*DCOS(PHI) + X0*X0)      
      SELECT CASE(MODEL)
C----------------------------
C
C----------------------------
      CASE(0)
        POT3D = -3.3392239215D0*DEXP(-R*R)
C----------------------------
C  model 1: Yukawa potential
C----------------------------
      CASE(1)
      R=RAD
      POT3D=DSIN(TET)
C----------------------------
C  model 2: Axial Yukawa potential
C----------------------------
      CASE(2)
      X=DSIN(TET)*DCOS(PHI)
      Y=DSIN(TET)*DSIN(PHI)
      Z=DCOS(TET)
      TX=.5D0
      TY=.5D0
      R=RAD*DSQRT(X**2/TX+Y**2/TY+Z**2)
      POT3D=-DEXP(-0.1D0*R)/R
C----------------------------
C  model 3: Axial Yukawa potential
C----------------------------
      CASE(3)
      X=DSIN(TET)*DCOS(PHI)
      Y=DSIN(TET)*DSIN(PHI)
      Z=DCOS(TET)
      TX=2.D0
      TY=2.D0
      R=RAD*DSQRT(X**2/TX+Y**2/TY+Z**2)
      POT3D=-DEXP(-0.1D0*R)/R
C----------------------------
C  model 4: Anizotropic Yukawa potential
C----------------------------
      CASE(4)
      X=DSIN(TET)*DCOS(PHI)
      Y=DSIN(TET)*DSIN(PHI)
      Z=DCOS(TET)
      TX=2.D0
      TY=0.5D0
      R=RAD*DSQRT(X**2/TX+Y**2/TY+Z**2)
      POT3D=-DEXP(-0.1D0*R)/R
C-----------------------------
C  model 5: Coulomb potential
C-----------------------------
      CASE(5)
      POT3D=-1.D0/RAD
C----------------------------
C  model 6: Anizotropic harmonic oscillator
C----------------------------
      CASE(6)
      X=0.9D0*RAD*DSIN(TET)*DCOS(PHI)
      Y=1.D0*RAD*DSIN(TET)*DSIN(PHI)
      Z=1.1D0*RAD*DCOS(TET)
      POT3D=0.5D0*(X**2+Y**2+Z**2-23.D0)
C      
C----------------------------
C  model 300: Angle test
C----------------------------
      CASE(300)
      X = DCOS(TET)
c      POT3D = DSQRT(1-X**2)*POLJ(1.D0,1.D0,0.D0,X)*DSIN(PHI)/
c     &	 DSQRT(4*DATAN(1.D0))
      POT3D = DSIN(PHI)*DSIN(TET)*DSQRT(.75D0)
C----------------------------
C  model 301: Angle test
C----------------------------
      CASE(301)
      POT3D = DCOS(TET)*DSQRT(3.D0/DATAN(1.D0))/4.D0
      
C----------------------------
C  model 302: Angle test
C----------------------------
      CASE(302)
      POT3D = DCOS(PHI)*DSQRT(2.D0)
C----------------------------
C  model 303: Angle test
C----------------------------
      CASE(303)
      POT3D = DEXP(DCOS(TET) + DSIN(2.D0*PHI))
C----------------------------
C  model 200: Rectangular potential
C----------------------------
      CASE(200)
      IF (RAD.LT.1.0000001D0) THEN
	POT3D = -1.125D2
      ELSE
	POT3D = 0.D0
      ENDIF
      END SELECT
      RETURN
      END
C=======================================================================
