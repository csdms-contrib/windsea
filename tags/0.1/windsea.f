      PROGRAM WINDSEA
C-----------------------------------------------------------------------
C This subroutine computes the deep water significant wave height and
C period at each point under a hurricane, as determined by empirical
C functions from CERC, U.S. Army Corps of Engineers.
C Input Variables:
C       Umax = maximum gradient windspeed 10 m above water (m/s)
C       R    = radius of maximum wind (km)
C       deltap = pressure difference between eye and ambient (mm Hg)
C       Vsubf = forward speed of hurricane (m/s)
C       im = maximum number of nodes in x direction (postive east)
C       jm = maximum number of nodes in y direction (positive north)
C       dx = space step in x direction (m)
C       dy = space step in y direction (m)
C       xo = x location of eye (m)
C       yo = y location of eye (m)
C       sdir = storm direction (degrees counterclockwise from east)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      Integer m,n,mm,nn,i,j,im,jm,i3,i4
      Parameter (m=12,n=7)
      Real x1a(m),x2a(n),ya(m,n),x1,x2,y2a(m,n),y,dx,dy,xo,yo,sdir,pdir
      Real alpha,Hzero,H(100,100),Ts(100,100),R,deltap,Vsubf,Usubr,g
      Real Tso,umax,pi,rad
C-----------------------------------------------------------------------
C     The following data statements and the spline subroutines
C     allow a representation of Fig. 3-43 in the CIRC Manual.
C     x1a is the x axis,
C     x2a is the y axis, and ya are the isoline values of the top half
C     of that figure.
C-----------------------------------------------------------------------
      Data x1a/-3,-2,-1,0,1,2,3,4,5,6,7,8/
      Data x2a/0,1,2,3,4,5,6/
      Data ya/.27,.43,.62,.48,1.0,.84,.72,.60,.51,.41,.30,.25,
     #        .27,.41,.59,.80,.87,.77,.68,.57,.49,.40,.30,.25,
     #        .25,.38,.50,.65,.70,.68,.59,.52,.45,.35,.28,.24,
     #        .23,.33,.40,.49,.55,.55,.52,.46,.40,.32,.28,.22,
     #        .21,.25,.30,.38,.44,.46,.43,.40,.35,.28,.23,.20,
     #        .20,.23,.25,.32,.35,.37,.35,.31,.28,.25,.20,.20,
     #        .20,.20,.22,.25,.26,.26,.26,.25,.20,.20,.20,.20/
C-----------------------------------------------------------------------
C     The following data statements are set some constant parameters
C-----------------------------------------------------------------------
      Data alpha,g,pi/1.0,9.8,3.1459/
      MM = 12
      NN = 7
C-----------------------------------------------------------------------
C     Open the usual files
C-----------------------------------------------------------------------
      OPEN(10, FILE= './windsea.in')
      OPEN(16, FILE= './HSIG.DATA')
      OPEN(17, FILE= './TSIG.DATA')
      OPEN(18, FILE= './WAVEDIR.DATA')
C-----------------------------------------------------------------------
C     Read input data
C-----------------------------------------------------------------------
      READ (10, *) UMAX, R, DELTAP, VSUBF, IM, JM, DX, DY, XO, YO, SDIR
C-----------------------------------------------------------------------
C     Calculate the deepwater significant wave height and period
C     under a hurricane at the point of maximum wind
C-----------------------------------------------------------------------
      R = R/1000.0
      USUBR = 0.865*UMAX + 0.5*VSUBF
      HZERO = 5.03*EXP(R*DELTAP/4700.0)*(1.0+0.29*ALPHA*VSUBF/SQRT(USUBR
     #   ))
      TSO = 8.6*EXP(R*DELTAP/9400.0)*(1.0+0.145*ALPHA*VSUBF/SQRT(USUBR))
      R = R*1000.0
      SDIR = SDIR*PI/180.0
C-----------------------------------------------------------------------
C     Now compute the remainder of the wave field. First precompute
C     the auxillary 2nd derivitive table for the splines, then
C     sweep the grid, redefining each x,y point in r/R space. The
C     storm is at xo,yo.
C-----------------------------------------------------------------------
      CALL SPLIE2 (X1A, X2A, YA, MM, NN, Y2A)
      DO J = 1, JM
         DO I = 1, IM
            X1 = I*DX
            X2 = J*DY
C-----------------------------------------------------------------------
C     In preparation for using Shore Protection Manual, Fig. 3-43 
C     (Anonymous,1988),locate the x,y point with respect to the
C     direction of storm
C-----------------------------------------------------------------------
            IF (X2 .GT. YO) THEN
               IF (X1 .GT. XO) THEN
                  PDIR = ATAN((X2-YO)/(X1-XO))
               ELSE IF (X1 .LT. XO) THEN
                  PDIR = PI - ATAN((X2-YO)/ABS(X1-XO))
               ELSE
                  PDIR = PI/2.0
               ENDIF
            ELSE IF (X2 .LT. YO) THEN
               IF (X1 .GT. XO) THEN
                  PDIR = 2.0*PI - ATAN(ABS(X2-YO)/(X1-XO))
               ELSE IF (X1 .LT. XO) THEN
                  PDIR = ATAN(ABS(X2-YO)/ABS(X1-XO)) + PI
               ELSE
                  PDIR = 1.5*PI
               ENDIF
            ELSE
               IF (X1 .GT. XO) THEN
                  PDIR = 0.0
               ELSE
                  PDIR = PI
               ENDIF
            ENDIF
            RAD = SQRT((X1-XO)**2+(X2-YO)**2)/R
            PDIR = PDIR - SDIR + PI/2.0
            X1 = RAD*COS(PDIR)
            X2 = ABS(RAD*SIN(PDIR))
            IF (X1.LE.8.0 .AND. X1.GE.(-3.0) .AND. X2.GE.0.0 .AND. X2
     #         .LE.7.0) THEN
C-----------------------------------------------------------------------
C     Compute the relative wave height
C-----------------------------------------------------------------------
               CALL SPLIN2 (X1A, X2A, YA, Y2A, MM, NN, X1, X2, Y)
               H(I,J) = HZERO*Y
               TS(I,J) = 12.1*SQRT(H(I,J)/G)
            ELSE
               H(I,J) = HZERO*0.2
               TS(I,J) = 12.1*SQRT(H(I,J)/G)
            ENDIF
         END DO
      END DO
C-----------------------------------------------------------------------
C     Write the results to files
C-----------------------------------------------------------------------
      I3 = 1
      I4 = 7
    3 CONTINUE
      IF (I3 .LE. IM) THEN
         I4 = MIN0(IM,I4)
         DO 4 J = 1, JM
            WRITE (16, 100) (H(I,J), I = I3, I4)
            WRITE (17, 100) (TS(I,J), I = I3, I4)
  100 FORMAT(7(1X,1PE9.2))
    4    CONTINUE
         I3 = I4 + 1
         I4 = I4 + 7
         GO TO 3
      ENDIF
      STOP 
      END
C-----------------------------------------------------------------------
C     The following four subroutines perform a bicubic spline
C     interpolation of the CIRC Fig. 3-43 isoline data. Given an
C     x,y point, they return an interpolated relative sig. wave height.
C     They are copied from Numerical Recipes (Press et al., 1986).
C-----------------------------------------------------------------------
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DO 3 J = 1, M
         DO 1 K = 1, N
            YTMP(K) = YA(J,K)
    1    CONTINUE
         CALL SPLINE (X2A, YTMP, N, 1.E30, 1.E30, Y2TMP)
         DO 2 K = 1, N
            Y2A(J,K) = Y2TMP(K)
    2    CONTINUE
    3 CONTINUE
      RETURN 
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),YYTMP(
     *NN)
      DO 2 J = 1, M
         DO 1 K = 1, N
            YTMP(K) = YA(J,K)
            Y2TMP(K) = Y2A(J,K)
    1    CONTINUE
         CALL SPLINT (X2A, YTMP, Y2TMP, N, X2, YYTMP(J))
    2 CONTINUE
      CALL SPLINE (X1A, YYTMP, M, 1.E30, 1.E30, Y2TMP)
      CALL SPLINT (X1A, YYTMP, Y2TMP, M, X1, Y)
      RETURN 
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1 .GT. .99E30) THEN
         Y2(1) = 0.
         U(1) = 0.
      ELSE
         Y2(1) = -0.5
         U(1) = (3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 1 I = 2, N - 1
         SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
         P = SIG*Y2(I-1) + 2.
         Y2(I) = (SIG-1.)/P
         U(I) = (6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I
     #      -1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    1 CONTINUE
      IF (YPN .GT. .99E30) THEN
         QN = 0.
         UN = 0.
      ELSE
         QN = 0.5
         UN = (3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 2 K = N - 1, 1, -1
         Y2(K) = Y2(K)*Y2(K+1) + U(K)
    2 CONTINUE
      RETURN 
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO = 1
      KHI = N
    1 CONTINUE
      IF (KHI - KLO .GT. 1) THEN
         K = (KHI+KLO)/2
         IF (XA(K) .GT. X) THEN
            KHI = K
         ELSE
            KLO = K
         ENDIF
         GO TO 1
      ENDIF
      H = XA(KHI) - XA(KLO)
      IF (H .EQ. 0.) THEN
         PAUSE 'Bad XA input.'
      ENDIF
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*
     #   H**2/6.
      RETURN 
      END
