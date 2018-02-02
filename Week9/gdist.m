function [srkm]  = gdist(lat1,lon1,lat2,lon2 )
%
%   
%         geodesic range
%              J.Clynch 9-98 from FTN routines
%                            which are based on NOAA M10 program
%    
%                           specificed by lat, lon (gd-lat)
%C
%C        GEODESIC INVERSE COMPUTATIONS
%C        INPUT LAT LON ,  ENDS OF LINE,
%C
%C        FROM NOAA routines in M10, VERSION 3
%C
%C      INPUT:
%C
%C             lat1    LATITUDE OF BEGINING IN DEG
%C             lon1    LONGITUDE "    "       " 
%C             lat2    LATITUDE OF END OF LINE IN DEG
%C             lon2    LONGITUDE "    "      "
%C
%C     OUTPUT:
%C
%C             srkm    DISTANCE IN KM ALONG GEODESIC
%C

%         Defaults to WGS 84
%                     Change via global variables if needed
%
%
%      ==============================================================
%                       ++ GLOBAL VARIABLES  ++
%
%       A_EARTH     if this and F_EARTH are present, defines ellipsoid
%       F_EARTH     "   "    "  A_EARTH  "    "         "      "     
%
%       FLAG_GDIST  0 converged, 1 error
%       ITER_GDIST  number of iterations used
%
     global   A_EARTH F_EARTH FLAG_GDIST ITER_GDIST
%
%    ================================================================
%   
%C
%C   -----------------------------------------------------------------
%C  
%C     HELMERT RAINSFORD INVERSE PROBLEM 
%C     DO NOT USE FOR MERIDIONAL ARCS AND BE CAREFUL ON EQUATOR
%C
%C  ------------------------------------------------------------------
%C
%C     AZIMUTHS FROM NORTH AND LONG POSITIVE EAST 
%C
%C -------------------------------------------------------------------
%C
%C
%C
%C  -------------------------------
      TOL    = 5.0D-15;
      flag   =  0;
%C
%C------------------------------------------------------------------
%C
      PI        =   pi;
      true      = 1;
      false     = 0;
%
      FLAT1     =  lat1;
      FLAT2     =  lat2;
      FLON1     =  lon1;
      FLON2     =  lon2;
%C
%C    FILL WITH WGS84 VALUES INITIALLY ON GEODETIC/GEOPHYSICAL
%C
%
%
%    -------------------------------------------------------
%
%
      RE      =   6378.137D00 ;
      FETH    =   0.335281066474D-2;
%
%              use global values if exist for ellipsoid
%
      tstell1 =   exist('A_EARTH');
      tstell2 =   exist('F_EARTH');
      if tstell1==1 &  tstell2==1, 
         szAE = size(A_EARTH);
         tstell3=max(szAE);
         if tstell3 > 0, 
            RE     = A_EARTH;
            FETH   = F_EARTH;
         end;
      end;
%
      EETH    =   8.18191908426D-2;
      ESETH   =   6.69437999013D-3;
%C
%C
      PIGDES    =   PI;
      FGDES     =   FETH;
      AGDES     =   1000.0D0*RE;
%C
%C
      PI       =   PIGDES;
      A        =   AGDES;
      F        =   FGDES;
%C
      TWOPI    =   2.0D0*PI;
      DTR      =   PI/180.0D0;
%C
      P1       =   DTR*FLAT1;
      E1       =   DTR*FLON1;
      P2       =   DTR*FLAT2;
      E2       =   DTR*FLON2;
%C
      ESQ      =   2.0D0*F-F*F;
%C  
      F0=(1.0D0-F);
      B=A*F0;
      EPSQ=ESQ/(1.0D0-ESQ);
      F2=F*F;     
      F3=F*F2;    
      F4=F*F3;   
%C
%C     TEST THE LONGITUDE DIFFERENCE
%C

      S=E2-E1;
      if abs(S) < TOL , E2=E2+TOL; end
%C
%C     THE LONGITUDE DIFFERENCE 
%C
      DLON=E2-E1;   
      AB=DLON;      
      KOUNT=0;    
%C
%C     THE REDUCED LATITUDES    
%C
      U1=F0*sin(P1)/cos(P1);
      U2=F0*sin(P2)/cos(P2);
%C
      U1=atan(U1);
      U2=atan(U2);
%C
      SU1=sin(U1);    
      CU1=cos(U1);    
%C
      SU2=sin(U2);
      CU2=cos(U2);
%C
%C     COUNTER FOR THE ITERATION OPERATION
%C
%
%1     CONTINUE

      done = false;

      while (~done)

%C
      KOUNT=KOUNT+1;     



%C
      CLON=cos(AB);   
      SLON=sin(AB);  
%C
      CSIG=SU1*SU2+CU1*CU2*CLON;
      SSIG=sqrt((SLON*CU2)^2+(SU2*CU1-SU1*CU2*CLON)^2);
%C
      SIG=atan2(SSIG,CSIG);
      SINALF=CU1*CU2*SLON/SSIG;
%C
      W =(1.0D0-SINALF*SINALF);
      T4=W*W;   
      T6=W*T4;   
%C
%C     THE COEFFICIENTS OF TYPE A      
%C
      AO=F-F2*(1.0D0+F+F2)*W/4.0D0+3.0D0*F3*(1.0D0+ ...
         9.0D0*F/4.0D0)*T4/16.0D0-25.0D0*F4*T6/128.0D0;
      A2=F2*(1.0D0+F+F2)*W/4.0D0-F3*(1.0D0+9.0D0*F/4.0D0)*T4/4.0D0+ ...
         75.0D0*F4*T6/256.0D0;
      A4=F3*(1.0D0+9.0D0*F/4.0D0)*T4/32.0D0-15.0D0*F4*T6/256.0D0;
      A6=5.0D0*F4*T6/768.0D0;
%C
%C     THE MULTIPLE ANGLE FUNCTIONS    
%C
      QO=0.0D0;
      if W>TOL , QO=-2.0D0*SU1*SU2/W; end
%C
      Q2=CSIG+QO;
      Q4=2.0D0*Q2*Q2-1.0D0;
      Q6=Q2*(4.0D0*Q2*Q2-3.0D0);
      R2=2.0D0*SSIG*CSIG;
      R3=SSIG*(3.0D0-4.0D0*SSIG*SSIG);
%C
%C     THE LONGITUDE DIFFERENCE 
%C
      S=SINALF*(AO*SIG+A2*SSIG*Q2+A4*R2*Q4+A6*R3*Q6);
      XZ=DLON+S;
%C
      XY=abs(XZ-AB);
      AB=DLON+S;
%C

%      if(XY.LT.0.5D-13) GOTO 4
%      if(KOUNT.LE.40)   GOTO 1

      if XY<0.5d-13 , done=true; end
      if KOUNT>30 , 
          done=true; 
          flag= 1;
      end

      end;

%C
%C     THE COEFFICIENTS OF TYPE B      
%C
%   4     CONTINUE
%C
      Z=EPSQ*W;
%C
      BO=1.0D0+Z*(1.0D0/4.0D0+Z*(-3.0D0/64.0D0+Z*(5.0D0/256.0D0- ...
         Z*175.0D0/16384.0D0)));
      B2=Z*(-1.0D0/4.0D0+Z*(1.0D0/16.0D0+Z*(-15.0D0/512.0D0+  ...
         Z*35.0D0/2048.0D0)));
      B4=Z*Z*(-1.0D0/128.0D0+Z*(3.0D0/512.0D0-Z*35.0D0/8192.0D0));
      B6=Z*Z*Z*(-1.0D0/1536.0D0+Z*5.0D0/6144.0D0);
%C
%C     THE DISTANCE IN METERS   
%C
      S=B*(BO*SIG+B2*SSIG*Q2+B4*R2*Q4+B6*R3*Q6);
%C
      SKM  =  0.001D0*S;
%
      srkm  =  SKM;
%C
%
%
      ITER_GDIST = KOUNT;
      FLAG_GDIST = flag;
%
%C
      return;

