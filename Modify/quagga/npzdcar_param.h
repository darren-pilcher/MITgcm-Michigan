c $Header: /u/gcmpack/MITgcm/pkg/npzdcar/npzdcar_param.h, 2007/12/12 
c
c     ==================================================================
c     HEADER npzdcar_param
c     ==================================================================
c
c     o Header file for the surface flux data. Used by the external
c       forcing package.
c
c     started: Christian Eckert eckert@mit.edu  30-Jun-1999
c
c     changed: Christian Eckert eckert@mit.edu  14-Jan-2000
c              - Restructured the original version in order to have a
c                better interface to the MITgcmUV.
c
c              Christian Eckert eckert@mit.edu  12-Feb-2000
c              - Changed some variables names (package prefix: npzdcar_)
c
c              Patrick Heimbach, heimbach@mit.edu  04-May-2000
c              - included npzdcar_iprec, npzdcar_yftype to enable easy
c                switch between 32bit/64 bit data format
c
c              Patrick Heimbach, heimbach@mit.edu  01-May-2001
c              - added obcs parameters
c
c     ==================================================================
c     HEADER npzdcar_param
c     ==================================================================
C     need to add other variables, only fice, so far (12/17/2007)

Cc     year in seconds
C      _RL     year2sec

c     Calendar data.
      _RL     npzdcar_repeatPeriod

C     npzdcar setup 
C     ice fraction
      integer ficestartdate1
      integer ficestartdate2
      _RL     ficestartdate
      _RL     ficeperiod
      _RL     ficeconst
      _RL     fice_npzdcarremo_intercept 
      _RL     fice_npzdcarremo_slope
      character*1 ficemask
      parameter(  ficemask = 's' )
C     iron atmosphere
c      integer featmosstartdate1
c      integer featmosstartdate2
c      _RL     featmosstartdate
c      _RL     featmosperiod
c      _RL     featmosconst
c      _RL     featmos_npzdcarremo_intercept 
c      _RL     featmos_npzdcarremo_slope
c      character*1 featmosmask
c      parameter(  featmosmask = 's' )
C     wind speed
      integer wspeedstartdate1
      integer wspeedstartdate2
      _RL     wspeedstartdate
      _RL     wspeedperiod
      _RL     wspeedconst
      _RL     wspeed_npzdcarremo_intercept 
      _RL     wspeed_npzdcarremo_slope
      character*1 wspeedmask
      parameter(  wspeedmask = 's' )
C     sur_par (PAR -- shortwave radiation)
      integer sur_parstartdate1
      integer sur_parstartdate2
      _RL     sur_parstartdate
      _RL     sur_parperiod
      _RL     sur_parconst
      _RL     sur_par_npzdcarremo_intercept 
      _RL     sur_par_npzdcarremo_slope
      character*1 sur_parmask
      parameter(  sur_parmask = 's' )
C     atmosp (atmospheric pressure)
      integer atmospstartdate1
      integer atmospstartdate2
      _RL     atmospstartdate
      _RL     atmospperiod
      _RL     atmospconst
      _RL     atmosp_npzdcarremo_intercept 
      _RL     atmosp_npzdcarremo_slope
      character*1 atmospmask
      parameter(  atmospmask = 's' )

C     Added by Darren Pilcher 10/2013

C     quagga (quagga mussel)
      integer quaggastartdate1
      integer quaggastartdate2
      _RL     quaggastartdate
      _RL     quaggaperiod
      _RL     quaggaconst
      _RL     quagga_npzdcarremo_intercept
      _RL     quagga_npzdcarremo_slope
      character*1 quaggamask
      parameter(  quaggamask = 's' )

c     File names.
      character*(128) ficefile
c      character*(128) featmosfile
      character*(128) wspeedfile
      character*(128) sur_parfile
      character*(128) atmospfile
      character*(128) quaggafile

C     useNpzdcarYearlyFields :: when set, automatically add extension
C                           _YEAR to input file names
C     twoDigitYear       :: when set, use 2-digit year extension YR
C                           instead of _YEAR for useNpzdcarYearlyFields
      logical useNpzdcarYearlyFields, NpzdcartwoDigitYear
      logical useNpzdcarCheckRange

      common /npzdcar_param_l/
     &                   useNpzdcarYearlyFields, NpzdcartwoDigitYear,
     &                   useNpzdcarCheckRange
      common /npzdcar_param_i/
     &                        ficestartdate1,    ficestartdate2,
c     &                        featmosstartdate1, featmosstartdate2,
     &                        wspeedstartdate1,  wspeedstartdate2,
     &                        sur_parstartdate1, sur_parstartdate2,
     &                        atmospstartdate1,  atmospstartdate2,
     &                        quaggastartdate1,  quaggastartdate2

      common /npzdcar_param_r/
C     &                          year2sec,            windstressmax,
C     &                          repeatPeriod,        fice_npzdcar_albedo,
C     &                          year2sec,
     &                    npzdcar_repeatPeriod,
     &                    ficeperiod,    ficestartdate,    ficeconst,
c     &                    featmosperiod, featmosstartdate, featmosconst,
     &                    wspeedperiod,  wspeedstartdate,  wspeedconst,
     &                    sur_parperiod, sur_parstartdate, sur_parconst,
     &                    atmospperiod,  atmospstartdate,  atmospconst,
     &                    quaggaperiod,  quaggastartdate,  quaggaconst

      common /npzdcar_param_trend_removal/
     &        fice_npzdcarremo_intercept,    fice_npzdcarremo_slope,
     &        wspeed_npzdcarremo_intercept,  wspeed_npzdcarremo_slope,
c     &        featmos_npzdcarremo_intercept, featmos_npzdcarremo_slope,
     &        sur_par_npzdcarremo_intercept, sur_par_npzdcarremo_slope,
     &        atmosp_npzdcarremo_intercept,  atmosp_npzdcarremo_slope,
     &        quagga_npzdcarremo_intercept,  quagga_npzdcarremo_slope

      common /npzdcar_param_c/
     &                        ficefile,
c     &                        featmosfile,
     &                        wspeedfile,
     &                        sur_parfile,
     &                        atmospfile,
     &                        quaggafile

c     file precision and field type

      common /npzdcar_param_type/ 
     &                        npzdcar_iprec,
     &                        npzdcar_yftype

      integer npzdcar_iprec
      character*(2) npzdcar_yftype

      _RL     npzdcar_inscal_fice
c      _RL     npzdcar_inscal_featmos
      _RL     npzdcar_inscal_wspeed
      _RL     npzdcar_inscal_sur_par
      _RL     npzdcar_inscal_atmosp 
      _RL     npzdcar_inscal_quagga     

      _RL     npzdcar_outscal_fice
c      _RL     npzdcar_outscal_featmos
      _RL     npzdcar_outscal_wspeed 
      _RL     npzdcar_outscal_sur_par
      _RL     npzdcar_outscal_atmosp 
      _RL     npzdcar_outscal_quagga

      common /npzdcar_param_scal/
     &        npzdcar_inscal_fice,    npzdcar_outscal_fice,
c     &        npzdcar_inscal_featmos, npzdcar_outscal_featmos,
     &        npzdcar_inscal_wspeed,  npzdcar_outscal_wspeed,
     &        npzdcar_inscal_sur_par, npzdcar_outscal_sur_par,
     &        npzdcar_inscal_atmosp,  npzdcar_outscal_atmosp, 
     &        npzdcar_inscal_quagga,  npzdcar_outscal_quagga

      character*(5) npzdcarexternalforcingversion
      character*(5) usenpzdcarscalendarversion
      parameter( npzdcarexternalforcingversion = '0.2.2' ) 
      parameter( usenpzdcarscalendarversion    = '0.2.0' )

      _RL npzdcar_one
      parameter( npzdcar_one  =  1.0d0)

C#ifndef USE_NPZDCAR_INTERPOLATION
c-- set dummy dimension 1
       integer MAX_LAT_INC
       parameter(MAX_LAT_INC = 1)
C#else
Cc for lat interpolation, arraysize currently set to 2176 max data values
C       integer MAX_LAT_INC
C       parameter(MAX_LAT_INC = 2176)
C      _RL fice_lon0, fice_lon_inc
C      _RL fice_lat0, fice_lat_inc(MAX_LAT_INC)
CC      INTEGER apressure_nlon,apressure_nlat
C
C      common /npzdcar_interpolation/
C     & fice_lon0, fice_lon_inc,
C     & fice_lat0, fice_lat_inc,
C     & fice_nlon, fice_nlat
C#endif
