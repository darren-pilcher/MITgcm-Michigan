c
#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_NPZDCAR

c===============================================================================
C===============================================================================
CStartofinterface
      SUBROUTINE NPZDCAR_INIT_FIXED(myThid)                        
C     =============== Global data ==========================================
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS.h"
#include "NPZDCAR.h"
#include "GCHEM.h"

        INTEGER myThid
C============== Local variables ============================================
      INTEGER i,j,k,bi,bj
      INTEGER tmp
      INTEGER prec
      CHARACTER*(MAX_LEN_MBUF) fn
C     /--------------------------------------------------------------\
C     |  initialise common block biochemical parameters               |
C     \--------------------------------------------------------------/

c      k0=        4.d-2                   !atten coefficient water(m^-1)
c Lake Superior should have an attenuation coefficient about 0.1
c Guilford has attenuation coefficient of 0.15 
c    tried 0.08 as well
       k0=        13.d-2                   !atten coefficient L Sup water (1/m)
c      k0=        4.d-2
c      kc=        5.d-2                    !atten coefficient chl (m2 (mg chl)-1)
       kc=        149.d-4
      phygrow=   1.4D0 /(1.D0*86400.d0)   !max dino growth rate (s-1)
      phygrow2=  1.6D0 /(1.D0*86400.d0)   !max diatom growth rate (s-1)
c      phygrow = 6.D0/(86400.d0)
c      phygrow2= 8.D0/(86400.d0)
c      phygrow = 8.D0/(86400.d0)
c      phygrow2= 10.D0/(86400.d0)
      zoograze=  1.D0 /(1.D0*86400.d0)   !max zoo grazing rate (s-1)
      mort=      1.D0 /(15.D0*86400.d0)  !dino mortality rate (s-1)
      mort2=     1.D0 /(15.D0*86400.d0)  !diatom mortality rate (s-1)
      mortz=     1.D0 /(15.D0*86400.d0)  !zoo mortality rate (s-1)
c     mortz=     0.2 /86400.d0 *16.d0    !for quadratic
      ilim=      25.d0                   !dino light 1/2 saturation (W/m2)
      ilim2=     15.d0                   !diatom light 1/2 saturation (W/m2)
c orig      ilim2=     12.d0                   !diatom light 1/2 saturation (W/m2)
c  original numbers:
c      nlim=      0.05                    !dino po4 1/2 saturation (uM)
c      nlim2 =    0.1                     !diatom po4 1/2 saturation (uM)
c tfunc1 nlim = 0.04, nlim2 = 0.07
c hope to reduce magnitude of P drawdown
c want the seasonal cycle to be a little less controlled by P
c need to increase growth rates some to get growth
       nlim=      0.05                    !dino po4 1/2 saturation (uM)
       nlim2 =    0.1                   !diatom po4 1/2saturation(uM)
       slim2=     5.d0                    !diatom sio2 1/2 saturation (uM)
c      flim=      1d-5                    !dino fe 1/2 saturation (uM)
c      flim2=     5d-5                    !diatom fe 1/2 saturation (uM)
      zplat=     9.d-1                    !palatability of dino
      zplat2=     0.7
c      zplat2=    4.d-1                   !palatability of diatoms
c MIGHT NEED TO CHANGE THIS TO GET ZOO !
c      plim=      1.2d0 /16.d0            !zoo phy 1/2 saturation (uM)
      plim=      1.5d0 /16.d0            !zoo phy 1/2 saturation (uM)
      reminn=    1.D0 /(30.D0*86400.d0)  !remin rate of DOP/POP(s-1) 
      remins=    1.d0 /(360.D0*86400.d0) !remin rate of DOS/POS (s-1)
c      reminf=    1.D0 /(70.D0*86400.d0)  !remin rate of DOF/POF (s-1)
      reminc=    1.D0 /(30.D0*86400.d0)  !remin rate of DOP/POP (s-1)
c      dissc=     1.d0 /(360.D0*86400.d0) !dissolution rate PIC (s-1)
c      remint=     1.d0/(86400.d0*30*365.25)
c      reminr=     1.d0/(365.25*180*86400.d0)
      refrac=     0.05
      slabfrac=   0.80
      labfrac=   0.15
      fexc =     0.13                    ! additional 13% of NPP is excreted 
c
      gampn=     0.3d0                   !zoo dino assimilation coeff
      gampn2=    0.3d0                   !zoo diatom assimilation coeff
      donfracmn1 = 0.50d0                 !fraction of mortality to DOP/POP
      donfracmn2 = 0.30d0
      donfracmn = 0.90d0
      donfracgn= 0.99d0                   !fraction of grazing to DOP/POP
c      donfracmn1 = 0.5d0                 !fraction of mortality to DOP/POP
c      donfracmn2 = 0.2d0
c      donfracmn = 0.2d0
c      donfracgn= 0.5d0                   !fraction of grazing to DOP/POP
c      donfracmf1= 0.5d0                  !fraction of mortality to DOF/POF
c      donfracmf2= 0.2d0
c      donfracmf=0.2d0
c      donfracgf= 0.5d0                   !fraction of grazing to DOF/POF
      donfracmc1= 0.99d0                  !fraction of mortality to DOF/POF
      donfracmc2= 0.90d0
      donfracmc=0.90d0
      donfracgc= 0.99d0                   !fraction of grazing to DOF/POF
c      donfracmc1= 0.5d0                  !fraction of mortality to DOF/POF
c      donfracmc2= 0.2d0
c      donfracmc=0.2d0
c      donfracgc= 0.5d0                   !fraction of grazing to DOF/POF
      Rsip2=     25.d0                   !ratio Si to P in diatoms
c      Rfep=      1.0d-4                  !ratio Fe to P in dino
c      Rfep2=     1.0d-3                  !ratio Fe to P in diatoms
c orgi      Rcp=      120                      !ratio C to P in dino
c orig      Rcp2=      120                     !ratio C to P in diatoms
       Rcp=      200                      !ratio C to P in dino
       Rcp2=      200                     !ratio C to P in diatoms
       Rop=       170                     !ratio O to P for plankton
       Rnp=        16                     !ratio N to P in ocean             
c      alpfe=     0.025                   !solubility of Fe dust
c      scav=      0.4/(360.d0*86400.d0)   !iron chem scavenging rate (s-1)
c      freefemax = .3e-3                  ! max free fe
c      ligand_tot=1.d-3                   !total ligand (uM)
c      ligand_stab=1.d5                   !ligand stability rate ratio
c      rainratio=0.07                     !rain ratio
      o2crit= 4.0d-6                     !critical O2 concentration
c how much P is SRP vs. DOP?
      fracsrp = 1                        ! input files are SRP only
      fracdop = 3                        ! TP is about 30% SRP
c minimum populations
c      phymin=    1.d-5                   !seed population of phy (uM)
      phymin=    1.d-3                   !seed population of phy (uM)
      grazlim=   5.d-4                   !grazlim (uM)
      zoomin=    1.d-3                   !egg population of zoo (uM)
c      zoominfe=  zoomin*Rfep2
      zoominc=   zoomin*Rcp2
cQQ make sure wn_sink*dt<delz
c Chai and Urban indicated particulate sinking rate of 2.3 m/day
c Noel said 0.5 -2 reasonable
c Added diatom sinking DJP 06/2014
      wphy2_sink = 5.d-1/(86400.d0) !sinking rate of Phy2 (m/s)
      wn_sink = 5.d-1/(86400.d0) !sinking rate of POP (m/s)
      wc_sink = 5.d-1/(86400.d0) !sinking rate of POC (m/s)
c      wn_sink= 2900.d0/(360.d0*86400.d0) !sinking rate of POP (m/s)
      ws_sink = 5.d-1/(86400.d0) !sinking rate of POS (m/s)
c      wf_sink= 2900.d0/(360.d0*86400.d0) !sinking rate of POF (m/s)
c      wc_sink= 2900.d0/(360.d0*86400.d0) !sinking rate of POC (m/s)
c      wpic_sink= 4000.d0/(360.d0*86400.d0) !sinking rate of PIC (m/s)
c for chl
c assuming based on C:P which is different in lake --
c adjust low and high values accordingly 
c August 28, 2008 Val Bennington
c      chlcmax = 0.2193 
      chlpmax = 28
      chlpmin = 20   ! better results vs. EPA to use Chl:P then Chl:C

      chlcmax = 0.3173
c      chlcmin = 0.0756
      chlcmin = 0.0959
      istar = 70                           ! uEin/m2/d
c for temperature function of growth rate
c Try using a version of Moore's for Lake Superior
c See npzd_forcing.F
c function goes from 0.73 to 2.1 from 0-25 degrees
c growth rate is 1.1&=*max rate at 10 degrees
c tfunc1 has base of -3400, norm of 280.15
c tfunc2 had base of -2600, norm of 275.15
c tfunc1 and tfunc2 are for Moore's tempfunc
c       temp_coeff_base = -3400
c       temp_coeff_norm = 280.15
c (Eppley) for temperature function of growth rate
      thetam = 1.0521
      thetar = 1.06
      temp_coeff_base=-3400              ! unitless, Moore
      temp_coeff_norm= 280.15             ! tweaked by Darren 06-19-12
c      temp_coeff_base=1.0368              ! unitless, Eppley
c      temp_coeff_norm=1.4                ! normalization, -2:35 mean=1
c      temp_coeff_norm= 3.6296             ! by steph previously  
c   temp_coeff_norm set to 3.0000 by val for first test 6-18-07
c temperature function for grazing
c       temp_zoo_base=1.066
c        temp_zoo_norm=2.3
ccccccccccccccccccccccccccccccccc
c parameters for air-sea parts
ccccccccccccccccccccccccccccccccc
C set up coefficients for DIC chemistry
C define Schmidt no. coefficients for CO2
c      sca1 = 2073.1
c      sca2 = -125.62
c      sca3 =    3.6276
c      sca4 =   -0.043219
ccccc freshwater values:
      sca1 = 1911.1
      sca2 = -118.11
      sca3 =    3.4527
      sca4 =   -0.041320
C define Schmidt no. coefficients for O2
C based on Keeling et al [GBC, 12, 141, (1998)]
c      sox1 = 1638.0
c      sox2 = -81.83
c      sox3 =    1.483
c      sox4 =   -0.008004
cccccc freshwater values from Wanninkhof (1992)
       sox1 = 1800.6
       sox2 = -120.10
       sox3 =    3.7818
       sox4 =  - 0.047608
C coefficients for determining saturation O2
c these don't need to change for freshwater
c --does hold with such low salinity values
c 0<=T<=35.5, 0<=S<=40
      oA0=  2.00907
      oA1=  3.22014
      oA2=  4.05010
      oA3=  4.94457
      oA4= -2.56847E-1
      oA5=  3.88767
      oB0= -6.24523E-3
      oB1= -7.37614E-3
      oB2= -1.03410E-2
      oB3= -8.17083E-3
      oC0= -4.88682E-7
c  garcia and gordon (1992) limnol. oceanogr.

C set carbon mol/m3 <---> mol/kg conversion factor
      permil = 1.0 / 1000

#ifdef ALLOW_MNC
      IF ( useMNC ) THEN
        CALL NPZDCAR_MNC_INIT( myThid )
      ENDIF
#endif /* ALLOW_MNC */


        RETURN
        END
C============================================================================
#endif
#endif

