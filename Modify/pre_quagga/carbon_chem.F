#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: CALC_PCO2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  pH and pCO2 solver for freshwater
cc  Val Bennington 01-2008
cc  all equations from CO2sys freshwater option 8
cc enter nutrients in mol/m3, T in Celsius
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE CALC_PCO2_APPROX_LAKE(
     I                       t,diclocal,ta,
     I                       k0local,k1local,k2local,
     I                       kwlocal,fugflocal,
     U                       phlocal,pCO2local)


      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "NPZDCAR.h"
c
c      ROUTINE ARGUMENTS
c      diclocal = total inorganic carbon (mol/m3)
c      ta = total alkalinity (mol eq/m3)
c      t = temp (C)
        _RL t, ta
        _RL pCO2local, diclocal, phlocal
        _RL k1local, k2local
        _RL fugflocal
        _RL kwlocal, k0local
        _RL kflocal, ftlocal
        _RL kslocal, stlocal
        
c  ta2 is alkalinity in umol/kg
c        LOCAL VARIABLES
        _RL phguess
        _RL cag
        _RL hguess
        _RL resid
        _RL denom
        _RL slope, deltaph
        _RL fCO2
        _RL calk
        _RL tc, OH        
        _RL lg10
C-------------------------------------------------------------------
c lake mean surface density is 1000 kg/m3
c make units mol/kg
c units converted back at end of program
        ta = ta*permil
        diclocal = diclocal*permil
cccccccccccccccccccccccccccccccccccccccccccccccccccc
        phguess = phlocal
        hguess = 10**(-phguess)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        OH = kwlocal/hguess
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  estimate carbonate alkalinitiy
cc    should be in umol/kg
        tc = diclocal
        denom = hguess*hguess+ k1local*hguess
     &          +k1local*k2local
        calk = tc*k1local*(hguess+2*k2local)/denom
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        tc = calk * (hguess*hguess + k1local*hguess
c     &       +k1local*k2local)
c     &       /(k1local*( hguess + 2*k2local))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        fCO2 = tc*hguess*hguess
     &        /(hguess*hguess + k1local*hguess
     &        + k1local*k2local)
     &        /k0local
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   use total carbon to guess at new pH:
        resid = ta - calk - OH + hguess
         lg10 = 2.30258509299
        slope = lg10 * (tc * k1local * hguess
     &          * (hguess * hguess + k1local * k2local + 
     &          + 4 * hguess * k2local)
     &          / denom / denom + OH + hguess)
        deltaph = resid/slope
        phlocal = phguess+deltaph
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        pCO2local = fCO2/fugflocal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
cc put back dic and alk in mol/m3
        diclocal = diclocal/permil
        ta = ta/permil
        return 
        end

c=================================================================
c *******************************************************************
c=================================================================
CStartOfInterFace
      SUBROUTINE CARBON_COEFFS_LAKE(
     I                   ttemp,
     I                   bi,bj,iMin,iMax,jMin,jMax)
C
C     /==========================================================\
C     | SUBROUTINE CARBON_COEFFS                                 |
C     | determine coefficients for surface carbon chemistry      |
C     | adapted from CO2_SYS:  SUBROUTINE CO2CALC                |
C     | val Bennington dec 2007                                  |
c     |                                                          |
C     \==========================================================/
C INPUT
C       ta  = total alkalinity (eq/m^3) 
C       pt  = inorganic phosphate (mol/^3) = 0
C       sit = inorganic silicate (mol/^3) = 0
C       t   = temperature (degrees C)
C       s   = salinity (PSU) = 0
C OUTPUT
C IMPORTANT: Some words about units - 
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c       all input concentrations (ta, pt, and st) should be 
c       given in mol/m^3; output arguments in mol/m^3.
C--------------------------------------------------------------------------

        IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "NPZDCAR.h"
C     == Routine arguments ==
C ttemp and stemp are local theta and salt arrays
C dont really need to pass T and S in, could use theta, salt in
C common block in DYNVARS.h, but this way keeps subroutine more
C general
        _RL  ttemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
        _RL  stemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
         INTEGER bi,bj,iMin,iMax,jMin,jMax
CEndOfInterface
C    LOCAL VARIABLES
       _RL  lnkw
       _RL  lnk0
       _RL  lnk1
       _RL  lnk2
       _RL  tk, t, logtk 
       _RL  P1atm, RT, Rgasconstant
       _RL  delta, B, ttt
cccccccccccccccccccccccccccccccccccccccccccccccccccc
        INTEGER i
        INTEGER j
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=imin,imax
         do j=jmin,jmax
          if (hFacC(i,j,1,bi,bj).gt.0.d0) then
           t = ttemp(i,j,1,bi,bj)
           tk = 273.15+t
           ttt = tk/100
           logtk = log(tk)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
           lnk0 = -60.2409 + 93.4517/ttt + 23.3585*log(ttt)
           ak0(i,j,bi,bj) = exp(lnk0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           lnkw = 148.9802 - 13847.26/tk - 23.6521 * logtk    
           akw(i,j,bi,bj) = exp(lnkw)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
           lnk1 = 290.9097 - 14554.21/tk - 45.0575 * logtk
           ak1(i,j,bi,bj) = exp(lnk1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
           lnk2 = 207.6548 - 11843.79/tk - 33.6485 * logtk
           ak2(i,j,bi,bj) = exp(lnk2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculate fugacity factor
           Rgasconstant = 83.1451
           P1atm = 1.01325
           RT = Rgasconstant*tk
           delta = 57.7 - 0.118*tk
           B = -1636.75 + 12.0408*tk - 0.0327957*tk*tk
           B = B + 3.16528*(10**-5)*tk*tk*tk
           fugf(i,j,bi,bj) = exp( (B+2*delta) *P1atm/RT )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculate ff - weiss and price correction for atmos
c   water vapor pressure and non-ideality of gags
           ff(i,j,bi,bj) = exp(-162.8301 + 218.2968/ttt +
     &         90.9241*log(ttt) - 1.47696*((ttt*ttt)) )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          else
            akf(i,j,bi,bj)=0.d0
            ft(i,j,bi,bj)=0.d0
            aks(i,j,bi,bj)=0.d0
            st(i,j,bi,bj)=0.d0
            ff(i,j,bi,bj)= 0.d0
            fugf(i,j,bi,bj)=0.d0
            ak0(i,j,bi,bj)= 0.d0
            ak1(i,j,bi,bj)= 0.d0
            ak2(i,j,bi,bj)= 0.d0
            akw(i,j,bi,bj)= 0.d0
          endif
         end do
        end do

        return
        end
           

c=================================================================
           
