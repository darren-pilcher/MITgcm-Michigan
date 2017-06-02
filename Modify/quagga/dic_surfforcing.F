#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: DIC_SURFFORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_SURFFORCING( PTR_DIC, PTR_ALK,
     O           GDC, 
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

C !DESCRIPTION:
C  Calculate the carbon air-sea flux terms              
C  following external_forcing_dic.F (OCMIP run) from Mick   
C  modified for lake chemistry by Val Bennington 1-18-08
C  model does not trace SIL or FE for lake chemistry         

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS.h"
#include "NPZDCAR.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
c  PTR_DIC              :: DIC tracer field
      INTEGER myIter, myThid
      _RL myTime
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
c      _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
c      _RL  PTR_SIL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

      INTEGER iMin,iMax,jMin,jMax, bi, bj

C !OUTPUT PARAMETERS: ===================================================
c GDC                   :: tendency due to air-sea exchange
      _RL  GDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef ALLOW_PTRACERS

C !LOCAL VARIABLES: ====================================================
       INTEGER I,J, kLev
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
      _RL SchmidtNoDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pCO2sat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pisvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C local variables for carbon chem
      _RL surfdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c      _RL surfphos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c      _RL surfsi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C
C variables for reading CO2 input files
      _RL tmp
      _RL aWght, bWght
      _RL co2atmos(1000)
      INTEGER it, ntim
CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      kLev=1

c default
       if (gchem_int1.eq.0) then
C PRE-INDUSTRIAL STEADY STATE pCO2 = 278.0 ppmv
         DO j=1,sNy
          DO i=1,sNx
c          AtmospCO2(i,j,bi,bj)=278.0d-6
             AtmospCO2(i,j,bi,bj)=385.0d-6
          ENDDO
         ENDDO
       endif

c set atmos co2 from file
        if (gchem_int1.eq.1) then
          DO j=1,sNy
           DO i=1,sNx
            AtmospCO2(i,j,bi,bj)=gchem_rl1
           ENDDO
          ENDDO
        endif

c read from ascii file co2atmos.dat
        if (gchem_int1.eq.2) then
c            print*,'myIter',myIter
            if (myIter.eq.nIter0+1) then
            OPEN(28,FILE='co2atmos.dat',STATUS='old')
            do it=1,gchem_int2
             READ(28,*) co2atmos(it)
c             print*,'co2atmos',co2atmos(it)
            enddo
          endif
          ntim=int((myIter-gchem_int3)/gchem_int4)+1
          aWght=0.5+float(myIter-gchem_int3)/float(gchem_int4)-
     &                      float(ntim-1)
          if (aWght.gt.1.d0) then
            ntim=ntim+1
            aWght=aWght-1.d0
          endif
          bWght=1.d0-aWght
          tmp=co2atmos(ntim)*bWght+co2atmos(ntim+1)*aWght
          DO j=1,sNy
           DO i=1,sNx
             if (j.lt.2.d0) then
                if (i.lt.2.d0) then
c                  print*,'atmospco2 ',tmp
                endif
             endif
            AtmospCO2(i,j,bi,bj)=tmp
           ENDDO
          ENDDO
        endif

c        print*,'AtmospCO2(10,10)',AtmospCO2(10,10,bi,bj)

C =================================================================
C determine inorganic carbon chem coefficients, convert to mol/m3
        DO j=1,sNy
         DO i=1,sNx
             surfdic(i,j) = PTR_DIC(i,j,klev)*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
             surfalk(i,j) = PTR_ALK(i,j,klev)*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
c             surfphos(i,j)  = PTR_PO4(i,j,klev)*1e-3
c     &                          * maskC(i,j,kLev,bi,bj)
c             surfsi(i,j)   = PTR_SIL(i,j,klev) *1e-3
c     &                          * maskC(i,j,kLev,bi,bj)
          ENDDO
         ENDDO

         CALL CARBON_COEFFS_LAKE(
     I                       theta,
     I                       bi,bj,iMin,iMax,jMin,jMax)
C====================================================================

c pCO2 solver...
C$TAF LOOP = parallel
       DO j=1,sNy
C$TAF LOOP = parallel
        DO i=1,sNx

          IF(maskC(i,j,kLev,bi,bj) .NE. 0.)THEN
          CALL CALC_PCO2_APPROX_LAKE(
     I        theta(i,j,kLev,bi,bj),surfdic(i,j),surfalk(i,j),
     I        ak0(i,j,bi,bj),ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I        akw(i,j,bi,bj),fugf(i,j,bi,bj),
     U        pH(i,j,bi,bj),pCO2(i,j,bi,bj) )
              

c            CALL CALC_PCO2_APPROX_LAKE(
c     I        theta(i,j,kLev,bi,bj),salt(i,j,kLev,bi,bj),
c     I        surfdic(i,j),surfphos(i,j),
c     I        surfsi(i,j),surfalk(i,j),
c     I        ak1(i,j,bi,bj),ak2(i,j,bi,bj),
c     I        ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj),
c     I        aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
c     I        aksi(i,j,bi,bj),akf(i,j,bi,bj),ff(i,j,bi,bj),
c     I        bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj),
c     U        pH(i,j,bi,bj),pCO2(i,j,bi,bj) )
          ELSE
             pCO2(i,j,bi,bj)=0. _d 0
             pH(i,j,bi,bj)=0. _d 0
          END IF
        ENDDO
       ENDDO
c             print*,'pH(10,10)',pH(10,10,1,1)
c             print*,'dic(10,10)',surfdic(10,10)
c             print*,'alk(10,10)',surfalk(10,10)
c             print*,'pCO2(10,10)',pCO2(10,10,1,1)
c program outputs pCO2 in atmospheres!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       DO j=1,sNy
        DO i=1,sNx
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            IF (maskC(i,j,kLev,bi,bj).NE.0.) THEN
C calculate SCHMIDT NO. for CO2
              SchmidtNoDIC(i,j) = 
     &            sca1 
     &          + sca2 * theta(i,j,kLev,bi,bj)
     &          + sca3 * theta(i,j,kLev,bi,bj)*theta(i,j,kLev,bi,bj)  
     &          + sca4 * theta(i,j,kLev,bi,bj)*theta(i,j,kLev,bi,bj) 
     &                *theta(i,j,kLev,bi,bj)

C Determine surface flux (FDIC)
C first correct pCO2at for surface atmos pressure
              pCO2sat(i,j) = 
     &          AtmosP(i,j,bi,bj)*AtmospCO2(i,j,bi,bj)
c              print*,'pCO2sat(i,j)',pCO2sat(i,j)
c find exchange coefficient
c  account for schmidt number and and varible piston velocity
c  divide by 3.6*10^5 converts cm/hr to meters/sec
              pisvel(i,j)=0.260*wspeed(i,j,bi,bj)**2/3.6d5
              Kwexch(i,j) =
     &             pisvel(i,j)
     &             / sqrt(SchmidtNoDIC(i,j)/600.0)
ccccccccccccccc  normalize to 600 for freshwater
c     &             pisvel(i,j)
c     &             / sqrt(SchmidtNoDIC(i,j)/660.0)
c OR use a constant  coeff
c             Kwexch(i,j) = 5e-5
c ice influence
              Kwexch(i,j)  =(1.d0-Fice(i,j,bi,bj))*Kwexch(i,j)


C Calculate flux in terms of DIC units using K0, solubility
C Flux = Vp * ([CO2sat] - [CO2])
C CO2sat = K0*pCO2atmos*P/P0
C Converting pCO2 to [CO2] using ff, as in CALC_PCO2
c put pCO2 (lake) into atm from uatm:
c              print*,'ak0 ',ak0(i,j,bi,bj)
c              print*,'pCO2sat ',pCO2sat(i,j)
c              print*,'ff ',ff(i,j,bi,bj)
c              print*,'fugf',fugf(i,j,bi,bj)
c              print*,'pCO2 ',pCO2(i,j,bi,bj)
c              print*,'kwexch ',Kwexch(i,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    corrected Nov. 2010 by Val Bennington
c    previously, atmospheric correction for water vapor pressure
c    was being applied to ocean, not atmos pCO2
c    here corrected, with fugacity factor applied to ocean and 
c    ff applied to atmosphere
              FluxCO2(i,j,bi,bj) = 
     &         maskC(i,j,kLev,bi,bj)*Kwexch(i,j)*( 
     &         ff(i,j,bi,bj)*pCO2sat(i,j) - 
     &         ak0(i,j,bi,bj)*fugf(i,j,bi,bj)*pCO2(i,j,bi,bj) 
     &         )
c              print*,'FluxCO2',FluxCO2(i,j,bi,bj) 
            ELSE
               FluxCO2(i,j,bi,bj) = 0.
            ENDIF
c               print*,'FluxCO2 ',FluxCO2(i,j,bi,bj)
C convert flux (mol kg-1 m s-1) to (mol m-2 s-1)
            FluxCO2(i,j,bi,bj) = FluxCO2(i,j,bi,bj)/permil
c convert flux (mol m-2 s-1) to (mmol m-2 s-1)
            FluxCO2(i,j,bi,bj) = FluxCO2(i,j,bi,bj)*1e3

            IF (maskC(i,j,kLev,bi,bj).NE.0.) THEN
c calculate virtual flux
c EminusPforV = dS/dt*(1/Sglob)
C NOTE: Be very careful with signs here!
C Positive EminusPforV => loss of water to atmos and increase
C in salinity. Thus, also increase in other surface tracers
C (i.e. positive virtual flux into surface layer)
C ...so here, VirtualFLux = dC/dt!
c              VirtualFlux(i,j,bi,bj)=
c     &            gsm_DIC*surfaceForcingS(i,j,bi,bj)/gsm_s
c OR
c let virtual flux be zero
              VirtualFlux(i,j,bi,bj)=0.d0
c
            ELSE
              VirtualFlux(i,j,bi,bj)=0. _d 0
            ENDIF
          ENDDO
         ENDDO
c         print*,'CO2(10,10) ',FluxCO2(10,10,1,1)
c         print*,'pCO2sat(10,10)',pCO2sat(10,10)
c         print*,'pCO2(10,10)',pCO2(10,10,1,1)
c         print*,'Kw(10,10)',Kwexch(10,10)
c         print*,'ak0(10,10)',ak0(10,10,1,1)
c         print*,'ff(10,10)',ff(10,10,1,1)
C update tendency      
         DO j=1,sNy
          DO i=1,sNx
           GDC(i,j)= hFacC(i,j,kLev,bi,bj)*recip_drF(kLev)*(
     &                    FluxCO2(i,j,bi,bj) + VirtualFlux(i,j,bi,bj)
     &                                              )
          ENDDO
         ENDDO

#endif
        RETURN
        END
