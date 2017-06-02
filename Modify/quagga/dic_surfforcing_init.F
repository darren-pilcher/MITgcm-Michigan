#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: DIC_SURFFORCING_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_SURFFORCING_INIT(
     I          myThid)

C !DESCRIPTION:
C  Calculate first guess of pH                            

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
#include "NPZDCAR_LOAD.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER  myThid

#ifdef ALLOW_PTRACERS

C !LOCAL VARIABLES: ====================================================
       INTEGER I,J, kLev, it
       INTEGER intime0,intime1
       _RL otime
       _RL aWght,bWght,rdt
       INTEGER nForcingPeriods,Imytm,Ifprd,Ifcyc,Iftm
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
C local variables for carbon chem
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      _RL surfalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c      _RL surfphos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c      _RL surfsi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      kLev=1

ccccccccccccccccccccccccccccccccccccccccc
C =================================================================
      DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)

        jMin=1
        jMax=sNy
        iMin=1
        iMax=sNx

C determine inorganic carbon chem coefficients
        DO j=1,sNy
         DO i=1,sNx

cQQQQ check ptracer numbers
c for lake chemistry, alkalinity is input into pCO2 calculator as umol/kg
c conversion from mol/kg to umol/kg done in carbon_chem.F
c here convert to mol/m3 from mmol/m3 --- surfalk
c for pH solver, all need to be in mol/m3 
c rest is done in carbon_chem.F
             surfalk(i,j) = PTRACER(i,j,klev,bi,bj,14)*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
c             surfphos(i,j)  = PTRACER(i,j,klev,bi,bj,1)*1e-3
c     &                          * maskC(i,j,kLev,bi,bj)
             surfdic(i,j) = PTRACER(i,j,klev,bi,bj,9)*1e-3
     &                          * maskC(i,j,kLev,bi,bj)
c             surfsi(i,j) = PTRACER(i,j,klev,bi,bj,7)*1e-3
c     &                          * maskC(i,j,kLev,bi,bj)
          ENDDO
         ENDDO

         CALL CARBON_COEFFS_LAKE(
     I                       theta,
     I                       bi,bj,iMin,iMax,jMin,jMax)
C====================================================================

c set number of iterations for [H+] solvers
C set guess of pH for first step here

        DO j=1,sNy
         DO i=1,sNx
                  pH(i,j,bi,bj) = 8.0
          ENDDO
         ENDDO
c         print*,'QQ: pCO2 approximation method'
c first approxmation
C$TAF LOOP = parallel
        DO j=1,sNy
C$TAF LOOP = parallel
        DO i=1,sNx
          IF(HFacC(i,j,kLev,bi,bj) .gt.0.d0)THEN
c         IF(maskC(i,j,kLev,bi,bj) .NE. 0.)THEN
C$TAF init dic_surf = static, 10
           do it=1,100
c          do it=1,10
C$TAF STORE pH(i,j,bi,bj), PTR_CO2(i,j,kLev)           = dic_surf
C$TAF STORE surfalk(i,j), surfphos(i,j), surfsi(i,j)   = dic_surf
            CALL CALC_PCO2_APPROX_LAKE(
     I         theta(i,j,kLev,bi,bj),surfdic(i,j),surfalk(i,j),
     I         ak0(i,j,bi,bj),ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I         akw(i,j,bi,bj),fugf(i,j,bi,bj),
     U         pH(i,j,bi,bj),pCO2(i,j,bi,bj) )
            enddo
         ENDIF
        ENDDO
        ENDDO

        ENDDO
        ENDDO
c        print*,'QQ first guess pH', pH(10,10,1,1), theta(10,10,1,1,1),
c     &         salt(10,10,1,1,1),
c     &        surfdic(10,10),
c     &        surfalk(10,10),
c     &        pCO2(10,10,1,1)
#endif
        RETURN
        END
