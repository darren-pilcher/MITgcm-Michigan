C $Header: /u/gcmpack/MITgcm/model/src/packages_readparms.F,v 1.45 2006/02/07 11:47:48 mlosch Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CStartOfInterface
CBOP
C     !ROUTINE: PACKAGES_READPARMS
C     !INTERFACE:
      SUBROUTINE PACKAGES_READPARMS(myThid)

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE PACKAGES_READPARMS                             
C     | o Read runtime package configuration parameters           
C     *==========================================================*
C     |   Packages can have internal runtime configuration        
C     |   parameters. A package provides a routine                
C     |   ${PKGNAME}_READ_PARMS to read these parameters. In      
C     |   general this routine should read parameters from a file 
C     |   called data.${pkgname}.                                 
C     |   This routine (S/R PACKAGES_READPARMS) calls per-package 
C     |   parameter reading routines.                             
C     *==========================================================*
C     \ev

C     !CALLING SEQUENCE:
C     PACKAGES_READPARMS
C       |
C       |-- MNC_READPARMS + MNC_INIT + MNC_CW_INIT
C       |
C       |-- RUNCLOCK_READPARMS
C       |
C       |-- GMREDI_READPARMS
C       |
C       |-- KPP_READPARMS
C       |
C       |-- PP81_READPARMS
C       |
C       |-- MY82_READPARMS
C       |
C       |-- GGL90_READPARMS
C       |
C       |-- OPPS_READPARMS
C       |
C       |-- SBO_READPARMS
C       |
C       |-- SEAICE_READPARMS
C       |
C       |-- SHELFICE_READPARMS
C       |
C       |-- SHAP_FILT_READPARMS
C       |
C       |-- ZONAL_FILT_READPARMS
C       |
C       |-- OBCS_READPARMS
C       |
C       |-- OFFLINE_READPARMS
C       |
C       |-- BULKF_READPARMS
C       |
C       |-- THSICE_READPARMS
C       |
C       |-- GCHEM_READPARMS
C       |
C       |-- OPTIM_READPARMS
C       |
C       |-- CTRL_READPARMS
C       |
C       |-- COST_READPARMS
C       |
C       |-- GRDCHK_READPARMS
C       |
C       |-- CAL_READPARMS
C       |
C       |-- EXF_READPARMS
C       |
C       |-- EBM_READPARMS
C       |
C       |-- ECCO_READPARMS
C       |
C       |-- CPL_READPARMS
C       |
C       |-- FIZHI_READPARMS
C       |
C       |-- CPL_READPARMS
C       |
C       |-- DIAGNOSTICS_READPARMS


C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     == Routine arguments ==
      INTEGER myThid
CEOP

#ifdef ALLOW_MNC
C     Read run-time flags and initialize look-up tables for MNC
      IF (useMNC) THEN
        CALL MNC_READPARMS( myThid )
        CALL MNC_INIT(myThid)
        CALL MNC_CW_INIT(sNx,sNy,OLx,OLy,nSx,nSy,nPx,nPy,
     &       Nr,myThid)
CEH3       IF ( mnc_echo_gvtypes ) THEN
CEH3       CALL MNC_CW_DUMP( myThid )
CEH3       ENDIF
      ENDIF
#endif

#ifdef ALLOW_RUNCLOCK
C--   Initialize RUNCLOCK parameters
      IF (useRunClock) CALL RUNCLOCK_READPARMS( myThid )
#endif

#ifdef ALLOW_GMREDI
C--   Initialize GM/Redi parameters
      IF (useGMRedi) CALL GMREDI_READPARMS( myThid )
#endif

#ifdef ALLOW_KPP
C--   Initialize KPP parameters
      IF (useKPP) CALL KPP_READPARMS( myThid )
#endif

#ifdef ALLOW_PP81
C--   Initialize PP81 parameters
      IF (usePP81) CALL PP81_READPARMS( myThid )
#endif

#ifdef ALLOW_MY82
C--   Initialize MY82 parameters
      IF (useMY82) CALL MY82_READPARMS( myThid )
#endif

#ifdef ALLOW_GGL90
C--   Initialize GGL90 parameters
      IF (useGGL90) CALL GGL90_READPARMS( myThid )
#endif

#ifdef ALLOW_OPPS
C--   Initialize OPPS parameters
      IF (useOPPS) CALL OPPS_READPARMS( myThid )
#endif

#ifdef ALLOW_SBO
C--   Initialize SBO parameters
      IF (useSBO) CALL SBO_READPARMS( myThid )
#endif

#ifdef ALLOW_SEAICE
C--   Initialize SEAICE parameters
      IF (useSEAICE) CALL SEAICE_READPARMS( myThid )
#endif

#ifdef ALLOW_SHELFICE
C--   Initialize SHELFICE parameters
      IF (useShelfIce) CALL SHELFICE_READPARMS( myThid )
#endif

#ifdef ALLOW_SHAP_FILT 
      IF (useSHAP_FILT) CALL SHAP_FILT_READPARMS( myThid )
#endif

#ifdef ALLOW_ZONAL_FILT 
      IF (useZONAL_FILT) CALL ZONAL_FILT_READPARMS( myThid )
#endif

#ifdef ALLOW_OBCS
      IF (useOBCS) CALL OBCS_READPARMS( myThid )
#endif

#ifdef ALLOW_OFFLINE
      call OFFLINE_READPARMS( myThid )
#endif

#ifdef ALLOW_BULK_FORCE
C -- initialize Bulkf parameters
      IF (useBulkForce) CALL BULKF_READPARMS( myThid )
#endif

#ifdef ALLOW_THSICE
C -- initialize Therm.Sea-Ice parameters
      IF (useThSIce) CALL THSICE_READPARMS( myThid )
#endif

#ifdef ALLOW_PTRACERS
      IF (usePTRACERS) CALL PTRACERS_READPARMS( myThid )
#endif

#ifdef ALLOW_GCHEM
      IF (useGCHEM) CALL GCHEM_READPARMS( myThid )
#endif

#ifdef ALLOW_RBCS
      IF (useRBCS) CALL RBCS_READPARMS( myThid )
#endif

#ifdef ALLOW_CTRL
C--   Initialise the optim. parameters.
      call optim_readparms( myThid )
C--   Initialise the control parameters
      call ctrl_readparms( myThid )
#endif

#ifdef ALLOW_COST
C--   Initialise the cost parameters
      call cost_readparms( mythid )
#endif

#ifdef ALLOW_GRDCHK
C--   Initialise gradient check parameters
      if (useGrdchk) call grdchk_readparms( myThid )
#endif

#ifdef ALLOW_CAL
C--   Initialise the calendar package.
ceh3 should have an IF ( useCAL ) THEN
      call cal_readparms( myThid )
#endif

C--   Custom routine to set forcing fields.
#ifdef ALLOW_EXF
ceh3 needs an IF ( useEXF ) THEN
      call exf_readparms  ( mythid )
#endif

C--- added start
#ifdef ALLOW_LAKEICE ! how open lakeice
C--   Initialize LAKEICE parameters
      IF (useLAKEICE) CALL LAKEICE_READPARMS( myThid )
#endif
C--- added start
c #ifdef ALLOW_RIVERS ! how open allow rivers
C--   Initialize RIVERS parameters
c      IF (useRIVERS) CALL RIVERS_READPARMS( myThid )
c       print*,'should have read river params by now'
c #endif

#if( defined (ALLOW_NPZDCAR) && defined (ALLOW_GCHEM) ) 
      CALL NPZDCAR_READPARMS( myThid )
#endif
C--- added end

#ifdef ALLOW_EBM
      IF (useEBM) call ebm_readparms( mythid )
#endif

#ifdef ALLOW_ECCO
c--   Initialize ECCO package parameters
      call ecco_readparms ( myThid )
#endif

#ifdef COMPONENT_MODULE
C--   set Coupling parameters
      IF ( useCoupler ) CALL CPL_READPARMS( myThid )
#endif /* COMPONENT_MODULE */

CAMM
#ifdef ALLOW_FIZHI
C--   set model date and time and other fizhi parameters
      IF ( usefizhi ) call fizhi_readparms(myThid )
#endif /* ALLOW_FIZHI */
CAMM

#ifdef ALLOW_DIAGNOSTICS
C--   set diagnostics choices
      IF ( useDiagnostics ) CALL DIAGNOSTICS_READPARMS( myThid )
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_MATRIX
C--    Initialize matrix package parameters
       IF (useMATRIX) CALL MATRIX_READPARMS ( myThid )
#endif

      RETURN
      END
