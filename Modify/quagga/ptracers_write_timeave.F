C $Header: /u/gcmpack/MITgcm/pkg/ptracers/ptracers_write_timeave.F,v 1.7 2005/09/29 12:19:52 edhill Exp $
C $Name:  $

#include "PTRACERS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 1
C     !ROUTINE: PTRACERS_WRITE_TIMEAVE

C     !INTERFACE:
      SUBROUTINE PTRACERS_WRITE_TIMEAVE(myTime, myIter, myThid)

C     !DESCRIPTION:
C     At the end of average period, write the time-average
C     state-variables on file; then reset for next period

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS.h"
#include "PTRACERS_STATV.h"
C     Functions
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      INTEGER  IO_ERRCOUNT
      EXTERNAL IO_ERRCOUNT

C     !INPUT PARAMETERS:
C     myThid  ::  Thread number for this instance of the routine.
C     myIter  ::  Iteration number
C     myTime  ::  Current time of simulation ( s )
      INTEGER myThid
      INTEGER myIter
      _RL     myTime
CEOP

#if ( defined ALLOW_PTRACERS && defined ALLOW_TIMEAVE )

C     !LOCAL VARIABLES:
C     suff    ::  Hold suffix part of a filename
C     TimeAve ::  total time over average
      INTEGER prevPrec
      CHARACTER*(MAX_LEN_FNAM) suff
      _RL DDTT
      INTEGER bi,bj,k, iTracer
      INTEGER beginIOErrCount
      INTEGER endIOErrCount
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_MNC
      INTEGER ii,jj, ilnb
      CHARACTER*(MAX_LEN_FNAM) name
C     Functions
      INTEGER ILNBLNK
      CHARACTER*(1) pf
#endif /*  ALLOW_MNC  */

C     Final Time Averages and Dump Files if needed
      IF ( DIFFERENT_MULTIPLE( PTRACERS_taveFreq, myTime, deltaTClock )
     &   ) THEN
        
        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
C           Normalize by integrated time
            DO iTracer=1,PTRACERS_numInUse
C              CALL TIMEAVE_NORMALIZ(
C     &             ptracerFluxtave(1-Olx,1-Oly,  1,1,iTracer),
C     &             ptracer_full, 1,bi,bj,myThid)
              CALL TIMEAVE_NORMALIZ(
     &             ptracertave    (1-Olx,1-Oly,1,1,1,iTracer),
     &             ptracer_half,Nr,bi,bj,myThid)
            ENDDO
          ENDDO
        ENDDO

C       Write to files
        _BARRIER
        _BEGIN_MASTER( myThid )

#ifdef ALLOW_MNC
        IF ( PTRACERS_timeave_mnc ) THEN
          IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
            pf(1:1) = 'D'
          ELSE
            pf(1:1) = 'R'
          ENDIF
          CALL MNC_CW_SET_UDIM('ptr_tave', -1, myThid)
          CALL MNC_CW_RL_W_S('D','ptr_tave',0,0,'T',myTime,myThid)
          CALL MNC_CW_SET_UDIM('ptr_tave', 0, myThid)
          CALL MNC_CW_I_W_S('I','ptr_tave',0,0,'iter',myIter,myThid)

C          CALL MNC_CW_SET_UDIM('ptr_flux_tave', -1, myThid)
C          CALL MNC_CW_RL_W_S('D','ptr_flux_tave',0,0,'T',myTime,
C     &         myThid)
C          CALL MNC_CW_SET_UDIM('ptr_flux_tave', 0, myThid)
C          CALL MNC_CW_I_W_S('I','ptr_flux_tave',0,0,'iter',myIter,
C     &         myThid)

          DO ii = 1,PTRACERS_numInUse
            CALL MNC_CW_RL_W(pf,'ptr_tave',0,0,PTRACERS_names(ii),
     &           ptracertave(1-OLx,1-OLy,1,1,1,ii),myThid)
            DO jj = 1,MAX_LEN_FNAM
              name(jj:jj) = ' '
            ENDDO
            ilnb = ILNBLNK(PTRACERS_names(ii))
            WRITE(name,'(a,a)') 'surf_',PTRACERS_names(ii)(1:ilnb)
C            CALL MNC_CW_RL_W(pf,'ptr_flux_tave',0,0,name,
C     &           ptracerFluxtave(1-OLx,1-OLy,1,1,ii),myThid)
          ENDDO
        ENDIF
#endif /*  ALLOW_MNC  */

        IF ( PTRACERS_timeave_mdsio ) THEN
        
C         Set IO "context" for writing state
#ifdef USE_DFILE
          CALL DFILE_SET_RW
          CALL DFILE_SET_CONT_ON_ERROR
#endif
C         Read IO error counter
          beginIOErrCount = IO_ERRCOUNT(myThid)
c         prevPrec = writeBinaryPrec
c         writeBinaryPrec = precFloat32
          
          do iTracer=1,PTRACERS_numInUse
C            WRITE(suff,'(A7,I2.2,A1,I10.10)') 'PtrFluxtave',
C     &           iTracer,'.',myIter
C            CALL WRITE_FLD_XY_RL (suff,' ',
C     &           ptracerFluxtave(1-Olx,1-Oly,  1,1,iTracer),
C     &           myIter,myThid)
            WRITE(suff,'(A7,I2.2,A1,I10.10)') 'PTRtave',
     &           iTracer,'.',myIter
            CALL WRITE_FLD_XYZ_RL(suff,' ',
     &           ptracertave    (1-Olx,1-Oly,1,1,1,iTracer),
     &           myIter,myThid)
          enddo

c         writeBinaryPrec = prevPrec
C         Reread IO error counter
          endIOErrCount = IO_ERRCOUNT(myThid)
          
C         Check for IO errors
          IF ( endIOErrCount .NE. beginIOErrCount ) THEN
            WRITE(msgBuf,'(A)')  'S/R PTRACERS_WRITE_TIMEAVE'
            CALL PRINT_ERROR( msgBuf, 1 )
            WRITE(msgBuf,'(A)')  'Error writing out data'
            CALL PRINT_ERROR( msgBuf, 1 )
            WRITE(msgBuf,'(A,I10)') 'Timestep ',myIter
            CALL PRINT_ERROR( msgBuf, 1 )
          ELSE
            WRITE(msgBuf,'(A,I10)')  
     &           '// PTRACER time-average data written, t-step', myIter
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit, 
     &           SQUEEZE_RIGHT, 1 )
            WRITE(msgBuf,'(A)')  ' '
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit, 
     &           SQUEEZE_RIGHT, 1 )
          ENDIF
        
        ENDIF

        _END_MASTER( myThid )
        _BARRIER
        
        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)

C           Like before the 1rst iteration,
C           ==> call TIMEAVE_STATVARS with myIter=nIter0 :
C           1) Reset the averages to zero ;
C           2) Start to cumulate state-variables with Half time step.

            CALL PTRACERS_STATVARS(myTime, nIter0, bi, bj, myThid)

          ENDDO
        ENDDO

      ENDIF
 
#endif /* ALLOW_PTRACERS and ALLOW_TIMEAVE */

      RETURN
      END
