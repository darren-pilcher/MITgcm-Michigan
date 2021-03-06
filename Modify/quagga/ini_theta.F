C $Header: /u/gcmpack/MITgcm/model/src/ini_theta.F,v 1.23 2005/11/08 23:01:10 cnh Exp $
C $Name:  $

#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_THETA
C     !INTERFACE:
      SUBROUTINE INI_THETA( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_THETA                                      
C     | o Set model initial temperature field.                    
C     *==========================================================*
C     | There are several options for setting the initial         
C     | temperature file                                          
C     |  1. Inline code                                           
C     |  2. Vertical profile ( uniform T in X and Y )             
C     |  3. Three-dimensional data from a file. For example from  
C     |     Levitus or from a checkpoint file from a previous     
C     |     integration.                                          
C     | In addition to setting the temperature field we also      
C     | set the initial temperature tendency term here.           
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid -  Number of this instance of INI_THETA
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj  - Loop counters
C     I,J,K
      INTEGER bi, bj
      INTEGER I, J, K, localWarnings
      _RL     Tfreezing
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C--   Initialise temperature field to the vertical reference profile
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K=1,Nr
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           theta(I,J,K,bi,bj) = tRef(K)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _BARRIER

      IF ( hydrogThetaFile .NE. ' ' ) THEN
C      _BEGIN_MASTER( myThid )
       CALL READ_FLD_XYZ_RL( hydrogThetaFile, ' ', theta, 0, myThid )
C      _END_MASTER(myThid)
       _EXCH_XYZ_R8(theta,myThid)
      ENDIF

C--   Apply mask and test consistency
      localWarnings=0
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K=1,Nr
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           IF (maskC(I,J,K,bi,bj).EQ.0.) theta(I,J,K,bi,bj) = 0.
          ENDDO
         ENDDO
         IF ( tRef(k).NE.0. ) THEN
          DO J=1,sNy
           DO I=1,sNx
            IF (  maskC(I,J,K,bi,bj).NE.0.
     &      .AND. theta(I,J,K,bi,bj).EQ.0. ) THEN
              localWarnings=localWarnings+1
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      IF (localWarnings.NE.0) THEN
       WRITE(msgBuf,'(A,A)')
     &  'S/R INI_THETA: theta = 0 identically. If this is intentional',
     &  'you will need to edit ini_theta.F to avoid this safety check'
       CALL PRINT_ERROR( msgBuf , myThid)
       STOP 'ABNORMAL END: S/R INI_THETA'
      ENDIF

C--   Check that there are no values of temperature below freezing point.
      Tfreezing=0. _d 0
C      Tfreezing=-1.9 _d 0
      IF ( allowFreezing ) THEN
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO K=1,Nr
          DO J=1-Oly,sNy+Oly
           DO I=1-Olx,sNx+Olx
            IF (theta(I,J,k,bi,bj) .LT. Tfreezing) THEN
               theta(I,J,K,bi,bj) = Tfreezing
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      CALL PLOT_FIELD_XYZRL( theta, 'Initial Temperature' ,
     &                       Nr, 1, myThid )

      RETURN
      END
