C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM.h,v 1.8 2005/08/18 18:23:02 stephd Exp $
C $Name:  $

#ifdef ALLOW_GCHEM

CBOP
C    !ROUTINE: GCHEM.h
C    !INTERFACE:
 
C    !DESCRIPTION:
C Contains tracer parameters and input files for chemical tracers.
C These can be read in from data.gchem
C
C  nsubtime   : : number of chemistry timesteps per deltaTtracer
C                 (default 1) 
C  WindFile   : : file name of wind speeds that may be needed for
C                 biogeochemical experiments
C  AtmospFile : : file name of atmospheric pressure that may be needed for
C                 biogeochemical experiments
C  IceFile    : : file name of seaice fraction that may be needed for
C                 biogeochemical experiments
C  IronFile   : : file name of aeolian iron flux that may be needed for
C                 biogeochemical experiments
C  SilicaFile : : file name of surface silica that may be needed for
C                 biogeochemical experiments

C  
      INTEGER nsubtime
      CHARACTER*(MAX_LEN_FNAM) WindFile
      CHARACTER*(MAX_LEN_FNAM) AtmospFile
      CHARACTER*(MAX_LEN_FNAM) IceFile
c      CHARACTER*(MAX_LEN_FNAM) IronFile
c      CHARACTER*(MAX_LEN_FNAM) SilicaFile
      CHARACTER*(MAX_LEN_FNAM) SurParFile 
      CHARACTER*(MAX_LEN_FNAM) Filename1
      CHARACTER*(MAX_LEN_FNAM) Filename2
      CHARACTER*(MAX_LEN_FNAM) Filename3
      CHARACTER*(MAX_LEN_FNAM) Filename4
      CHARACTER*(MAX_LEN_FNAM) Filename5
      INTEGER gchem_int1
      INTEGER gchem_int2
      INTEGER gchem_int3
      INTEGER gchem_int4
      INTEGER gchem_int5
      _RL     gchem_rl1
      _RL     gchem_rl2
      _RL     gchem_rl3
      _RL     gchem_rl4
      _RL     gchem_rl5

      COMMON /GCHEM_PARAMS/
     &                   WindFile,
     &                   AtmospFile,
     &                   IceFile,
c     &                   IronFile,
c     &                   SilicaFile,
     &                   SurParFile,
     &                   Filename1,
     &                   Filename2,
     &                   Filename3,
     &                   Filename4,
     &                   Filename5,
     &                   nsubtime,
     &           gchem_int1, gchem_int2, gchem_int3,
     &           gchem_int4, gchem_int5,
     &           gchem_rl1, gchem_rl2, gchem_rl3,
     &           gchem_rl4, gchem_rl5
CEOP

#endif /* ALLOW_GCHEM */
