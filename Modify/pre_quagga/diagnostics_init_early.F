C $Header: /u/gcmpack/MITgcm/pkg/diagnostics/diagnostics_init_early.F,v 1.11 2005/07/09 17:19:48 jmc Exp $
C $Name:  $

#include "DIAG_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: DIAGNOSTICS_INIT_EARLY

C     !INTERFACE:
      SUBROUTINE DIAGNOSTICS_INIT_EARLY( myThid )

C     !DESCRIPTION:
C     Initialize available diagnostics list: set the following attributes:
C     name (=cdiag), parsing code (=gdiag), units (=udiag), and title (=tdiag)
C     Note: 1) diagnostics defined here are not presently filled. To use
C           one of them, one just needs to add a call to S/R DIAGNOSTICS_FILL
C           with the diagnostic name 'SDIAG...' or 'UDIAG...'
C           2) GDIAG is defined as character*16 and can be to character*1
C           parse(16) with the following codes currently defined:

C     \begin{center}
C       \begin{tabular}[h]{|c|c|}\hline
C         \textbf{Positions}  &  \textbf{Characters}
C         &  \textbf{Meanings} \\\hline
C         parse(1)  &  S  &  scalar \\
C                   &  U  &  vector component in X direction \\
C                   &  V  &  vector component in Y direction \\
C                   &  W  &  vector component in vertical direction \\
C         parse(2)  &  U  &  C-grid U-Point  \\
C                   &  V  &  C-grid V-Point  \\
C                   &  M  &  C-grid Mass Point  \\
C                   &  Z  &  C-grid Corner Point  \\
C         parse(3)  &  --  &  --Unused--  \\
C         parse(4)  &  P  &  positive definite  \\
C         parse(5 ) &  C  &  with counter array  \\
C                   &  D  &  disable an array for output  \\
C         parse(6--8) & '123'  &  3-digit mate number \\
C         parse(9)  &  U  &  model-level plus 1/2  \\
C                   &  M  &  model-level middle  \\
C                   &  L  &  model-level minus 1/2  \\
C         parse(10) &  0  &  levels = 0  \\
C                   &  1  &  levels = 1  \\
C                   &  R  &  levels = Nr  \\
C                   &  L  &  levels = MAX(Nr,NrPhys)  \\
C                   &  M  &  levels = MAX(Nr,NrPhys) - 1  \\
C                   &  G  &  levels = Ground_level Number \\
C                   &  I  &  levels = sea-Ice_level Number \\
C       \end{tabular}
C     \end{center}

C     !USES:
      implicit none
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS.h"

C     !INPUT PARAMETERS:
      integer myThid
CEOP

      integer n

      do n=1,ndiagMax
        kdiag(n) = 0
        cdiag(n) = '        '
        gdiag(n) = '                '
        udiag(n) = '                '
        tdiag(n) = ' '
      enddo

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     For each output variable,
C     specify Name (cdiag, 8c), Descriptions (tdiag, *c), Units (udiag, 16c)
C         and parsing code (location on C grid, 2D/3D, ...) (gdiag, 16c)
C----------------------------------------------------------------------
      cdiag( 1)  = 'SDIAG1  '
      tdiag( 1)  = 'User-Defined Surface   Diagnostic #1'
      udiag( 1)  = 'user-defined    '
      gdiag( 1)  = 'SM      L1      '

      cdiag( 2)  = 'SDIAG2  '
      tdiag( 2)  = 'User-Defined Surface   Diagnostic #2'
      udiag( 2)  = 'user-defined    '
      gdiag( 2)  = 'SM      L1      '

      cdiag( 3) = 'SDIAG3  '
      tdiag( 3) = 'User-Defined Surface   Diagnostic #3    '
      udiag( 3) = 'user-defined    '
      gdiag( 3) = 'SM      L1      '

      cdiag( 4) = 'TNDyTr12        '
      tdiag( 4) = 'ALK tendency due to merid advect       '
      udiag( 4) = 'mmol/m3/sec     '
      gdiag( 4) = 'SM      ML      '

      cdiag( 5) = 'TNDxTr12  '
      tdiag( 5) = 'alk tendency due to zonal advect       '
      udiag( 5) = 'mmol/m3/se      '
      gdiag( 5) = 'SM      ML      '

      cdiag( 6) = 'TNDyTr7        '
      tdiag( 6) = 'DIC tendency due to merid advect        '
      udiag( 6) = 'mmol/m3/sec     '
      gdiag( 6) = 'SM      ML      '

      cdiag( 7) = 'TNDxTr7'
      tdiag( 7) = 'DIC tendency due to zonal advect        '
      udiag( 7) = 'mmol/m3/s       '
      gdiag( 7) = 'SM      ML      '

      cdiag( 8) = 'TNDyTr9 '
      tdiag( 8) = 'DOC tendency from merid advect     '
      udiag( 8) = 'mmol/m3/s     '
      gdiag( 8) = 'SM      ML      '

      cdiag( 9) = 'TNDxTr9 '
      tdiag( 9) = 'DOC tendency from zonal advect        '
      udiag( 9) = 'mmol/m3/s     '
      gdiag( 9) = 'SM      ML      '

      cdiag(10) = 'TNRTr12  '
      tdiag(10) = 'tendency of alk due to vert advect      '
      udiag(10) = 'concentration/s '
      gdiag(10) = 'SM      ML      '

      cdiag(11)  = 'TNDxTr10  '
      tdiag(11)  = 'tendency of poc due to zonal advect    '
      udiag(11)  = 'mmol/m3/s  '
      gdiag(11)  = 'SM      ML      '

      cdiag(12)  = 'TNDyTr10  '
      tdiag(12)  = 'tendency of poc due to merid advect   '
      udiag(12)  = 'mmol/m3/sec'
      gdiag(12)  = 'SM      ML      '

      cdiag(13) = 'TNRTr7    '
      tdiag(13) = 'tendency of DIC due to vert advect             '
      udiag(13) = 'mmol/m3/sec     '
      gdiag(13) = 'SM      ML      '

      cdiag(14) = 'REMDIC    '
      tdiag(14) = 'tendency due to remineralization        '
      udiag(14) = 'mmol DIC /m3/s  '
      gdiag(14) = 'SM      ML      '

      cdiag(15) = 'BIODIC1   '
      tdiag(15) = 'tendency due to small phyto growth      '
      udiag(15) = 'mmol DIC/m3/s   '
      gdiag(15) = 'SM      ML      '

      cdiag(16) = 'BIODIC2   '
      tdiag(16) = 'tendency due to diatom growth           '
      udiag(16) = 'mmol /m3 /s     '
      gdiag(16) = 'SM      ML      '

      cdiag(17) = 'PICALK    '
      tdiag(17) = 'tendency of alk due to pic tendency     '
      udiag(17) = 'mmol equiv/m3/s '
      gdiag(17) = 'SM      ML      '

      cdiag(18) = 'PHY1ALK   '
      tdiag(18) = 'tendency of alk due to phy1 uptake      '
      udiag(18) = 'mmol/m3/s       '
      gdiag(18) = 'SM      ML      '

      cdiag(19) = 'PHY2ALK   '
      tdiag(19) = 'tendency of alk due to phy2 uptake      '
      udiag(19) = 'mmol/m3/s       '
      gdiag(19) = 'SM      ML      '

      cdiag(20) = 'REMALK    '
      tdiag(20) = 'tendency of alk due to remin of DON      '
      udiag(20) = 'mmol/m3/sec     '
      gdiag(20) = 'SM      ML      '

      cdiag(21) = 'docp1m'
      tdiag(21) = 'tendency of doc due to mortality of phy1   '
      udiag(21) = 'mmol/m3/sec     '
      gdiag(21) = 'SM      ML      '

      cdiag(22) = 'doczmort'
      tdiag(22) = 'tendency of doc due to mortality of zoo   '
      udiag(22) = 'mmol/m3/sec     '
      gdiag(22) = 'SM      ML      '
     
      cdiag(23) = 'docslop'
      tdiag(23) = 'tendency of doc from sloppy feeding   '
      udiag(23) = 'mmol/m3/sec     '
      gdiag(23) = 'SM      ML      '
      
      cdiag(24) = 'docrem'
      tdiag(24) = 'tendency of doc from remineralization   '
      udiag(24) = 'mmol/m3/sec     '
      gdiag(24) = 'SM      ML      '
      
      cdiag(25) = 'pocp1m'
      tdiag(25) = 'tendency of poc due to mortality of phy1   '
      udiag(25) = 'mmol/m3/sec     '
      gdiag(25) = 'SM      ML      '
      
      cdiag(26) = 'poczmort'
      tdiag(26) = 'tendency of poc due to mortality of zoo   '
      udiag(26) = 'mmol/m3/sec     '
      gdiag(26) = 'SM      ML      '
      
      cdiag(27) = 'pocslop'
      tdiag(27) = 'tendency of poc from sloppy feeding   '
      udiag(27) = 'mmol/m3/sec     '
      gdiag(27) = 'SM      ML      '

      cdiag(28) = 'pocrem'
      tdiag(28) = 'tendency of poc due to remin   '
      udiag(28) = 'mmol/m3/sec     '
      gdiag(28) = 'SM      ML      '

      cdiag(29) = 'pocsink'
      tdiag(29) = 'tendency of poc from sinking in/out of layer  '
      udiag(29) = 'mmol/m3/sec     '
      gdiag(29) = 'SM      ML      '

      cdiag(30) = 'SDIAGC  '
      tdiag(30) = 'User-Defined Counted Surface Diagnostic         '
      udiag(30) = 'user-defined    '
      gdiag(30) = 'SM  C022L1      '

      cdiag(31) = 'VALK    '
      tdiag(31) = 'Tendency of ALK due to air-sea        '
      udiag(31) = 'mmol eq / m3 / s'
      gdiag(31) = 'SM      L1      '

      cdiag(32) = 'TNRTr9  '
      tdiag(32) = 'tendency of doc due to vert advect      '
      udiag(32) = 'concentration/s '
      gdiag(32) = 'SM      ML      '

      cdiag(33) = 'TNRTr10  '
      tdiag(33) = 'tendency of poc due to vert advect      '
      udiag(33) = 'concentration/s '
      gdiag(33) = 'SM      ML      '

      cdiag(34) = 'docp2m'
      tdiag(34) = 'tendency of doc due to mortality of phy2   '
      udiag(34) = 'mmol/m3/sec     '
      gdiag(34) = 'SM      ML      '

      cdiag(35) = 'pocp2m'
      tdiag(35) = 'tendency of poc due to mortality of phy2   '
      udiag(35) = 'mmol/m3/sec     '
      gdiag(35) = 'SM      ML      '
c-    set the total number of available diagnostics
      ndiagt = 35

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      RETURN
      END
