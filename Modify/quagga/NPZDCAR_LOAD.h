#include "GCHEM_OPTIONS.h"

      COMMON/ npzdcar_load /
c     &    fice0, fice1, featmos0, featmos1,
     &    fice0, fice1, 
     &    wspeed0, wspeed1, sur_par0, sur_par1,
     &    atmosp0, atmosp1,
     &    quagga0, quagga1
      _RS fice0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS fice1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RS featmos0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RS featmos1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wspeed0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wspeed1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS sur_par0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS sur_par1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS atmosp0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS atmosp1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS quagga0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS quagga1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
