#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_NPZDCAR

c=============================================================
c subroutine NPZDCAR_forcing
c step forward bio-chemical tracers in time
C==============================================================
c ******** Modified by Darren 3/2012 to remove river effects *********
         SUBROUTINE NPZDCAR_Forcing(Nut,Phy,Phy2,Zoo,DON,PON,
     &                              Sil,POS,DIC,ZoC,DOC,POC,Oxy,
     &                              Alk,
     &                       bi,bj,imin,imax,jmin,jmax,
     &                        myIter,myTime,myThid)
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS.h"
#include "NPZDCAR.h"
#include "GCHEM.h"
#ifdef PAR_FROMSW
#include "FFIELDS.h"
#endif
c #ifdef ALLOW_RIVERS
c #include "rivers_param.h"
c #include "rivers.h"
c #include "rivers_fields.h"
c #include "RIVERS_OPTIONS.h"
c #endif

C     === Global variables ===
c tracers
      _RL  Nut(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Phy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Phy2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Zoo(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  DON(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  PON(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Sil(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  POS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  FeT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  ZFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  DOF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  POF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  ZoC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  DOC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  POC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  PIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Oxy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  Alk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  TDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  RDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
       INTEGER myIter
       _RL myTime
       INTEGER myThid
       INTEGER bi,bj,imin,imax,jmin,jmax


C============== Local variables ============================================
c tracer tendencies
      _RL  GNut(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GPhy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GPhy2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GZoo(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GDON(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GPON(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GSil(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GPOS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  GFeT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  GZFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  GDOF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  GPOF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GZoC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GDOC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GPOC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  GPIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GOxy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  GAlk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy) 
c      _RL  GTDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy) 
c      _RL  GRDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy) 
c#ifdef READ_PAR
c      _RL  GPAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy) 
c      _RL  growpar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c#endif
c
c      _RL  freefe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  SURC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  SURO(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  SURA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  CHL(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  parsfc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c #ifdef ALLOW_RIVERS
c      _RL  GRALK(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  GRDOC(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  GRDIC(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  rflow(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  ralk(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  rdic(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  rdoc(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c
c      _RL  rivflow(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c      _RL  rivalk(1-Olx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c #endif

#ifdef ALLOW_DIAGNOSTICS
c     IF ( useDiagnostics ) THEN
      CHARACTER*8 diagName
      _RL  valk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  remalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  phy1alk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  phy2alk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      _RL  picalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  biodic1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  remdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  biodic2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  docp2m(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  docp1m(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  doczmort(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  docslop(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  docrem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  pocp1m(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  pocp2m(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  poczmort(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  pocsink(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  pocrem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
      _RL  pocslop(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx, nSy)
c      ENDIF
#endif

c light variables
      _RL  sfac(1-OLy:sNy+OLy)
      _RL  atten, lite
      _RL  atten0
c
      _RL  mu, mu2, graz, graz2
      _RL  nu, py, py2, zo, dn
      _RL  pyup, py2up 
      _RL  pym, py2m, zom, zfm, zcm
      _RL  pn,  pnup
      _RL  zeep, zp
      _RL  si, ps, psup
      _RL  fe, zf, df,pf, pfup
      _RL  ca, zc, dc, pc, pcup, pcc, pccup
      _RL  ox, al
c      _RL tc, rfc 
c      _RL  freefu
      _RL  tempfunc
      _RL  temp_zoofunc
      _RL  reminn_new, reminc_new
c      _RL  reminr_new
      _RL  newtime
        INTEGER I,J,K,it, ktmp
c
c
       DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        surc(i,j,bi,bj)=0.d0
        suro(i,j,bi,bj)=0.d0
        sura(i,j,bi,bj)=0.d0
c #ifdef ALLOW_RIVERS
c        GRDOC(i,j,bi,bj)=0.d0
c        GRALK(i,j,bi,bj)=0.d0
c        GRDIC(i,j,bi,bj)=0.d0
c #endif
c        do k=1,Nr
c           freefe(i,j,k,bi,bj)=0.d0
c        enddo
       ENDDO
       ENDDO
c
c bio-chemical time loop
c--------------------------------------------------
       do it=1,nsubtime
c -------------------------------------------------
         newtime=myTime-dTtracerLev(1)+
     &           float(it)*dTtracerLev(1)/float(nsubtime)
c        print*,'it  ',it,newtime,nsubtime,myTime

c carbon air-sea interaction
       CALL DIC_SURFFORCING( DIC, ALK,
     &                     SURC,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

c      print*,'forc, co2 flux',SURC(5,3,1,1)
c
c alkalinity air-sea interaction
       CALL ALK_SURFFORCING( ALK, SURA,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

c      print*,'forc, alk flux',SURA(10,10,1,1)
c
c carbon air-sea interaction
       CALL O2_SURFFORCING( Oxy, SURO,
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)

c       print*,'forc, o2 flux',SURO(10,10,1,1)
c #ifdef ALLOW_RIVERS
c         DO j=1-OLy,sNy+OLy
c         DO i=1-OLx,sNx+OLx
c             rflow(i,j,bi,bj) = rivflow(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c             ralk(i,j,bi,bj) = rivalk(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c             rdic(i,j,bi,bj) = rivdic(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c             rdoc(i,j,bi,bj) = rivdoc(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c             print*,'rivflow ',rivflow(i,j,bi,bj)
c             print*,'rflow ',rflow(i,j,bi,bj)
c             print*,'ralk ',ralk(i,j,bi,bj)
c         ENDDO
c         ENDDO
c              CALL RIVERS_DOCFORCE( rflow, rdoc, 
c     &                            GRDOC,
c     &             bi,bj,imin,imax,jmin,jmax,
c     &             myIter,myTime,myThid)
c              CALL RIVERS_ALKFORCE( rflow, ralk, ALK,
c     &                            GRALK,
c     &             bi,bj,imin,imax,jmin,jmax,
c     &             myIter,myTime,myThid)
c              CALL RIVERS_DICFORCE( rflow, rdic, DIC,
c     &                            GRDIC,
c     &             bi,bj,imin,imax,jmin,jmax,
c     &             myIter,myTime,myThid)
c        print*,'river dic effect ',GRDIC(5,3,1,1)
c #endif
c
c find free iron
c         call fe_chem(bi,bj,iMin,iMax,jMin,jMax, FeT, freefe,
c     &                myIter, mythid)
c find incidence light
#ifndef  READ_PAR
#ifdef PAR_FROMSW
         DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           sur_par(i,j,bi,bj)=-0.4*Qsw(i,j,bi,bj)*maskC(i,j,1,bi,bj)
         ENDDO
         ENDDO
#else
         DO j=1-OLy,sNy+OLy
          sfac(j)=0.d0
         ENDDO
         call insol(newTime,sfac,bj)
         DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
c           sur_par(i,j,bi,bj)=sfac(j)*maskC(i,j,1,bi,bj)
           sur_par(i,j,bi,bj)=0.45*sfac(j)*maskC(i,j,1,bi,bj)
         ENDDO
         ENDDO
#endif
#endif
ccccccccccccccccc Val April 2009 cccccccccccccccccccccccccccccc
#ifdef READ_PAR
c   USE THIS PAR FOR CHLOROPHYLL
c   USE instantaneous par for growth rates
c         print*,' reading PAR file'
         DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           parsfc(i,j,bi,bj)=0.45*sur_par(i,j,bi,bj)*maskC(i,j,1,bi,bj)
c           growpar(i,j,bi,bj)=-0.45*Qsw(i,j,bi,bj)*maskC(i,j,1,bi,bj)
         ENDDO
         ENDDO
c         print*,' surface_PAR ',parsfc(10,10,1,1)
#endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccc Spatial Loop cccccccccccccccccccccccccccc
         DO j=1,sNy 
         DO i=1,sNx
          atten=0.d0
c for each layer find
          do k= 1, NR   
           if (hFacC(i,j,k,bi,bj).gt.0.d0) then
c temperature dependepence of remin
            reminn_new=reminn
            reminc_new=reminc
c            remint_new=remint
c            reminr_new=reminr
c            print*,('reminn',reminn_new)
c            print*,('reminc',reminc_new)
c
c make sure we only deal with positive definite
c numbers
             nu  = max(Nut(i,j,k,bi,bj),0.D0)
             py =  max(Phy(i,j,k,bi,bj),0.D0)
             py2 = max(Phy2(i,j,k,bi,bj),0.D0)
             zo  = max(Zoo(i,j,k,bi,bj),0.D0)
             dn  = max(DON(i,j,k,bi,bj),0.D0)
             if (k.eq.1) then
              pyup = 0.d0
              py2up = 0.d0
             else
              pyup = max(Phy(i,j,k-1,bi,bj),0.D0)
              py2up = max(Phy2(i,j,k-1,bi,bj),0.D0)
             endif
             si = max(Sil(i,j,k,bi,bj),0.d0)
c             fe = max(FeT(i,j,k,bi,bj),0.d0)
c             zf = max(ZFe(i,j,k,bi,bj),0.d0)
c             df = max(DOF(i,j,k,bi,bj),0.d0)
c             freefu = max(freefe(i,j,k,bi,bj),0.d0)
             ca = max(DIC(i,j,k,bi,bj),0.d0)
             zc = max(ZoC(i,j,k,bi,bj),0.d0)
             dc = max(DOC(i,j,k,bi,bj),0.d0)
c             rfc = max(RDC(i,j,k,bi,bj),0.d0)
             ox = max(Oxy(i,j,k,bi,bj),0.d0)
             al = max(Alk(i,j,k,bi,bj),0.d0)
c             tc = max(TDC(i,j,k,bi,bj),0.d0)
c             if (k.lt.5) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             if (k.lt.30) then
c changed oct 2009 - assuming meant top 50m in N. Atlantic
c layer thicknesses different in superior model - 10 layers top 50m
c changed 03/2012 for 29 layers in lake Michigan -Darren 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               pym = max((py-phymin),0.D0)
               py2m =max((py2-phymin),0.D0)
               zom = max((zo-zoomin),0.D0)
c               zfm = max((zf-zoominfe),0.D0)
               zcm = max((zc-zoominc),0.D0)
             else
               pym = py
               py2m = py2
               zom = zo
c               zfm = zf
               zcm = zc
             endif
c if explicit sinking of particulate matter
             pn  = max(PON(i,j,k,bi,bj),0.D0)
             ps  = max(POS(i,j,k,bi,bj),0.D0)
c             pf  = max(POF(i,j,k,bi,bj),0.D0)
             pc  = max(POC(i,j,k,bi,bj),0.D0)
c             pcc  = max(PIC(i,j,k,bi,bj),0.D0) 
             if (k.eq.1) then
               pnup =0.d0
               psup = 0.d0
c               pfup = 0.d0
               pcup = 0.d0
c               pccup = 0.d0
             else
               pnup = max(PON(i,j,k-1,bi,bj),0.D0)
               psup = max(POS(i,j,k-1,bi,bj),0.d0)
c               pfup = max(POF(i,j,k-1,bi,bj),0.d0)
               pcup = max(POC(i,j,k-1,bi,bj),0.d0)
c               pccup = max(PIC(i,j,k-1,bi,bj),0.d0)
             endif
c for start, provide value to par (need to make better)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Val Apr 2009
#ifdef READ_PAR
            if (myIter.eq.nIter0) then
             if (k.eq.1) atten0= 0.D0
             atten0= atten0 +( k0 +kc*((py+py2)*16) )*5.D-1*(drF(k))
             if (k.gt.1) atten0 = atten0
     &            +( k0+kc*(CHL(i,j,k-1)) )*5.D-1*(drF(k-1))
             PAR(i,j,k,bi,bj)=
     &          parsfc(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten0)
c             GPAR(i,j,k,bi,bj) = 
c     &           growpar(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten0)
            endif
#else
            if (myIter.eq.nIter0) then
             if (k.eq.1) atten0= 0.D0
             atten0= atten0 +( k0 +kc*((py+py2)*16) )*5.D-1*(drF(k))
             if (k.gt.1) atten0 = atten0
     &            +( k0+kc*(CHL(i,j,k-1)) )*5.D-1*(drF(k-1))
             PAR(i,j,k,bi,bj)=
     &          sur_par(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten0)
            endif
#endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c USING 24HR RUNNING MEAN PAR FOR THIS !!! 
c REDFIELD RATIO C:P DIFFERENT IN LAKE !!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c Nutrient stress no effect on Chlorophyll here:
c Val Apr 2009
            CHL(i,j,k) = ( ( chlpmax - (chlpmax-chlpmin) *
     &            min(PAR(i,j,k,bi,bj)/istar, 1.d0)) ) *
     &            (py + py2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c cholorphyll concentrations
c            CHL(i,j,k)=(chlpmax-(chlpmax-chlpmin)*
c     &           min(PAR(i,j,k,bi,bj)/istar, 1.d0) )*
c     &           (py*(nu/(nu+nlim)) +
c     &            py2*(nu/(nu+nlim2)) ) 
c
c
c updated by val bennington June 2008
c calculates chlorophyll from biomass only
c since light is zero at night, were getting incorrect
c daily chlorophyll values as previously written
c could also use 24-hr average light to calc chlorophyll
c conversion from steph dutkiewicz for N. Atlantic is:
c (phy1+phy2)*1.5*16
c lake has different C:P ratio --> 200 instead of 120
c conversion to CHL scaled by this:
c            CHL(i,j,k) = (py+py2)*16*200/120*1.5
c
c cholorphyll concentrations
c            CHL(i,j,k)=(chlpmax-(chlpmax-chlpmin)*
c     &           min(PAR(i,j,k,bi,bj)/istar, 1.d0) )*
c     &           (py*(nu/(nu+nlim)) +
c     &            py2*(nu/(nu+nlim2)) ) 
c            CHL(i,j,k)=(chlpmax-(chlpmax-chlpmin)*
c     &           min(PAR(i,j,k,bi,bj)/istar, 1.d0) )*
c     &           (py*min(nu/(nu+nlim),fe/(fe+flim)) +
c     &            py2*min(fe/(fe+flim2),
c     &                      min(nu/(nu+nlim2),si/(si+slim2)) )  )
c attenuate light
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Val Apr 2009 cccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef READ_PAR
             if (k.eq.1) atten= 0.D0
             atten= atten +( k0 +kc*(CHL(i,j,k)) )*5.D-1*(drF(k))      
             if (k.gt.1) atten = atten                 
     &            +( k0+kc*(CHL(i,j,k-1)) )*5.D-1*(drF(k-1))     
             PAR(i,j,k,bi,bj)=
     &          parsfc(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten)
c             GPAR(i,j,k,bi,bj)=
c     &          growpar(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten)
             lite = PAR(i,j,k,bi,bj)
#else
             if (k.eq.1) atten= 0.D0
             atten= atten +( k0 +kc*(CHL(i,j,k)) )*5.D-1*(drF(k))      
             if (k.gt.1) atten = atten                 
     &            +( k0+kc*(CHL(i,j,k-1)) )*5.D-1*(drF(k-1))     
             PAR(i,j,k,bi,bj)=
     &          sur_par(i,j,bi,bj)*(1.d0-fice(i,j,bi,bj))*exp(-atten)
             lite=PAR(i,j,k,bi,bj)
#endif
c             print*,'PAR(10,10,k,1,1) ',PAR(10,10,k,1,1)
c phytoplankton growth rates
c             mu = phygrow*min(lite/(ilim+lite),nu/(nu+nlim))
             mu= phygrow*lite/(ilim+lite)   
             mu= mu*(nu/(nu+nlim))
c             mu= mu* min(nu/(nu+nlim),fe/(fe+flim))
c             mu2 = phygrow2*min(lite/(ilim2+lite),nu/(nu+nlim2))
             mu2= phygrow2*lite/(ilim2+lite)
c             mu2= mu2*(nu/(nu+nlim2))
             mu2= mu2*min(nu/(nu+nlim2),si/(si+slim2))
c take into account ice coverage
c this is a double-counting of the effect of ice!
c already took ice into account with the light
c             mu= (1.d0 -fice(i,j,bi,bj))*mu
c             mu2= (1.d0 -fice(i,j,bi,bj))*mu2
c temperature function
c             tempfunc=1.d0/temp_coeff_norm * 
c     &          temp_coeff_base**theta(i,j,k,bi,bj)
c June 20, 2009
c Use a rendition of Moore's function (same shape, easier to tune)
c to increase growth in April without slowing growth in August:
             if (theta(i,j,k,bi,bj).gt.10.d0) then
             temp_coeff_base = -3400
             temp_coeff_norm = 280.15
             else
             temp_coeff_base = temp_coeff_base
             temp_coeff_norm = temp_coeff_norm
             endif 
             tempfunc = exp(temp_coeff_base*((1/(theta(i,j,k,bi,bj)+273.15)-
     &          (1/temp_coeff_norm))))
             mu = mu * tempfunc
             mu2 = mu2 * tempfunc
             reminn_new=reminn_new*tempfunc
             reminc_new=reminc_new*tempfunc
c  semi-labile pool not affected by temperature in remin?
c             remint_new=remint_new*tempfunc
c zooplankton grazing rates
c trying to increase wintering population
c make it so that zooplankton can only "see" phy above
c the seed population since mortality only acts on phy 
c above the seed population:
c             zeep= zplat*py + zplat2*py2
             zeep = zplat*pym + zplat2*py2m
             zp=max(zeep-grazlim,0.d0)
c             temp_zoofunc=1/temp_zoo_norm*
c    &                    temp_zoo_base**theta(i,j,k,bi,bj)
             if (zeep.gt.0.d0) then
               graz= zoograze*zplat*py/zeep*
     &                          (zp)/(zp+plim)
               graz=graz*tempfunc
               graz2= zoograze*zplat2*py2/zeep*
     &                          (zp)/(zp+plim)
               graz2=graz2*tempfunc
             else
               graz=0.d0
               graz2=0.d0
             endif
c calculate quagga mussel grazing 
c added by Darren Pilcher 11/13
c quaggatempfunc is based on a best fit line to Diggins et al 2001 temperature 
c growth rates.  The if loop is to check whether the next depth layer is beyond 
c lake bottom.  This way, quaggatempfunc will only be calculated if the current 
c depth layer is the bottom layer.  Otherwise, quaggatempfunc will be equal
c to 0, and all quagga functions will be equal to 0
c removed temperature dependence for now DJP 02/2014

             if (k.lt.29.d0) then 
             if (hFacC(i,j,(k+1.d0),bi,bj).eq.0.d0) then
c             quaggatempfunc = -2.d-12 * (theta(i,j,k,bi,bj) **2) + 8.d-11 
c     &              * theta(i,j,k,bi,bj) + 2.d-10
             quaggatempfunc = 1.d0
             else
             quaggatempfunc = 0.d0
             endif
             else
c             quaggatempfunc = -2.d-12 * (theta(i,j,k,bi,bj) **2) + 8.d-11
c     &              * theta(i,j,k,bi,bj) + 2.d-10
             quaggatempfunc = 1.d0
             endif    
             quaggagraze = quagga(i,j,bi,bj) * quaggatempfunc * qplat
             quaggagraze2= quagga(i,j,bi,bj) * quaggatempfunc * qplat2
c             print*,' quaggatempfunc value',quaggatempfunc
c             print*,' quagga population',quagga(i,j,bi,bj)
c             if (quaggagraze.gt.0.d0) then
c             print*,' quagga pop',quagga(i,j,bi,bj)
c             print*,' quaggatempfunc',quaggatempfunc 
c             print*,' qplat',qplat
c             print*,' quagga graze',quaggagraze
c             endif
c
c calculate tendency terms
             GNut(i,j,k,bi,bj)=
     &              -mu*py -mu2*py2 +reminn_new*dn
             GPhy(i,j,k,bi,bj)=
     &              +mu*py -mort*pym -graz*zo
     &              - quaggagraze*pym
             GPhy2(i,j,k,bi,bj)=
     &              +mu2*py2 -mort2*py2m -graz2*zo
     &              - quaggagraze2*py2m
             GZoo(i,j,k,bi,bj)=
     &              +gampn*graz*zo +gampn2*graz2*zo
     &              -mortz*zom
             GDON(i,j,k,bi,bj)=
     &              +donfracmn1*mort*pym
     &              +donfracmn2*mort2*py2m
     &              +donfracmn*mortz*zom
     &              +donfracgn*
     &               ( (1.d0-gampn)*graz*zo 
     &                +(1.d0-gampn2)*graz2*zo )
     &              -reminn_new*dn
     &              +donfracq*quaggagraze*pym
     &              +donfracq*quaggagraze2*py2m
c
             GSil(i,j,k,bi,bj)=
     &              -Rsip2*mu2*py2 
c
c             GFeT(i,j,k,bi,bj)=
c     &              -Rfep*mu*py -Rfep2*mu2*py2 
c     &              +reminf*df -scav*freefu
c             if (k.eq.1)  GFeT(i,j,k,bi,bj)=
c     &              GFeT(i,j,k,bi,bj) 
c     &              +alpfe*inputFe(i,j,bi,bj)/drF(1)
c             GZFe(i,j,k,bi,bj)= 
c     &              +Rfep*gampn*graz*zo + Rfep2*gampn2*graz2*zo
c     &              -mortz*zfm
c             GDOF(i,j,k,bi,bj)=
c     &              +Rfep*( donfracmf1*mort*pym
c     &                     +donfracgf*(1.d0-gampn)*graz*zo )
c     &              +Rfep2*( donfracmf2*mort2*py2m
c     &                      +donfracgf*(1.d0-gampn2)*graz2*zo )
c     &              +donfracmf*mortz*zfm
c     &              -reminf*df
c
#ifdef ALLOW_DIAGNOSTICS
c        IF ( useDiagnostics ) THEN
c            biodic1(i,j,k,bi,bj) = -Rcp*(1.d0+rainratio)*mu*py 
            biodic1(i,j,k,bi,bj) = -Rcp*(1.d0)*mu*py 
            biodic2(i,j,k,bi,bj) = -Rcp2*mu2*py2
            remdic(i,j,k,bi,bj) = reminc_new*dc
            remalk(i,j,k,bi,bj) = reminn_new*dn*Rnp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c added oct 2009 for carbon pathways in superior:
c docphy1mort is doc tendency due to mortality of small phyto
c docphy2mort is doc tendency due to mortality of large phyto
c doczoomort is doc tendency due to zoo mortality
c docsloppy is doc tendency due to sloppy grazing
c docremloss is loss of doc due to remineralization
            docp1m(i,j,k,bi,bj) = Rcp*donfracmc1*mort*pym
            docp2m(i,j,k,bi,bj) = Rcp2*donfracmc2*mort2*py2m
            doczmort(i,j,k,bi,bj) = donfracmc*mortz*zcm
            docslop(i,j,k,bi,bj) = Rcp*donfracgc*
     &                             ( (1.d0-gampn)*graz*zo
     &                             +(1.d0-gampn2)*graz2*zo )
            docrem(i,j,k,bi,bj) = -reminn_new*dc
            pocp1m(i,j,k,bi,bj) = (1.D0-donfracmc1)*Rcp*mort*pym
            pocp2m(i,j,k,bi,bj) = (1.D0-donfracmc2)*Rcp2*mort2*py2m
            poczmort(i,j,k,bi,bj) = (1.D0-donfracmc)*mortz*zcm
            pocslop(i,j,k,bi,bj) = (1.D0-donfracgc)*
     &                             ( Rcp*(1.d0-gampn)*graz*zo
     &                             + Rcp2*(1.d0-gampn2)*graz2*zo )
            ktmp=min(nR,k+1)
            if (hFacC(i,j,ktmp,bi,bj).eq.0.d0.or.k.eq.Nr) then
            pocsink(i,j,k,bi,bj) = wc_sink*(pcup)/
     &                      (drf(k)*hFacC(i,j,k,bi,bj)) 
            else
            pocsink(i,j,k,bi,bj)=wc_sink*(pcup-pc)/(drf(k)
     &                               *hFacC(i,j,k,bi,bj))
            endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ENDIF
#endif

c#ifdef ALLOW_DIAGNOSTICS
c        IF ( useDiagnostics ) THEN
c            remalk(i,j,k,bi,bj) = reminn_new*dn*Rnp
c            phy1alk(i,j,k,bi,bj) = biodic1(i,j,k,bi,bj)*Rnp/Rcp
c            phy2alk(i,j,k,bi,bj) = biodic2(i,j,k,bi,bj)*Rnp/Rcp2
c        ENDIF
c#endif
c
            GDIC(i,j,k,bi,bj)=
c     &              -Rcp*(1.d0+rainratio)*mu*py -Rcp2*mu2*py2
     &              -Rcp*(1.d0)*mu*py -Rcp2*mu2*py2
     &              +reminc_new*dc
c     &              +remint_new*tc
c     &              +reminr_new*rfc 
            if (k.eq.1)  then
c            print*,'before air-lake ',GDIC(i,j,k,bi,bj)
            GDIC(i,j,k,bi,bj)=
     &              GDIC(i,j,k,bi,bj)
     &              +SURC(i,j,bi,bj)
c            print*,'air-lake ',SURC(i,j,bi,bj)
            endif
c new --- terrestrial carbon:
c             GTDC(i,j,k,bi,bj)=
c     &              -remint_new*tc
c             GRDC(i,j,k,bi,bj)=
c     &              -reminr_new*rfc
c     &              -reminc_new*tc*labilefrac
cccccccccccccccccccccccccccccccccc
             GZoC(i,j,k,bi,bj)=
     &              +Rcp*gampn*graz*zo + Rcp2*gampn2*graz2*zo
     &              -mortz*zcm
             GDOC(i,j,k,bi,bj)=
     &              +Rcp*( donfracmc1*mort*pym
     &                     +donfracgc*(1.d0-gampn)*graz*zo )
     &              +Rcp2*( donfracmc2*mort2*py2m
     &                      +donfracgc*(1.d0-gampn2)*graz2*zo )
     &              +donfracmc*mortz*zcm
     &              -reminc_new*dc
     &              +Rcp*donfracq*quaggagraze*pym
     &              +Rcp2*donfracq*quaggagraze2*py2m
c             GPIC(i,j,k,bi,bj)=Rcp*rainratio*mu*py
             GOxy(i,j,k,bi,bj)=Rop*(mu*py +mu2*py2
     &               -reminn_new*dn - reminn_new*pn)
             if (k.eq.1)  GOxy(i,j,k,bi,bj)=
     &              GOxy(i,j,k,bi,bj)
     &              +SURO(i,j,bi,bj) 
c Added by Val Bennington July 2008 
c Allows river inputs to affect DOC concentrations 
c See added pkg: rivers
c #ifdef ALLOW_RIVERS
c              if (k.eq.1) then
c              GDOC(i,j,k,bi,bj) = GDOC(i,j,k,bi,bj)
c     &           +GRDOC(i,j,bi,bj)*labfrac
c              GTDC(i,j,k,bi,bj) = GTDC(i,j,k,bi,bj)
c     &          +GRDOC(i,j,bi,bj)*slabfrac
c              GRDC(i,j,k,bi,bj) = GRDC(i,j,k,bi,bj)
c     &          +GRDOC(i,j,bi,bj)*refrac
c              endif
c #endif
c
c particulate matter for each layer
             GPON(i,j,k,bi,bj)=
     &            +(1.D0 -donfracmn1)*mort*pym
     &            +(1.D0 -donfracmn2)*mort2*py2m
     &            +(1.D0 -donfracmn)*mortz*zom
     &            +(1.D0 -donfracgn)*
     &              ( (1.d0 -gampn)*graz*zo 
     &               +(1.d0 -gampn2)*graz2*zo)
     &            +(1.D0 -donfracq)*quaggagraze*pym
     &            +(1.D0 -donfracq)*quaggagraze2*py2m            
             GPOS(i,j,k,bi,bj)=
     &            +Rsip2*
     &             ( (1.D0)*(mort2*py2m)
     &              +(1.D0)*(graz2*zo) )
c             GPOF(i,j,k,bi,bj)=
c     &            (1.D0 -donfracmf1)*Rfep*mort*pym
c     &           +(1.D0 -donfracmf2)*Rfep2*mort2*py2m
c     &           +(1.D0 -donfracmf)*mortz*zfm
c     &           +(1.D0 -donfracgf)*
c     &              ( Rfep*(1.d0-gampn)*graz*zo
c     &               +Rfep2*(1.d0-gampn2)*graz2*zo )
             GPOC(i,j,k,bi,bj)=
     &            (1.D0 -donfracmc1)*Rcp*mort*pym
     &           +(1.D0 -donfracmc2)*Rcp2*mort2*py2m
     &           +(1.D0 -donfracmc)*mortz*zcm
     &           +(1.D0 -donfracgc)*
     &              ( Rcp*(1.d0-gampn)*graz*zo
     &               +Rcp2*(1.d0-gampn2)*graz2*zo )
     &           +(1.D0 -donfracq)*Rcp*quaggagraze*pym
     &           +(1.D0 -donfracq)*Rcp2*quaggagraze2*py2m
c for implicit sinking
c pom out of bottom is instantly remineralized
             ktmp=min(nR,k+1)
             if (hFacC(i,j,ktmp,bi,bj).eq.0.d0.or.k.eq.Nr) then
               GPON(i,j,k,bi,bj)=GPON(i,j,k,bi,bj)
     &            +wn_sink*(pnup)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPOS(i,j,k,bi,bj)=GPOS(i,j,k,bi,bj)
     &            +ws_sink*(psup)/(drf(k)*hFacC(i,j,k,bi,bj))
c               GPOF(i,j,k,bi,bj)=GPOF(i,j,k,bi,bj)
c     &            +wf_sink*(pfup)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPOC(i,j,k,bi,bj)=GPOC(i,j,k,bi,bj)
     &            +wc_sink*(pcup)/(drf(k)*hFacC(i,j,k,bi,bj))
c               GPIC(i,j,k,bi,bj)=GPIC(i,j,k,bi,bj)
c     &            +wpic_sink*(pccup)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPHY2(i,j,k,bi,bj)=GPHY2(i,j,k,bi,bj)
     &            +wphy2_sink*(py2up)/(drf(k)*hFacC(i,j,k,bi,bj))
             else 
               GPON(i,j,k,bi,bj)=GPON(i,j,k,bi,bj)
     &            +wn_sink*(pnup-pn)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPOS(i,j,k,bi,bj)=GPOS(i,j,k,bi,bj)
     &            +ws_sink*(psup-ps)/(drf(k)*hFacC(i,j,k,bi,bj))
c               GPOF(i,j,k,bi,bj)=GPOF(i,j,k,bi,bj)
c     &            +wf_sink*(pfup-pf)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPOC(i,j,k,bi,bj)=GPOC(i,j,k,bi,bj)
     &            +wc_sink*(pcup-pc)/(drf(k)*hFacC(i,j,k,bi,bj))
c               GPIC(i,j,k,bi,bj)=GPIC(i,j,k,bi,bj)
c     &            +wpic_sink*(pccup-pcc)/(drf(k)*hFacC(i,j,k,bi,bj))
               GPHY2(i,j,k,bi,bj)=GPHY2(i,j,k,bi,bj)
     &            +wphy2_sink*(py2up-py2m)/(drf(k)*hFacC(i,j,k,bi,bj))
             endif
c remineralization of pom
             GPON(i,j,k,bi,bj)=GPON(i,j,k,bi,bj) -reminn_new*pn
             GPOS(i,j,k,bi,bj)=GPOS(i,j,k,bi,bj) -remins*ps
c             GPOF(i,j,k,bi,bj)=GPOF(i,j,k,bi,bj) -reminf*pf
             GPOC(i,j,k,bi,bj)=GPOC(i,j,k,bi,bj) -reminc_new*pc
c             GPIC(i,j,k,bi,bj)=GPIC(i,j,k,bi,bj) -dissc*pcc
c
             GNUT(i,j,k,bi,bj)=GNUT(i,j,k,bi,bj) +reminn_new*pn
             GSIL(i,j,k,bi,bj)=GSIL(i,j,k,bi,bj) +remins*ps
c             GFeT(i,j,k,bi,bj)=GFeT(i,j,k,bi,bj) +reminf*pf
             GDIC(i,j,k,bi,bj)=GDIC(i,j,k,bi,bj) +reminc_new*pc
c     &                                           +dissc*pcc
             remdic(i,j,k,bi,bj) = remdic(i,j,k,bi,bj) 
     &                             +reminc_new*pc
c +dissc*pcc
c             picalk(i,j,k,bi,bj) = -2.d0*GPIC(i,j,k,bi,bj)

             GALK(i,j,k,bi,bj)=-Rnp*GNut(i,j,k,bi,bj)
c     &                        -2.d0*GPIC(i,j,k,bi,bj) 
           
             if (k.eq.1)  GALK(i,j,k,bi,bj)=
     &              GALK(i,j,k,bi,bj)
     &              +SURA(i,j,bi,bj)
c #ifdef ALLOW_RIVERS
c              if (k.eq.1) then
c              print*,'GALK ',GALK(i,j,k,bi,bj)
c              GALK(i,j,k,bi,bj) =
c     &          GALK(i,j,k,bi,bj)+GRALK(i,j,bi,bj)
c         print*,'GDIC before river ',GDIC(i,j,k,bi,bj)
c              GDIC(i,j,k,bi,bj) =
c     &          GDIC(i,j,k,bi,bj)+GRDIC(i,j,bi,bj)
c              print*,'GRDIC+GDIC ',GDIC(i,j,k,bi,bj)
c              endif
c #endif
           
#ifdef ALLOW_DIAGNOSTICS
c            IF ( useDiagnostics ) THEN
             pocrem(i,j,k,bi,bj)  = -reminc_new*pc
             if (k.eq.1)  valk(i,j,bi,bj) = SURA(i,j,bi,bj)
c            ENDIF
#endif
     


#ifdef ALLOW_TIMEAVE
c save averages
             Phygrow1ave(i,j,k,bi,bj)=Phygrow1ave(i,j,k,bi,bj)+
     &                             mu*py*dTtracerLev(k)
     &                       /float(nsubtime)
             Phygrow2ave(i,j,k,bi,bj)=Phygrow2ave(i,j,k,bi,bj)+
     &                             mu2*py2*dTtracerLev(k)
     &                       /float(nsubtime)
             Zoograzave(i,j,k,bi,bj)=Zoograzave(i,j,k,bi,bj)+
     &                      (gampn*graz*zo +gampn2*graz2*zo)*
     &                       dTtracerLev(k)/float(nsubtime)
             Parave(i,j,k,bi,bj)=Parave(i,j,k,bi,bj)+
     &                           PAR(i,j,k,bi,bj)*
     &                           dTtracerLev(k)/float(nsubtime)
             Chlave(i,j,k,bi,bj)=Chlave(i,j,k,bi,bj)+
     &                           CHL(i,j,k)
     &                          * dTtracerLev(k)/float(nsubtime)
             Quaggagrazave(i,j,k,bi,bj)=Quaggagrazave(i,j,k,bi,bj)+
     &                      (quaggagraze+quaggagraze2)*
     &                      dTtracerLev(k)/float(nsubtime)
             if (k.eq.1) then
               FluxCO2ave(i,j,bi,bj)= FluxCO2ave(i,j,bi,bj)+
     &               FluxCO2(i,j,bi,bj)*dTtracerLev(k)/float(nsubtime)
               FluxO2ave(i,j,bi,bj)= FluxO2ave(i,j,bi,bj)+
     &               FluxO2(i,j,bi,bj)*dTtracerLev(k)/float(nsubtime)
               pHave(i,j,bi,bj)=pHave(i,j,bi,bj)+
     &               pH(i,j,bi,bj)*dTtracerLev(k)/float(nsubtime)
               pCO2ave(i,j,bi,bj)=pCO2ave(i,j,bi,bj)+
     &               pCO2(i,j,bi,bj)*dTtracerLev(k)/float(nsubtime)
               virtCO2fluxave(i,j,bi,bj)=virtCO2fluxave(i,j,bi,bj)+
     &           VirtualFlux(i,j,bi,bj)*dTtracerLev(k)/float(nsubtime)
              endif
#endif
#ifdef ALLOW_DIAGNOSTICS
                 IF ( useDiagnostics ) THEN
                   diagName = 'VALK    '
                   if (k.eq.1) CALL DIAGNOSTICS_FILL(valk,diagName,
     &                                         k,1, 2,bi,bj,myThid)
                 ENDIF
#endif           
         endif        

            
c end if hFac>0
         enddo
c end layer loop
c
c step forward in time
        do  k=1,Nr
         Nut(i,j,k,bi,bj) = Nut(i,j,k,bi,bj) +
     &                      GNut(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Phy(i,j,k,bi,bj) = Phy(i,j,k,bi,bj)+
     &                       GPhy(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Phy2(i,j,k,bi,bj) = Phy2(i,j,k,bi,bj)+
     &                       GPhy2(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Zoo(i,j,k,bi,bj) = Zoo(i,j,k,bi,bj) +
     &                      GZoo(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         DON(i,j,k,bi,bj) = DON(i,j,k,bi,bj) +
     &                      GDON(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Sil(i,j,k,bi,bj) = Sil(i,j,k,bi,bj) +
     &                      GSil(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
c         FeT(i,j,k,bi,bj) = FeT(i,j,k,bi,bj) +
c     &                      GFeT(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)
cQQQ - add 31/08/03
c        FeT(i,j,k,bi,bj) = min(FeT(i,j,k,bi,bj),ligand_tot)
cQQQ - end add 
c         Zfe(i,j,k,bi,bj) = Zfe(i,j,k,bi,bj)+
c     &                      GZfe(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)
c         DOF(i,j,k,bi,bj) = DOF(i,j,k,bi,bj) +
c     &                      GDOF(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)

         DIC(i,j,k,bi,bj) = DIC(i,j,k,bi,bj) +
     &                      GDIC(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         ZoC(i,j,k,bi,bj) = ZoC(i,j,k,bi,bj) +
     &                      GZoC(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         DOC(i,j,k,bi,bj) = DOC(i,j,k,bi,bj) +
     &                      GDOC(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Oxy(i,j,k,bi,bj) = Oxy(i,j,k,bi,bj) +
     &                      GOxy(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         Alk(i,j,k,bi,bj) = Alk(i,j,k,bi,bj) +
     &                      GAlk(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
c for explicit sinking particulate matter
         PON(i,j,k,bi,bj) = PON(i,j,k,bi,bj) +
     &                      GPON(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
         POS(i,j,k,bi,bj) = POS(i,j,k,bi,bj) +
     &                      GPOS(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
c         POF(i,j,k,bi,bj) = POF(i,j,k,bi,bj) +
c     &                      GPOF(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)
         POC(i,j,k,bi,bj) = POC(i,j,k,bi,bj) +
     &                      GPOC(i,j,k,bi,bj)*dTtracerLev(k)
     &                       /float(nsubtime)
c         TDC(i,j,k,bi,bj) = TDC(i,j,k,bi,bj) +
c     &                      GTDC(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)
c         RDC(i,j,k,bi,bj) = RDC(i,j,k,bi,bj) +
c     &                      GRDC(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)
c         PIC(i,j,k,bi,bj) = PIC(i,j,k,bi,bj) +
c     &                      GPIC(i,j,k,bi,bj)*dTtracerLev(k)
c     &                       /float(nsubtime)


        enddo
       ENDDO
       ENDDO
c         print*,'GALK(10,10) ',GALK(10,10,1,1,1)

#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         diagName = 'BIODIC1 '
         CALL DIAGNOSTICS_FILL(biodic1,diagName, 0,23, 2,
     &                           bi, bj, myThid)
         diagName = 'REMDIC  '
         CALL DIAGNOSTICS_FILL(remdic,diagName, 0,23, 2,
     &                           bi,bj, myThid)
         diagName = 'BIODIC2 '
         CALL DIAGNOSTICS_FILL(biodic2,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'docp1m'
         CALL DIAGNOSTICS_FILL(docp1m,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'docp2m'
         CALL DIAGNOSTICS_FILL(docp2m,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'doczmort'
         CALL DIAGNOSTICS_FILL(doczmort,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'docslop '
         CALL DIAGNOSTICS_FILL(docslop,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'docrem  '
         CALL DIAGNOSTICS_FILL(docrem,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'pocp1m'
         CALL DIAGNOSTICS_FILL(pocp1m,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'pocp2m'
         CALL DIAGNOSTICS_FILL(pocp2m,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'poczmort'
         CALL DIAGNOSTICS_FILL(poczmort,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'pocslop '
         CALL DIAGNOSTICS_FILL(pocslop,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'pocsink '
         CALL DIAGNOSTICS_FILL(pocsink,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'pocrem  '
         CALL DIAGNOSTICS_FILL(pocrem,diagName, 0,23, 2,
     &                          bi,bj, myThid)
        ENDIF
#endif



#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
c         diagName = 'PICALK  '
c         CALL DIAGNOSTICS_FILL(picalk,diagName, 0,23, 2,
c     &                          bi,bj, myThid)
         diagName = 'PHY1ALK '
         CALL DIAGNOSTICS_FILL(phy1alk,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'PHY2ALK '
         CALL DIAGNOSTICS_FILL(phy2alk,diagName, 0,23, 2,
     &                          bi,bj, myThid)
         diagName = 'REMALK  '
         CALL DIAGNOSTICS_FILL(remalk,diagName, 0,23, 2,
     &                          bi,bj, myThid)
        ENDIF
#endif




c find free iron
c         call fe_chem(bi,bj,iMin,iMax,jMin,jMax, FeT, freefe,
c     &                myIter, mythid)

c
#ifdef ALLOW_TIMEAVE
c save averages
       do k=1,nR
         npzd_timeave(k,bi,bj)=npzd_timeave(k,bi,bj)
     &                         +dTtracerLev(k)
     &                       /float(nsubtime)
       enddo
#endif
c
c -----------------------------------------------------
       ENDDO
c -----------------------------------------------------
c end of bio-chemical time loop
c 
        RETURN
        END
#endif  /*NPZDCAR*/
#endif  /*ALLOW_PTRACERS*/

C============================================================================
