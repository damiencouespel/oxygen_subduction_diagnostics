function ldfslp, prd, pn2, nmln, mld, umask, vmask, omlmask
@common
;%       !!----------------------------------------------------------------------
;%       !!
;%       !! ** Purpose :   Compute the slopes of neutral surface (slope of isopycnal
;%       !!              surfaces referenced locally) (ln_traldf_iso=T).
;%       !!
;%       !! ** Method  :   The slope in the i-direction is computed at U- and
;%       !!      W-points (uslp, wslpi) and the slope in the j-direction is
;%       !!      computed at V- and W-points (vslp, wslpj).
;%       !!      They are bounded by 1/100 over the whole ocean, and within the
;%       !!      surface layer they are bounded by the distance to the surface
;%       !!      ( slope<= depth/l  where l is the length scale of horizontal
;%       !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
;%       !!      of 10cm/s)
;%       !!        A horizontal shapiro filter is applied to the slopes
;%       !!        ln_sco=T, s-coordinate, add to the previously computed slopes
;%       !!      the slope of the model level surface.
;%       !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
;%       !!      [slopes already set to zero at level 1, and to zero or the ocean
;%       !!      bottom slope (ln_sco=T) at level jpk in inildf]
;%       !!
;%       !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
;%       !!               of now neutral surfaces at u-, w- and v- w-points, resp.
;%       !!----------------------------------------------------------------------
;%       USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
;%       USE oce     , ONLY:   zgru => ua       , zww => va   ! (ua,va) used as workspace
;%       USE oce     , ONLY:   zgrv => ta       , zwz => sa   ! (ta,sa) used as workspace
;%       USE wrk_nemo, ONLY:   zdzr => wrk_3d_1               ! 3D workspace
;%       !!
;%       INTEGER , INTENT(in)                   ::   kt    ! ocean time-step index
;%       REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   prd   ! in situ density
;%       REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   pn2   ! Brunt-Vaisala frequency (locally ref.)
;%       !!
;%       INTEGER  ::   ji , jj , jk    ! dummy loop indices
;%       INTEGER  ::   ii0, ii1, iku   ! temporary integer
;%       INTEGER  ::   ij0, ij1, ikv   ! temporary integer
;%       REAL(wp) ::   zeps, zm1_g, zm1_2g, z1_16     ! local scalars
;%       REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
;%       REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
;%       REAL(wp) ::   zck, zfk,      zbw             !   -      -
fse3u=e3t
fse3v=e3t
fse3w=e3w
fsdept=gdept
fsdepw=gdepw

wslpi=fltarr(jpi, jpj, jpk)
wslpj=fltarr(jpi, jpj, jpk)
hmlpt = fltarr(jpi, jpj)
FOR  i = 0, jpi-1 DO BEGIN 
  FOR  j = 0, jpj-1 DO BEGIN 
    hmlpt(i, j) = gdept(i, j, nmln(i, j)) ;attention indice nmln au passage a idl
  ENDFOR 
ENDFOR 

hmlp=mld
jpjm1=jpj-1
jpim1=jpi-1
jpkm1=jpk-1
grav=9.80664999999999942   ; %gravité  m/s^2

zeps   =  1.e-20   ;       % !==   Local constant initialization   ==!
z1_16  =  1./16    
zm1_g  = -1./grav  
zm1_2g = -0.5/grav              

zww = fltarr(jpi, jpj, jpk) 
zwz = fltarr(jpi, jpj, jpk) 
zgru = fltarr(jpi, jpj, jpk) 
zgrv = fltarr(jpi, jpj, jpk) 
      
print,  'Calcule la pente des isoneutres'
FOR  jk = 0, jpk-1 DO BEGIN      ;%==   i- & j-gradient of density   ==!
  FOR  jj = 0, jpjm1-1 DO BEGIN 
    FOR  ji = 0, jpim1-1 DO BEGIN                                             ;% vector opt.
      zgru(ji, jj, jk) = umask[ji, jj, jk]*(prd(ji+1, jj, jk)-prd(ji, jj, jk)) ;
      zgrv(ji, jj, jk) = vmask[ji, jj, jk]*(prd(ji, jj+1, jk)-prd(ji, jj, jk)) ;
    ENDFOR 
  ENDFOR 
ENDFOR 
      
zdzr = fltarr(jpi, jpj, jpk) ;   Local vertical density gradient at T-point   == !   (evaluated from N^2)      

FOR  jk = 1, jpkm1-1 DO BEGIN 
;%          !                                ! zdzr = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
;%          !                                !   trick: tmask(ik  )  = 0   =>   all pn2   = 0   =>   zdzr = 0
;%          !                                !    else  tmask(ik+1)  = 0   =>   pn2(ik+1) = 0   =>   zdzr divides by 1
;%          !                                !          umask(ik+1) /= 0   =>   all pn2  /= 0   =>   zdzr divides by 2
;%          !                                ! NB: 1/(tmask+1) = (1-.5*tmask)  substitute a / by a *  ==> faster
  zdzr(*, *, jk) = zm1_g*(prd(*, *, jk)+1)*(pn2(*, *, jk)+pn2(*, *, jk+1))*(1-0.5*tmask(*, *, jk+1))
ENDFOR


;%       !
;%       !                          !==   Slopes just below the mixed layer   ==!
slpml=ldf_slp_mxl( prd,pn2,zgru,zgrv,zdzr,nmln, umask, vmask);        %! output: uslpml, vslpml, wslpiml, wslpjml
uslpml = slpml(*, *, 0)
vslpml = slpml(*, *, 1)
wslpiml = slpml(*, *, 2)
wslpjml = slpml(*, *, 3)


;%             ! I.  slopes at u and v point      | uslp = d/di( prd ) / d/dz( prd )
;%       ! ===========================      | vslp = d/dj( prd ) / d/dz( prd )
;%       !
FOR  jk = 1, jpkm1-1 DO BEGIN                             ;!* Slopes at u and v points
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1 DO BEGIN   ; vector opt.
;%                !                                      ! horizontal and vertical density gradient at u- and v-points
      zau = zgru(ji, jj, jk) / e1u(ji, jj)                                    
      zav = zgrv(ji, jj, jk) / e2v(ji, jj)
      zbu = 0.5 * ( zdzr(ji, jj, jk) + zdzr(ji+1, jj, jk) )
      zbv = 0.5 * ( zdzr(ji, jj, jk) + zdzr(ji, jj+1, jk) )
;%                !                                      ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
;%                !                                      ! + kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
      zbu = min([zbu, -100*abs(zau), -7e+3/fse3u(ji, jj, jk)*abs(zau)])         ;
      zbv = min([zbv, -100*abs(zav), -7e+3/fse3v(ji, jj, jk)*abs(zav)])        ;
 ;%              !                                      ! uslp and vslp output in zwz and zww, resp.
      zfi = max([omlmask(ji, jj, jk), omlmask(ji+1, jj, jk)])
      zfj = max([omlmask(ji, jj, jk), omlmask(ji, jj+1, jk)])
      zwz(ji, jj, jk) = ((1-zfi)*zau/(zbu-zeps)+zfi*uslpml(ji, jj)*0.5*(fsdept(ji+1, jj, jk)+fsdept(ji, jj, jk)-fse3u(ji, jj, 1))/max([hmlpt(ji, jj), hmlpt(ji+1, jj), 5]))*umask[ji, jj, jk]
      zww(ji, jj, jk) = ((1-zfj)*zav/(zbv-zeps)+zfj*vslpml(ji, jj)*0.5*(fsdept(ji, jj+1, jk)+fsdept(ji, jj, jk)-fse3v(ji, jj, 1))/max([hmlpt(ji, jj), hmlpt(ji, jj+1), 5]))*vmask[ji, jj, jk]

    ENDFOR
  ENDFOR
ENDFOR
; %     !                                            !* horizontal Shapiro filter
uslp = fltarr(jpi, jpj, jpk)
vslp = fltarr(jpi, jpj, jpk)
FOR  jk = 1, jpkm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN 
;   FOR  jj = 2:max(1,jpj-3):jpjm1 %                     ! rows jj=2 and =jpjm1 only
    jj = 2-1
    uslp(ji, jj, jk) = z1_16*(zwz(ji-1, jj-1, jk)+ zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+4*zwz(ji, jj, jk))
    vslp(ji, jj, jk) = z1_16*(zww(ji-1, jj-1, jk)+ zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+4*zww(ji, jj, jk))
    jj = jpjm1-1
    uslp(ji, jj, jk) = z1_16*(zwz(ji-1, jj-1, jk)+ zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+4*zwz(ji, jj, jk))
    vslp(ji, jj, jk) = z1_16*(zww(ji-1, jj-1, jk)+ zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+4*zww(ji, jj, jk))
  ENDFOR 


  FOR  jj = 2, jpj-3 DO BEGIN   ;      %     ! other rows
    FOR  ji = 2, jpim1-1 DO BEGIN ;%! vector opt.
      uslp(ji, jj, jk) = z1_16*(zwz(ji-1, jj- 1, jk)+zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+ 4*zwz(ji, jj, jk)) 
      vslp(ji, jj, jk) = z1_16*(zww(ji-1, jj- 1, jk)+zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+ 4*zww(ji, jj, jk)) 
    ENDFOR 
  ENDFOR 
; %        !                                        !* decrease along coastal boundaries
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN  ; vector opt.
      uslp(ji, jj, jk) = uslp(ji, jj, jk)*(umask[ji, jj+1, jk]+umask[ji, jj-1, jk])*0.5*(umask[ji, jj, jk]+umask[ji, jj, jk+1] )*0.5 ;
      vslp(ji, jj, jk) = vslp(ji, jj, jk)*(vmask[ji+1, jj, jk]+vmask[ji-1, jj, jk])*0.5*(vmask[ji, jj, jk]+vmask[ji, jj, jk+1] )*0.5 ;
    ENDFOR 
  ENDFOR 
ENDFOR 
 
;%       ! II.  slopes at w point           | wslpi = mij( d/di( prd ) / d/dz( prd )
;%       ! ===========================      | wslpj = mij( d/dj( prd ) / d/dz( prd )
;%       !
FOR  jk = 1, jpkm1-1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1 DO BEGIN ;%! vector opt.
;%               !                                  !* Local vertical density gradient evaluated from N^2
      zbw = zm1_2g*pn2(ji, jj, jk)*(prd(ji, jj, jk)+prd(ji, jj, jk-1)+2)
;%                !                                  !* Slopes at w point
;%                !                                        ! i- & j-gradient of density at w-points
      zci = max([umask[ji-1, jj, jk]+umask[ji, jj, jk]+umask[ji-1, jj, jk-1]+umask[ji, jj, jk-1], zeps])*e1t(ji, jj)
      zcj = max([vmask[ji, jj-1, jk]+vmask[ji, jj, jk-1]+vmask[ji, jj-1, jk-1]+vmask[ji, jj, jk], zeps])*e2t(ji, jj)
      zai = (zgru(ji-1, jj, jk)+zgru(ji, jj, jk-1)+zgru(ji-1, jj, jk-1)+zgru(ji, jj, jk))/zci*tmask(ji, jj, jk)
      zaj = (zgrv(ji, jj-1, jk)+zgrv(ji, jj, jk-1)+zgrv(ji, jj-1, jk-1)+zgrv(ji, jj, jk))/zcj*tmask(ji, jj, jk)
;%                !                                        ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
;%                !                                        ! + kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
      zbi = min([zbw, -100*abs(zai), -7e+3/fse3w(ji, jj, jk)*abs(zai)])
      zbj = min([zbw, -100*abs(zaj), -7e+3/fse3w(ji, jj, jk)*abs(zaj)])
;%               !                                        ! wslpi and wslpj with ML flattening (output in zwz and zww, resp.)
      zfk = max([omlmask(ji, jj, jk), omlmask(ji, jj, jk-1)] ) ; %   ! zfk=1 in the ML otherwise zfk=0
      zck = fsdepw(ji, jj, jk)/max([hmlp(ji, jj), 10])
      zwz(ji, jj, jk) = (zai/(zbi-zeps)*(1-zfk)+zck*wslpiml(ji, jj)*zfk)*tmask(ji, jj, jk)
      zww(ji, jj, jk) = (zaj/(zbj-zeps)*(1-zfk)+zck*wslpjml(ji, jj)*zfk)*tmask(ji, jj, jk)
    ENDFOR 
  ENDFOR 
ENDFOR 

FOR  jk = 1, jpkm1-1 DO BEGIN 
;         FOR  jj=2:max([1 jpj-3]):jpjm1 %                       ! rows jj=2 and =jpjm1 only
  FOR  ji = 1, jpim1-1 DO BEGIN 
    jj = 2- 1
    wslpi(ji, jj, jk) = (zwz(ji-1, jj-1, jk)+zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+4*zwz(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
    wslpj(ji, jj, jk) = (zww(ji-1, jj-1, jk)+zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+4*zww(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
    jj = jpjm1- 1
    wslpi(ji, jj, jk) = (zwz(ji-1, jj-1, jk)+zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+4*zwz(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
    wslpj(ji, jj, jk) = (zww(ji-1, jj-1, jk)+zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+4*zww(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
  ENDFOR
  FOR  jj = 2, jpj-3 DO BEGIN                              ;! other rows
    FOR  ji = 1, jpim1-1  DO BEGIN  ;! vector opt.
      wslpi(ji, jj, jk) = (zwz(ji-1, jj-1, jk)+zwz(ji+1, jj-1, jk)+zwz(ji-1, jj+1, jk)+zwz(ji+1, jj+1, jk)+2*(zwz(ji, jj-1, jk)+zwz(ji-1, jj, jk)+zwz(ji+1, jj, jk)+zwz(ji, jj+1, jk))+4*zwz(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
      wslpj(ji, jj, jk) = (zww(ji-1, jj-1, jk)+zww(ji+1, jj-1, jk)+zww(ji-1, jj+1, jk)+zww(ji+1, jj+1, jk)+2*(zww(ji, jj-1, jk)+zww(ji-1, jj, jk)+zww(ji+1, jj, jk)+zww(ji, jj+1, jk))+4*zww(ji, jj, jk))*z1_16*tmask(ji, jj, jk)
    ENDFOR 
  ENDFOR 
;%         !                                        !* decrease along coastal boundaries
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN ;! vector opt.
      zck = (umask[ji, jj, jk]+umask[ji-1, jj, jk])*(vmask[ji, jj, jk]+vmask[ji, jj-1, jk])*0.25 
      wslpi(ji, jj, jk) = wslpi(ji, jj, jk)*zck                                                  
      wslpj(ji, jj, jk) = wslpj(ji, jj, jk)*zck                                                  
    ENDFOR
  ENDFOR
ENDFOR
slp = fltarr(jpi, jpj, jpk, 4)

slp(*, *, *, 0) = wslpi
slp(*, *, *, 1) = wslpj
slp(*, *, *, 2) = uslp
slp(*, *, *, 3) = vslp

ind = where(finite(slp) EQ 0)
IF ind(0) NE -1 THEN slp(ind) = 0.

return, slp

END
