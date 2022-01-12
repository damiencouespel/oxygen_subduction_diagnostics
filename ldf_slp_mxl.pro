FUNCTION ldf_slp_mxl, prd, pn2, p_gru, p_grv, p_dzr, nmln, umask, vmask
@common
;function[uslpml,vslpml,wslpiml,wslpjml,omlmask]=ldf_slp_mxl( prd, pn2, p_gru, p_grv, p_dzr, nmln )
;global e1u e2v e3u e3v umask vmask e1t e2t tmask e3w
;%       !!----------------------------------------------------------------------
;%       !!                  ***  ROUTINE ldf_slp_mxl  ***
;%       !!
;%       !! ** Purpose :   Compute the slopes of iso-neutral surface just below
;%       !!              the mixed layer.
;%       !!
;%       !! ** Method  :   The slope in the i-direction is computed at u- & w-points
;%       !!              (uslpml, wslpiml) and the slope in the j-direction is computed
;%       !!              at v- and w-points (vslpml, wslpjml) with the same bounds as
;%       !!              in ldf_slp.
;%       !!
;%       !! ** Action  :   uslpml, wslpiml :  i- &  j-slopes of neutral surfaces
;%       !!                vslpml, wslpjml    just below the mixed layer
;%       !!                omlmask         :  mixed layer mask
;%       !!----------------------------------------------------------------------
;%       REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   prd            ! in situ density
;%       REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pn2            ! Brunt-Vaisala frequency (locally ref.)
;%       REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_gru, p_grv   ! i- & j-gradient of density (u- & v-pts)
;%       REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_dzr          ! z-gradient of density      (T-point)
;%       !!
;%       INTEGER  ::   ji , jj , jk         ! dummy loop indices
;%       INTEGER  ::   iku, ikv, ik, ikm1   ! local integers
;%       REAL(wp) ::   zeps, zm1_g, zm1_2g            ! local scalars
;%       REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
;%       REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
;%       REAL(wp) ::   zck, zfk,      zbw             !   -      -
;%       !!----------------------------------------------------------------------

print, 'calcule la pente des isoneutres au niveau de la ml'
jpjm1=jpj-1
jpim1=jpi-1
jpkm1=jpk-1
fse3u=e3t
fse3v=e3t
fse3w=e3w

zeps   =  1.e-20     ; ! = =   Local constant initialization   = = !
grav = 9.80664999999999942      ; %gravité  m/s^2
zm1_g  = -1.0/ grav
zm1_2g = -0.5/ grav 

uslpml = fltarr(jpi, jpj)
vslpml = fltarr(jpi, jpj)
wslpiml = fltarr(jpi, jpj)
wslpjml = fltarr(jpi, jpj)

uslpml (0, *) = 0
uslpml (jpi-1, *) = 0
vslpml (0, *) = 0
vslpml (jpi-1, *) = 0
wslpiml(0, *) = 0
wslpiml(jpi-1, *) = 0
wslpjml(0, *) = 0
wslpjml(jpi-1, *) = 0



;%       ! Slopes of isopycnal surfaces just before bottom of mixed layer
;%       ! --------------------------------------------------------------
;%       ! The slope are computed as in the 3D case.
;%       ! A key point here is the definition of the mixed layer at u- and v-points.
;%       ! It is assumed to be the maximum of the two neighbouring T-point mixed layer depth.
;%       ! Otherwise, a n2 value inside the mixed layer can be involved in the computation
;%       ! of the slope, resulting in a too steep diagnosed slope and thus a spurious eddy
;%       ! induce velocity field near the base of the mixed layer.
;%       !-----------------------------------------------------------------------
;%       !

FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN 

;%             !                    !==   Slope at u- & v-points just below the Mixed Layer   ==!
;%             !
;%             !                          !- vertical density gradient for u- and v-slopes (from dzr at T-point)
    iku = min([  max([0, nmln(ji, jj),  nmln(ji+1, jj)]),   jpkm1-1  ])       ;%   
    ikv = min([  max([0, nmln(ji, jj),  nmln(ji, jj+1)]),   jpkm1-1  ])         ;   
    zbu = 0.5* ( p_dzr(ji, jj, iku) + p_dzr(ji+1, jj, iku) )            ;
    zbv = 0.5* ( p_dzr(ji, jj, ikv) + p_dzr(ji, jj+1, ikv) )            ;
;          %  !                          !- horizontal density gradient at u- & v-points
    zau = p_gru(ji, jj, iku) / e1u(ji, jj)         ;
    zav = p_grv(ji, jj, ikv) / e2v(ji, jj)         ;
;%            !                          !- bound the slopes: abs(zw.)<= 1/100 and zb..<0
;%             !                                kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
    zbu = min([  zbu,  -100*abs(zau),  -7e+3/fse3u(ji, jj, iku)*abs(zau)])
    zbv = min([  zbv,  -100*abs(zav),  -7e+3/fse3v(ji, jj, ikv)*abs(zav)])
;%            !                          !- Slope at u- & v-points (uslpml, vslpml)
    uslpml(ji, jj) = zau / ( zbu - zeps ) * umask[ji, jj, iku]
    vslpml(ji, jj) = zav / ( zbv - zeps ) * vmask[ji, jj, ikv]
;%             !
;%             !                    !==   i- & j-slopes at w-points just below the Mixed Layer   ==!
;%             !
    ik   = min([ nmln(ji, jj)+1, jpk-1] )
    ikm1 = max([ 0, ik-1] )
;%            !                          !- vertical density gradient for w-slope (from N^2)
    zbw = zm1_2g * pn2 (ji, jj, ik) * ( prd (ji, jj, ik) + prd (ji, jj, ikm1) + 2 )
;%            !                          !- horizontal density i- & j-gradient at w-points
    zci = max([umask[ji-1, jj, ik]+umask[ji, jj, ik]+ umask[ji-1, jj, ikm1]+umask[ji, jj, ikm1], zeps]) * e1t(ji, jj)
    zcj = max([vmask[ji, jj-1, ik]+vmask[ji, jj, ik]+ vmask[ji, jj-1, ikm1]+vmask[ji, jj, ikm1], zeps]) * e2t(ji, jj)
    zai = (p_gru(ji-1, jj, ik)+p_gru(ji, jj, ik)+p_gru(ji-1, jj, ikm1)+p_gru(ji, jj, ikm1))/zci*tmask(ji, jj, ik)
    zaj = (p_grv(ji, jj-1, ik)+p_grv(ji, jj, ik)+p_grv(ji, jj-1, ikm1)+p_grv(ji, jj, ikm1))/zcj*tmask(ji, jj, ik)
;%             !                          !- bound the slopes: abs(zw.)<= 1/100 and zb..<0.
;%             !                             kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
    zbi = min([zbw,  -100*abs(zai), -7.e+3/fse3w(ji, jj, ik)*abs(zai) ])         ;
    zbj = min([zbw,  -100*abs(zaj), -7.e+3/fse3w(ji, jj, ik)*abs(zaj) ])         ;
;%           !                          !- i- & j-slope at w-points (wslpiml, wslpjml)
    wslpiml(ji, jj) = zai / ( zbi - zeps ) * tmask (ji, jj, ik)        ;
    wslpjml(ji, jj) = zaj / ( zbj - zeps ) * tmask (ji, jj, ik)        ;
  ENDFOR
ENDFOR

slpml = fltarr(jpi, jpj, 4)
slpml(*, *, 0) = uslpml
slpml(*, *, 1) = vslpml
slpml(*, *, 2) = wslpiml
slpml(*, *, 3) = wslpjml


return, slpml

END
