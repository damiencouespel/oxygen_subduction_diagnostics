FUNCTION trsubgm, mld,mldu, mldv, trmldu, trmldv, trmld, prd,pn2,un,vn,wn,wslpi,wslpj,uslp,vslp, umask,  vmask, mbathy, ff, nmlnu,  nmlnv
@common

;compute total velocity including eiv:
tra = traeiv(prd,pn2,un,vn,wn,wslpi,wslpj,uslp,vslp, umask,  vmask, mbathy, ff)

;remove eulerien velocity and keep GM velocity only:
uneiv = tra(*, *, *, 0)-un
vneiv = tra(*, *, *, 1)-vn
wneiv = tra(*, *, *, 2)-wn

wneivmld = fltarr(jpi, jpj)
uneivmld = fltarr(jpi, jpj)
vneivmld = fltarr(jpi, jpj)

FOR  j = 0, jpj-1 DO BEGIN  
  FOR  i = 0, jpi-1 DO BEGIN  
    wneivmld(i, j) = interpol(wneiv(i, j, *), gdepw(i, j, *), mld(i, j))      ; %wneiv at MLD base
    uneivmld(i, j) = interpol(uneiv(i, j, *), gdept(i, j, *), mldu(i, j))     ; %uneiv at MLD base
    vneivmld(i, j) = interpol(vneiv(i, j, *), gdept(i, j, *), mldv(i, j))     ; %vneiv at MLD base
  ENDFOR
ENDFOR

tot = trsubhadv(mld,uneivmld,vneivmld,trmldu, trmldv, umask, vmask, nmlnu,  nmlnv) + trsubwadv(wneivmld, trmld)

return,  tot

END 

