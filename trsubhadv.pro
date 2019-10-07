FUNCTION trsubhadv,  mld, umld, vmld,trmldu, trmldv, umask_i,  vmask_i,  nmlnu,  nmlnv
@common
suadv = fltarr(jpi, jpj)
svadv = fltarr(jpi, jpj)

FOR i = 0,  jpi-2 DO BEGIN
  FOR j = 0,  jpj-1 DO BEGIN 
    suadv(i, j) = -umld(i, j)*trmldu(i, j)* e2u(i, j) * (mld(i+1, j)-mld(i, j)) * umask_i[i, j, nmlnu(i, j)]
  ENDFOR
ENDFOR

FOR i = 0,  jpi-1 DO BEGIN
  FOR j = 0,  jpj-2 DO BEGIN 
    svadv(i, j) = -vmld(i, j)*trmldv(i, j)* e1v(i, j) * (mld(i, j+1)-mld(i, j)) * vmask_i[i, j, nmlnv(i, j)]
  ENDFOR
ENDFOR


sadv = suadv+svadv

return, sadv
END

