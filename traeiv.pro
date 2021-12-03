FUNCTION traeiv,  prd, pn2, un, vn, wn, wslpi,wslpj,uslp,vslp, umask,  vmask, mbathy, ff
@common
jpjm1=jpj-1
jpim1=jpi-1
jpkm1=jpk-1

zu_eiv = fltarr(jpi,jpj)       
zv_eiv = fltarr(jpi,jpj)     
zw_eiv = fltarr(jpi,jpj)     
puneiv = fltarr(jpi,jpj, jpk)      
pvneiv=fltarr(jpi,jpj,jpk)
pwneiv=fltarr(jpi,jpj,jpk)

pun = fltarr(jpi, jpj, jpk)
pvn = fltarr(jpi, jpj, jpk)
pwn = fltarr(jpi, jpj, jpk)

FOR  jk = 0, jpk-1 DO BEGIN 
  pun(*, *, jk) = e2u * e3t(*, *, jk) *un(*, *, jk)     ;           %    ! eulerian transport only
  pvn(*, *, jk) = e1v * e3t(*, *, jk) *vn(*, *, jk)     ;           %    ! eulerian transport only
  pwn(*, *, jk) = e1t * e2t *wn(*, *, jk)               ;           %    ! eulerian transport only
ENDFOR 

print,  'calcul des coeff de diffusion eiv'
aei = ldfeiv(pn2, wslpi, wslpj, mbathy, ff)
fsaeiu = aei(*, *, *, 0)
fsaeiv = aei(*, *, *, 1)
fsaeiw = aei(*, *, *, 2)
 
 FOR  jk = 0, jpkm1-1 DO BEGIN                           
   FOR  jj = 0, jpjm1-1 DO BEGIN 
     FOR  ji = 0, jpim1-1  DO BEGIN 
       zuwk = ( wslpi(ji, jj, jk  ) + wslpi(ji, jj, jk+1  ) ) * fsaeiu(ji, jj, jk  ) * umask[ji, jj, jk  ] ;
       zuwk1 = ( wslpi(ji, jj, jk+1) + wslpi(ji, jj, jk+1) ) * fsaeiu(ji, jj, jk+1) * umask[ji,jj,jk+1] ;
       zvwk = ( wslpj(ji, jj, jk  ) + wslpj(ji,jj+1,jk  ) ) * fsaeiv(ji, jj, jk  ) * vmask[ji, jj, jk  ]  ;
       zvwk1 = ( wslpj(ji,jj,jk+1) + wslpj(ji, jj+1, jk+1) ) * fsaeiv(ji,jj,jk+1) * vmask[ji,jj,jk+1] ;
       
       zu_eiv(ji, jj) = 0.5 * umask[ji, jj, jk] * ( zuwk - zuwk1 )      
       zv_eiv(ji, jj) = 0.5 * vmask[ji, jj, jk] * ( zvwk - zvwk1 )      
       
       puneiv(ji, jj, jk) = pun(ji, jj, jk) + e2u(ji, jj) * zu_eiv(ji, jj)   
       pvneiv(ji, jj, jk) = pvn(ji, jj, jk) + e1v(ji, jj) * zv_eiv(ji, jj)   
     ENDFOR
   ENDFOR 
   
   IF ( jk GE 1 ) THEN BEGIN 
     FOR  jj = 1, jpjm1-1 DO BEGIN 
       FOR  ji = 1, jpim1-1  DO BEGIN
         zuwi  = ( wslpi(ji, jj, jk)+wslpi(ji-1, jj, jk) ) * fsaeiu(ji-1, jj, jk) * e2u( ji-1, jj) * umask[ji-1, jj, jk] ;
         zuwi1 = ( wslpi(ji, jj, jk)+wslpi(ji+1, jj, jk) ) * fsaeiu(ji, jj, jk) * e2u(ji, jj ) * umask[ji, jj, jk]   
         zvwj  = ( wslpj(ji, jj, jk)+wslpj(ji,jj-1,jk) ) * fsaeiv(ji,jj-1,jk) * e1v(ji, jj-1) * vmask[ji,jj-1,jk] 
         zvwj1 = ( wslpj(ji, jj, jk)+wslpj(ji, jj+1, jk) ) * fsaeiv(ji, jj, jk) * e1v( ji, jj) * vmask[ji, jj, jk]   
         
         zw_eiv(ji, jj) = - 0.5 * tmask(ji, jj, jk) * ( zuwi1 - zuwi + zvwj1 - zvwj ) 


         pwneiv(ji, jj, jk) = pwn(ji, jj, jk) + zw_eiv(ji, jj) ;
       ENDFOR 
     ENDFOR 
   ENDIF 
 ENDFOR 

 FOR  jk = 0, jpk-1 DO BEGIN 
   puneiv(*, *, jk) = 1./(e2u * e3t(*, *, jk)) *puneiv(*, *, jk)    ;           %    ! eulerian transport only
   pvneiv(*, *, jk) = 1./(e1v * e3t(*, *, jk)) *pvneiv(*, *, jk)    ;           %    ! eulerian transport only
   pwneiv(*, *, jk) = 1./(e1t * e2t)*pwneiv(*, *, jk)               ;           %    ! eulerian transport only
 ENDFOR 

tra = fltarr(jpi, jpj, jpk, 3)
tra(*, *, *, 0) = puneiv
tra(*, *, *, 1) = pvneiv
tra(*, *, *, 2) = pwneiv

return, tra

END 



