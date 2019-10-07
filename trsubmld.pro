FUNCTION trsubmld,  mld, mldm1, mldp1, tr, trm1
@common
;computing d(MLD)/dt and center it in time
; Subduciton is equal to -(trcmld2-trcmld1)=-tr(n)*(mld(n+1)-mld(n))
  trcmld1 = fltarr(jpi, jpj)    ;  =mld(n)*tr(n)
  trcmld2 = fltarr(jpi, jpj)    ;  =mld(n+1)*tr(n)
  mld1 = (mld+mldm1)*.5         ;  1er du mois en cours
  mld2 = (mld+mldp1)*.5         ;  1er du mois suivant
  trcentre = (tr+trm1)*.5       ;

; Mixed layer content: all levels above MLD + prorata of the level in
; which the ml is located
  kmldinf1 = fltarr(jpi, jpj)
  kmldinf2 = fltarr(jpi, jpj)
     
  FOR  i = 0, jpi-1 DO BEGIN 
    FOR  j = 0, jpj-1 DO BEGIN 
      a = where (gdepw(i, j, *) LE mld1(i, j))
      kmldinf1(i, j) = max(a)
       
      FOR  k = 0, kmldinf1(i, j)-1 DO BEGIN 
        trcmld1(i, j) = trcmld1(i, j)+e1t(i, j)*e2t(i, j)*e3t(i, j, k)*trcentre(i, j, k)*tmask(i, j, k)
      ENDFOR
      
      b = where (gdepw(i, j, *) LE mld2(i, j))
      kmldinf2(i, j) = max(b)
      
      FOR  k = 0, kmldinf2(i, j)-1 DO BEGIN 
        trcmld2(i, j) = trcmld2(i, j)+e1t(i, j)*e2t(i, j)*e3t(i, j, k)*trcentre(i, j, k)*tmask(i, j, k)
      ENDFOR      
      
      trcmld1(i, j) =   trcmld1(i, j)  +e1t(i, j)*e2t(i, j)*(mld1(i, j)-gdepw(i, j, kmldinf1(i, j)))* trcentre(i, j, kmldinf1(i, j))*tmask(i, j, kmldinf1(i, j))
      trcmld2(i, j) =   trcmld2(i, j)  +e1t(i, j)*e2t(i, j)*(mld2(i, j)-gdepw(i, j, kmldinf2(i, j)))* trcentre(i, j, kmldinf2(i, j))*tmask(i, j, kmldinf2(i, j))        
      
    ENDFOR 
  ENDFOR 

zz = trcmld1-trcmld2
return, zz

END 
