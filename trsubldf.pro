FUNCTION trsubldf, ptb, fsah, mld,uslp,vslp,wslpi, wslpj, umask,  vmask, nmln, nmlnu, nmlnv
@common
jpjm1=jpj-1
jpim1=jpi-1
jpkm1=jpk-1
fse3u=e3t
fse3v=e3t
pahtb0=0

; difusion coef
fsahtu=fsah
fsahtv=fsah
fsahtw=fsah
           
zdit =fltarr(jpi, jpj, jpk) 
zdjt = fltarr(jpi, jpj, jpk) 
zdkt = fltarr(jpi, jpj)
zftu =fltarr(jpi, jpj, jpk) 
zftv =fltarr(jpi, jpj, jpk) 
zxsub=fltarr(jpi, jpj)
zysub=fltarr(jpi, jpj)
zzsub=fltarr(jpi, jpj)
      
;%          !
;%          !!----------------------------------------------------------------------
;%          !!   I - masked horizontal derivative
;%          !!----------------------------------------------------------------------
         
       

;%       ! Horizontal tracer gradient
FOR  jk = 0, jpkm1-1 DO BEGIN 
  FOR  jj = 0, jpjm1-1 DO BEGIN 
    FOR  ji = 0, jpim1-1 DO BEGIN 
      zdit(ji, jj, jk) = ( ptb(ji+1, jj, jk) - ptb(ji, jj, jk) ) * umask[ji, jj, jk]
      zdjt(ji, jj, jk) = ( ptb(ji, jj+1, jk) - ptb(ji, jj, jk) ) * vmask[ji, jj, jk]
    ENDFOR 
  ENDFOR
ENDFOR



;%          !!----------------------------------------------------------------------
;%          !!   II - horizontal trend  (full)
;%          !!----------------------------------------------------------------------
;% 
;%          !                                                ! ===============
FOR  jk = 0, jpkm1-1 DO BEGIN    

;           % ! Horizontal slab
;           % !                                             %! ===============
;%             ! 1. Vertical tracer gradient at level jk and jk+1
;%             ! ------------------------------------------------
;%             ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

  zdk1t = ( ptb(*, *, jk) - ptb(*, *, jk+1) )* tmask(*, *, jk+1)
            
  IF ( jk EQ  0 ) THEN  zdkt = zdk1t ELSE zdkt = ( ptb(*, *, jk-1) - ptb(*, *, jk) ) * tmask(*, *, jk)

        

;            ! 2. Horizontal fluxes
;             ! --------------------
  FOR  jj = 0, jpjm1-1 DO BEGIN 
    FOR  ji = 0, jpim1-1 DO BEGIN    
      zabe1 = ( fsahtu) * e2u(ji, jj) * fse3u(ji, jj, jk) / e1u(ji, jj)         
      zabe2 = ( fsahtv) * e1v(ji, jj) * fse3v(ji, jj, jk) / e2v(ji, jj)         
      
      zmsku = 1./max([tmask(ji+1, jj, jk)+tmask(ji, jj, jk+1)+tmask(ji+1, jj, jk+1)+tmask(ji, jj, jk), 1.]) 
      zmskv = 1./max([tmask(ji, jj+1, jk)+tmask(ji, jj, jk+1)+tmask(ji, jj+1, jk+1)+tmask(ji, jj, jk), 1.]) 
      
      zcof1 = - fsahtu * e2u(ji, jj) * uslp(ji, jj, jk) * zmsku          
      zcof2 = - fsahtv * e1v(ji, jj) * vslp(ji, jj, jk) * zmskv          
    
      zftu(ji, jj, jk) = (zabe1 * zdit(ji, jj, jk)+ zcof1 * (  zdkt(ji+1, jj) + zdk1t(ji, jj)+ zdk1t(ji+1, jj) + zdkt(ji, jj)  )  ) * umask[ji, jj, jk]             
      zftv(ji, jj, jk) = (zabe2 * zdjt(ji, jj, jk)+ zcof2 * (  zdkt(ji, jj+1) + zdk1t(ji, jj)+ zdk1t(ji, jj+1) + zdkt(ji, jj)  )  ) * vmask[ji, jj, jk]             
    ENDFOR
  ENDFOR
ENDFOR

diffumld = fltarr(jpi, jpj)
diffvmld = fltarr(jpi, jpj)
       
FOR  ji = 1, jpim1-1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    diffumld(ji, jj) = interpol(zftu(ji, jj, *)/fse3u(ji, jj, *), gdept(ji, jj, *),  mld(ji, jj))*umask[ji, jj, nmlnu(ji, jj)]
    diffvmld(ji, jj) = interpol(zftv(ji, jj, *)/fse3u(ji, jj, *), gdept(ji, jj, *),  mld(ji, jj))*vmask[ji, jj, nmlnv(ji, jj)]
    IF umask[ji, jj, nmlnu(ji, jj)+1] EQ 0 THEN diffumld(ji, jj) = zftu(ji, jj, nmlnu(ji, jj))/fse3u(ji, jj, *)*umask[ji, jj, nmlnu(ji, jj)] 
    IF vmask[ji, jj, nmlnv(ji, jj)+1] EQ 0 THEN diffvmld(ji, jj) = zftv(ji, jj, nmlnv(ji, jj))/fse3u(ji, jj, *)*vmask[ji, jj, nmlnv(ji, jj)]
    zxsub(ji, jj) = -diffumld(ji, jj)*(mld(ji+1, jj)-mld(ji, jj))                                
    zysub(ji, jj) = -diffvmld(ji, jj)*(mld(ji, jj+1)-mld(ji, jj))                                
  ENDFOR
ENDFOR


;%          !!----------------------------------------------------------------------
;%          !!   III - vertical trend of T & S (extra diagonal terms only)
;%          !!----------------------------------------------------------------------
;% 
;%          ! Local constant initialization
;%          ! -----------------------------
ztfw = fltarr(jpi, jpj, jpk)     



;%          ! interior (2=<jk=<jpk-1)
FOR  jk = 1, jpkm1-1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN ; vector opt.

      zcoef0 = - fsahtw* tmask(ji, jj, jk) 
                  
      zmsku = 1./max([umask[ji, jj, jk-1]+umask[ji-1, jj, jk]+umask[ji-1, jj, jk-1]+umask[ji, jj, jk], 1.])
      zmskv = 1./max([vmask[ji, jj, jk-1]+vmask[ji, jj-1, jk]+vmask[ji, jj-1, jk-1]+vmask[ji, jj, jk], 1.])
      
      zcoef3 = zcoef0 * e2t(ji, jj)* zmsku * wslpi(ji,jj,jk)
      zcoef4 = zcoef0 * e1t(ji, jj)* zmskv * wslpj(ji,jj,jk)
                  
      ztfw(ji, jj, jk) = zcoef3*(zdit(ji, jj, jk-1)+zdit(ji-1, jj, jk)+zdit(ji-1, jj, jk-1)+zdit(ji, jj, jk)) +zcoef4*(zdjt(ji, jj, jk-1)+zdjt(ji, jj-1, jk)+zdjt(ji, jj-1, jk-1)+zdjt(ji, jj, jk))
    ENDFOR 
  ENDFOR 
ENDFOR 
FOR  ji= 2, jpi-1 DO BEGIN 
    FOR  jj= 2, jpj-1 DO BEGIN 
      zzsub(ji, jj)=interpol(ztfw(ji, jj, *), gdepw(ji, jj, *),mld(ji, jj))*tmask[ji, jj, nmln(ji, jj)]
      IF tmask(ji, jj, nmln(ji, jj)+1) EQ 0 THEN zzsub(ji, jj) = ztfw(ji, jj, nmln(ji, jj))
    ENDFOR 
  ENDFOR 

zsub = zxsub+zysub+zzsub
return, zsub

END 


