FUNCTION trsubzdf, p2dt, ptb,wslpi,wslpj,mldlev, fsavs, aht0
@common
;%       !!----------------------------------------------------------------------
;%       !!                  ***  ROUTINE tra_zdf_imp  ***
;%       !!
;%       !! ** Purpose :   Compute the after tracer through a implicit computation
;%;       !!     of the vertical tracer diffusion (including the vertical component
;%;       !!     of lateral mixing (only for 2nd order operator, for fourth order
;%       !!     it is already computed and add to the general trend in traldf)
;%       !!
;%       !! ** Method  :  The vertical diffusion of the tracer t  is given by:
;%;       !!                  difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
;%       !!      It is computed using a backward time scheme (t=ta).
;%       !!      If lk_zdfddm=T, use avs for salinity or for passive tracers
;%       !!      Surface and bottom boundary conditions: no diffusive flux on
;%       !!      both tracers (bottom, applied through the masked field avt).
;%       !!      If iso-neutral mixing, add to avt the contribution due to lateral mixing.
;%       !!
;%       !! ** Action  : - pta  becomes the after tracer
;%       !!---------------------------------------------------------------------
      
;%       !
;%       INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
;%       CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
;%       INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
;%       REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt     ! vertical profile of tracer time-step
;%       REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb      ! before and now tracer fields
;%       REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta      ! tracer trend
;%       !
;     
; %     !!---------------------------------------------------------------------

;ptb = ptb*tmask

fse3t_n = e3t                 ;
fse3w = e3w                   ;

jpjm1 = jpj-1
jpim1=jpi-1
jpkm1=jpk-1

fsahtw = aht0*(fltarr(jpi, jpj, jpk) +1.)      
pta = fltarr(jpi, jpj, jpk)              
sub=fltarr(jpi,jpj)

zwd=fltarr(jpi,jpj,jpk)
zwi=zwd
zws=zwd
zwt=zwd
;%                                                     ! ============= !
;%          !                                            ! ============= !
;%          !
;%          !  Matrix construction
;%          ! --------------------
;%          ! Build matrix if temperature or salinity (only in double diffusion case) or first passive tracer
;%          !
;%          
            
zwt(*, *, 1:jpk-1) = fsavs(*, *, 1:jpk-1) 
                    
;%          %   ! isoneutral diffusion: add the contribution
;%             
FOR  jk = 1, jpkm1-1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1 DO BEGIN  
      zwt(ji,jj,jk) = zwt(ji,jj,jk)+fsahtw(ji,jj,jk)*(wslpi(ji,jj,jk)*wslpi(ji,jj,jk)+wslpj(ji,jj,jk)*wslpj(ji,jj,jk)) ;
    ENDFOR
  ENDFOR
ENDFOR

;           % ! Diagonal, lower (i), upper (s)  (including the bottom boundary condition since avt is masked)
FOR  jk = 0, jpkm1-1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN 
      ze3ta =  1.                                                                                   ; 
      ze3tn =  fse3t_n(ji, jj, jk)                                                                 ; %  ! now   scale factor at T-point
      zwi(ji, jj, jk) = - p2dt * zwt(ji, jj, jk  ) / ( ze3tn * fse3w(ji, jj, jk  ) )               ;
      zws(ji, jj, jk) = - p2dt * zwt(ji, jj, jk+1) / ( ze3tn * fse3w(ji, jj, jk+1) )               ;
      zwd(ji, jj, jk) = (ze3ta - zwi(ji, jj, jk) - zws(ji, jj, jk))                                ;
    ENDFOR
  ENDFOR
ENDFOR
;%             !
;%             !! Matrix inversion from the first level
;%             !!----------------------------------------------------------------------
;%             !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
;%             !
;%;;             !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
;%             !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
;%             !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
;%             !        (        ...               )( ...  ) ( ...  )
;%             !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
;%             !
;%             !   m is decomposed in the product of an upper and lower triangular matrix.
;%             !   The 3 diagonal terms are in 3d arrays: zwd, zws, zwi.
;%             !   Suffices i,s and d indicate "inferior" (below diagonal), diagonal
;%             !   and "superior" (above diagonal) components of the tridiagonal system.
;%             !   The solution will be in the 4d array pta.
;%             !   The 3d array zwt is used as a work space array.
;%             !   En route to the solution pta is used a to evaluate the rhs and then
;%;             !   used as a work space array: its value is modified.
;%;             !
;%             ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
;%             ! done once for all passive tracers (so included in the IF instruction)
FOR   jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN 
    zwt(ji, jj, 0) = zwd(ji, jj, 0) 
  ENDFOR
ENDFOR
FOR  jk = 1, jpkm1-1  DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN 
      zwt(ji, jj, jk) = (zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk-1) / zwt(ji, jj, jk-1)) ;
    ENDFOR
  ENDFOR
ENDFOR
;; print,  'profil', pta(180, 50, *) ;DC


;%          !
;%          ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
 
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN 
    
    ze3tb = 1.                                                                    ;
    ze3tn = 1.                                                                    ;
    pta(ji, jj, 0) = ze3tb * ptb(ji, jj, 0) + p2dt * ze3tn * pta(ji, jj, 0)      ;
  ENDFOR 
ENDFOR

FOR  jk = 1, jpkm1-1  DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN 
      ze3tb = 1.                                                                                              ;
      ze3tn = 1.                                                                                              ;
      zrhs = ze3tb * ptb(ji, jj, jk) + p2dt * ze3tn * pta(ji, jj, jk)                                        ; %  ! zrhs=right hand side
      pta(ji, jj, jk) = (zrhs - zwi(ji, jj, jk) / zwt(ji, jj, jk-1) * pta(ji, jj, jk-1))   ;
    ENDFOR
  ENDFOR
ENDFOR


 
;%          ! third recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk   (result is the after tracer)
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN 
    pta(ji, jj, jpkm1-1) = pta(ji, jj, jpkm1-1) / zwt(ji, jj, jpkm1-1) * tmask(ji, jj, jpkm1-1) 
  ENDFOR 
ENDFOR 



FOR  jk = jpk-3, 0, -1 DO BEGIN 
  FOR  jj = 1, jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1  DO BEGIN 
      pta(ji, jj, jk) = (pta(ji, jj, jk)-zws(ji, jj, jk)*pta(ji, jj, jk+1))/ zwt(ji, jj, jk) * tmask(ji, jj, jk) 
    ENDFOR 
  ENDFOR 
ENDFOR 


;%          !                                            ! ================= !
;                                       %  !  end tracer loop  !
;%       !                                               ! ================= !
;%       !
;    % tendance entre ptb et pta
dpt = ptb-pta                   ;
FOR  jj = 0, jpj-1 DO BEGIN 
  FOR  ji = 0, jpi-1  DO BEGIN 
    FOR  jk = 0,  mldlev(ji, jj) DO BEGIN 
      sub(ji, jj) = sub(ji, jj)+dpt(ji, jj, jk)*e2t(ji, jj)*e1t(ji, jj)*e3t(ji, jj, jk)*tmask(ji, jj, jk)
    ENDFOR
  ENDFOR 
ENDFOR 

return, sub
                
END 
