FUNCTION ldfeiv, rn2b,wslpi,wslpj, mbathy, ff
@common

print, 'calcule les coefficients de mélange'

omega=0.729211508304606178E-04 ; %s-1
rad=0.174532925199432955E-01 ; %deg to rad

;%           ! 0. Local initialization
;%       ! -----------------------

jpjm1=jpj-1
jpim1=jpi-1
jpkm1=jpk-1

aeiw=fltarr(jpi,jpj)
fse3w=e3w
zn   =  fltarr(jpi,jpj)                                                              
zhw  =  fltarr(jpi,jpj)                                                            
zah   =  fltarr(jpi,jpj)                                                         
zross =  fltarr(jpi,jpj)                                                                 



FOR  jk = 0, jpk-1 DO BEGIN 
  FOR  jj = 1,  jpjm1-1 DO BEGIN 
    FOR  ji = 1, jpim1-1 DO BEGIN 
;                   ! Take the max of N^2 and zero then take the vertical sum
;                   ! of the square root of the resulting N^2 ( required to compute
;                   ! internal Rossby radius Ro = .5 * sum_jpk(N) / f
      zn2 = max([rn2b(ji, jj, jk), 0])
      zn(ji, jj) = zn(ji, jj)+sqrt(zn2)*fse3w(ji, jj, jk)
;                   ! Compute elements required for the inverse time scale of baroclinic
;                   ! eddies using the isopycnal slopes calculated in ldfslp.F :
;                   ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
     ze3w = fse3w(ji, jj, jk)*tmask(ji, jj, jk)
      zah(ji, jj) = zah(ji, jj)+zn2*(wslpi(ji, jj, jk)*wslpi(ji, jj, jk)+wslpj(ji, jj, jk)*wslpj(ji, jj, jk))*ze3w
      zhw(ji, jj) = zhw(ji, jj)+ze3w
    ENDFOR
  ENDFOR
ENDFOR
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1  DO BEGIN                                  ;% ! vector opt.
    zfw = max([abs(2*omega*sin(rad*gphit(ji, jj))), 1e-10] )
;            % Rossby radius at w-point taken < 40km and  > 2km
    zross(ji, jj) = max( [min( [.4*zn(ji, jj)/zfw, 40e3 ]), 2e3] )
;            % Compute aeiw by multiplying Ro^2 and T^-1
    aeiw(ji, jj) = zross(ji, jj)*zross(ji, jj)*sqrt(zah(ji, jj)/zhw(ji, jj))*tmask(ji, jj, 1)
  ENDFOR 
ENDFOR 

;       %   ! ORCA R2
mbkt = fltarr(jpi, jpj)
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1  DO BEGIN ;% ! vector opt.
;  %             ! Take the minimum between aeiw and 1000 m2/s over shelves (depth shallower than 650 m)
    mbkt(ji, jj) = max([mbathy(ji, jj), 1])
    IF ( mbkt(ji, jj) LE  20 )  THEN  aeiw(ji, jj) = min( [aeiw(ji, jj), 1000 ])
  ENDFOR
ENDFOR
      

; %     ! Decrease the coefficient in the tropics (20N-20S)
zf20 = 2 * omega * sin( rad * 20)
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN  ;% ! vector opt.
    aeiw(ji, jj) = min([1, abs(ff(ji, jj)/zf20)])*aeiw(ji, jj)
  ENDFOR 
ENDFOR 
aeiw2d=aeiw
  
; %     ! Average the diffusive coefficient at u- v- points
aeiv2d = fltarr(jpi, jpj)
aeiu2d = fltarr(jpi, jpj)
FOR  jj = 1, jpjm1-1 DO BEGIN 
  FOR  ji = 1, jpim1-1 DO BEGIN                         ;   %! vector opt.
    aeiu2d(ji, jj) = 0.5* ( aeiw(ji, jj) + aeiw(ji+1, jj) ) ;
    aeiv2d(ji, jj) = 0.5* ( aeiw(ji, jj) + aeiw(ji, jj+1) ) ;
  ENDFOR
ENDFOR


print,  'traldf2d, les coefficients sont les mêmes sur toute la colonne deau'
aeiu = fltarr(jpi, jpj, jpk)
aeiv = fltarr(jpi, jpj, jpk)
aeiw = fltarr(jpi, jpj, jpk)

for jk = 0, jpk-1 DO BEGIN 
  aeiu(*, *, jk) = aeiu2d
  aeiv(*, *, jk) = aeiv2d
  aeiw(*, *, jk) = aeiw2d
ENDFOR 

aei = fltarr(jpi, jpj, jpk, 3)
aei(*, *, *, 0)= aeiu
aei(*, *, *, 1)= aeiv
aei(*, *, *, 2)= aeiw
return,  aei
end
