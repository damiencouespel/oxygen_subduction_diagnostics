;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;
;;;;;;;;;;        BILAN_OFFLINE         ;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Purpose 
; -------
; compute and save the offline computation of the tracers subduction

; Usage
; -----
; fdir='/data/daclod/DATA_GREENLAND02SV/
; files='piControl2_19900101_20991231_1M'
; nametr='O2'
; - EX1: Sub O2, varying mld
; savefile='sub_offline_' + nametr + '_' + files + '.sav'
; bilan_offline, fdir, files, SAVEFILE=savefile, TRC=nametr
; - EX2: Sub O2sat, varying mld
; savefile='/data/daclod/MY_IDL/DATA_TMP/sub_offline_O2sat_' + files + '.sav'
; bilan_offline, fdir, files, SAVEFILE=savefile, TRC=nametr, /O2SAT
; - EX3: Sub O2, fixed depth
; savefile='sub_offline_' + nametr + '_' + files + '_dep200m.sav'
; bilan_offline, fdir, files, SAVEFILE=savefile, TRC=nametr, DEP=200
 
; WARNINGS
; --------
; You might want to change some parameters set in the code such as,
; number of time step, diagnostics to be computed... Look for the tag !CHECK!

; Parameters
; ----------
; fdir : string
;    input files directory.
; files : string
;    common files prefix.

; Keywords
; --------
; trc: string
;    Tracer name.
; savefile: string
;    Output file name.
; filec: string
;    File prefix for the tracers. If not given, files is used.
; shortfile: integer (1, 2 or 3)
;    1: files for the dynamic have only 12 months.
;    2: files for tracer have only 12 months. 
;    3: files for mld have only 12 months.
; o2sat
;    If activated O2 concentration is computed from temperature and
;    salinity folowing Garcia and Gordon (1992).
; filemld : string
;    File prefix for the mixed layer depth.
; mld_cte
;    If activated, the mixed layer depth is kept constant. file_mld
;    has to be given.
; dep : float
;    Subduction computed accross a fixed depth dep.

; History
; -------
; Adapted from the original code of M. Lévy


;----------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRO bilan_offline, fdir, files, FILEML = fileml, FILEC = filec, SAVEFILE = savefile, SHORTFILE = shortfile, O2SAT = o2sat, MLD_CTE = mld_cte, TRC = trc, DEP = dep
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;----------------------------------------------------------------------------

@common

;---------------------------------------------------------------------------
;%%%%%%%%%%% INPUT FILE AND VARIABLE NAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

; !CHECK!
meshmask = '/data/daclod/mesh_mask.nc'
namemld  = 'somxl010'   ; mixed layer depth         
nameu    = 'vozocrtx'   ; zonal speed                 
namev    = 'vomecrty'   ; meridional speed            
namew    = 'vovecrtz'   ; vertical speed              
namekz   = 'votkeavt'   ; diffusion coef              
nametemp = 'votemper'   ; temperature                  
namesal  = 'vosaline'   ; salinity                     
; !CHECK!

nametr   = TRC 

fileu    = fdir + files + '_'  + nameu    + '.nc' ; zonal speed       
filev    = fdir + files + '_'  + namev    + '.nc' ; meridional speed  
filew    = fdir + files + '_'  + namew    + '.nc' ; vertical speed    
filekz   = fdir + files + '_'  + namekz   + '.nc' ; diffusion coef    
filetemp = fdir + files + '_'  + nametemp + '.nc' ; temperature       
filesal  = fdir + files + '_'  + namesal  + '.nc' ; salinity          
; mixed layer depth
IF KEYWORD_SET(fileml) THEN BEGIN
   filemld = fdir + fileml + '_' + namemld + '.nc' 
ENDIF ELSE BEGIN
   IF KEYWORD_SET(mld_cte) THEN BEGIN
      print, "If MLD_CTE keyword cte, FILEML has to be set too"
      RETALL
   ENDIF ELSE  filemld = fdir + files + '_' + namemld   + '.nc' 
ENDELSE 
; tracer
IF KEYWORD_SET(filec) THEN BEGIN
   IF KEYWORD_SET(o2sat) THEN BEGIN
      filetempO2sat  = fdir + filec + '_' + nametemp + '.nc'
      filesaliO2sat  = fdir + filec + '_' + namesal  + '.nc'
   ENDIF ELSE filetr = fdir + filec + '_' + nametr + '.nc' 
ENDIF ELSE BEGIN
   IF KEYWORD_SET(o2sat) THEN BEGIN
      filetempO2sat  = fdir + files + '_' + nametemp + '.nc'
      filesaliO2sat  = fdir + files + '_' + namesal  + '.nc'
   ENDIF ELSE filetr = fdir + files + '_' + nametr   + '.nc' 
ENDELSE 


;---------------------------------------------------------------------------
;%%%%%%%%%%% PARAMETRES %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

@init_run_clim_dc

; !CHECK!
rdt  = 30 * 86400 ; time step
; rdt has to be define depending on the frequency of the
; inputs. This is because the entrainment/detrainment term is computed
; as the difference between two time steps.
; if it is monthly inputs, rdt=30*24*3600 and unit will be ???/month/mesh
; if it is daily inputs, rdt=24*3600 and unit will be ???/month/mesh
aht0 = 2000       ; isoneutral mixing coef [m2/s]
rau0 = 1035.      ; reference density
avt0 = 1.e-4      ; vertical diffusion coef [m2/s]
; !CHECK!

;---------------------------------------------------------------------------
;%%%%%%%%%% OPTIONS DE CALCULS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

; !CHECK!
nstart = 0    ; first time step
nend   = 1319 ; last time step
; !CHECK!
nstep  = 1 + nend - nstart ; time step number

; If shortfile keyword is set, shorter files will be readen
; cyclically
; ncycl=1, cyclic file, ncycl=0, non cyclic file
IF KEYWORD_SET(shortfile) THEN BEGIN 
   IF shortfile EQ 1 THEN BEGIN ; mld and dyn are short files
      ncycldyn = 1
      ncyclmld = 1
      ncycltr  = 0
      print, 'ncycldyn : ', ncycldyn
      print, 'ncyclmld : ', ncyclmld
      print, 'ncycltr  : ', ncycltr
   ENDIF ELSE IF shortfile EQ 2 THEN BEGIN ; tracer is a shortfile
      ncycldyn = 0
      ncyclmld = 0
      ncycltr  = 1
      print, 'ncycldyn : ', ncycldyn
      print, 'ncyclmld : ', ncyclmld
      print, 'ncycltr  : ', ncycltr
   ENDIF ELSE IF shortfile EQ 3 THEN BEGIN ; mld only is a shortfile
      ncycldyn = 0
      ncyclmld = 1
      ncycltr  = 0
      print, 'ncycldyn : ', ncycldyn
      print, 'ncyclmld : ', ncyclmld
      print, 'ncycltr  : ', ncycltr
   ENDIF ELSE BEGIN
      print, "The SHORTFILE keyword can only be 1, 2 or 3"
      RETALL
   ENDELSE
ENDIF ELSE BEGIN ; default value
   ncycldyn = 0
   ncyclmld = 0
   ncycltr  = 0
   print, 'ncycldyn : ', ncycldyn
   print, 'ncycldyn : ', ncycldyn
   print, 'ncycltr  : ', ncycltr
ENDELSE 

; !CHECK!
ninterpolkz = 0 ; 0= KZ cte,  1=KZ is readen in a file and interpolated at MLD
ISLOPE = 1 ; islope=0, isopycnals forced to 0
; Chosen diagnostics
ISUBMLD = 1 ; subduction by entrainment/detrainment
ISUBH   = 1 ; subduction by lateral advection
ISUBW   = 1 ; subduction by vertical advection
ISUBGM  = 1 ; subdcution by eddies (Gent McWilliams param)
ISUBLDF = 1 ; subdcution by lateral diffusion along isopycnals
ISUBZDF = 1 ; subdcution by vertical diffusion
; !CHECK!

;---------------------------------------------------------------------------
;%%%%%%%%%%%% Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

mbathy = ncdf_lec(meshmask, var = 'mbathy')
ff     = ncdf_lec(meshmask, var = 'ff')
e3u    = ncdf_lec(meshmask, var = 'e3u')
e3v    = ncdf_lec(meshmask, var = 'e3v')
e3t    = ncdf_lec(meshmask, var = 'e3t')
e3w    = ncdf_lec(meshmask, var = 'e3w')
gdepv  = ncdf_lec(meshmask, var = 'gdepv')
gdepu  = ncdf_lec(meshmask, var = 'gdepu')
gdept  = ncdf_lec(meshmask, var = 'gdept')
gdepw  = ncdf_lec(meshmask, var = 'gdepw')

;---------------------------------------------------------------------------
;%%%%%%%%%%%% Subduction terms %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

submld = fltarr(jpi, jpj, nstep)
subw   = fltarr(jpi, jpj, nstep)
subh   = fltarr(jpi, jpj, nstep)
subgm  = fltarr(jpi, jpj, nstep)
subzdf = fltarr(jpi, jpj, nstep)
subldf = fltarr(jpi, jpj, nstep)

;---------------------------------------------------------------------------
;%%%%%%%%%% U and V masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

umask = umask()      
umask[0                , *    , *] = 0.
umask[jpi-1            , *    , *] = 0.
umask[*                , jpj-1, *] = 0.
umask[((jpi/2)+1):jpi-1, jpj-2, *] = 0.
umask[1                , jpj-2, *] = 0.

vmask = vmask()
vmask[0    , *    , *] = 0.
vmask[jpi-1, *    , *] = 0.
vmask[*    , jpj-1, *] = 0.
vmask[*    , jpj-2, *] = 0.


;---------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%% TEMPORAL LOOP %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

print, 'TEMPORAL LOOP START'

FOR n = nstart, nend DO BEGIN

; Set the time step for the dynamics and tracer which can be different
; if the keyword shortfile is activated
    IF KEYWORD_SET(shortfile) THEN BEGIN
       IF ( KEYWORD_SET(filec) NE 1 ) AND ( KEYWORD_SET(fileml) NE 1 ) THEN BEGIN
          print, "There is no meaning to have keyword SHORTFILE set and FILEC or FILEML not set"
          RETALL
       ENDIF
       IF shortfile EQ 1 THEN BEGIN ; mld and dyn 1 year of monthly inputs
          nmld      = n MOD 12
          nstartmld = 0
          nendmld   = 11
          ndyn      = n MOD 12
          nstartdyn = 0
          nenddyn   = 11
          ntr       = n
          nstarttr  = nstart
          nendtr    = nend
       ENDIF ELSE IF shortfile EQ 2 THEN BEGIN ; tracer 1 year of monthly inputs
          nmld      = n
          nstartmld = nstart
          nendmld   = nend
          ndyn      = n
          nstartdyn = nstart
          nenddyn   = nend
          ntr       = n MOD 12
          nstarttr  = 0
          nendtr    = 11
       ENDIF ELSE IF shortfile EQ 3 THEN BEGIN ; mld 1 year of monthly inpputs
          nmld      = n MOD 12 
          nstartmld = 0
          nendmld   = 11
          ndyn      = n
          nstartdyn = nstart
          nenddyn   = nend
          ntr       = n
          nstarttr  = nstart
          nendtr    = nend
       ENDIF ELSE BEGIN
          print, "The SHORTFILE keyword can only be 1, 2 or 3"
          RETALL
       ENDELSE
    ENDIF ELSE BEGIN
       ndyn      = n
       nmld      = n
       ntr       = n
       nstartdyn = nstart
       nenddyn   = nend
       nstartmld = nstart
       nendmld   = nend
       nstarttr  = nstart
       nendtr    = nend
    ENDELSE
    print, '___________________________________'
    print, 'dynamic time step : ', ndyn
    print, 'mld time step     : ', nmld
    print ,'tracer time step  : ', ntr


;---------------------------------------------------------------------------
;%%%%%%%%%%%% DATA PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

;___________________________________________________________________________
;%%%%%%%% READ TRACER n-1, n, n+1

    IF NOT KEYWORD_SET(o2sat) THEN BEGIN 
       ; First case, o2sat not activated, tracer is read

       print, 'Read tracer : ', filetr
       tr = read_ncdf(nametr, ntr, /timestep, file = filetr,  /nostruct)
       tr = orca_repli(tr) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
       IF ncycltr EQ 1 THEN BEGIN 
          ; cyclical case, first time step = last time step
          IF ntr  NE nstarttr THEN  trm1 = read_ncdf(nametr, ntr-1, /timestep, file = filetr,  /nostruct) ELSE trm1 = read_ncdf(nametr, nendtr  , /timestep, file = filetr,  /nostruct)
          IF ntr  NE nendtr   THEN  trp1 = read_ncdf(nametr, ntr+1, /timestep, file = filetr,  /nostruct) ELSE trp1 = read_ncdf(nametr, nstarttr, /timestep, file = filetr,  /nostruct) 
       ENDIF ELSE BEGIN 
          ; non cylical case, first and last time step are duplicated
          IF ntr  NE nstarttr THEN  trm1 = read_ncdf(nametr, ntr-1, /timestep, file = filetr,  /nostruct) ELSE trm1 = read_ncdf(nametr, nstarttr, /timestep, file = filetr,  /nostruct)
          IF ntr  NE nendtr   THEN  trp1 = read_ncdf(nametr, ntr+1, /timestep, file = filetr,  /nostruct) ELSE trp1 = read_ncdf(nametr, nendtr  , /timestep, file = filetr,  /nostruct) 
       ENDELSE 

    ENDIF ELSE BEGIN 
        ; Second case, o2sat activated, O2
        ; concentration computed from
        ; temperature and salinity (Garcia and
        ; Gordon 1992)
      
        print, 'O2sat computation start', ' : ', filetempO2sat
        tempO2sat = read_ncdf(nametemp, ntr, /timestep, file = filetempO2sat,  /nostruct)
        tempO2sat = orca_repli(tempO2sat) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
        saliO2sat = read_ncdf(namesal , ntr, /timestep, file = filesaliO2sat,  /nostruct)
        saliO2sat = orca_repli(saliO2sat) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
        ts        = alog((298.15 - tempO2sat)/(273.15+tempO2sat))
        ; coef de Garcia and Gordon 1992 for oxygen
        a0 = 2.00907
        a1 = 3.22014
        a2 = 4.05010
        a3 = 4.94457
        a4 = -0.256847
        a5 = 3.88767
        b0 = -6.24523e-3
        b1 = -7.37614e-3
        b2 = -1.03410e-2
        b3 = -8.17083e-3
        c0 = -4.88682e-7
        logsolub = a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4 + a5*ts^5 + saliO2sat*(b0 + b1*ts + b2*ts^2 + b3*ts^3) + c0*saliO2sat^2 ; unit cm3/dm3
        tr = exp(logsolub)/1000./22.3916 ; tracer concentration in kmol/m3

        IF ncycltr EQ 1 THEN BEGIN
           ; case: cyclical, first time step = last time step

           IF ntr NE nstarttr  THEN BEGIN 
              tempO2satm1 = read_ncdf(nametemp, ntr-1, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satm1 = orca_repli(tempO2satm1)
              saliO2satm1 = read_ncdf(namesal , ntr-1, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satm1 = orca_repli(saliO2satm1)
              tsm1        = alog((298.15 - tempO2satm1)/(273.15+tempO2satm1))
              logsolubm1  = a0 + a1*tsm1 + a2*tsm1^2 + a3*tsm1^3 + a4*tsm1^4 + a5*tsm1^5 + saliO2satm1*(b0 + b1*tsm1 + b2*tsm1^2 + b3*tsm1^3) + c0*saliO2satm1^2
              trm1        = exp(logsolubm1)/1000./22.3916 ; kmol/m3
           ENDIF ELSE BEGIN               
              tempO2satm1 = read_ncdf(nametemp, nendtr, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satm1 = orca_repli(tempO2satm1)
              saliO2satm1 = read_ncdf(namesal , nendtr, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satm1 = orca_repli(saliO2satm1)
              tsm1        = alog((298.15 - tempO2satm1)/(273.15+tempO2satm1))
              logsolubm1  = a0 + a1*tsm1 + a2*tsm1^2 + a3*tsm1^3 + a4*tsm1^4 + a5*tsm1^5 + saliO2satm1*(b0 + b1*tsm1 + b2*tsm1^2 + b3*tsm1^3) + c0*saliO2satm1^2
              trm1        = exp(logsolubm1)/1000./22.3916 ; kmol/m3
           ENDELSE

           IF ntr NE nendtr THEN BEGIN 
              tempO2satp1 = read_ncdf(nametemp, ntr+1, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satp1 = orca_repli(tempO2satp1)
              saliO2satp1 = read_ncdf(namesal , ntr+1, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satp1 = orca_repli(saliO2satp1)
              tsp1        = alog((298.15 - tempO2satp1)/(273.15+tempO2satp1))
              logsolubp1  = a0 + a1*tsp1 + a2*tsp1^2 + a3*tsp1^3 + a4*tsp1^4 + a5*tsp1^5 + saliO2satp1*(b0 + b1*tsp1 + b2*tsp1^2 + b3*tsp1^3) + c0*saliO2satp1^2
              trp1        = exp(logsolubp1)/1000./22.3916 ; kmol/m3
           ENDIF ELSE BEGIN                      
              tempO2satp1 = read_ncdf(nametemp, nstarttr, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satp1 = orca_repli(tempO2satp1)
              saliO2satp1 = read_ncdf(namesal , nstarttr, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satp1 = orca_repli(saliO2satp1)
              tsp1        = alog((298.15 - tempO2satp1)/(273.15+tempO2satp1))
              logsolubp1  = a0 + a1*tsp1 + a2*tsp1^2 + a3*tsp1^3 + a4*tsp1^4 + a5*tsp1^5 + saliO2satp1*(b0 + b1*tsp1 + b2*tsp1^2 + b3*tsp1^3) + c0*saliO2satp1^2
              trp1        = exp(logsolubp1)/1000./22.3916 ; kmol/m3
           ENDELSE

        ENDIF ELSE BEGIN
           ; case: non cylical, first and last time step are duplicated

           IF ntr NE nstarttr THEN BEGIN 
              tempO2satm1 = read_ncdf(nametemp, ntr-1, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satm1 = orca_repli(tempO2satm1)
              saliO2satm1 = read_ncdf(namesal , ntr-1, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satm1 = orca_repli(saliO2satm1)
              tsm1        = alog((298.15 - tempO2satm1)/(273.15+tempO2satm1))
              logsolubm1  = a0 + a1*tsm1 + a2*tsm1^2 + a3*tsm1^3 + a4*tsm1^4 + a5*tsm1^5 + saliO2satm1*(b0 + b1*tsm1 + b2*tsm1^2 + b3*tsm1^3) + c0*saliO2satm1^2
              trm1        = exp(logsolubm1)/1000./22.3916 ; kmol/m3
           ENDIF ELSE BEGIN
              tempO2satm1 = read_ncdf(nametemp, nstarttr, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satm1 = orca_repli(tempO2satm1)
              saliO2satm1 = read_ncdf(namesal , nstarttr, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satm1 = orca_repli(saliO2satm1)
              tsm1        = alog((298.15 - tempO2satm1)/(273.15+tempO2satm1))
              logsolubm1  = a0 + a1*tsm1 + a2*tsm1^2 + a3*tsm1^3 + a4*tsm1^4 + a5*tsm1^5 + saliO2satm1*(b0 + b1*tsm1 + b2*tsm1^2 + b3*tsm1^3) + c0*saliO2satm1^2
              trm1        = exp(logsolubm1)/1000./22.3916 ; kmol/m3
           ENDELSE

           IF ntr NE nendtr THEN BEGIN
              tempO2satp1 = read_ncdf(nametemp, ntr+1, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satp1 = orca_repli(tempO2satp1)
              saliO2satp1 = read_ncdf(namesal , ntr+1, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satp1 = orca_repli(saliO2satp1)
              tsp1        = alog((298.15 - tempO2satp1)/(273.15+tempO2satp1))
              logsolubp1  = a0 + a1*tsp1 + a2*tsp1^2 + a3*tsp1^3 + a4*tsp1^4 + a5*tsp1^5 + saliO2satp1*(b0 + b1*tsp1 + b2*tsp1^2 + b3*tsp1^3) + c0*saliO2satp1^2
              trp1        = exp(logsolubp1)/1000./22.3916 ; kmol/m3
           ENDIF ELSE BEGIN
              tempO2satp1 = read_ncdf(nametemp, nendtr, /timestep, file = filetempO2sat,  /nostruct)
              tempO2satp1 = orca_repli(tempO2satp1)
              saliO2satp1 = read_ncdf(namesal , nendtr, /timestep, file = filesaliO2sat,  /nostruct)
              saliO2satp1 = orca_repli(saliO2satp1)
              tsp1        = alog((298.15 - tempO2satp1)/(273.15+tempO2satp1))
              logsolubp1  = a0 + a1*tsp1 + a2*tsp1^2 + a3*tsp1^3 + a4*tsp1^4 + a5*tsp1^5 + saliO2satp1*(b0 + b1*tsp1 + b2*tsp1^2 + b3*tsp1^3) + c0*saliO2satp1^2
              trp1        = exp(logsolubp1)/1000./22.3916 ; kmol/m3
           ENDELSE

        ENDELSE                 

     ENDELSE                    

;___________________________________________________________________________
;%%%%%%%% READ MLD n-1, n, n+1

    print, 'Read mld : ', filemld

    IF KEYWORD_SET(mld_cte) THEN BEGIN
       ; case: constant mld 

       print, 'mld is constant over time, nmld = 0'
       nmld = 0
       mld = read_ncdf(namemld, nmld, /timestep, file = filemld,  /nostruct)
       mld = orca_repli(mld) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
       mldm1 = mld
       mldp1 = mld

    ENDIF ELSE BEGIN
       ; case: varying mld 

       mld = read_ncdf(namemld, nmld, /timestep, file = filemld,  /nostruct)
       mld = orca_repli(mld) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)

       IF ncyclmld EQ 1       THEN BEGIN 
          ; case: cyclical, first time step = last time step
          IF nmld NE nstartmld THEN  mldm1 = orca_repli(read_ncdf(namemld, nmld-1, /timestep, file = filemld,  /nostruct)) ELSE mldm1 = orca_repli(read_ncdf(namemld, nendmld  , /timestep, file = filemld,  /nostruct) )
          IF nmld NE nendmld   THEN  mldp1 = orca_repli(read_ncdf(namemld, nmld+1, /timestep, file = filemld,  /nostruct)) ELSE mldp1 = orca_repli(read_ncdf(namemld, nstartmld, /timestep, file = filemld,  /nostruct) )

       ENDIF ELSE BEGIN        
          ; case: non cylical, first and last time step are duplicated
          IF nmld NE nstartmld THEN  mldm1 = orca_repli(read_ncdf(namemld, nmld-1, /timestep, file = filemld,  /nostruct)) ELSE mldm1 = orca_repli(read_ncdf(namemld, nstartmld, /timestep, file = filemld,  /nostruct) )
          IF nmld NE nendmld   THEN  mldp1 = orca_repli(read_ncdf(namemld, nmld+1, /timestep, file = filemld,  /nostruct)) ELSE mldp1 = orca_repli(read_ncdf(namemld, nendmld  , /timestep, file = filemld,  /nostruct) )

       ENDELSE  
       
       IF KEYWORD_SET(DEP) THEN BEGIN
          ; case: fixed depth, no more mld
          print, 'Subduction accros a fixed depth: ', dep, 'm'
          indice = where(mld LT 1.e+19)
          mld(where(mld LT 1.e+19))     = dep
          mldm1(where(mldm1 LT 1.e+19)) = dep
          mldp1(where(mldp1 LT 1.e+19)) = dep
       ENDIF

    ENDELSE
    ; DC2

;___________________________________________________________________________
;%%%%%%%% MLD INDICES

; vertical index at mixed layer base
    nmln = intarr(jpi, jpj)  
    FOR  j = 0, jpj-1 DO BEGIN 
        FOR  i = 0, jpi-1 DO BEGIN 
            delta = abs( gdepw(i, j, *)-mld(i, j)*tmask(i, j, 0))
            ind = where (delta EQ min(delta))  
            ind = max(ind)
            IF tmask(i, j, ind) EQ 0 THEN ind = max([0, ind-1])  
            nmln(i, j) = ind   
        ENDFOR
    ENDFOR 

;!==   surface mixed layer mask   !
    omlmask = fltarr(jpi, jpj, jpk)
    FOR  jk = 0, jpk-1 DO BEGIN ;! = 1 inside the mixed layer, = 0 otherwise
        FOR  jj = 0, jpj-1 DO BEGIN 
            FOR  ji = 0, jpi-1 DO BEGIN 
                ik = nmln(ji, jj) - 1     
                IF ( jk lE ik ) THEN  omlmask(ji, jj, jk) = 1 ELSE omlmask(ji, jj, jk) = 0 
            ENDFOR
        ENDFOR
    ENDFOR

; MLD at U and V grid points
    mldu = 0.5* (mld +shift(mld, -1, 0))*umask
    mldv = 0.5* (mld+ shift(mld, 0, -1))*vmask

; vertical index at the base of the mixed layer in U and V grid faces
; (nmlnu, nmlnv)
    nmlnu = intarr(jpi, jpj)  
    nmlnv = intarr(jpi, jpj)  
    FOR  j = 0, jpj-1 DO BEGIN 
        FOR  i = 0, jpi-1 DO BEGIN 
            delta = abs( gdepw(i, j, *)-mldu(i, j)*umask[i, j, 0]) 
            ind = where (delta EQ min(delta))  
            ind = max(ind)
            IF umask[i, j, ind] EQ 0 THEN ind = max([0, ind-1])  
            nmlnu(i, j) = ind   

            delta = abs( gdepw(i, j, *)-mldv(i, j)*vmask[i, j, 0]) 
            ind = where (delta EQ min(delta))  
            ind = max(ind)
            IF vmask[i, j, ind] EQ 0 THEN ind = max([0, ind-1])  
            nmlnv(i, j) = ind   
        ENDFOR
    ENDFOR 

;___________________________________________________________________________
;%%%%%%%% INTERPOLATION OF TR AT MLD

    print,  'interpolation of tracer concentration at the mixed layer base'
    trmld = fltarr(jpi, jpj)

    FOR j = 0, jpj-1 DO BEGIN 
        FOR i = 0, jpi-1 DO BEGIN 
            trmld(i, j) = interpol(tr(i, j, *), gdept(i, j, *), mld(i, j)) 
            IF tmask(i, j, nmln(i, j)+1) EQ 0 THEN trmld(i, j) = tr(i, j, nmln(i, j)) 
        ENDFOR
    ENDFOR

; trmld at U and V grid point
    trmldu = 0.5*(trmld+shift(trmld, -1, 0))*umask  
    trmldv = 0.5*(trmld+shift(trmld, 0, -1))*vmask

;___________________________________________________________________________
;%%%%%%%%%% VERTICAL SPEED

    wmld = fltarr(jpi, jpj)
    print, 'Read vertical speed : ', filew
    wn = read_ncdf(namew, ndyn, /timestep, file = filew,  /nostruct)
    wn = orca_repli(wn)  ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
    print,  'interpolation of W at the mixed layer base'
    FOR  j = 0, jpj-1 DO BEGIN  
       FOR  i = 0, jpi-1 DO BEGIN  
          wmld(i, j) = interpol(wn(i, j, *), gdepw(i, j, *), mld(i, j))
          IF tmask(i, j, nmln(i, j)+1) EQ 0 THEN wmld(i, j) = wn(i, j, nmln(i, j))
       ENDFOR 
    ENDFOR 


;___________________________________________________________________________
;%%%%%%%%%% VERTICAL DIFFUSION COEF KZ

    IF  (ninterpolkz EQ 0) THEN BEGIN
       ; case: kz constant
       print, 'kz = ', avt0
       kz = avt0*(fltarr(jpi, jpj, jpk)+1.)*tmask
    ENDIF ELSE BEGIN
       ; case: kz read
       print, 'Read kz : ', filekz 
       kz = read_ncdf(namekz, ndyn, /timestep, file = filekz,  /nostruct)
    ENDELSE
    
    kz = orca_repli(kz)

;___________________________________________________________________________
;%%%%%%%%%% HORIZONTAL SPEED

    umld = fltarr(jpi, jpj)
    vmld = fltarr(jpi, jpj)
    print, 'Read zonal speed: ', fileu
    un = read_ncdf(nameu, ndyn, /timestep, file = fileu,  /nostruct)
    un = orca_repli(un) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
    print, 'Read meridional speed: ', filev
    vn = read_ncdf(namev, ndyn, /timestep, file = filev,  /nostruct)
    vn = orca_repli(vn) ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
    print,  'interpolation of U and V at the mixed layer base'
    FOR  j = 0, jpj-1 DO BEGIN  
       FOR  i = 0, jpi-1 DO BEGIN  
          umld(i, j) = interpol(un(i, j, *), gdept(i, j, *), mldu(i, j))*umask[i, j, nmlnu(i, j)]
          vmld(i, j) = interpol(vn(i, j, *), gdept(i, j, *), mldv(i, j))*vmask[i, j, nmlnv(i, j)]
          IF umask[i, j, nmlnu(i, j)+1] EQ 0 THEN umld(i, j) = un(i, j, nmlnu(i, j))*umask[i, j, nmlnu(i, j)] 
          IF vmask[i, j, nmlnv(i, j)+1] EQ 0 THEN vmld(i, j) = vn(i, j, nmlnv(i, j))*vmask[i, j, nmlnv(i, j)] 
       ENDFOR 
    ENDFOR 
    
;___________________________________________________________________________
;%%%%%%%%%%% READ TEMPERATURE AND SALINITY, COMPUTE DENSITY AND ISO SLOPES

    IF ISLOPE NE 0 THEN BEGIN 
       ; case: isopycnal slope computed
       IF  ISUBLDF or ISUBZDF or ISUBGM THEN BEGIN       
          
          print, 'Read temperature: ', filetemp
          tn = read_ncdf(nametemp, ndyn, /timestep, file = filetemp,  /nostruct)
          tn = orca_repli(tn)   ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
          print, 'Read salinity: ', filesal
          sn = read_ncdf(namesal , ndyn, /timestep, file = filesal ,  /nostruct)
          sn = orca_repli(sn)   ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
          print, 'Compute density'
          prd = rhoeos(tn, sn, rau0)
          print, 'Compute Brunt-Vaisala frequency'
          pn2 = bvf(tn, sn)
          
          print, 'Compute isopycnal slopes'
          slp = ldfslp(prd, pn2, nmln, mld, umask, vmask, omlmask)
          wslpi = slp(*, *, *, 0)
          wslpj = slp(*, *, *, 1)
          uslp = slp(*, *, *, 2)
          vslp = slp(*, *, *, 3)
       ENDIF 
    ENDIF ELSE BEGIN 
       slp = fltarr(jpi, jpj, jpk, 4)
       wslpi = slp(*, *, *, 0)
       wslpj = slp(*, *, *, 1)
       uslp = slp(*, *, *, 2)
       vslp = slp(*, *, *, 3)
    ENDELSE 
    
    IF ISLOPE EQ 0 AND ISUBGM THEN BEGIN 
       ; case: isopycnal slopes forced to zero
       print, 'Read temperature: ', filetemp
       tn = read_ncdf(nametemp, ndyn, /timestep, file = filetemp,  /nostruct)
       tn = orca_repli(tn)      ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
       print, 'Read salinity: ', filesal
       sn = read_ncdf(namesal , ndyn, /timestep, file = filesal ,  /nostruct)
       sn = orca_repli(sn)      ; reconstruction of cyclical conditions on x, y (specific to ORCA grid)
       print, 'Compute density'
       prd = rhoeos(tn, sn, rau0)
       print, 'Compute Brunt-Vaisala frequency'
       pn2 = bvf(tn, sn)
    ENDIF 
    
;---------------------------------------------------------------------------
;%%%%%%%%%%%% SUBDUCTION TERMS COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

;%%%%%%%%%%% 1-ML term:
    IF  ISUBMLD THEN BEGIN 
        print, 'Calcul de la subduction par variation de MLD'
        submld(*, *, n - nstart) = trsubmld(mld, mldm1, mldp1, tr, trm1)
    ENDIF 

;%%%%%%%%%%% 2-W term:
    IF  ISUBW THEN BEGIN 
        print, 'Calcul de la subduction par W'
        subw(*, *, n - nstart) = rdt*trsubwadv(wmld, trmld)
    ENDIF 

;%%%%%%%%%%% 3-U and V terms
    IF  ISUBH THEN BEGIN 
        print, 'Calcul de la subduction par U et V'
        subh(*, *, n - nstart) = rdt*trsubhadv(mld, umld, vmld, trmldu, trmldv, umask,  vmask,  nmlnu,  nmlnv)
    ENDIF 

;%%%%%%%%%%% 4-Gent et Mc Williams advection - 
    IF  ISUBGM THEN BEGIN 
        print, 'Calcul de la subduction par gent et mc williams'
        subgm(*, *, n - nstart) =  rdt* trsubgm(mld, mldu, mldv, trmldu, trmldv, trmld, prd, pn2, un, vn, wn, wslpi, wslpj, uslp, vslp, umask,  vmask, mbathy, ff, nmlnu,  nmlnv)
    ENDIF 

;%%%%%%%%%%% 5-Vertical diffusion term
    IF  ISUBZDF THEN BEGIN 
        print, 'Calcul de la subduction par diffusion verticale'
        subzdf(*, *, n - nstart) = trsubzdf(rdt, tr, wslpi, wslpj, nmln, kz, aht0)
        ; !CHECKED! subzdf(*, *, n - nstart) = trsubzdf_v2(rdt, tr, wslpi, wslpj, nmln, kz, aht0) ; vertical contribution of iso-neutral mixing removed, but be careful it is not included in subldf. If you want to get the value of the vertical contribution of isopycnal mixing, you have to compute subzdf twice : one with and one without.
    ENDIF 

;%%%%%%%%%%% 6-Isopycnal diffusion term
    IF  ISUBLDF THEN  BEGIN 
        print, 'calcul de la subduction par diffusion laterale'
        subldf(*, *, n - nstart) = rdt*trsubldf(tr, aht0, mld, uslp, vslp, wslpi, wslpj, umask, vmask, nmln, nmlnu, nmlnv)
    ENDIF 

ENDFOR

;---------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%% FIN DE BOUCLE TEMPPORELLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------


; unit: kmol/month/mesh 
; note that the unit is partialy define by the frequency of the input
; if it is monthly inputs, unit is ???/month/mesh
; if it is daily inputs, unit is ???/month/mesh
IF KEYWORD_SET(savefile) THEN BEGIN
   ; !CHECK!
    save,  subw, subh,  subgm, subzdf, subldf, submld, file = savefile 
   ; !CHECK!
    print, 'file saved : ', savefile
ENDIF

END 
