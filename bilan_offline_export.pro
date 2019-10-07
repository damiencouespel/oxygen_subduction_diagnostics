;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;           
;;;;;;;;;;     BILAN_OFFLINE_EXPORT     ;;;;;;;;;; 
;;;;;;;;;;                              ;;;;;;;;;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Purpose 
; -------
; compute and save the offline computation organic matter export
; accross the mixed layer depth. Note it only compute the sinking flux
; of particulate organic carbon.

; Usage
; -----
; fdir='/data/daclod/DATA_GREENLAND02SV/
; files='piControl2_19900101_20991231_1M'
; - EX1: export accross varying mld
; savefile='export_offline_' + files + '.sav'
; bilan_offline_export, fdir, files, SAVEFILE=savefile
; - EX2: export accross 200m depth
; savefile='export_offline_' + files + '_dep200m.sav'
; bilan_offline_export, fdir, files, SAVEFILE=savefile, DEP=200
 
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
; savefile: string
;    Output file name.
; fileml : string
;    File prefix for the mixed layer depth.
; shortfile: integer (1, 2 or 3)
;    1: files for the mld have only 12 months.
;    2: files for tracer have only 12 months. 
; dep : float
;    Export computed accross a fixed depth dep.
; mld_cte
;    If activated, the mixed layer depth is kept constant. file_mld
;    has to be given.

PRO bilan_offline_export, fdir, files, FILEML = fileml, SHORTFILE = shortfile, DEP = dep, SAVEFILE = savefile, MLD_CTE = mld_cte

@common
@init_run_clim_dc

;---------------------------------------------------------------------------
;%%%%%%%%%%% INPUT FILE AND VARIABLE NAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

; !CHECK!
nameMLD='somxl010' ; mixed layer depth
namePOC='POC'      ; small particles
nameGOC='GOC'      ; large particles
; !CHECK!

filePOC = fdir + files + '_' + namePOC + '.nc' ; small particles
fileGOC = fdir + files + '_' + nameGOC + '.nc' ; large particles
; mixed layer depth
IF KEYWORD_SET(fileml) THEN BEGIN
   fileMLD = fdir + fileml + '_' + namemld + '.nc' 
ENDIF ELSE BEGIN
   fileMLD = fdir + files + '_' + namemld   + '.nc' 
ENDELSE 

meshmask= '/data/daclod/mesh_mask.nc'

;---------------------------------------------------------------------------
;%%%%%%%%%%% PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

; Nombre de pas de temps dans les fichiers d'entrée
; !CHECK!
nstart = 0    ; first time step
nend   = 1319 ; last time step
; !CHECK!
nstep  = 1 + nend - nstart ; time step number

wGOC = -30.   ; GOC sinking speed at mixed layer base in [m/day]
wPOC = -2.    ; POC sinking speed at mixed layer base in [m/day]
rdt = 30.     ; conversion from 1/day en 1/month

;---------------------------------------------------------------------------
;%%%%%%%%%%%% Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

e3u = ncdf_lec(meshmask, var = 'e3u')
e3v = ncdf_lec(meshmask, var = 'e3v')
e3t = ncdf_lec(meshmask, var = 'e3t')
e3w = ncdf_lec(meshmask, var = 'e3w')
gdepv = ncdf_lec(meshmask, var = 'gdepv')
gdepu = ncdf_lec(meshmask, var = 'gdepu')
gdept = ncdf_lec(meshmask, var = 'gdept')
gdepw = ncdf_lec(meshmask, var = 'gdepw')

;---------------------------------------------------------------------------
;%%%%%%%%%%%% Export term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------
export = fltarr(jpi,jpj,nstep)

;---------------------------------------------------------------------------
;%%%%%%%%%%%% TEMPORAL LOOP %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;---------------------------------------------------------------------------

print, 'TEMPORAL LOOP START'

FOR n = nstart, nend DO BEGIN

; Set the time step for the dynamics and tracer which can be different
; if the keyword shortfile is activated
   IF KEYWORD_SET(shortfile) THEN BEGIN
      IF KEYWORD_SET(fileml) NE 1 THEN BEGIN
         print, "There is no meaning to have keyword SHORTFILE set and FILEML not set"
         RETALL
      ENDIF
      IF shortfile EQ 1 THEN BEGIN ; MLD: 1 year with 12 monthly inputs
         nml      = n MOD 12
         ntr      = n
         nstartml = 0
         nendml   = 11
         nstarttr = nstart
         nendtr   = nend
      ENDIF ELSE IF shortfile EQ 2 THEN BEGIN ; POC adn GOC: 1 year of monthly inputs
         nml      = n
         ntr      = n MOD 12
         nstartml = nstart
         nendml   = nend
         nstarttr = 0
         nendtr   = 11
      ENDIF ELSE BEGIN
         print, "The SHORTFILE keyword can only be 1 or 2"
         RETALL
      ENDELSE
   ENDIF ELSE BEGIN
       nml      = n
       ntr      = n
       nstartml = nstart
       nendml   = nend
       nstarttr = nstart
       nendtr   = nend
    ENDELSE
    print, '___________________________________'
    print, 'mld time step       : ', nml
    print ,'POC, GOC time step  : ', ntr

;___________________________________________________________________________
;%%%%%%%% READ INPUT

    print, 'Read input'

    print, '>>> POC'
    POC = read_ncdf(namePOC, ntr, /timestep, file=filePOC, /nostruct)
    POC = orca_repli(POC)

    print, '>>> GOC'
    GOC = read_ncdf(nameGOC, ntr, /timestep, file=fileGOC, /nostruct)
    GOC = orca_repli(GOC)
    
    print, '>>> MLD'
    IF KEYWORD_SET(mld_cte) THEN BEGIN
       print, 'mld is constant over time, nml = 0'
       nml = 0
       MLD = read_ncdf(nameMLD, nml, /timestep, file=fileMLD, /nostruct)
    ENDIF ELSE BEGIN
       MLD = read_ncdf(nameMLD, nml, /timestep, file=fileMLD, /nostruct)
    ENDELSE
    MLD = orca_repli(MLD)

    IF KEYWORD_SET(DEP) THEN BEGIN
       print, '>>> suppression mld, export a  ', dep, 'm'
       mld(where(mld LT 1.e+19))     = dep
    ENDIF


;___________________________________________________________________________
;%%%%%%%% MLD INDEX

    nmln = intarr(jpi, jpj)  
    FOR  jj = 0, jpj-1 DO BEGIN 
        FOR  ii = 0, jpi-1 DO BEGIN 
            delta = abs( gdepw(ii, jj, *)-mld(ii, jj)*tmask(ii, jj, 0))
            ind = where (delta EQ min(delta))  
            ind = max(ind)
            IF tmask(ii, jj, ind) EQ 0 THEN ind = max([0, ind-1])  
            nmln(ii, jj) = ind
        ENDFOR
    ENDFOR 

;______________________________________________________________
;%%%%%%%% INTERPOLATION AT MLD

    pocmld = fltarr(jpi,jpj)
    gocmld = fltarr(jpi,jpj)
    FOR ii = 0, jpi-1 DO BEGIN
        FOR jj = 0, jpj-1 DO BEGIN
            pocmld(ii, jj) = interpol(poc(ii, jj, *), gdept(ii, jj, *), mld(ii, jj))
            gocmld(ii, jj) = interpol(goc(ii, jj, *), gdept(ii, jj, *), mld(ii, jj))
            IF tmask(ii, jj, nmln(ii, jj)+1) EQ 0 THEN BEGIN
                pocmld(ii, jj) = poc(ii, jj, nmln(ii, jj))
                gocmld(ii, jj) = goc(ii, jj, nmln(ii, jj))
            END
        END
    END

;______________________________________________________________
;%%%%%%%% COMPUTE EXPORT
    print, 'compute export'

    subpoc=fltarr(jpi,jpj)
    subgoc=fltarr(jpi,jpj)
    surf=e1t*e2t

    subpoc=rdt*surf*pocmld*wpoc
    subgoc=rdt*surf*gocmld*wgoc
    export(*,*,n)=subpoc+subgoc

END

;conversion from kmol/month/mesh to Gmol/month/mesh
export = export*1.e-6

IF KEYWORD_SET(savefile) THEN BEGIN
    save, export, file = savefile 
    print, 'file saved : ', savefile
ENDIF

END
