;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;;; FIG_BILAN_GLOBAL_FLX_VS_TIME ;;;;;;;;;;   
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Purpose 
; -------
; Plot timeseries of global subduction and oxygen consumption in the control
; simulation and the cumulative change in the global warming
; simulation 

;==================================================
; PARAM
;==================================================

namefigbase = 'fig_bilan_global_flx_vs_time' ; base of the file name to save the figures 
                                                                            
fdir = '/data/daclod/MY_IDL/DATA_TMP/' ; directory of the files to restore  

;__________________________________________________
; name of the files to restore for the climate change simulation

; O2 subduction
fTOp     = 'sub_offline_O2_v3.rcp85.1_19900101_20991231_1M.sav'
fzdfTOp  = 'sub_offline_O2_v3.rcp85.1_19900101_20991231_1M_subzdf.sav' ; when not commented, do not forget to uncomment the restore command below
fzdf2TOp = 'sub_offline_v10_O2_v3.rcp85.1_19900101_20991231_1M_subzdf.sav' ; when not commented, do not forget to uncomment the restore command below (ZDF sans diff iso-neutral)
; Saturated O2 subduction
faTOsatp = 'sub_offline_v3.rcp85.1_19900101_19991231_1M_ymonmean_O2sat_v3.rcp85.1_19900101_20991231_1M.sav'
; Organic matter (POC and GOC) sinking
fEXPp    = 'export_offline_v3.rcp85.1_19900101_20991231_1M.sav'
; Organic matter (DOC, POC, GOC, ZOO, ZOO2, PHYT) and NH4 subduction
fTDOCp   = 'sub_offline_v8_DOC_v3.rcp85.1_19900101_20991231_1M.sav'
fTPOCp   = 'sub_offline_v9_POC_v3.rcp85.1_19900101_20991231_1M.sav'
fTGOCp   = 'sub_offline_v9_GOC_v3.rcp85.1_19900101_20991231_1M.sav'
fTZOOp   = 'sub_offline_v9_ZOO_v3.rcp85.1_19900101_20991231_1M.sav'
fTZOO2p  = 'sub_offline_v9_ZOO2_v3.rcp85.1_19900101_20991231_1M.sav'
fTPHYTp  = 'sub_offline_v9_PHYT_v3.rcp85.1_19900101_20991231_1M.sav'
fTNH4p   = 'sub_offline_v9_NH4_v3.rcp85.1_19900101_20991231_1M.sav'

;__________________________________________________
; name of the files to restore for the control simulation

; O2 subduction
fTOc     = 'sub_offline_O2_piControl2_19900101_20991231_1M.sav'
fzdfTOc  = 'sub_offline_O2_piControl2_19900101_20991231_1M_subzdf.sav' ; when not commented, do not forget to uncomment the restore command below
fzdf2TOc = 'sub_offline_v10_O2_piControl2_19900101_20991231_1M_subzdf.sav' ; when not commented, do not forget to uncomment the restore command below (ZDF sans diff iso-neutral)
; Saturated O2 subduction
faTOsatc = 'sub_offline_piControl2_19900101_19991231_1M_ymonmean_O2sat_piControl2_19900101_20991231_1M.sav'
; Organic matter (POC and GOC) sinking
fEXPc    = 'export_offline_piControl2_19900101_20991231_1M.sav'
; Organic matter (DOC, POC, GOC, ZOO, ZOO2, PHYT) and NH4 subduction
fTDOCc   = 'sub_offline_v8_DOC_piControl2_19900101_20991231_1M.sav'
fTPOCc   = 'sub_offline_v9_POC_piControl2_19900101_20991231_1M.sav'
fTGOCc   = 'sub_offline_v9_GOC_piControl2_19900101_20991231_1M.sav'
fTZOOc   = 'sub_offline_v9_ZOO_piControl2_19900101_20991231_1M.sav'
fTZOO2c  = 'sub_offline_v9_ZOO2_piControl2_19900101_20991231_1M.sav'
fTPHYTc  = 'sub_offline_v9_PHYT_piControl2_19900101_20991231_1M.sav'
fTNH4c   = 'sub_offline_v9_NH4_piControl2_19900101_20991231_1M.sav'

; conversion factor of oxygen subduction
zfactsub   =               1000. * 1.e-15  ; kmolO2 -> PmolO2
; conversion factor of organic matter subduction into O2 consumption
zfactbio   = -160./122.  * 1000. * 1.e-15  ; kmolC -> PmolO2  
; conversion factor of organic matter sinking into O2 consumption 
zfactexp   =  160./122.  * 1.e9 * 1.e-15   ; GmolC -> PmolO2
; conversion factor of ammonium subduction into O2 consumption (nitrification)
zfactnit   = - 40./122.  * 1000. * 1.e-15 ; kmolC -> PmolO2  nitrification amonium
; 160/122: the ratio Carbon/Oxygen take into account remineralization
; (X % of 132/122) and nitification (1-X % of 40/122)


;__________________________________________________
; name of the files for oxygen and saturated oxygen data

fdirO2 = '/data/daclod/DATA_CMIP5/' ; directory of oxygen data
; oxygen data
fO2C = 'piControl2_19900101_20991231_1Y_O2.nc' ; in control simulation
fO2P = 'v3.rcp85.1_19900101_20991231_1Y_O2.nc' ; in climate change simulation
; data for o2sat
fTC = 'piControl2_19900101_20991231_1Y_votemper.nc' ; temperature in the control simulation
fTP = 'v3.rcp85.1_19900101_20991231_1Y_votemper.nc' ; temperature in the climate change simulation
fSC = 'piControl2_19900101_20991231_1Y_vosaline.nc' ; salinity in the control simulation       
fSP = 'v3.rcp85.1_19900101_20991231_1Y_vosaline.nc' ; salinity in the climate change simulation
; mixed layer depth data
fMLDP = 'v3.rcp85.1_19900101_20991231_1Y_somxl010.nc' ; in the control simulation       
fMLDC = 'piControl2_19900101_20991231_1Y_somxl010.nc' ; in the climate change simulation
; first and last date to read
start1 = 19900101
end1 = 20991231
zfactoxy = 1000. * 1.e-15  ; to convert from kmol/m3 to PmolO2/m3


;==================================================
; PROCESS SUBDUCTION AND EXPORT DATA
;==================================================

@common
@init_run_clim_dc
surf = e1t*e2t ; mesh surface

;__________________________________________________
; Oxygen subduction 

; in the control simulation
subdu = restore_subdu(fdir, fTOC, 1e13, FILEZDF1 = fzdfTOc, FILEZDF2 = fzdf2TOc) ; varying mld, correction diffusion isoneutral
; unit : kmolO2 / mesh / month
TOc   = bilan_global_subdu_vs_time(subdu, zfactsub)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTOP, 1e13, FILEZDF1 = fzdfTOp, FILEZDF2 = fzdf2TOp) ; varying mld, correction diffusion isoneutral
; unit : kmolO2 / mesh / month
TOp   = bilan_global_subdu_vs_time(subdu, zfactsub)
; unit : PmolO2

;__________________________________________________
; Saturated oxygen subduction 

; in the control simulation
subdu = restore_subdu(fdir, faTOsatc, 1e13)
; unit : kmolO2 / mesh / month
aTOsatc   = bilan_global_subdu_vs_time(subdu, zfactsub)
; unit : PmolO2

;in the climate change simulation
subdu = restore_subdu(fdir, faTOsatp, 1e13)
; unit : kmolO2 / mesh / month
aTOsatp   = bilan_global_subdu_vs_time(subdu, zfactsub)
; unit : PmolO2

;__________________________________________________
; DOC subduction 

; in the control simulation 
subdu = restore_subdu(fdir, fTDOCc, 1e12)
; unit : kmolO2 / mesh / month
TDOCc   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTDOCp, 1e12)
; unit : kmolO2 / mesh / month
TDOCp   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; POC subduction 

; in the control simulation
subdu = restore_subdu(fdir, fTPOCc, 1e10)
; unit : kmolO2 / mesh / month
TPOCc   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTPOCp, 1e10)
; unit : kmolO2 / mesh / month
TPOCp   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; GOC subduction

; in the control simulation
subdu = restore_subdu(fdir, fTGOCc, 1e9)
; unit : kmolO2 / mesh / month
TGOCc   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTGOCp, 1e9)
; unit : kmolO2 / mesh / month
TGOCp   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; ZOO subduction 

; in the control simulation
subdu = restore_subdu(fdir, fTZOOc, 1e10)
; unit : kmolO2 / mesh / month
TZOOc   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTZOOp, 1e10)
; unit : kmolO2 / mesh / month
TZOOp   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; ZOO2 subduction

; in the control simulation
subdu = restore_subdu(fdir, fTZOO2c, 1e10)
; unit : kmolO2 / mesh / month
TZOO2c   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTZOO2p, 1e10)
; unit : kmolO2 / mesh / month
TZOO2p   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; Total phytoplankton subduction 

; in the control simulation
subdu = restore_subdu(fdir, fTPHYTc, 1e11)
; unit : kmolO2 / mesh / month
TPHYTc   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTPHYTp, 1e11)
; unit : kmolO2 / mesh / month
TPHYTp   = bilan_global_subdu_vs_time(subdu, zfactbio)
; unit : PmolO2

;__________________________________________________
; NH4 subduction

; in the control simulation
subdu = restore_subdu(fdir, fTNH4c, 1e11)
; unit : kmolO2 / mesh / month
TNH4c   = bilan_global_subdu_vs_time(subdu, zfactnit)
; unit : PmolO2

; in the climate change simulation
subdu = restore_subdu(fdir, fTNH4p, 1e11)
; unit : kmolO2 / mesh / month
TNH4p   = bilan_global_subdu_vs_time(subdu, zfactnit)
; unit : PmolO2

;__________________________________________________
; POC and GOC sinking 

; restore 
restorefile = fdir + fexpC
print, 'file restored : ', restorefile
restore, file = restorefile, /verbose
; remove outliers
export( where( abs(export) GT 1e15 ) ) = 0.
; to remove the overlap in the indian ocean DC2
; and north pole
export   = orca2_3d_overlap_remove(export)
; unit : GmolC / mesh / month
; annual integral
jpt = 12
zwork = FLTARR(jpi, jpj, 110)
FOR yy  = 0, 109 DO BEGIN
   zwork(*, *, yy) = grossemoyenne(export(*, *, yy*12 : yy*12 + 11)  , 't', /INT, /NAN) * zfactexp / surf
ENDFOR
; unit : GmolC * zfactexp / m2 / y
; remove overlap on the north pole and east-west
zwork = orca2_3d_overlap_remove(zwork)
; horizontal integral
jpt = 110
expC = grossemoyenne(zwork, 'xy', /NAN, /INT)
; unit : GmolC * zfactexp / y
  
; restore
restorefile = fdir + fexpP
print, 'file restored : ', restorefile
restore, file = restorefile, /verbose
; remove outliers
export( where( abs(export) GT 1e15 ) ) = 0.
; to remove the overlap in the indian ocean DC2
; and north pole
export   = orca2_3d_overlap_remove(export)
; unit : GmolC / mesh / month
; annual integral
jpt = 12
zwork = FLTARR(jpi, jpj, 110)
FOR yy  = 0, 109 DO BEGIN
   zwork(*, *, yy) = grossemoyenne(export(*, *, yy*12 : yy*12 + 11)  , 't', /INT, /NAN) * zfactexp / surf
ENDFOR
; unit : GmolC * zfactexp / m2 / y
; remove overlap on the north pole and east-west
zwork = orca2_3d_overlap_remove(zwork)
; horizontal integral
jpt = 110
expP = grossemoyenne(zwork, 'xy', /NAN, /INT)
; unit : GmolC * zfactexp / y

;__________________________________________________
; SUM

; O2 consumpution due to organic matter subduction in the control simulation
subbioC    = TDOCc_dyn + TPOCc_dyn + TGOCc_dyn + TZOOc_dyn + TZOO2c_dyn + TPHYTc_dyn + TNH4c_dyn
; O2 consumption due to total export in the control simulation
bioC       = subbioC + expC
; total = O2 subduction + O2 consumption in the control simulation
totC       = bioC + TOC_dyn

; O2 consumpution due to organic matter subduction in the climate
; change simulation
subbioP    = TDOCp_dyn + TPOCp_dyn + TGOCp_dyn + TZOOp_dyn + TZOO2p_dyn + TPHYTp_dyn + TNH4p_dyn
; O2 consumption due to total export in the climate change simulation
bioP       = subbioP + expP
; total = O2 subduction + O2 consumption in the climate change simulation
totP       = bioP + TOp_dyn

;__________________________________________________
; Change = climate change - control

dTO_dyn = total( TOp.dyn - TOc.dyn, /cumulative)  
dTO_had = total( TOp.had - TOc.had, /cumulative)  
dTO_zad = total( TOp.zad - TOc.zad, /cumulative)  
dTO_zdf = total( TOp.zdf - TOc.zdf, /cumulative)
dTO_mld = total( TOp.mld - TOc.mld, /cumulative)  
dTO_ldf = total( TOp.ldf - TOc.ldf, /cumulative)  
dTO_gm  = total( TOp.gm  - TOc.gm , /cumulative)  
dTO_adv = total( TOp.adv - TOc.adv, /cumulative)  
dTO_zmx = total( TOp.zmx - TOc.zmx, /cumulative)  
dTO_edd = total( TOp.edd - TOc.edd, /cumulative)  

daTOsat_dyn = total( aTOsatp.dyn - aTOsatc.dyn, /cumulative)  
daTOsat_had = total( aTOsatp.had - aTOsatc.had, /cumulative)  
daTOsat_zad = total( aTOsatp.zad - aTOsatc.zad, /cumulative)  
daTOsat_zdf = total( aTOsatp.zdf - aTOsatc.zdf, /cumulative)  
daTOsat_mld = total( aTOsatp.mld - aTOsatc.mld, /cumulative)  
daTOsat_ldf = total( aTOsatp.ldf - aTOsatc.ldf, /cumulative)  
daTOsat_gm  = total( aTOsatp.gm  - aTOsatc.gm , /cumulative)  
daTOsat_adv = total( aTOsatp.adv - aTOsatc.adv, /cumulative)  
daTOsat_zmx = total( aTOsatp.zmx - aTOsatc.zmx, /cumulative)  
daTOsat_edd = total( aTOsatp.edd - aTOsatc.edd, /cumulative)  

dTAOU_dyn = daTOsat_dyn - dTO_dyn
dTAOU_had = daTOsat_had - dTO_had
dTAOU_zad = daTOsat_zad - dTO_zad
dTAOU_zdf = daTOsat_zdf - dTO_zdf
dTAOU_mld = daTOsat_mld - dTO_mld
dTAOU_ldf = daTOsat_ldf - dTO_ldf
dTAOU_gm  = daTOsat_gm  - dTO_gm 
dTAOU_adv = daTOsat_adv - dTO_adv
dTAOU_zmx = daTOsat_zmx - dTO_zmx
dTAOU_edd = daTOsat_edd - dTO_edd

dTDOC_dyn  = total( TDOCp.dyn  - TDOCc.dyn, /cumulative)  
dTPOC_dyn  = total( TPOCp.dyn  - TPOCc.dyn, /cumulative)  
dTGOC_dyn  = total( TGOCp.dyn  - TGOCc.dyn, /cumulative)  
dTZOO_dyn  = total( TZOOp.dyn  - TZOOc.dyn, /cumulative)  
dTZOO2_dyn = total( TZOO2p.dyn - TZOO2c.dyn, /cumulative)  
dTPHYT_dyn = total( TPHYTp.dyn - TPHYTc.dyn, /cumulative)  
dTNH4_dyn  = total( TNH4p.dyn  - TNH4c.dyn, /cumulative)  

dsubbio = total( subbioP - subbioC, /cumulative)
dexp    = total( expP - expC, /cumulative)
dbio = dsubbio + dexp
dtot = dbio + dTO_dyn

;=================================================
; PROCESS DATA O2 and O2SAT
;=================================================

@common
@init_run_clim_dc
surf = e1t*e2t ; mesh surface

;__________________________________________________
; Read data

oxy_P = read_ncdf('O2'      , start1, end1, file=fdirO2 + fO2P, /nostru) * zfactoxy
oxy_C = read_ncdf('O2'      , start1, end1, file=fdirO2 + fO2C, /nostru) * zfactoxy

tem_P = read_ncdf('votemper', start1, end1, file=fdirO2 + fTP , /nostru) 
tem_C = read_ncdf('votemper', start1, end1, file=fdirO2 + fTC , /nostru)
sal_P = read_ncdf('vosaline', start1, end1, file=fdirO2 + fSP , /nostru)
sal_C = read_ncdf('vosaline', start1, end1, file=fdirO2 + fSC , /nostru)

; MLD
IF KEYWORD_SET(dep) THEN BEGIN
   mldP = FLTARR(jpi, jpj, jpt) + dep
   mldC = FLTARR(jpi, jpj, jpt) + dep
ENDIF ELSE BEGIN
   mldP = read_ncdf('somxl010', start1, end1, file=fdirO2 + fMLDP, /nostruct)
   mldC = read_ncdf('somxl010', start1, end1, file=fdirO2 + fMLDC, /nostruct)
ENDELSE

;__________________________________________________
; Mesh volume

tvolume = fltarr(jpi, jpj, jpk)
FOR ii = 0, jpi-1 DO BEGIN
   FOR jj = 0, jpj-1 DO BEGIN
      FOR kk = 0, jpk-1 DO BEGIN
         tvolume(ii, jj, kk) = surf(ii, jj) * e3t(kk)
      ENDFOR 
   ENDFOR
ENDFOR
tvolume = tvolume * tmask

; __________________
; Saturated oxygen from Garcia and Gordon 1992

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
; control simulation
ts = alog( (298.15 - temC) / (273.15+temC) )
logsolub = a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4 + a5*ts^5 + salC*(b0 + b1*ts + b2*ts^2 + b3*ts^3) + c0*salC^2 ; [cm3/dm3]
o2sC = exp(logsolub)/1000./22.3916 *zfactoxy ; [kmol/m3 * zfactoxy]
; climate change simulation
ts = alog( (298.15 - temP) / (273.15+temP) )
logsolub = a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4 + a5*ts^5 + salP*(b0 + b1*ts + b2*ts^2 + b3*ts^3) + c0*salP^2 ; [cm3/dm3]
o2sP = exp(logsolub)/1000./22.3916 * zfactoxy ; [kmol/m3 * zfactoxy]

;__________________________________________________
; Oxygen content in each mesh

jpt = 110
Qoxy_P = FLTARR(jpi, jpj, jpk, jpt) ; oxygen content in the climate change simulation
Qoxy_C = FLTARR(jpi, jpj, jpk, jpt) ; oxygen content in the control simulation        
Qo2sC = FLTARR(jpi, jpj, jpk, jpt)  ; saturated oxygen content in the climate change simulation
Qo2sP = FLTARR(jpi, jpj, jpk, jpt)  ; saturated oxygen content in the control simulation       
FOR tt = 0, jpt-1 DO BEGIN
   QoxyP(*, *, *, tt) = oxyP(*, *, *, tt) * tvolume
   QoxyC(*, *, *, tt) = oxyC(*, *, *, tt) * tvolume
   Qo2sP(*, *, *, tt) = o2sP(*, *, *, tt) * tvolume
   Qo2sC(*, *, *, tt) = o2sC(*, *, *, tt) * tvolume
END

;__________________________________________________
; Oxygen content under the mixed layer depth

; initialization 
oxyP_umld   = fltarr(jpi, jpj, jpt)
oxyC_umld   = fltarr(jpi, jpj, jpt)
o2sP_umld   = fltarr(jpi, jpj, jpt)
o2sC_umld   = fltarr(jpi, jpj, jpt)

jpt = 110
FOR tt = 0, jpt-1 DO BEGIN 
   FOR ii = 0, jpi-1 DO BEGIN 
      FOR jj = 0, jpj-1 DO BEGIN
         zmldP = mldP(ii, jj, tt)
         zmldC = mldC(ii, jj, tt)
         ; index where mld is the deepest
         indP = where( gdepw GE zmldP ) 
         indC = where( gdepw GE zmldC ) 
         ; Oxygen content under MLD = all
         ; levels below MLD + prorata of the
         ; level in which the ml is located     
         zoxy_P   = total( QoxyP(ii, jj, indP, tt), /NAN )
         zoxy_C   = total( QoxyC(ii, jj, indC, tt), /NAN )
         zo2s_P   = total( Qo2sP(ii, jj, indP, tt), /NAN )
         zo2s_C   = total( Qo2sC(ii, jj, indC, tt), /NAN )
         ; control simulation
         indm1 = indC(0) - 1
         IF (indm1 GE 0) THEN BEGIN 
            frac = ( gdepw(indC(0)) - zmldC ) / e3t(indm1)
            zoxy_C = zoxy_C + QoxyC(ii, jj, indm1, tt) * frac
            zo2s_C = zo2s_C + Qo2sC(ii, jj, indm1, tt) * frac
         ENDIF    
         ; in climate change simulation
         indm1 = indP(0) - 1
         IF (indm1 GE 0) THEN BEGIN 
            frac = ( gdepw(indP(0)) - zmldP ) / e3t(indm1)
            zoxy_P = zoxy_P + QoxyP(ii, jj, indm1, tt) * frac
            zo2s_P = zo2s_P + Qo2sP(ii, jj, indm1, tt) * frac
         ENDIF         
         ; convert from 1/mesh to 1/m2
         oxyP_umld(ii, jj, tt)   = zoxy_P / surf(ii, jj)
         oxyC_umld(ii, jj, tt)   = zoxy_C / surf(ii, jj)
         o2sP_umld(ii, jj, tt)   = zo2s_P / surf(ii, jj)
         o2sC_umld(ii, jj, tt)   = zo2s_C / surf(ii, jj)
         ; unit : kmolO2 * zfactoxy / m2
      END 
   END 
END

; remove overlap on the north pole and east-west
oxyP_umld = orca2_3d_overlap_remove(oxyP_umld)
oxyC_umld = orca2_3d_overlap_remove(oxyC_umld)
o2sP_umld = orca2_3d_overlap_remove(o2sP_umld)
o2sc_umld = orca2_3d_overlap_remove(o2sC_umld)
; horizontal integral
jpt = 110
toxyP = grossemoyenne(oxyP_umld, 'xy', /NAN, /INT)
toxyC = grossemoyenne(oxyC_umld, 'xy', /NAN, /INT)
to2sP = grossemoyenne(o2sP_umld, 'xy', /NAN, /INT)
to2sC = grossemoyenne(o2sC_umld, 'xy', /NAN, /INT)
; unit : kmolO2 * zfactoxy / y

;__________________________________________________
; Change

; drift in the control simulation 
doxyC = ( mean(toxyC[100:109]) - mean(toxyC[0:9]) )
; Change = trend in climate change simulation - drift in the control simulation
doxy = ( toxyP - mean(toxyP[0:9]) ) - ( toxyC - mean(toxyC[0:9]) )
do2s = ( to2sP - mean(to2sP[0:9]) ) - ( to2sC - mean(to2sC[0:9]) )

stop

;=================================================
; PRINT
;=================================================

print, 'doxyC  ', doxyc
print, 'tot    ', total(totC    , /NAN)  
print, 'TOc_dyn', total(TOc.dyn , /NAN)
print, 'bio    ', total(bioC    , /NAN)
print, 'subbio ', total(subbioC , /NAN)
print, 'expC   ', total(expC    , /NAN)
print, 'TOc_adv', total(TOc.adv , /NAN)
print, 'TOc_zmx', total(TOc.zmx , /NAN)
print, 'TOc_edd', total(TOc.edd , /NAN)
print, 'TOc_zad', total(TOc.zad , /NAN)
print, 'TOc_had', total(TOc.had , /NAN)
print, 'TOc_mld', total(TOc.mld , /NAN)
print, 'TOc_zdf', total(TOc.zdf , /NAN)
print, 'TOc_ldf', total(TOc.ldf , /NAN)
print, 'TOc_gm ', total(TOc.gm  , /NAN)

;=================================================
; PLOT
;=================================================


stop

; ________________________________________
; ARTICLE

w      = WINDOW(/BUFFER, DIMENSIONS = [300, 600]) ; dimension = [width, height]
x1        = 1990 + findgen(110)
props    = {BUFFER:1, OVERPLOT:1, THICK:1., FONT_SIZE:8, $
            XRANGE:[min(x1), max(x1)], XTITLE:'Time [years]', XTHICK:0.5,  $
            YRANGE: [-12, 10], YTICKINTERVAL:2, $
            YTITLE:'Deoxygenation     [PmolO2]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}

tmp1 = dTO_adv + dTO_gm + dTO_mld
tmp2 = dTO_ldf + dTO_zdf
p1    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 1], TITLE = 'Change in O2 fluxes')
p1_1   = plot(x1, dtot     , NAME = 'Total'       , COLOR = 'black'      , _EXTRA = props)  
p1_2   = plot(x1, doxy     , NAME = 'O2 ss MLD'   , COLOR = 'gray'       , _EXTRA = props)
p1_3   = plot(x1, dbio     , NAME = 'Respiration' , COLOR = 'green'      , _EXTRA = props)
p1_4   = plot(x1, tmp1     , NAME = 'Kinematic'   , COLOR = 'sky blue'   , _EXTRA = props)  
p1_5   = plot(x1, tmp2     , NAME = 'Diffusive '  , COLOR = 'tomato'     , _EXTRA = props)  

tmp1 = daTOsat_adv + daTOsat_gm + daTOsat_mld
tmp2 = daTOsat_ldf + daTOsat_zdf
p2    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 2], TITLE = 'Change in Osat fluxes')
p2_1   = plot(x1, daTOsat_dyn , NAME = 'Total'         , COLOR = 'black'   , LINESTYLE = '__', _EXTRA = props)
p2_2   = plot(x1, do2s        , NAME = 'Osat ss MLD'   , COLOR = 'gray'    , LINESTYLE = '__', _EXTRA = props)
p2_3   = plot(x1, tmp1        , NAME = 'Kinematic'     , COLOR = 'sky blue', LINESTYLE = '__', _EXTRA = props)
p2_4   = plot(x1, tmp2        , NAME = 'Diffusive'     , COLOR = 'tomato'  , LINESTYLE = '__', _EXTRA = props)

tmp1 = -(dTAOU_adv + dTAOU_gm + dTAOU_mld)
tmp2 = -(dTAOU_ldf + dTAOU_zdf)
tmp3 = tmp1 + tmp2 + dbio
tmp4 = -(do2s - doxy)
p3    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 3], TITLE = 'Change in residual fluxes')
p3_1   = plot(x1, tmp3     , NAME = 'Total'       , COLOR = 'black'   , LINESTYLE = '-.', _EXTRA = props)  
p3_2   = plot(x1, tmp4     , NAME = '-AOU ss MLD' , COLOR = 'gray'    , LINESTYLE = '-.', _EXTRA = props)
p3_3   = plot(x1, dbio     , NAME = 'Respiration' , COLOR = 'green'   , _EXTRA = props)
p3_4   = plot(x1, tmp1     , NAME = 'Kinematic'   , COLOR = 'sky blue', LINESTYLE = '-.', _EXTRA = props)  
p3_5   = plot(x1, tmp2     , NAME = 'Diffusive'   , COLOR = 'tomato'  , LINESTYLE = '-.', _EXTRA = props)

props_leg = {AUTO_TEXT_COLOR:1, FONT_SIZE:8}
leg = LEGEND(TARGET=[p1_1, p1_2, p1_3, p1_4, p1_5]            , POSITION = [0.7, 0.8], _EXTRA = props_leg)

ax          = p1.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       
ax          = p2.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       
ax          = p3.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

filename = '/data/daclod/FIG/' + namefigbase + '.eps'
p1.save, filename

stop

; ________________________________________
; ARTICLE BARPLOT

Sk = dTO_adv + dTO_gm + dTO_mld
Sd = dTO_ldf + dTO_zdf
Sko2s = daTOsat_adv + daTOsat_gm + daTOsat_mld
Sdo2s = daTOsat_ldf + daTOsat_zdf
Skaou = -(dTAOU_adv + dTAOU_gm + dTAOU_mld)
Sdaou = -(dTAOU_ldf + dTAOU_zdf)
daou = -(daTOsat_dyn - dtot)

tmp1 = [daTOsat_dyn[-1]]
tmp2 = [daou   [-1]]
tmp3 = [Sko2s  [-1]]
tmp4 = [Sdo2s  [-1]]
tmp5 = [Skaou  [-1]]
tmp6 = [Sdaou  [-1]]
tmp7 = [dbio   [-1]]

props    = {BUFFER:1, OVERPLOT:1, FONT_SIZE:8, NBARS:4, $
            LINESTYLE:'-', THICK:1., $
            XRANGE:[-1, 1], $
            XTHICK:0,  XTICKLEN:0, $
            YRANGE: [-14, 10], YTICKINTERVAL:2, $
            YTITLE:'Deoxygenation     [PmolO2]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}

w      = WINDOW(/BUFFER, DIMENSIONS = [300, 200]) ; dimension = [width, height]

p0 = plot(0.*findgen(2), /CURRENT, /BUFFER, LINESTYLE = 'none')
b01 = barplot(tmp1       , INDEX = 0                      , FILL_COLOR = 'slate gray', NAME = 'O2sat'       , _EXTRA = props)
b02 = barplot(tmp1 + tmp2, INDEX = 0, BOTTOM_VALUES = tmp1, FILL_COLOR = 'black'     , NAME = 'O2'          , _EXTRA = props)
b03 = barplot(tmp3       , INDEX = 1                      , FILL_COLOR = 'sky blue'  , NAME = 'Sk'          , _EXTRA = props)
b04 = barplot(tmp3 + tmp4, INDEX = 1, BOTTOM_VALUES = tmp3, FILL_COLOR = 'tomato'    , NAME = 'Sd'          , _EXTRA = props)
b05 = barplot(tmp5       , INDEX = 2                      , FILL_COLOR = 'sky blue'  , NAME = 'Sk'          , _EXTRA = props)
b06 = barplot(tmp5 + tmp6, INDEX = 2, BOTTOM_VALUES = tmp5, FILL_COLOR = 'tomato'    , NAME = 'Sd'          , _EXTRA = props)
b07 = barplot(tmp7       , INDEX = 3                      , FILL_COLOR = 'green'     , NAME = 'Respiration' , _EXTRA = props)

ax                   = p0.AXES
p0.XTICKVALUES       = [-0.375, -0.125, 0.125, 0.375]
p0.XTICKNAME         = ['O2', 'Temp', 'Circ', 'Bio']
ax[2].HIDE  = 1
ax[3].MAJOR = 0
ax[3].MINOR = 0

filename = '/data/daclod/FIG/' + namefigbase + '_barplot.eps'
p0.save, filename




stop

; ________________________________________
; TOUT
w      = WINDOW(/BUFFER, DIMENSIONS = [1000, 600]) ; dimension = [width, height]
x1        = 1990 + findgen(110)
props    = {BUFFER:1, OVERPLOT:1, LINESTYLE:'-', THICK:2.5, FONT_SIZE:6, $
            XRANGE:[min(x1), max(x1)+2], XTITLE:'Time [years]', XTHICK:0.5,  $
            YRANGE: [-14, 10], YTICKINTERVAL:2, $
            YTITLE:'Deoxygenation     [PmolO2]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}

p1    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 1])
p1_1   = plot(x1, dTO_dyn     , NAME = 'Ventilation'   , COLOR = 'royal blue'     ,  _EXTRA = props)
p1_2   = plot(x1, dbio        , NAME = 'Respiration'   , COLOR = 'green'          ,  _EXTRA = props)  
p1_3   = plot(x1, dtot        , NAME = 'Ventil + respi', COLOR = 'black'          ,  _EXTRA = props)  
p1_4   = plot(x1, doxy        , NAME = 'O2 ss MLD'     , COLOR = 'slate gray'     ,  _EXTRA = props)
p1_5   = plot(x1, daTOsat_dyn , NAME = 'Ventil Osat'   , COLOR = 'royal blue'     , THICK = 2.5, LINESTYLE = '--', /BUFFER, /OVERPLOT)
p1_6   = plot(x1, do2s        , NAME = 'Osat ss MLD'   , COLOR = 'slate gray'     , THICK = 2.5, LINESTYLE = '--', /BUFFER, /OVERPLOT)

p4    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 4])
p4_1   = plot(x1, dTO_had   , NAME = 'HAD'   , COLOR = 'sky blue'   , _EXTRA = props)  
p4_2   = plot(x1, dTO_zad   , NAME = 'ZAD'   , COLOR = 'dark violet', _EXTRA = props)
p4_3   = plot(x1, dTO_zdf   , NAME = 'ZDF'   , COLOR = 'tomato'     , _EXTRA = props)
p4_4   = plot(x1, dTO_mld   , NAME = 'MLD'   , COLOR = 'dark red'   , _EXTRA = props)
p4_5   = plot(x1, dTO_ldf   , NAME = 'LDF'   , COLOR = 'gold'       , _EXTRA = props)
p4_6   = plot(x1, dTO_gm    , NAME = 'GM'    , COLOR = 'deep pink'  , _EXTRA = props)

p5    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 5])
p5_1   = plot(x1, daTOsat_had   , NAME = 'HAD'   , COLOR = 'sky blue'   , _EXTRA = props)  
p5_2   = plot(x1, daTOsat_zad   , NAME = 'ZAD'   , COLOR = 'dark violet', _EXTRA = props)
p5_3   = plot(x1, daTOsat_zdf   , NAME = 'ZDF'   , COLOR = 'tomato'     , _EXTRA = props)
p5_4   = plot(x1, daTOsat_mld   , NAME = 'MLD'   , COLOR = 'dark red'   , _EXTRA = props)
p5_5   = plot(x1, daTOsat_ldf   , NAME = 'LDF'   , COLOR = 'gold'       , _EXTRA = props)
p5_6   = plot(x1, daTOsat_gm    , NAME = 'GM'    , COLOR = 'deep pink'  , _EXTRA = props)

p6    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 6])
p6_1   = plot(x1, dTAOU_had   , NAME = 'HAD'   , COLOR = 'sky blue'   , _EXTRA = props)  
p6_2   = plot(x1, dTAOU_zad   , NAME = 'ZAD'   , COLOR = 'dark violet', _EXTRA = props)
p6_3   = plot(x1, dTAOU_zdf   , NAME = 'ZDF'   , COLOR = 'tomato'     , _EXTRA = props)
p6_4   = plot(x1, dTAOU_mld   , NAME = 'MLD'   , COLOR = 'dark red'   , _EXTRA = props)
p6_5   = plot(x1, dTAOU_ldf   , NAME = 'LDF'   , COLOR = 'gold'       , _EXTRA = props)
p6_6   = plot(x1, dTAOU_gm    , NAME = 'GM'    , COLOR = 'deep pink'  , _EXTRA = props)

tmp1 = dTO_adv + dTO_gm
tmp2 = tmp1 + dTO_mld
tmp3 = dTO_ldf + dTO_zdf
p7    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 7])
p7_1   = plot(x1, dTO_adv   , NAME = 'ADV'           , COLOR = 'blue'       , _EXTRA = props)  
p7_2   = plot(x1, tmp1      , NAME = 'ADV + GM'      , COLOR = 'tan'        , _EXTRA = props)
p7_3   = plot(x1, tmp2      , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse' , _EXTRA = props)
p7_4   = plot(x1, tmp3      , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown', _EXTRA = props)

tmp1 = daTOsat_adv + daTOsat_gm
tmp2 = tmp1 + daTOsat_mld
tmp3 = daTOsat_ldf + daTOsat_zdf
p8    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 8])
p8_1   = plot(x1, daTOsat_adv, NAME = 'ADV'           , COLOR = 'blue'       , _EXTRA = props)  
p8_2   = plot(x1, tmp1       , NAME = 'ADV + GM'      , COLOR = 'tan'        , _EXTRA = props)
p8_3   = plot(x1, tmp2       , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse' , _EXTRA = props)
p8_4   = plot(x1, tmp3       , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown', _EXTRA = props)

tmp1 = dTAOU_adv + dTAOU_gm
tmp2 = tmp1 + dTAOU_mld
tmp3 = dTAOU_ldf + dTAOU_zdf
p9    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [3, 3, 9])
p9_1   = plot(x1, dTAOU_adv , NAME = 'ADV'           , COLOR = 'blue'       , _EXTRA = props)  
p9_2   = plot(x1, tmp1      , NAME = 'ADV + GM'      , COLOR = 'tan'        , _EXTRA = props)
p9_3   = plot(x1, tmp2      , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse' , _EXTRA = props)
p9_4   = plot(x1, tmp3      , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown', _EXTRA = props)

props_leg = {AUTO_TEXT_COLOR:1, FONT_SIZE:6}
leg = LEGEND(TARGET=[p1_1, p1_2, p1_3, p1_4]            , POSITION = [0.5, 0.9], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p4_1, p4_2, p4_3, p4_4, p4_5, p4_6], POSITION = [0.7, 0.9], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p7_1, p7_2, p7_3, p7_4]            , POSITION = [0.9, 0.9], _EXTRA = props_leg)

;; p1.TITLE  = 'Ventilation, respiration and delta O2'
ax          = p1.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

p4.TITLE  = 'O2'
p5.TITLE  = 'Osat'
p6.TITLE  = 'AOU'

ax          = p4.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0

ax          = p5.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

ax          = p6.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

ax          = p7.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

ax          = p8.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

ax          = p9.axes 
ax[2].MAJOR = 0       
ax[2].MINOR = 0       
ax[3].MAJOR = 0       
ax[3].MINOR = 0       

filename = '/data/daclod/FIG/' + namefigbase + '_tout.eps'
p1.save, filename

stop
END

