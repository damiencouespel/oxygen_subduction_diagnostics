;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;; FIG_BILAN_OFFLINE_ZONAL_MEAN ;;;;;;;;;;   
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Purpose 
; -------
; Plot meridional integrals of subduction and oxygen consumption in the control
; simulation and the cumulative change in the global warming
; simulation 

;==================================================
; PARAM
;==================================================

namefigbase = 'fig_bilan_offline_zonal_mean' ; base of the file name to save the figures 
                                                                            
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
zfactsub   =               1000. ; kmolO2 -> molO2
; conversion factor of organic matter subduction into O2 consumption
zfactbio   = -160./122.  * 1000. ; kmolC -> molO2  
; conversion factor of organic matter sinking into O2 consumption 
zfactexp   =  160./122.  * 1.e9  ; GmolC -> molO2
; conversion factor of ammonium subduction into O2 consumption (nitrification)
zfactnit   = - 40./122.  * 1000. ; kmolC -> molO2  nitrification amonium
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
zfactoxy = 1.e-6 ; to convert from kmol/m3 to Gmol/m3

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
   ; mldP = moyenne(mldP, 't', /nan)
   mldC = read_ncdf('somxl010', start1, end1, file=fdirO2 + fMLDC, /nostruct)
   ; mldC = moyenne(mldC, 't', /nan)
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

;__________________________________________________
; Oxygen content in each mesh

jpt = 110
Qoxy_P = FLTARR(jpi, jpj, jpk, jpt) & $ ; oxygen content in the climate change simulation
Qoxy_C = FLTARR(jpi, jpj, jpk, jpt) & $ ; oxygen content in the control simulation
FOR tt = 0, jpt-1 DO BEGIN & $ 
   Qoxy_P(*, *, *, tt) = oxy_P(*, *, *, tt) * tvolume & $
   Qoxy_C(*, *, *, tt) = oxy_C(*, *, *, tt) * tvolume & $
END

;__________________________________________________
; Oxygen content under the mixed layer depth

; initialization 
oxy_P_umld   = fltarr(jpi, jpj, jpt)
oxy_C_umld   = fltarr(jpi, jpj, jpt)

jpt = 110
FOR tt = 0, jpt-1 DO BEGIN & $
   FOR ii = 0, jpi-1 DO BEGIN & $
      FOR jj = 0, jpj-1 DO BEGIN & $
   
         zmldP = mldP(ii, jj, tt)
         zmldC = mldC(ii, jj, tt)
         
         ; index where mld is the deepest
         indP = where( gdepw GT zmldP )  & $
         indC = where( gdepw GT zmldC )  & $

         ; Oxygen content under MLD = all
         ; levels below MLD + prorata of the
         ; level in which the ml is located     

         zoxy_P   = total( Qoxy_P(ii, jj, indP, tt), /nan ) & $
         zoxy_C   = total( Qoxy_C(ii, jj, indC, tt), /nan ) & $
         
         indm1 = indP(0) - 1 & $
         IF (indm1 GT 0) THEN BEGIN & $
            frac = ( gdepw(indP(0)) - zmldP ) / e3t(indm1) & $
            zoxy_P   = zoxy_P + Qoxy_P(ii, jj, indm1, tt) * frac & $
         ENDIF           & $
         
         indm1 = indC(0) - 1 & $
         IF (indm1 GT 0) THEN BEGIN & $
            frac = ( gdepw(indC(0)) - zmldC ) / e3t(indm1) & $
            zoxy_C   = zoxy_C + Qoxy_C(ii, jj, indm1, tt) * frac & $
         ENDIF           & $
         
         ; convert from 1/mesh to 1/m2
         oxy_P_umld(ii, jj, tt) = zoxy_P / surf(ii, jj) & $
         oxy_C_umld(ii, jj, tt) = zoxy_C / surf(ii, jj) & $

      END & $
   END & $
END

;__________________________________________________
; Zonal mean

; to remove the overlap in the indian ocean dc2
; and north pole
oxy_P_umld = orca2_3d_overlap_remove(oxy_P_umld)
oxy_C_umld = orca2_3d_overlap_remove(oxy_C_umld)

oxy_P_umld = grossemoyenne(oxy_P_umld, 'x', /NAN, /INT)
oxy_C_umld = grossemoyenne(oxy_C_umld, 'x', /NAN, /INT)
o2s_P_umld = grossemoyenne(o2s_P_umld, 'x', /NAN, /INT)
o2s_C_umld = grossemoyenne(o2s_C_umld, 'x', /NAN, /INT)

;__________________________________________________
; Change

; drift in the control simulation
doxy_ctl =  mean(oxy_C_umld[*, 100:109], dimension = 2) - mean(oxy_C_umld[*, 0:9], dimension = 2) 
; Change = trend in climate change simulation - drift in the control simulation
doxy = ( mean(oxy_P_umld[*, 100:109], dimension = 2)  - mean(oxy_P_umld[*, 0:9], dimension = 2)) - ( mean(oxy_C_umld[*, 100:109], dimension = 2) - mean(oxy_C_umld[*, 0:9], dimension = 2) )

; __________________
; Meridional running mean

doxy_ctl_smooth = FLTARR(149)
FOR ind = 0, 138 DO doxy_Ctl_smooth(ind+5) = MEAN(doxy_ctl[ind:ind+9])
doxy_smooth = FLTARR(149)
FOR ind = 0, 138 DO doxy_smooth(ind+5) = MEAN(doxy[ind:ind+9])

Z_doxy = doxy
Z_doxy_smooth = doxy_smooth

;==================================================
; PROCESS SUBDUCTION AND EXPORT DATA
;==================================================

@common
@init_run_clim_dc
surf = e1t*e2t ; mesh surface

;__________________________________________________
; Oxygen subduction in the control simulation

subdu = restore_subdu(fdir, fTOC, 1e13, FILEZDF1 = fzdfTOc, FILEZDF2 = fzdf2TOc) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TOc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactsub / surf ; lateral advection                               
TOc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactsub / surf ; vertical advection                              
TOc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactsub / surf ; entrainment                                                        
TOc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactsub / surf ; vertical diffusion                                                 
TOc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactsub / surf ; advection by bolus velocities (Gent McWilliams)                     
TOc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactsub / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TOc_edd = TOc_ldf + TOc_gm  ; eddies          
TOc_zmx = TOc_zdf + TOc_mld ; vertical mixing 
TOc_adv = TOc_had + TOc_zad ; advection       
TOc_dyn = TOc_adv + TOc_zmx + TOc_edd

;__________________________________________________
; Saturated oxygen subduction in the control simulation

subdu = restore_subdu(fdir, faTOsatc, 1e13)
; unit : kmolO2 / mesh / month

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
aTOsatc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactsub / surf ; lateral advection                               
aTOsatc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactsub / surf ; vertical advection                              
aTOsatc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactsub / surf ; entrainment                                                        
aTOsatc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactsub / surf ; vertical diffusion                                                 
aTOsatc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactsub / surf ; advection by bolus velocities (Gent McWilliams)                     
aTOsatc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactsub / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

aTOsatc_edd = aTOsatc_ldf + aTOsatc_gm  ; eddies          
aTOsatc_zmx = aTOsatc_zdf + aTOsatc_mld ; vertical mixing 
aTOsatc_adv = aTOsatc_had + aTOsatc_zad ; advection       
aTOsatc_dyn = aTOsatc_adv + aTOsatc_zmx + aTOsatc_edd

;__________________________________________________
; AOU subduction = sub O2sat - sub O2 in the control simulation

TAOUc_had = aTOsatc_had - TOc_had  ; lateral advection                               
TAOUc_zad = aTOsatc_zad - TOc_zad  ; vertical advection                              
TAOUc_mld = aTOsatc_mld - TOc_mld  ; entrainment                                     
TAOUc_zdf = aTOsatc_zdf - TOc_zdf  ; vertical diffusion                              
TAOUc_gm  = aTOsatc_gm  - TOc_gm   ; advection by bolus velocities (Gent McWilliams) 
TAOUc_ldf = aTOsatc_ldf - TOc_ldf  ; isopycnal diffusion                             

TAOUc_edd = TAOUc_ldf + TAOUc_gm  ; eddies          
TAOUc_zmx = TAOUc_zdf + TAOUc_mld ; vertical mixing 
TAOUc_adv = TAOUc_had + TAOUc_zad ; advection       
TAOUc_dyn = TAOUc_adv + TAOUc_zmx + TAOUc_edd


;__________________________________________________
; DOC subduction in the control simulation

subdu = restore_subdu(fdir, fTDOCc, 1e12)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TDOCc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio / surf ; lateral advection                               
TDOCc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio / surf ; vertical advection                              
TDOCc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio / surf ; entrainment                                                        
TDOCc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio / surf ; vertical diffusion                                                 
TDOCc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio / surf ; advection by bolus velocities (Gent McWilliams)                     
TDOCc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TDOCc_edd = TDOCc_ldf + TDOCc_gm  ; eddies          
TDOCc_zmx = TDOCc_zdf + TDOCc_mld ; vertical mixing 
TDOCc_adv = TDOCc_had + TDOCc_zad ; advection       
TDOCc_dyn = TDOCc_adv + TDOCc_zmx + TDOCc_edd

;__________________________________________________
; POC subduction in the control simulation

subdu = restore_subdu(fdir, fTPOCc, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TPOCc_had = grossemoyenne(subh  , 't', /INT, /NAN)  * zfactbio / surf ; lateral advection                               
TPOCc_zad = grossemoyenne(subw  , 't', /INT, /NAN)  * zfactbio / surf ; vertical advection                              
TPOCc_mld = grossemoyenne(submld, 't', /INT, /NAN)  * zfactbio / surf ; entrainment                                                        
TPOCc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN)  * zfactbio / surf ; vertical diffusion                                                 
TPOCc_gm  = grossemoyenne(subgm , 't', /INT, /NAN)  * zfactbio / surf ; advection by bolus velocities (Gent McWilliams)                     
TPOCc_ldf = grossemoyenne(subldf, 't', /INT, /NAN)  * zfactbio / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TPOCc_edd = TPOCc_ldf + TPOCc_gm  ; eddies          
TPOCc_zmx = TPOCc_zdf + TPOCc_mld ; vertical mixing 
TPOCc_adv = TPOCc_had + TPOCc_zad ; advection       
TPOCc_dyn = TPOCc_adv + TPOCc_zmx + TPOCc_edd

;__________________________________________________
; GOC subduction in the control simulation

subdu = restore_subdu(fdir, fTGOCc, 1e9)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TGOCc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TGOCc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TGOCc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TGOCc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TGOCc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TGOCc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TGOCc_edd = TGOCc_ldf + TGOCc_gm  ; eddies          
TGOCc_zmx = TGOCc_zdf + TGOCc_mld ; vertical mixing 
TGOCc_adv = TGOCc_had + TGOCc_zad ; advection       
TGOCc_dyn = TGOCc_adv + TGOCc_zmx + TGOCc_edd

;__________________________________________________
; ZOO subduction in the control simulation

subdu = restore_subdu(fdir, fTZOOc, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TZOOc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TZOOc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TZOOc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TZOOc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TZOOc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TZOOc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TZOOc_edd = TZOOc_ldf + TZOOc_gm  ; eddies          
TZOOc_zmx = TZOOc_zdf + TZOOc_mld ; vertical mixing 
TZOOc_adv = TZOOc_had + TZOOc_zad ; advection       
TZOOc_dyn = TZOOc_adv + TZOOc_zmx + TZOOc_edd

;__________________________________________________
; ZOO2 subduction in the control simulation

subdu = restore_subdu(fdir, fTZOO2c, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TZOO2c_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TZOO2c_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TZOO2c_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TZOO2c_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TZOO2c_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TZOO2c_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TZOO2c_edd = TZOO2c_ldf + TZOO2c_gm  ; eddies          
TZOO2c_zmx = TZOO2c_zdf + TZOO2c_mld ; vertical mixing 
TZOO2c_adv = TZOO2c_had + TZOO2c_zad ; advection       
TZOO2c_dyn = TZOO2c_adv + TZOO2c_zmx + TZOO2c_edd

;__________________________________________________
; Total phytoplankton subduction in the control simulation

subdu = restore_subdu(fdir, fTPHYTc, 1e11)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TPHYTc_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TPHYTc_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TPHYTc_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TPHYTc_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TPHYTc_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TPHYTc_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TPHYTc_edd = TPHYTc_ldf + TPHYTc_gm  ; eddies          
TPHYTc_zmx = TPHYTc_zdf + TPHYTc_mld ; vertical mixing 
TPHYTc_adv = TPHYTc_had + TPHYTc_zad ; advection       
TPHYTc_dyn = TPHYTc_adv + TPHYTc_zmx + TPHYTc_edd

;__________________________________________________
; NH4 subduction in the control simulation

subdu = restore_subdu(fdir, fTNH4c, 1e11)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TNH4c_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactnit  / surf ; lateral advection                               
TNH4c_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactnit  / surf ; vertical advection                              
TNH4c_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactnit  / surf ; entrainment                                                        
TNH4c_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactnit  / surf ; vertical diffusion                                                 
TNH4c_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactnit  / surf ; advection by bolus velocities (Gent McWilliams)                     
TNH4c_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactnit  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TNH4c_edd = TNH4c_ldf + TNH4c_gm  ; eddies          
TNH4c_zmx = TNH4c_zdf + TNH4c_mld ; vertical mixing 
TNH4c_adv = TNH4c_had + TNH4c_zad ; advection       
TNH4c_dyn = TNH4c_adv + TNH4c_zmx + TNH4c_edd

;__________________________________________________
; POC and GOC sinking in the control simulation

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

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
expC = grossemoyenne(export, 't', /NAN, /INT) * zfactexp / surf
; unit : PmolO2 / m2 / 110year

;__________________________________________________
; SUM

; O2 consumpution due to organic matter subduction in the control simulation
subbioC    = TDOCc_dyn + TPOCc_dyn + TGOCc_dyn + TZOOc_dyn + TZOO2c_dyn + TPHYTc_dyn + TNH4c_dyn
; O2 consumption due to total export in the control simulation
bioC       = subbioC + expC
; total = O2 subduction + O2 consumption in the control simulation
totC       = bioC + TOC_dyn


;__________________________________________________
; Oxygen subduction in the climate change simulation

subdu = restore_subdu(fdir, fTOP, 1e13, FILEZDF1 = fzdfTOp, FILEZDF2 = fzdf2TOp) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TOp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactsub / surf ; lateral advection                               
TOp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactsub / surf ; vertical advection                              
TOp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactsub / surf ; entrainment                                                        
TOp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactsub / surf ; vertical diffusion                                                 
TOp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactsub / surf ; advection by bolus velocities (Gent McWilliams)                     
TOp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactsub / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TOp_edd = TOp_ldf + TOp_gm  ; eddies          
TOp_zmx = TOp_zdf + TOp_mld ; vertical mixing 
TOp_adv = TOp_had + TOp_zad ; advection       
TOp_dyn = TOp_adv + TOp_zmx + TOp_edd

;__________________________________________________
; Saturated oxygen subduction in the climate change simulation

subdu = restore_subdu(fdir, faTOsatp, 1e13)
; unit : kmolO2 / mesh / month

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
aTOsatp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactsub / surf ; lateral advection                               
aTOsatp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactsub / surf ; vertical advection                              
aTOsatp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactsub / surf ; entrainment                                                        
aTOsatp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactsub / surf ; vertical diffusion                                                 
aTOsatp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactsub / surf ; advection by bolus velocities (Gent McWilliams)                     
aTOsatp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactsub / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

aTOsatp_edd = aTOsatp_ldf + aTOsatp_gm  ; eddies          
aTOsatp_zmx = aTOsatp_zdf + aTOsatp_mld ; vertical mixing 
aTOsatp_adv = aTOsatp_had + aTOsatp_zad ; advection       
aTOsatp_dyn = aTOsatp_adv + aTOsatp_zmx + aTOsatp_edd

;__________________________________________________
; AOU subduction = sub O2sat - sub O2 in the climate change simulation

TAOUp_had = aTOsatp_had - TOp_had  ; lateral advection                               
TAOUp_zad = aTOsatp_zad - TOp_zad  ; vertical advection                              
TAOUp_mld = aTOsatp_mld - TOp_mld  ; entrainment                                     
TAOUp_zdf = aTOsatp_zdf - TOp_zdf  ; vertical diffusion                              
TAOUp_gm  = aTOsatp_gm  - TOp_gm   ; advection by bolus velocities (Gent McWilliams) 
TAOUp_ldf = aTOsatp_ldf - TOp_ldf  ; isopycnal diffusion                             

TAOUp_edd = TAOUp_ldf + TAOUp_gm  ; eddies          
TAOUp_zmx = TAOUp_zdf + TAOUp_mld ; vertical mixing 
TAOUp_adv = TAOUp_had + TAOUp_zad ; advection       
TAOUp_dyn = TAOUp_adv + TAOUp_zmx + TAOUp_edd

;__________________________________________________
; DOC subduction in the climate change simulation

subdu = restore_subdu(fdir, fTDOCp, 1e12)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TDOCp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio / surf ; lateral advection                               
TDOCp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio / surf ; vertical advection                              
TDOCp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio / surf ; entrainment                                                        
TDOCp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio / surf ; vertical diffusion                                                 
TDOCp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio / surf ; advection by bolus velocities (Gent McWilliams)                     
TDOCp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TDOCp_edd = TDOCp_ldf + TDOCp_gm  ; eddies          
TDOCp_zmx = TDOCp_zdf + TDOCp_mld ; vertical mixing 
TDOCp_adv = TDOCp_had + TDOCp_zad ; advection       
TDOCp_dyn = TDOCp_adv + TDOCp_zmx + TDOCp_edd

;__________________________________________________
; POC subduction in the climate change simulation

subdu = restore_subdu(fdir, fTPOCp, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TPOCp_had = grossemoyenne(subh  , 't', /INT, /NAN)  * zfactbio / surf ; lateral advection                               
TPOCp_zad = grossemoyenne(subw  , 't', /INT, /NAN)  * zfactbio / surf ; vertical advection                              
TPOCp_mld = grossemoyenne(submld, 't', /INT, /NAN)  * zfactbio / surf ; entrainment                                                        
TPOCp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN)  * zfactbio / surf ; vertical diffusion                                                 
TPOCp_gm  = grossemoyenne(subgm , 't', /INT, /NAN)  * zfactbio / surf ; advection by bolus velocities (Gent McWilliams)                     
TPOCp_ldf = grossemoyenne(subldf, 't', /INT, /NAN)  * zfactbio / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TPOCp_edd = TPOCp_ldf + TPOCp_gm  ; eddies          
TPOCp_zmx = TPOCp_zdf + TPOCp_mld ; vertical mixing 
TPOCp_adv = TPOCp_had + TPOCp_zad ; advection       
TPOCp_dyn = TPOCp_adv + TPOCp_zmx + TPOCp_edd

;__________________________________________________
; GOC subduction in the climate change simulation

subdu = restore_subdu(fdir, fTGOCp, 1e9)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TGOCp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TGOCp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TGOCp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TGOCp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TGOCp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TGOCp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TGOCp_edd = TGOCp_ldf + TGOCp_gm  ; eddies          
TGOCp_zmx = TGOCp_zdf + TGOCp_mld ; vertical mixing 
TGOCp_adv = TGOCp_had + TGOCp_zad ; advection       
TGOCp_dyn = TGOCp_adv + TGOCp_zmx + TGOCp_edd

;__________________________________________________
; ZOO subduction in the climate change simulation

subdu = restore_subdu(fdir, fTZOOp, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TZOOp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TZOOp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TZOOp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TZOOp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TZOOp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TZOOp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TZOOp_edd = TZOOp_ldf + TZOOp_gm  ; eddies          
TZOOp_zmx = TZOOp_zdf + TZOOp_mld ; vertical mixing 
TZOOp_adv = TZOOp_had + TZOOp_zad ; advection       
TZOOp_dyn = TZOOp_adv + TZOOp_zmx + TZOOp_edd

;__________________________________________________
; ZOO2 subduction in the climate change simulation

subdu = restore_subdu(fdir, fTZOO2p, 1e10)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TZOO2p_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TZOO2p_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TZOO2p_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TZOO2p_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TZOO2p_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TZOO2p_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TZOO2p_edd = TZOO2p_ldf + TZOO2p_gm  ; eddies          
TZOO2p_zmx = TZOO2p_zdf + TZOO2p_mld ; vertical mixing 
TZOO2p_adv = TZOO2p_had + TZOO2p_zad ; advection       
TZOO2p_dyn = TZOO2p_adv + TZOO2p_zmx + TZOO2p_edd

;__________________________________________________
; Total phytoplankton subduction in the climate change simulation

subdu = restore_subdu(fdir, fTPHYTp, 1e11)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TPHYTp_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio  / surf ; lateral advection                               
TPHYTp_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio  / surf ; vertical advection                              
TPHYTp_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio  / surf ; entrainment                                                        
TPHYTp_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio  / surf ; vertical diffusion                                                 
TPHYTp_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio  / surf ; advection by bolus velocities (Gent McWilliams)                     
TPHYTp_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TPHYTp_edd = TPHYTp_ldf + TPHYTp_gm  ; eddies          
TPHYTp_zmx = TPHYTp_zdf + TPHYTp_mld ; vertical mixing 
TPHYTp_adv = TPHYTp_had + TPHYTp_zad ; advection       
TPHYTp_dyn = TPHYTp_adv + TPHYTp_zmx + TPHYTp_edd

;__________________________________________________
; NH4 subduction in the climate change simulation

subdu = restore_subdu(fdir, fTNH4p, 1e11)

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolC / mesh / month

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
TNH4p_had = grossemoyenne(subh  , 't', /INT, /NAN) * zfactnit  / surf ; lateral advection                               
TNH4p_zad = grossemoyenne(subw  , 't', /INT, /NAN) * zfactnit  / surf ; vertical advection                              
TNH4p_mld = grossemoyenne(submld, 't', /INT, /NAN) * zfactnit  / surf ; entrainment                                                        
TNH4p_zdf = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactnit  / surf ; vertical diffusion                                                 
TNH4p_gm  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactnit  / surf ; advection by bolus velocities (Gent McWilliams)                     
TNH4p_ldf = grossemoyenne(subldf, 't', /INT, /NAN) * zfactnit  / surf ; isopycnal diffusion                                                
; unit : PmolO2 / m2 / 110year

TNH4p_edd = TNH4p_ldf + TNH4p_gm  ; eddies          
TNH4p_zmx = TNH4p_zdf + TNH4p_mld ; vertical mixing 
TNH4p_adv = TNH4p_had + TNH4p_zad ; advection       
TNH4p_dyn = TNH4p_adv + TNH4p_zmx + TNH4p_edd

;__________________________________________________
; POC and GOC sinking in the climate change simulation

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

; TEMPORAL INTEGRAL
jpt = 1320
time = findgen(jpt)
expP = grossemoyenne(export, 't', /NAN, /INT) * zfactexp / surf
; unit : PmolO2 / m2 / 110year

;__________________________________________________
; SUM

; O2 consumpution due to organic matter subduction in the climate
; change simulation
subbioP    = TDOCp_dyn + TPOCp_dyn + TGOCp_dyn + TZOOp_dyn + TZOO2p_dyn + TPHYTp_dyn + TNH4p_dyn
; O2 consumption due to total export in the climate change simulation
bioP       = subbioP + expP
; total = O2 subduction + O2 consumption in the climate change simulation
totP       = bioP + TOp_dyn

;==================================================
; ZONAL INTEGRAL
;==================================================

xint_TOc_dyn = moyenne(TOc_dyn, 'x', /INT, /NAN)
xint_TOc_edd = moyenne(TOc_edd, 'x', /INT, /NAN)
xint_TOc_zmx = moyenne(TOc_zmx, 'x', /INT, /NAN)
xint_TOc_adv = moyenne(TOc_adv, 'x', /INT, /NAN)
xint_TOc_had = moyenne(TOc_had, 'x', /INT, /NAN)
xint_TOc_zad = moyenne(TOc_zad, 'x', /INT, /NAN)
xint_TOc_mld = moyenne(TOc_mld, 'x', /INT, /NAN)
xint_TOc_zdf = moyenne(TOc_zdf, 'x', /INT, /NAN)
xint_TOc_ldf = moyenne(TOc_ldf, 'x', /INT, /NAN)
xint_TOc_gm  = moyenne(TOc_gm , 'x', /INT, /NAN)

xint_TOp_dyn = moyenne(TOp_dyn, 'x', /INT, /NAN)
xint_TOp_edd = moyenne(TOp_edd, 'x', /INT, /NAN)
xint_TOp_zmx = moyenne(TOp_zmx, 'x', /INT, /NAN)
xint_TOp_adv = moyenne(TOp_adv, 'x', /INT, /NAN)
xint_TOp_had = moyenne(TOp_had, 'x', /INT, /NAN)
xint_TOp_zad = moyenne(TOp_zad, 'x', /INT, /NAN)
xint_TOp_mld = moyenne(TOp_mld, 'x', /INT, /NAN)
xint_TOp_zdf = moyenne(TOp_zdf, 'x', /INT, /NAN)
xint_TOp_ldf = moyenne(TOp_ldf, 'x', /INT, /NAN)
xint_TOp_gm  = moyenne(TOp_gm , 'x', /INT, /NAN)

xint_aTOsatc_dyn = moyenne(aTOsatc_dyn, 'x', /INT, /NAN)
xint_aTOsatc_edd = moyenne(aTOsatc_edd, 'x', /INT, /NAN)
xint_aTOsatc_zmx = moyenne(aTOsatc_zmx, 'x', /INT, /NAN)
xint_aTOsatc_adv = moyenne(aTOsatc_adv, 'x', /INT, /NAN)
xint_aTOsatc_had = moyenne(aTOsatc_had, 'x', /INT, /NAN)
xint_aTOsatc_zad = moyenne(aTOsatc_zad, 'x', /INT, /NAN)
xint_aTOsatc_mld = moyenne(aTOsatc_mld, 'x', /INT, /NAN)
xint_aTOsatc_zdf = moyenne(aTOsatc_zdf, 'x', /INT, /NAN)
xint_aTOsatc_ldf = moyenne(aTOsatc_ldf, 'x', /INT, /NAN)
xint_aTOsatc_gm  = moyenne(aTOsatc_gm , 'x', /INT, /NAN)

xint_aTOsatp_dyn = moyenne(aTOsatp_dyn, 'x', /INT, /NAN)
xint_aTOsatp_edd = moyenne(aTOsatp_edd, 'x', /INT, /NAN)
xint_aTOsatp_zmx = moyenne(aTOsatp_zmx, 'x', /INT, /NAN)
xint_aTOsatp_adv = moyenne(aTOsatp_adv, 'x', /INT, /NAN)
xint_aTOsatp_had = moyenne(aTOsatp_had, 'x', /INT, /NAN)
xint_aTOsatp_zad = moyenne(aTOsatp_zad, 'x', /INT, /NAN)
xint_aTOsatp_mld = moyenne(aTOsatp_mld, 'x', /INT, /NAN)
xint_aTOsatp_zdf = moyenne(aTOsatp_zdf, 'x', /INT, /NAN)
xint_aTOsatp_ldf = moyenne(aTOsatp_ldf, 'x', /INT, /NAN)
xint_aTOsatp_gm  = moyenne(aTOsatp_gm , 'x', /INT, /NAN)

xint_TAOUc_dyn = moyenne(TAOUc_dyn, 'x', /INT, /NAN)
xint_TAOUc_edd = moyenne(TAOUc_edd, 'x', /INT, /NAN)
xint_TAOUc_zmx = moyenne(TAOUc_zmx, 'x', /INT, /NAN)
xint_TAOUc_adv = moyenne(TAOUc_adv, 'x', /INT, /NAN)
xint_TAOUc_had = moyenne(TAOUc_had, 'x', /INT, /NAN)
xint_TAOUc_zad = moyenne(TAOUc_zad, 'x', /INT, /NAN)
xint_TAOUc_mld = moyenne(TAOUc_mld, 'x', /INT, /NAN)
xint_TAOUc_zdf = moyenne(TAOUc_zdf, 'x', /INT, /NAN)
xint_TAOUc_ldf = moyenne(TAOUc_ldf, 'x', /INT, /NAN)
xint_TAOUc_gm  = moyenne(TAOUc_gm , 'x', /INT, /NAN)

xint_TAOUp_dyn = moyenne(TAOUp_dyn, 'x', /INT, /NAN)
xint_TAOUp_edd = moyenne(TAOUp_edd, 'x', /INT, /NAN)
xint_TAOUp_zmx = moyenne(TAOUp_zmx, 'x', /INT, /NAN)
xint_TAOUp_adv = moyenne(TAOUp_adv, 'x', /INT, /NAN)
xint_TAOUp_had = moyenne(TAOUp_had, 'x', /INT, /NAN)
xint_TAOUp_zad = moyenne(TAOUp_zad, 'x', /INT, /NAN)
xint_TAOUp_mld = moyenne(TAOUp_mld, 'x', /INT, /NAN)
xint_TAOUp_zdf = moyenne(TAOUp_zdf, 'x', /INT, /NAN)
xint_TAOUp_ldf = moyenne(TAOUp_ldf, 'x', /INT, /NAN)
xint_TAOUp_gm  = moyenne(TAOUp_gm , 'x', /INT, /NAN)

xint_bioC = moyenne(bioC, 'x', /INT, /NAN)
xint_bioP = moyenne(bioP, 'x', /INT, /NAN)

xint_totC = moyenne(totC, 'x', /INT, /NAN)
xint_totP = moyenne(totP, 'x', /INT, /NAN)

;==================================================
; Change = climate change - control
;==================================================

Dxint_TO_dyn = xint_TOp_dyn - xint_TOc_dyn 
Dxint_TO_edd = xint_TOp_edd - xint_TOc_edd 
Dxint_TO_zmx = xint_TOp_zmx - xint_TOc_zmx 
Dxint_TO_adv = xint_TOp_adv - xint_TOc_adv 
Dxint_TO_had = xint_TOp_had - xint_TOc_had 
Dxint_TO_zad = xint_TOp_zad - xint_TOc_zad 
Dxint_TO_mld = xint_TOp_mld - xint_TOc_mld 
Dxint_TO_zdf = xint_TOp_zdf - xint_TOc_zdf 
Dxint_TO_ldf = xint_TOp_ldf - xint_TOc_ldf 
Dxint_TO_gm  = xint_TOp_gm  - xint_TOc_gm 

Dxint_aTOsat_dyn = xint_aTOsatp_dyn - xint_aTOsatc_dyn 
Dxint_aTOsat_edd = xint_aTOsatp_edd - xint_aTOsatc_edd 
Dxint_aTOsat_zmx = xint_aTOsatp_zmx - xint_aTOsatc_zmx 
Dxint_aTOsat_adv = xint_aTOsatp_adv - xint_aTOsatc_adv 
Dxint_aTOsat_had = xint_aTOsatp_had - xint_aTOsatc_had 
Dxint_aTOsat_zad = xint_aTOsatp_zad - xint_aTOsatc_zad 
Dxint_aTOsat_mld = xint_aTOsatp_mld - xint_aTOsatc_mld 
Dxint_aTOsat_zdf = xint_aTOsatp_zdf - xint_aTOsatc_zdf 
Dxint_aTOsat_ldf = xint_aTOsatp_ldf - xint_aTOsatc_ldf 
Dxint_aTOsat_gm  = xint_aTOsatp_gm  - xint_aTOsatc_gm 

Dxint_TAOU_dyn = xint_TAOUp_dyn - xint_TAOUc_dyn 
Dxint_TAOU_edd = xint_TAOUp_edd - xint_TAOUc_edd 
Dxint_TAOU_zmx = xint_TAOUp_zmx - xint_TAOUc_zmx 
Dxint_TAOU_adv = xint_TAOUp_adv - xint_TAOUc_adv 
Dxint_TAOU_had = xint_TAOUp_had - xint_TAOUc_had 
Dxint_TAOU_zad = xint_TAOUp_zad - xint_TAOUc_zad 
Dxint_TAOU_mld = xint_TAOUp_mld - xint_TAOUc_mld 
Dxint_TAOU_zdf = xint_TAOUp_zdf - xint_TAOUc_zdf 
Dxint_TAOU_ldf = xint_TAOUp_ldf - xint_TAOUc_ldf 
Dxint_TAOU_gm  = xint_TAOUp_gm  - xint_TAOUc_gm 
                                            
Dxint_bio = xint_bioP - xint_bioC
Dxint_tot = xint_totP - xint_totC


;==================================================
; SMOOTH : RUNNING MEAN
;==================================================

xint_totC_smooth = FLTARR(149)
xint_bioC_smooth = FLTARR(149)
xint_TOc_dyn_smooth = FLTARR(149)
xint_TOc_edd_smooth = FLTARR(149)
xint_TOc_zmx_smooth = FLTARR(149)
xint_TOc_adv_smooth = FLTARR(149)
xint_TOc_had_smooth = FLTARR(149)
xint_TOc_zad_smooth = FLTARR(149)
xint_TOc_mld_smooth = FLTARR(149)
xint_TOc_zdf_smooth = FLTARR(149)
xint_TOc_ldf_smooth = FLTARR(149)
xint_TOc_gm_smooth  = FLTARR(149)

xint_totP_smooth = FLTARR(149)
xint_bioP_smooth = FLTARR(149)
xint_TOp_dyn_smooth = FLTARR(149)
xint_TOp_edd_smooth = FLTARR(149)
xint_TOp_zmx_smooth = FLTARR(149)
xint_TOp_adv_smooth = FLTARR(149)
xint_TOp_had_smooth = FLTARR(149)
xint_TOp_zad_smooth = FLTARR(149)
xint_TOp_mld_smooth = FLTARR(149)
xint_TOp_zdf_smooth = FLTARR(149)
xint_TOp_ldf_smooth = FLTARR(149)
xint_TOp_gm_smooth  = FLTARR(149)

xint_aTOsatc_dyn_smooth = FLTARR(149)
xint_aTOsatc_edd_smooth = FLTARR(149)
xint_aTOsatc_zmx_smooth = FLTARR(149)
xint_aTOsatc_adv_smooth = FLTARR(149)
xint_aTOsatc_had_smooth = FLTARR(149)
xint_aTOsatc_zad_smooth = FLTARR(149)
xint_aTOsatc_mld_smooth = FLTARR(149)
xint_aTOsatc_zdf_smooth = FLTARR(149)
xint_aTOsatc_ldf_smooth = FLTARR(149)
xint_aTOsatc_gm_smooth  = FLTARR(149)

xint_TAOUc_dyn_smooth = FLTARR(149)
xint_TAOUc_edd_smooth = FLTARR(149)
xint_TAOUc_zmx_smooth = FLTARR(149)
xint_TAOUc_adv_smooth = FLTARR(149)
xint_TAOUc_had_smooth = FLTARR(149)
xint_TAOUc_zad_smooth = FLTARR(149)
xint_TAOUc_mld_smooth = FLTARR(149)
xint_TAOUc_zdf_smooth = FLTARR(149)
xint_TAOUc_ldf_smooth = FLTARR(149)
xint_TAOUc_gm_smooth  = FLTARR(149)

Dxint_tot_smooth = FLTARR(149)
Dxint_bio_smooth = FLTARR(149)
Dxint_TO_dyn_smooth = FLTARR(149)
Dxint_TO_edd_smooth = FLTARR(149)
Dxint_TO_zmx_smooth = FLTARR(149)
Dxint_TO_adv_smooth = FLTARR(149)
Dxint_TO_had_smooth = FLTARR(149)
Dxint_TO_zad_smooth = FLTARR(149)
Dxint_TO_mld_smooth = FLTARR(149)
Dxint_TO_zdf_smooth = FLTARR(149)
Dxint_TO_ldf_smooth = FLTARR(149)
Dxint_TO_gm_smooth  = FLTARR(149)

Dxint_aTOsat_dyn_smooth = FLTARR(149)
Dxint_aTOsat_edd_smooth = FLTARR(149)
Dxint_aTOsat_zmx_smooth = FLTARR(149)
Dxint_aTOsat_adv_smooth = FLTARR(149)
Dxint_aTOsat_had_smooth = FLTARR(149)
Dxint_aTOsat_zad_smooth = FLTARR(149)
Dxint_aTOsat_mld_smooth = FLTARR(149)
Dxint_aTOsat_zdf_smooth = FLTARR(149)
Dxint_aTOsat_ldf_smooth = FLTARR(149)
Dxint_aTOsat_gm_smooth  = FLTARR(149)

Dxint_TAOU_dyn_smooth = FLTARR(149)
Dxint_TAOU_edd_smooth = FLTARR(149)
Dxint_TAOU_zmx_smooth = FLTARR(149)
Dxint_TAOU_adv_smooth = FLTARR(149)
Dxint_TAOU_had_smooth = FLTARR(149)
Dxint_TAOU_zad_smooth = FLTARR(149)
Dxint_TAOU_mld_smooth = FLTARR(149)
Dxint_TAOU_zdf_smooth = FLTARR(149)
Dxint_TAOU_ldf_smooth = FLTARR(149)
Dxint_TAOU_gm_smooth  = FLTARR(149)

FOR ind = 0, 138 DO BEGIN 
   xint_totC_smooth(ind+5) = MEAN(xint_totC(ind:ind+9))
   xint_bioC_smooth(ind+5) = MEAN(xint_bioC(ind:ind+9))
   xint_TOc_dyn_smooth(ind+5) = MEAN(xint_TOc_dyn(ind:ind+9))
   xint_TOc_edd_smooth(ind+5) = MEAN(xint_TOc_edd(ind:ind+9))
   xint_TOc_zmx_smooth(ind+5) = MEAN(xint_TOc_zmx(ind:ind+9))
   xint_TOc_adv_smooth(ind+5) = MEAN(xint_TOc_adv(ind:ind+9))
   xint_TOc_had_smooth(ind+5) = MEAN(xint_TOc_had(ind:ind+9))
   xint_TOc_zad_smooth(ind+5) = MEAN(xint_TOc_zad(ind:ind+9))
   xint_TOc_mld_smooth(ind+5) = MEAN(xint_TOc_mld(ind:ind+9))
   xint_TOc_zdf_smooth(ind+5) = MEAN(xint_TOc_zdf(ind:ind+9))
   xint_TOc_ldf_smooth(ind+5) = MEAN(xint_TOc_ldf(ind:ind+9))
   xint_TOc_gm_smooth (ind+5) = MEAN(xint_TOc_gm (ind:ind+9))

   xint_totP_smooth(ind+5) = MEAN(xint_totP(ind:ind+9))
   xint_bioP_smooth(ind+5) = MEAN(xint_bioP(ind:ind+9))
   xint_TOp_dyn_smooth(ind+5) = MEAN(xint_TOp_dyn(ind:ind+9))
   xint_TOp_edd_smooth(ind+5) = MEAN(xint_TOp_edd(ind:ind+9))
   xint_TOp_zmx_smooth(ind+5) = MEAN(xint_TOp_zmx(ind:ind+9))
   xint_TOp_adv_smooth(ind+5) = MEAN(xint_TOp_adv(ind:ind+9))
   xint_TOp_had_smooth(ind+5) = MEAN(xint_TOp_had(ind:ind+9))
   xint_TOp_zad_smooth(ind+5) = MEAN(xint_TOp_zad(ind:ind+9))
   xint_TOp_mld_smooth(ind+5) = MEAN(xint_TOp_mld(ind:ind+9))
   xint_TOp_zdf_smooth(ind+5) = MEAN(xint_TOp_zdf(ind:ind+9))
   xint_TOp_ldf_smooth(ind+5) = MEAN(xint_TOp_ldf(ind:ind+9))
   xint_TOp_gm_smooth (ind+5) = MEAN(xint_TOp_gm (ind:ind+9))

   xint_aTOsatc_dyn_smooth(ind+5) = MEAN(xint_aTOsatc_dyn(ind:ind+9))
   xint_aTOsatc_edd_smooth(ind+5) = MEAN(xint_aTOsatc_edd(ind:ind+9))
   xint_aTOsatc_zmx_smooth(ind+5) = MEAN(xint_aTOsatc_zmx(ind:ind+9))
   xint_aTOsatc_adv_smooth(ind+5) = MEAN(xint_aTOsatc_adv(ind:ind+9))
   xint_aTOsatc_had_smooth(ind+5) = MEAN(xint_aTOsatc_had(ind:ind+9))
   xint_aTOsatc_zad_smooth(ind+5) = MEAN(xint_aTOsatc_zad(ind:ind+9))
   xint_aTOsatc_mld_smooth(ind+5) = MEAN(xint_aTOsatc_mld(ind:ind+9))
   xint_aTOsatc_zdf_smooth(ind+5) = MEAN(xint_aTOsatc_zdf(ind:ind+9))
   xint_aTOsatc_ldf_smooth(ind+5) = MEAN(xint_aTOsatc_ldf(ind:ind+9))
   xint_aTOsatc_gm_smooth (ind+5) = MEAN(xint_aTOsatc_gm (ind:ind+9))

   xint_TAOUc_dyn_smooth(ind+5) = MEAN(xint_TAOUc_dyn(ind:ind+9))
   xint_TAOUc_edd_smooth(ind+5) = MEAN(xint_TAOUc_edd(ind:ind+9))
   xint_TAOUc_zmx_smooth(ind+5) = MEAN(xint_TAOUc_zmx(ind:ind+9))
   xint_TAOUc_adv_smooth(ind+5) = MEAN(xint_TAOUc_adv(ind:ind+9))
   xint_TAOUc_had_smooth(ind+5) = MEAN(xint_TAOUc_had(ind:ind+9))
   xint_TAOUc_zad_smooth(ind+5) = MEAN(xint_TAOUc_zad(ind:ind+9))
   xint_TAOUc_mld_smooth(ind+5) = MEAN(xint_TAOUc_mld(ind:ind+9))
   xint_TAOUc_zdf_smooth(ind+5) = MEAN(xint_TAOUc_zdf(ind:ind+9))
   xint_TAOUc_ldf_smooth(ind+5) = MEAN(xint_TAOUc_ldf(ind:ind+9))
   xint_TAOUc_gm_smooth (ind+5) = MEAN(xint_TAOUc_gm (ind:ind+9))

   Dxint_tot_smooth(ind+5) = MEAN(Dxint_tot(ind:ind+9))
   Dxint_bio_smooth(ind+5) = MEAN(Dxint_bio(ind:ind+9))
   Dxint_TO_dyn_smooth(ind+5) = MEAN(Dxint_TO_dyn(ind:ind+9))
   Dxint_TO_edd_smooth(ind+5) = MEAN(Dxint_TO_edd(ind:ind+9))
   Dxint_TO_zmx_smooth(ind+5) = MEAN(Dxint_TO_zmx(ind:ind+9))
   Dxint_TO_adv_smooth(ind+5) = MEAN(Dxint_TO_adv(ind:ind+9))
   Dxint_TO_had_smooth(ind+5) = MEAN(Dxint_TO_had(ind:ind+9))
   Dxint_TO_zad_smooth(ind+5) = MEAN(Dxint_TO_zad(ind:ind+9))
   Dxint_TO_mld_smooth(ind+5) = MEAN(Dxint_TO_mld(ind:ind+9))
   Dxint_TO_zdf_smooth(ind+5) = MEAN(Dxint_TO_zdf(ind:ind+9))
   Dxint_TO_ldf_smooth(ind+5) = MEAN(Dxint_TO_ldf(ind:ind+9))
   Dxint_TO_gm_smooth (ind+5) = MEAN(Dxint_TO_gm (ind:ind+9))

   Dxint_aTOsat_dyn_smooth(ind+5) = MEAN(Dxint_aTOsat_dyn(ind:ind+9))
   Dxint_aTOsat_edd_smooth(ind+5) = MEAN(Dxint_aTOsat_edd(ind:ind+9))
   Dxint_aTOsat_zmx_smooth(ind+5) = MEAN(Dxint_aTOsat_zmx(ind:ind+9))
   Dxint_aTOsat_adv_smooth(ind+5) = MEAN(Dxint_aTOsat_adv(ind:ind+9))
   Dxint_aTOsat_had_smooth(ind+5) = MEAN(Dxint_aTOsat_had(ind:ind+9))
   Dxint_aTOsat_zad_smooth(ind+5) = MEAN(Dxint_aTOsat_zad(ind:ind+9))
   Dxint_aTOsat_mld_smooth(ind+5) = MEAN(Dxint_aTOsat_mld(ind:ind+9))
   Dxint_aTOsat_zdf_smooth(ind+5) = MEAN(Dxint_aTOsat_zdf(ind:ind+9))
   Dxint_aTOsat_ldf_smooth(ind+5) = MEAN(Dxint_aTOsat_ldf(ind:ind+9))
   Dxint_aTOsat_gm_smooth (ind+5) = MEAN(Dxint_aTOsat_gm (ind:ind+9))

   Dxint_TAOU_dyn_smooth(ind+5) = MEAN(Dxint_TAOU_dyn(ind:ind+9))
   Dxint_TAOU_edd_smooth(ind+5) = MEAN(Dxint_TAOU_edd(ind:ind+9))
   Dxint_TAOU_zmx_smooth(ind+5) = MEAN(Dxint_TAOU_zmx(ind:ind+9))
   Dxint_TAOU_adv_smooth(ind+5) = MEAN(Dxint_TAOU_adv(ind:ind+9))
   Dxint_TAOU_had_smooth(ind+5) = MEAN(Dxint_TAOU_had(ind:ind+9))
   Dxint_TAOU_zad_smooth(ind+5) = MEAN(Dxint_TAOU_zad(ind:ind+9))
   Dxint_TAOU_mld_smooth(ind+5) = MEAN(Dxint_TAOU_mld(ind:ind+9))
   Dxint_TAOU_zdf_smooth(ind+5) = MEAN(Dxint_TAOU_zdf(ind:ind+9))
   Dxint_TAOU_ldf_smooth(ind+5) = MEAN(Dxint_TAOU_ldf(ind:ind+9))
   Dxint_TAOU_gm_smooth (ind+5) = MEAN(Dxint_TAOU_gm (ind:ind+9))
ENDFOR

;==================================================
; PLOT
;==================================================

x1tmp     = moyenne(gphit, 'x')
x1 = x1tmp[0:147]
x1(0)  =  -80

props1_ctl    = {BUFFER:1, OVERPLOT:1, LINESTYLE:'-', THICK:1.7, FONT_SIZE:8, $
                XRANGE:[min(x1), max(x1)], XTITLE:'Latitude', XTHICK:0.5,  $
                YRANGE: [-70, 30], YTICKINTERVAL:10, $
                YTITLE:'Deoxygenation     [GmolO2/m]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}
props1_delta  = {BUFFER:1, OVERPLOT:1, LINESTYLE:'-', THICK:1.7, FONT_SIZE:8, $
                XRANGE:[min(x1), max(x1)], XTITLE:'Latitude', XTHICK:0.5,  $
                YRANGE: [-10, 10], YTICKINTERVAL:2, $
                YTITLE:'Deoxygenation     [GmolO2/m]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}

props2_ctl    = {BUFFER:1, OVERPLOT:1, LINESTYLE:'-', THICK:2.5, FONT_SIZE:8, $
                XRANGE:[min(x1), max(x1)], XTITLE:'Latitude', XTHICK:0.5,  $
                YRANGE: [-10, 10], YTICKINTERVAL:2, $
                YTITLE:'Deoxygenation     [GmolO2/m]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}
props2_delta  = {BUFFER:1, OVERPLOT:1, LINESTYLE:'-', THICK:2.5, FONT_SIZE:8, $
                XRANGE:[min(x1), max(x1)], XTITLE:'Latitude', XTHICK:0.5,  $
                YRANGE: [-2, 2], YTICKINTERVAL:.4, $
                YTITLE:'Deoxygenation     [GmolO2/m]     Oxygenation', YTHICK:0.5, YTICKLEN:1, YGRIDSTYLE:2, YMINOR:0}

stop

;; ____________________
;; ARTICLE


w      = WINDOW(/BUFFER, DIMENSIONS = [300, 600]) ; dimension = [width, height]


tmp1 = xint_TOc_adv_smooth + xint_TOc_gm_smooth + xint_TOc_mld_smooth
tmp2 = xint_TOc_ldf_smooth + xint_TOc_zdf_smooth
p1     = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 1], TITLE = 'CTL simulation')
p1_1   = plot(x1, xint_totC_smooth(0:147), NAME = 'Total'       , COLOR = 'black'   , _EXTRA = props1_ctl)
p1_2   = plot(x1, xint_bioC_smooth(0:147), NAME = 'Respiration' , COLOR = 'green'   , _EXTRA = props1_ctl)
p1_3   = plot(x1, tmp1            (0:147), NAME = 'Kinematic'   , COLOR = 'sky blue', _EXTRA = props1_ctl)  
p1_4   = plot(x1, tmp2            (0:147), NAME = 'Diffusive'   , COLOR = 'tomato'  , _EXTRA = props1_ctl)  

tmp1 = xint_TOp_adv_smooth + xint_TOp_gm_smooth + xint_TOp_mld_smooth
tmp2 = xint_TOp_ldf_smooth + xint_TOp_zdf_smooth
p2     = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 2], TITLE = 'CC simulation')
p2_1   = plot(x1, xint_totP_smooth(0:147), NAME = 'Total'       , COLOR = 'black'   , _EXTRA = props1_ctl)
p2_2   = plot(x1, xint_bioP_smooth(0:147), NAME = 'Respiration' , COLOR = 'green'   , _EXTRA = props1_ctl)
p2_3   = plot(x1, tmp1            (0:147), NAME = 'Kinematic'   , COLOR = 'sky blue', _EXTRA = props1_ctl)  
p2_4   = plot(x1, tmp2            (0:147), NAME = 'Diffusive'   , COLOR = 'tomato'  , _EXTRA = props1_ctl)  
  
tmp1 = Dxint_TO_adv_smooth + Dxint_TO_gm_smooth + Dxint_TO_mld_smooth
tmp2 = Dxint_TO_ldf_smooth + Dxint_TO_zdf_smooth

p3    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [1, 3, 3], TITLE = 'Change')
p3_1  = plot(x1, Dxint_tot_smooth(0:147)     , NAME = 'Total'         , COLOR = 'black'         , _EXTRA = props1_delta)
p3_3  = plot(x1, Dxint_bio_smooth(0:147)     , NAME = 'Respiration'   , COLOR = 'green'         , _EXTRA = props1_delta)
p3_4  = plot(x1, tmp1(0:147)                 , NAME = 'Kinematic'     , COLOR = 'sky blue'      , _EXTRA = props1_delta)
p3_5  = plot(x1, tmp2(0:147)                 , NAME = 'Diffusive'     , COLOR = 'tomato'        , _EXTRA = props1_delta)
p3_2  = plot(x1, Z_doxy          (0:147)     , NAME = 'O2 ss MLD'     , COLOR = 'slate gray'    , LINESTYLE = '--', THICK = 2.5, /BUFFER, /OVERPLOT)

props_leg = {AUTO_TEXT_COLOR:1, FONT_SIZE:8, SAMPLE_WIDTH:.01}
leg = LEGEND(TARGET=[p3_1, p3_2, p3_3, p3_4, p3_5]  , POSITION = [1, 0.5], _EXTRA = props_leg)

filename = '/data/daclod/FIG/' + namefigbase + '.eps'
p1.save, filename

;; ____________________
;; TOUT


w      = WINDOW(/BUFFER, DIMENSIONS = [800, 1000]) ; dimension = [width, height]


p1     = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 1], TITLE = 'Control simulation')
p1_1   = plot(x1, xint_TOc_dyn_smooth(0:147)     , NAME = 'Ventilation'   , COLOR = 'royal blue'    , _EXTRA = props1_ctl)
p1_2   = plot(x1,    xint_bioC_smooth(0:147)     , NAME = 'Respiration'   , COLOR = 'green'         , _EXTRA = props1_ctl)  
p1_3   = plot(x1,    xint_totC_smooth(0:147)     , NAME = 'Ventil + respi', COLOR = 'black'         , _EXTRA = props1_ctl)  
  
p2    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 2], TITLE = 'Change')
p2_1  = plot(x1, Dxint_TO_dyn_smooth(0:147)     , NAME = 'Ventilation'   , COLOR = 'royal blue'    , _EXTRA = props1_delta)
p2_2  = plot(x1,    Dxint_bio_smooth(0:147)     , NAME = 'Respiration'   , COLOR = 'green'         , _EXTRA = props1_delta)
p2_3  = plot(x1,    Dxint_tot_smooth(0:147)     , NAME = 'Ventil + Respi', COLOR = 'black'         , _EXTRA = props1_delta)
p2_4  = plot(x1,       Z_doxy(0:147)            , NAME = 'O2 ss MLD'  , COLOR = 'slate gray'       , _EXTRA = props1_delta)

p3    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 3])
p3_1  = plot(x1, xint_TOc_had_smooth(0:147)        , NAME = 'HAD'           , COLOR = 'sky blue'   , _EXTRA = props1_ctl)
p3_2  = plot(x1, xint_TOc_zad_smooth(0:147)        , NAME = 'ZAD'           , COLOR = 'dark violet', _EXTRA = props1_ctl)

p4    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 4])
p4_1  = plot(x1, Dxint_TO_had_smooth(0:147)        , NAME = 'HAD'           , COLOR = 'sky blue'   , _EXTRA = props1_delta)  
p4_2  = plot(x1, Dxint_TO_zad_smooth(0:147)        , NAME = 'ZAD'           , COLOR = 'dark violet', _EXTRA = props1_delta)  

p5    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 5])
p5_1   = plot(x1, xint_TOc_zdf_smooth(0:147)        , NAME = 'ZDF'           , COLOR = 'tomato'     ,_EXTRA = props2_ctl)  
p5_2   = plot(x1, xint_TOc_mld_smooth(0:147)        , NAME = 'MLD'           , COLOR = 'dark red'   ,_EXTRA = props2_ctl)  

p6    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 6])
p6_1  = plot(x1, Dxint_TO_zdf_smooth(0:147)        , NAME = 'ZDF'           , COLOR = 'tomato'     , _EXTRA = props2_delta)   
p6_2  = plot(x1, Dxint_TO_mld_smooth(0:147)        , NAME = 'MLD'           , COLOR = 'dark red'   , _EXTRA = props2_delta)  
                                                                                                                                
p7    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 7])
p7_1   = plot(x1, xint_TOc_ldf_smooth(0:147)        , NAME = 'LDF'           , COLOR = 'gold'      , _EXTRA = props2_ctl)  
p7_2   = plot(x1, xint_TOc_gm_smooth (0:147)        , NAME = 'GM'            , COLOR = 'deep pink' , _EXTRA = props2_ctl)  

p8    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 8])
p8_1  = plot(x1, Dxint_TO_ldf_smooth(0:147)        , NAME = 'LDF'           , COLOR = 'gold'       , _EXTRA = props2_delta)   
p8_2  = plot(x1, Dxint_TO_gm_smooth (0:147)        , NAME = 'GM'            , COLOR = 'deep pink'  , _EXTRA = props2_delta)  

tmp1 = xint_TOc_adv_smooth + xint_TOc_gm_smooth
tmp2 = tmp1 + xint_TOc_mld_smooth
tmp3 = xint_TOc_ldf_smooth + xint_TOc_zdf_smooth
p9    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 9])
p9_1  = plot(x1, xint_TOc_adv_smooth(0:147), NAME = 'ADV'           , COLOR = 'blue'        , _EXTRA = props1_ctl)   
p9_2  = plot(x1, tmp1(0:147)               , NAME = 'ADV + GM'      , COLOR = 'tan'         , _EXTRA = props1_ctl)  
p9_3  = plot(x1, tmp2(0:147)               , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse'  , _EXTRA = props1_ctl)  
p9_4  = plot(x1, tmp3(0:147)               , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown' , _EXTRA = props1_ctl)  

tmp1 = Dxint_TO_adv_smooth + Dxint_TO_gm_smooth
tmp2 = tmp1 + Dxint_TO_mld_smooth
tmp3 = Dxint_TO_ldf_smooth + Dxint_TO_zdf_smooth
p10    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 10])
p10_1  = plot(x1, Dxint_TO_adv_smooth(0:147) , NAME = 'ADV'           , COLOR = 'blue'       , _EXTRA = props1_delta)   
p10_2  = plot(x1, tmp1(0:147)                , NAME = 'ADV + GM'      , COLOR = 'tan'        , _EXTRA = props1_delta)  
p10_3  = plot(x1, tmp2(0:147)                , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse' , _EXTRA = props1_delta)  
p10_4  = plot(x1, tmp3(0:147)                , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown', _EXTRA = props1_delta)  

props_leg = {AUTO_TEXT_COLOR:1, FONT_SIZE:6, SAMPLE_WIDTH:.01}
leg = LEGEND(TARGET=[p2_1, p2_2, p2_3]  , POSITION = [0.7, 0.83], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p4_1, p4_2]        , POSITION = [0.6, 0.67], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p6_1, p6_2]        , POSITION = [0.6, 0.5], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p8_1, p8_2]        , POSITION = [0.6, 0.33], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p10_1, p10_2, p10_3, p10_4], POSITION = [0.7, 0.17], _EXTRA = props_leg)

filename = '/data/daclod/FIG/' + namefigbase + '_tout.eps'
p1.save, filename

stop
;; ____________________
;; DaTOsat, DTAOU

w      = WINDOW(/BUFFER, DIMENSIONS = [800, 1000]) ; dimension = [width, height]


p1     = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 1], TITLE = 'Osat')
p1_1   = plot(x1, Dxint_aTOsat_dyn_smooth(0:147)     , NAME = 'Ventilation'   , COLOR = 'royal blue'    , _EXTRA = props1_delta)
  
p2    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 2], TITLE = 'Residu')
p2_1  = plot(x1, Dxint_TAOU_dyn_smooth(0:147)        , NAME = 'Ventilation'   , COLOR = 'royal blue', _EXTRA = props1_delta)

p3    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 3])
p3_1  = plot(x1, Dxint_aTOsat_had_smooth(0:147)        , NAME = 'HAD'           , COLOR = 'sky blue'   , _EXTRA = props1_delta)
p3_2  = plot(x1, Dxint_aTOsat_zad_smooth(0:147)        , NAME = 'ZAD'           , COLOR = 'dark violet', _EXTRA = props1_delta)

p4    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 4])
p4_1  = plot(x1, Dxint_TAOU_had_smooth(0:147)        , NAME = 'HAD'           , COLOR = 'sky blue'   , _EXTRA = props1_delta)  
p4_2  = plot(x1, Dxint_TAOU_zad_smooth(0:147)        , NAME = 'ZAD'           , COLOR = 'dark violet', _EXTRA = props1_delta)  

p5    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 5])
p5_1   = plot(x1, Dxint_aTOsat_zdf_smooth(0:147)        , NAME = 'ZDF'           , COLOR = 'tomato'     ,_EXTRA = props2_delta)  
p5_2   = plot(x1, Dxint_aTOsat_mld_smooth(0:147)        , NAME = 'MLD'           , COLOR = 'dark red'   ,_EXTRA = props2_delta)  

p6    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 6])
p6_1  = plot(x1, Dxint_TAOU_zdf_smooth(0:147)        , NAME = 'ZDF'           , COLOR = 'tomato'     , _EXTRA = props2_delta)   
p6_2  = plot(x1, Dxint_TAOU_mld_smooth(0:147)        , NAME = 'MLD'           , COLOR = 'dark red'   , _EXTRA = props2_delta)  
                                                                                                                                
p7    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 7])
p7_1   = plot(x1, Dxint_aTOsat_ldf_smooth(0:147)        , NAME = 'LDF'           , COLOR = 'gold'      , _EXTRA = props2_delta)  
p7_2   = plot(x1, Dxint_aTOsat_gm_smooth (0:147)        , NAME = 'GM'            , COLOR = 'deep pink' , _EXTRA = props2_delta)  

p8    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 8])
p8_1  = plot(x1, Dxint_TAOU_ldf_smooth(0:147)        , NAME = 'LDF'           , COLOR = 'gold'       , _EXTRA = props2_delta)   
p8_2  = plot(x1, Dxint_TAOU_gm_smooth (0:147)        , NAME = 'GM'            , COLOR = 'deep pink'  , _EXTRA = props2_delta)  

tmp1 = Dxint_aTOsat_adv_smooth + Dxint_aTOsat_gm_smooth
tmp2 = tmp1 + Dxint_aTOsat_mld_smooth
tmp3 = Dxint_aTOsat_ldf_smooth + Dxint_aTOsat_zdf_smooth
p9    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 9])
p9_1  = plot(x1, Dxint_aTOsat_adv_smooth(0:147), NAME = 'ADV'           , COLOR = 'blue'        , _EXTRA = props1_delta)   
p9_2  = plot(x1, tmp1(0:147)                   , NAME = 'ADV + GM'      , COLOR = 'tan'         , _EXTRA = props1_delta)  
p9_3  = plot(x1, tmp2(0:147)                   , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse'  , _EXTRA = props1_delta)  
p9_4  = plot(x1, tmp3(0:147)                   , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown' , _EXTRA = props1_delta)  

tmp1 = Dxint_TAOU_adv_smooth + Dxint_TAOU_gm_smooth
tmp2 = tmp1 + Dxint_TAOU_mld_smooth
tmp3 = Dxint_TAOU_ldf_smooth + Dxint_TAOU_zdf_smooth
p10    = plot(x1, 0 * x1, /CURRENT, /BUFFER, THICK = 0.5, LINESTYLE = 0, LAYOUT = [2, 5, 10])
p10_1  = plot(x1, Dxint_TAOU_adv_smooth(0:147), NAME = 'ADV'           , COLOR = 'blue'       , _EXTRA = props1_delta)   
p10_2  = plot(x1, tmp1(0:147)                 , NAME = 'ADV + GM'      , COLOR = 'tan'        , _EXTRA = props1_delta)  
p10_3  = plot(x1, tmp2(0:147)                 , NAME = 'ADV + GM + MLD', COLOR = 'chartreuse' , _EXTRA = props1_delta)  
p10_4  = plot(x1, tmp3(0:147)                 , NAME = 'LDF + ZDF'     , COLOR = 'sandy brown', _EXTRA = props1_delta)  

props_leg = {AUTO_TEXT_COLOR:1, FONT_SIZE:6, SAMPLE_WIDTH:.01}
leg = LEGEND(TARGET=[p2_1]              , POSITION = [0.7, 0.83], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p4_1, p4_2]        , POSITION = [0.6, 0.67], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p6_1, p6_2]        , POSITION = [0.6, 0.5], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p8_1, p8_2]        , POSITION = [0.6, 0.33], _EXTRA = props_leg)
leg = LEGEND(TARGET=[p10_1, p10_2, p10_3, p10_4], POSITION = [0.7, 0.17], _EXTRA = props_leg)

filename = '/data/daclod/FIG/' + namefigbase + '_aTOsat_TAOU.eps'
p1.save, filename


stop

END
