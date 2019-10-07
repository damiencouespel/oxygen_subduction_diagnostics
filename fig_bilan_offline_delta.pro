;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;;;   FIG_BILAN_OFFLINE_DELTA    ;;;;;;;;;;   
;;;;;;;;;;                              ;;;;;;;;;;             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Purpose 
; -------
; Plot maps of cumulative changes of subduction and respiration

@init_run_clim_dc

;==================================================
; PARAM
;==================================================

namefigbase = 'fig_bilan_offline_delta'; base of the file name to save the figures 
                                                                            
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

;==================================================
; PROCESS DATA
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

; suppression des valeurs aberrantes
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

; suppression des valeurs aberrantes
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
; Change = climate change - control
;==================================================

dtot = totP - totC          ; total = O2 subduction + O2 consumption
dbio = bioP - bioC          ; total export (O2 consumption)
dsubbio = subbioP - subbioC ; organic matter subduction (O2 consumption)
dexp = expP - expC          ; organic matter sinking (O2 consumption)
dTO_dyn = TOp_dyn - TOc_dyn ; total subduction
dTO_adv = TOp_adv - TOc_adv ; advection
dTO_zmx = TOp_zmx - TOc_zmx ; vertical mixing
dTO_edd = TOp_edd - TOc_edd ; eddies
dTO_zdf = TOp_zdf - TOc_zdf ; vertical diffusion
dTO_mld = TOp_mld - TOc_mld ; entrainment
dTO_had = TOp_had - TOc_had ; lateral advection
dTO_zad = TOp_zad - TOc_zad ; vertical advection
dTO_ldf = TOp_ldf - TOc_ldf ; isopycnal diffusion
dTO_gm  = TOp_gm  - TOc_gm  ; advection by bolus velocities (Gent McWilliams)

;==================================================
; PLOT
;==================================================

@init_orca2_shift

unit2 = 'Deoxygenation                      molO2/m2                      Oxygenation'
props2 = {CB_TITLE:unit2, LCT:42, NOCONTOUR:1, REALCONT:1, SUBTITLE:'', FORMAT:'(e11.2)', INV:1}

stop 

saveplot = namefigbase + ''
tmp1 = dTO_adv + dTO_gm + dTO_mld
tmp2 = dTO_ldf + dTO_zdf
print, 'plot saved: ', saveplot
openps, FILE = saveplot, /LANDSCAPE
plt, shift(tmp1(1:180, 0:147), [30, 0]), -500, 500, small=[2, 2, 1], title = 'Kinematic flux change' , _extra = props2, win = 1, /LANDSCAPE, /inv
plt, shift(tmp2(1:180, 0:147), [30, 0]), -500, 500, small=[2, 2, 2], title = 'Diffusive flux change' , _extra = props2, /noerase, /inv
plt, shift(dbio(1:180, 0:147), [30, 0]), -500, 500, small=[2, 2, 3], title = 'Respiration change'    , _extra = props2, /noerase, /inv
closeps

saveplot = namefigbase + '_smooth'
tmp1 = dTO_adv + dTO_gm + dTO_mld
tmp2 = dTO_ldf + dTO_zdf
print, 'plot saved: ', saveplot
openps, FILE = saveplot, /LANDSCAPE
plt, shift(GAUSS_SMOOTH(tmp1(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 2, 1], title = 'Kinematic flux change' , _extra = props2, win = 1, /LANDSCAPE, /inv
plt, shift(GAUSS_SMOOTH(tmp2(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 2, 2], title = 'Diffusive flux change' , _extra = props2, /noerase, /inv
plt, shift(GAUSS_SMOOTH(dbio(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 2, 3], title = 'Respiration change'    , _extra = props2, /noerase, /inv
closeps

saveplot = namefigbase + '_dyn_bio'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(dTO_dyn(1:180, 0:147) , [30, 0]), -200, 200, small=[1, 2, 1], title = 'Ventil O2 change'   , _extra = props2, win = 1, /portrait
plt, shift(dbio(1:180, 0:147)    , [30, 0]), -200, 200, small=[1, 2, 2], title = 'Respiration change' , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_dynO2_exp_subbio'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape
plt, shift(dTO_dyn(1:180, 0:147), [30, 0]), -2000, 2000, small=[2, 2, 1], title = 'Ventil O2 change'  , _extra = props2, win = 1, /landscape
plt, shift(dsubbio(1:180, 0:147), [30, 0]), -2000, 2000, small=[2, 2, 2], title = 'Ventil OM change'  , _extra = props2, /noerase
plt, shift(dexp(1:180, 0:147)   , [30, 0]), -2000, 2000, small=[2, 2, 3], title = 'Export change'     , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_exp_subbio'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(dsubbio(1:180, 0:147), [30, 0]), -500, 500, small=[1, 2, 1], title = 'Ventil OM change'  , _extra = props2, win = 1, /portrait
plt, shift(dexp(1:180, 0:147)   , [30, 0]), -500, 500, small=[1, 2, 2], title = 'Export change'     , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_exp_subbio_smooth'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(GAUSS_SMOOTH(dsubbio(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[1, 2, 1], title = 'Ventil OM change'  , _extra = props2, win = 1, /portrait
plt, shift(GAUSS_SMOOTH(dexp(1:180, 0:147)   , /EDGE_WRAP), [30, 0]), -500, 500, small=[1, 2, 2], title = 'Export change'     , _extra = props2, /noerase
closeps

stop

;; DYN O2

saveplot = namefigbase + '_adv_zmx_edd_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape
plt, shift(dTO_dyn(1:180, 0:147) , [30, 0]), -2000, 2000, small=[2, 2, 1], title = 'Ventil O2 change'  , _extra = props2, win = 1, /landscape
plt, shift(dTO_adv(1:180, 0:147) , [30, 0]), -2000, 2000, small=[2, 2, 2], title = 'ADV O2 change'     , _extra = props2, /noerase
plt, shift(dTO_zmx(1:180, 0:147) , [30, 0]), -2000, 2000, small=[2, 2, 3], title = 'ZMIX O2 change'    , _extra = props2, /noerase
plt, shift(dTO_edd(1:180, 0:147) , [30, 0]), -2000, 2000, small=[2, 2, 4], title = 'EDD MIX O2 change' , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_zmx_zdf_mld_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape
plt, shift(dTO_zmx(1:180, 0:147), [30, 0]), -200, 200, small=[2, 2, 1], title = 'ZMIX O2 change'    , _extra = props2, win = 1, /landscape
plt, shift(dTO_zdf(1:180, 0:147), [30, 0]), -200, 200, small=[2, 2, 2], title = 'ZDF O2 change'    , _extra = props2, /noerase
plt, shift(dTO_mld(1:180, 0:147), [30, 0]), -200, 200, small=[2, 2, 3], title = 'MLD O2 change'    , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_zdf_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape
plt, shift(dTO_zdf(1:180, 0:147), [30, 0]), -500, 500, title = 'ZDF O2 change'    , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_adv_had_zad_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape
plt, shift(dTO_adv(1:180, 0:147), [30, 0]), -2000, 2000, small=[2, 2, 1], title = 'ADV O2 change'  , _extra = props2, win = 1, /landscape
plt, shift(dTO_had(1:180, 0:147), [30, 0]), -2000, 2000, small=[2, 2, 2], title = 'HAD O2 change'  , _extra = props2, /noerase
plt, shift(dTO_zad(1:180, 0:147), [30, 0]), -2000, 2000, small=[2, 2, 3], title = 'ZAD O2 change'  , _extra = props2, /noerase
closeps

saveplot = namefigbase + '_dyn_adv_had_zad_zdf_mld_edd_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(dTO_dyn(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 1], title = 'DYN O2'   , _extra = props2, win = 1, /portrait, /inv
plt, shift(dTO_adv(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 2], title = 'ADV O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_had(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 3], title = 'HAD O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_zad(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 4], title = 'ZAD O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_zdf(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 5], title = 'ZDF O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_mld(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 6], title = 'MLD O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_edd(1:180, 0:147), [30, 0]), -500, 500, small=[2, 4, 7], title = 'EDD O2'   , _extra = props2, /noerase, /inv
closeps

saveplot = namefigbase + '_had_zad_zdf_mld_ldf_gm_O2'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(dTO_had(1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 1], title = 'HAD O2'   , _extra = props2, win = 1, /portrait, /inv
plt, shift(dTO_zad(1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 2], title = 'ZAD O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_gm (1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 3], title = 'GM O2'    , _extra = props2, /noerase, /inv
plt, shift(dTO_mld(1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 4], title = 'MLD O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_zdf(1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 5], title = 'ZDF O2'   , _extra = props2, /noerase, /inv
plt, shift(dTO_ldf(1:180, 0:147), [30, 0]), -500, 500, small=[2, 3, 6], title = 'LDF O2'   , _extra = props2, /noerase, /inv
closeps

saveplot = namefigbase + '_had_zad_zdf_mld_ldf_gm_O2_smooth'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait
plt, shift(GAUSS_SMOOTH(dTO_had(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 1], title = 'HAD O2'   , _extra = props2, win = 1, /portrait, /inv
plt, shift(GAUSS_SMOOTH(dTO_zad(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 2], title = 'ZAD O2'   , _extra = props2, /noerase, /inv
plt, shift(GAUSS_SMOOTH(dTO_gm (1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 3], title = 'GM O2'    , _extra = props2, /noerase, /inv
plt, shift(GAUSS_SMOOTH(dTO_mld(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 4], title = 'MLD O2'   , _extra = props2, /noerase, /inv
plt, shift(GAUSS_SMOOTH(dTO_zdf(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 5], title = 'ZDF O2'   , _extra = props2, /noerase, /inv
plt, shift(GAUSS_SMOOTH(dTO_ldf(1:180, 0:147), /EDGE_WRAP), [30, 0]), -500, 500, small=[2, 3, 6], title = 'LDF O2'   , _extra = props2, /noerase, /inv
closeps

stop
END
