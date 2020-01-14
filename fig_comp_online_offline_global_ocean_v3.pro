; FIG_COMP_ONLINE_OFFLINE_GLOBAL_OCEAN
;------------------------------------------------------------
; Trace le bilan regional de subdu pour les calculs online et offline
; - subduction totale
; - par mécanismes : adv horizontale, verticale, zmix, eddy mixing
; - de l'export

; 4 JUNE 2019: DC4
; fig_comp_online_offline_global_ocean_v2.pro, modif pour ajouter
; calcul respiration
; 14 MAR 2019: DC3
; modification pour affiner le calcul
; 24 JAN 2019: DC2
; from fig_comp_online_offline_region_v2
; WED 15 OCT 2014: DC
; - ajout de la comp online-offline pour l'export
; MON 13 OCT 2014: DC

@common
@init_run_clim_dc
surf = e1t*e2t

;==================================================
; PARAM
;==================================================

namefigbase = 'fig_comp_online_offline_global_ocean_v3'

fdirOFF  = '/data/daclod/MY_IDL/DATA_TMP/'
fsubOFF     = 'sub_offline_orca2_clim_O2_ORCA2_CLIM_1m_00500101_00511231.sav'
fsubOFF_zdf = 'sub_offline_orca2_clim_v2_O2_ORCA2_CLIM_1m_00500101_00511231_subzdf.sav'
fexpOFF    = 'export_offline_orca2_clim_ORCA2_CLIM_1m_00500101_00511231.sav'
fsubtocOFF = 'sub_offline_orca2_clim_TOC_ORCA2_CLIM_1m_00500101_00511231.sav'
fsubdocOFF = 'sub_offline_orca2_clim_DOC_ORCA2_CLIM_1m_00500101_00511231.sav'
fsubnh4OFF = 'sub_offline_orca2_clim_NH4_ORCA2_CLIM_1m_00500101_00511231.sav'
;; fOFF_subw = 'sub_offline_orca2_clim_v3_O2_ORCA2_CLIM_1m_00500101_00511231_subwadv.sav'
fdirON    = '/data/daclod/DATA_ORCA2_CLIM/'
fsubON    = 'ORCA2_CLIM_1m_00500101_00511231_oxysub.nc'
fexpON    = 'ORCA2_CLIM_1m_00500101_00511231_tocsub.nc'
fsubtocON = 'ORCA2_CLIM_1m_00500101_00511231_tocsub.nc'
fsubdocON = 'ORCA2_CLIM_1m_00500101_00511231_docsub.nc'
fsubnh4ON = 'ORCA2_CLIM_1m_00500101_00511231_nh4sub.nc'

;==================================================
; DATA
;==================================================

; ________________________________________
; SUB OFFLINE

subdu = restore_subdu(fdirOFF, fsubOFF, 1e13) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / maille / month

; selection de 12 mois
subh   = subh   (*, *, 12:23) 
subw   = subw   (*, *, 12:23)
submld = submld (*, *, 12:23) 
subzdf = subzdf (*, *, 12:23)
subgm  = subgm  (*, *, 12:23)
subldf = subldf (*, *, 12:23)

; TEMPORAL INTEGRAL
zfactsub = 1000. ; kmolO2 -> molO2
jpt = 12
time = findgen(jpt)
hadOFF = grossemoyenne(subh  , 't', /INT, /NAN) * zfactsub / surf
zadOFF = grossemoyenne(subw  , 't', /INT, /NAN) * zfactsub / surf
mldOFF = grossemoyenne(submld, 't', /INT, /NAN) * zfactsub / surf                   
zdfOFF = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactsub / surf                   
ldfOFF = grossemoyenne(subldf, 't', /INT, /NAN) * zfactsub / surf                   
gmwOFF  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactsub / surf                    
; unit : molO2 / m2 / ans

; to remove the overlap in the indian ocean and north pole
hadOFF = orca2_2d_overlap_remove(hadOFF)
zadOFF = orca2_2d_overlap_remove(zadOFF)
zdfOFF = orca2_2d_overlap_remove(zdfOFF)
mldOFF = orca2_2d_overlap_remove(mldOFF)
ldfOFF = orca2_2d_overlap_remove(ldfOFF)
gmwOFF = orca2_2d_overlap_remove(gmwOFF)
; unit : PmolO2/m2/year

; horizontal integral
hadOFF_tot = moyenne(hadOFF, 'xy', /NAN, /INT)
zadOFF_tot = moyenne(zadOFF, 'xy', /NAN, /INT)
mldOFF_tot = moyenne(mldOFF, 'xy', /NAN, /INT)
zdfOFF_tot = moyenne(zdfOFF, 'xy', /NAN, /INT)
gmwOFF_tot = moyenne(gmwOFF, 'xy', /NAN, /INT)
ldfOFF_tot = moyenne(ldfOFF, 'xy', /NAN, /INT)

; totaux
kinOFF = hadOFF + zadOFF + mldOFF + gmwOFF
difOFF = zdfOFF + ldfOFF
dynOFF = kinOFF + difOFF

kinOFF_tot = hadOFF_tot + zadOFF_tot + mldOFF_tot + gmwOFF_tot
difOFF_tot = zdfOFF_tot + ldfOFF_tot
dynOFF_tot = kinOFF_tot + difOFF_tot

; ________________________________________
; SUB ONLINE


restorefile = fdirON + fsubON

;%%% Masques U et V (sur les bords uniquement)
uumask = fltarr(182, 149, 24)+1.
vvmask = fltarr(182, 149, 24)+1.

uumask(0, *, *) = 0.
uumask(181, *, *) = 0.
uumask(*, 148, *) = 0.
uumask(92:181, 147, *) = 0.
uumask(1, 147, *) = 0.

vvmask(0, *, *) = 0.
vvmask(181, *, *) = 0.
vvmask(*, 148, *) = 0.

ttmask =  fltarr(182, 149, 24)
for tt=0, 23 DO ttmask[*, *, tt] = tmask[*, *, 0]

;Advection by u,v,w
xad = read_ncdf('xad_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
xad(where (xad GE 1.e19)) = 0.
yad = read_ncdf('yad_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
yad(where(yad GE 1.e19)) = 0.
xad = xad*uumask
yad = yad*vvmask
zad = read_ncdf('zad_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask

;advection totale= advection by u,v,w + adv GM
xeiv = read_ncdf('xei_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
xeiv(where (xeiv GE 1.e19)) = 0.
yeiv = read_ncdf('yei_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
yeiv(where (yeiv GE 1.e19)) = 0.
xeiv = xeiv*uumask
yeiv = yeiv*vvmask
zeiv = read_ncdf('zei_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask
heiv=xeiv+yeiv
eiv=xeiv+yeiv+zeiv

;Diffusion laterale
;warning: signe - sur xlf et ylf pour corriger bug dans le code
xlf = read_ncdf('xlf_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
ylf = read_ncdf('ylf_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
zlf = read_ncdf('zlf_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)
xlf(where (xlf GE 1.e19)) = 0.
ylf(where (ylf GE 1.e19)) = 0.
lf = -xlf*uumask-ylf*vvmask+zlf*ttmask
;; hlf =  -xlf*uumask-ylf*vvmask
;; zlf =  zlf*ttmask

;entrainement de la MLD
mld = read_ncdf('mld_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion verticale
zdf = read_ncdf('zdf_sub_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;surface boundary condition: effet des precipitations/evaporation
sbc = read_ncdf('sbc_mld_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion dans bottom boundary layer
bbl = read_ncdf('bbl_mld_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement de la surface libre (compte en negatif car dans la subdu)
rof = read_ncdf('rof_mld_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement correction des valeurs negatives 
rad = read_ncdf('rad_mld_O2', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;Facteur de conversion pour avoir des unites en molO2/an and flux
;toward the subsurface positive 
zfact = 365. * 86400. * 1000.
xad = zfact * xad;(*, *, 0:11)
yad = zfact * yad;(*, *, 0:11)
zad = zfact * zad;(*, *, 0:11)
eiv = zfact * eiv;(*, *, 0:11)
lf  = zfact * lf ;(*, *, 0:11)
;; hlf = zfact * hlf;(*, *, 0:11)
;; zlf = zfact * zlf;(*, *, 0:11)
mld = zfact * mld;(*, *, 0:11)
zdf = zfact * zdf;(*, *, 0:11)
rof = zfact * rof;(*, *, 0:11)
sbc = zfact * sbc;(*, *, 0:11)
rad = zfact * rad;(*, *, 0:11)
bbl = zfact * bbl;(*, *, 0:11)

;quelques bilans:
hadON = (xad + yad)
zadON = zad
mldON = mld
zdfON = zdf
ldfON = lf
;; zlfON = zlf
;; hlfON = hlf
gmwON = (eiv - xad - yad - zad)
otherON = rof + sbc + rad + bbl

;; ; temporal integral
jpt=12
time = findgen(jpt)
hadON = grossemoyenne(hadON, 't', /NAN)
zadON = grossemoyenne(zadON, 't', /NAN)
mldON = grossemoyenne(mldON, 't', /NAN)
zdfON = grossemoyenne(zdfON, 't', /NAN)
ldfON = grossemoyenne(ldfON, 't', /NAN)
;; zlfON = grossemoyenne(zlfON, 't', /NAN)
;; hlfON = grossemoyenne(hlfON, 't', /NAN)
gmwON = grossemoyenne(gmwON, 't', /NAN)
otherON = grossemoyenne(otherON, 't', /INT, /NAN)

hadON = orca2_2d_overlap_remove(hadON) / surf
zadON = orca2_2d_overlap_remove(zadON) / surf
mldON = orca2_2d_overlap_remove(mldON) / surf
zdfON = orca2_2d_overlap_remove(zdfON) / surf
ldfON = orca2_2d_overlap_remove(ldfON) / surf
;; zlfON = orca2_2d_overlap_remove(zlfON) / surf
;; hlfON = orca2_2d_overlap_remove(hlfON) / surf
gmwON = orca2_2d_overlap_remove(gmwON) / surf
otherON = orca2_2d_overlap_remove(otherON) / surf
; unit : PmolO2/m2/an

;; ; horizontal integral
hadON_tot = moyenne(hadON, 'xy', /NAN, /INT)
zadON_tot = moyenne(zadON, 'xy', /NAN, /INT)
mldON_tot = moyenne(mldON, 'xy', /NAN, /INT)
zdfON_tot = moyenne(zdfON, 'xy', /NAN, /INT)
ldfON_tot = moyenne(ldfON, 'xy', /NAN, /INT)
;; zlfON_tot = moyenne(zlfON, 'xy', /NAN, /INT)
;; hlfON_tot = moyenne(hlfON, 'xy', /NAN, /INT)
gmwON_tot = moyenne(gmwON, 'xy', /NAN, /INT)

; totaux
kinON = hadON + zadON + mldON + gmwON
difON = zdfON + ldfON
dynON = kinON + difON

kinON_tot = hadON_tot + zadON_tot + mldON_tot + gmwON_tot
difON_tot = zdfON_tot + ldfON_tot
dynON_tot = kinON_tot + difON_tot

; ________________________________________
; SUB TOC OFFLINE

zfactbio   = -160./122.  * 1000. ; kmolC -> molO2  

subdu = restore_subdu(fdirOFF, fsubtocOFF, 1e9) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / maille / month

; selection de 12 mois
subh   = subh   (*, *, 12:23) 
subw   = subw   (*, *, 12:23)
submld = submld (*, *, 12:23) 
subzdf = subzdf (*, *, 12:23)
subgm  = subgm  (*, *, 12:23)
subldf = subldf (*, *, 12:23)

; TEMPORAL INTEGRAL
jpt = 12
time = findgen(jpt)
hadtocOFF = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio / surf
zadtocOFF = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio / surf
mldtocOFF = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio / surf                   
zdftocOFF = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio / surf                   
ldftocOFF = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio / surf                   
gmwtocOFF  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio / surf                    
; unit : molO2 / m2 / ans

; to remove the overlap in the indian ocean and north pole
hadtocOFF = orca2_2d_overlap_remove(hadtocOFF)
zadtocOFF = orca2_2d_overlap_remove(zadtocOFF)
zdftocOFF = orca2_2d_overlap_remove(zdftocOFF)
mldtocOFF = orca2_2d_overlap_remove(mldtocOFF)
ldftocOFF = orca2_2d_overlap_remove(ldftocOFF)
gmwtocOFF = orca2_2d_overlap_remove(gmwtocOFF)
; unit : PmolO2/m2/year

; horizontal integral
hadtocOFF_tot = moyenne(hadtocOFF, 'xy', /NAN, /INT)
zadtocOFF_tot = moyenne(zadtocOFF, 'xy', /NAN, /INT)
mldtocOFF_tot = moyenne(mldtocOFF, 'xy', /NAN, /INT)
zdftocOFF_tot = moyenne(zdftocOFF, 'xy', /NAN, /INT)
gmwtocOFF_tot = moyenne(gmwtocOFF, 'xy', /NAN, /INT)
ldftocOFF_tot = moyenne(ldftocOFF, 'xy', /NAN, /INT)

; totaux
kintocOFF = hadtocOFF + zadtocOFF + mldtocOFF + gmwtocOFF
diftocOFF = zdftocOFF + ldftocOFF
dyntocOFF = kintocOFF + diftocOFF

kintocOFF_tot = hadtocOFF_tot + zadtocOFF_tot + mldtocOFF_tot + gmwtocOFF_tot
diftocOFF_tot = zdftocOFF_tot + ldftocOFF_tot
dyntocOFF_tot = kintocOFF_tot + diftocOFF_tot

; ________________________________________
; SUB TOC ONLINE

restorefile = fdirON + fsubtocON

;%%% Masques U et V (sur les bords uniquement)
uumask = fltarr(182, 149, 24)+1.
vvmask = fltarr(182, 149, 24)+1.

uumask(0, *, *) = 0.
uumask(181, *, *) = 0.
uumask(*, 148, *) = 0.
uumask(92:181, 147, *) = 0.
uumask(1, 147, *) = 0.

vvmask(0, *, *) = 0.
vvmask(181, *, *) = 0.
vvmask(*, 148, *) = 0.

ttmask =  fltarr(182, 149, 24)
for tt=0, 23 DO ttmask[*, *, tt] = tmask[*, *, 0]

;Advection by u,v,w
xad = read_ncdf('xad_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xad(where (xad GE 1.e19)) = 0.
yad = read_ncdf('yad_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
yad(where(yad GE 1.e19)) = 0.
xad = xad*uumask
yad = yad*vvmask
zad = read_ncdf('zad_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask

;advection totale= advection by u,v,w + adv GM
xeiv = read_ncdf('xei_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xeiv(where (xeiv GE 1.e19)) = 0.
yeiv = read_ncdf('yei_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
yeiv(where (yeiv GE 1.e19)) = 0.
xeiv = xeiv*uumask
yeiv = yeiv*vvmask
zeiv = read_ncdf('zei_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask
heiv=xeiv+yeiv
eiv=xeiv+yeiv+zeiv

;Diffusion laterale
;warning: signe - sur xlf et ylf pour corriger bug dans le code
xlf = read_ncdf('xlf_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
ylf = read_ncdf('ylf_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
zlf = read_ncdf('zlf_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xlf(where (xlf GE 1.e19)) = 0.
ylf(where (ylf GE 1.e19)) = 0.
lf = -xlf*uumask-ylf*vvmask+zlf*ttmask
;; hlf =  -xlf*uumask-ylf*vvmask
;; zlf =  zlf*ttmask

;entrainement de la MLD
mld = read_ncdf('mld_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion verticale
zdf = read_ncdf('zdf_sub_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;surface boundary condition: effet des precipitations/evaporation
sbc = read_ncdf('sbc_mld_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion dans bottom boundary layer
bbl = read_ncdf('bbl_mld_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement de la surface libre (compte en negatif car dans la subdu)
rof = read_ncdf('rof_mld_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement correction des valeurs negatives 
rad = read_ncdf('rad_mld_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;Facteur de conversion pour avoir des unites en molO2/an 
zfactbio   = -160./122. ; molC -> molO2  
zfact = 365. * 86400. * 1000. * zfactbio
xad = zfact * xad;(*, *, 0:11)
yad = zfact * yad;(*, *, 0:11)
zad = zfact * zad;(*, *, 0:11)
eiv = zfact * eiv;(*, *, 0:11)
lf  = zfact * lf ;(*, *, 0:11)
;; hlf = zfact * hlf;(*, *, 0:11)
;; zlf = zfact * zlf;(*, *, 0:11)
mld = zfact * mld;(*, *, 0:11)
zdf = zfact * zdf;(*, *, 0:11)
rof = zfact * rof;(*, *, 0:11)
sbc = zfact * sbc;(*, *, 0:11)
rad = zfact * rad;(*, *, 0:11)
bbl = zfact * bbl;(*, *, 0:11)

;quelques bilans:
hadtocON = (xad + yad)
zadtocON = zad
mldtocON = mld
zdftocON = zdf
ldftocON = lf
;; zlftocON = zlf
;; hlftocON = hlf
gmwtocON = (eiv - xad - yad - zad)
othertocON = rof + sbc + rad + bbl

;; ; temporal integral
jpt=12
time = findgen(jpt)
hadtocON = grossemoyenne(hadtocON, 't', /NAN)
zadtocON = grossemoyenne(zadtocON, 't', /NAN)
mldtocON = grossemoyenne(mldtocON, 't', /NAN)
zdftocON = grossemoyenne(zdftocON, 't', /NAN)
ldftocON = grossemoyenne(ldftocON, 't', /NAN)
;; zlftocON = grossemoyenne(zlftocON, 't', /NAN)
;; hlftocON = grossemoyenne(hlftocON, 't', /NAN)
gmwtocON = grossemoyenne(gmwtocON, 't', /NAN)
othertocON = grossemoyenne(othertocON, 't', /INT, /NAN)

hadtocON = orca2_2d_overlap_remove(hadtocON) / surf
zadtocON = orca2_2d_overlap_remove(zadtocON) / surf
mldtocON = orca2_2d_overlap_remove(mldtocON) / surf
zdftocON = orca2_2d_overlap_remove(zdftocON) / surf
ldftocON = orca2_2d_overlap_remove(ldftocON) / surf
;; zlftocON = orca2_2d_overlap_remove(zlftocON) / surf
;; hlftocON = orca2_2d_overlap_remove(hlftocON) / surf
gmwtocON = orca2_2d_overlap_remove(gmwtocON) / surf
othertocON = orca2_2d_overlap_remove(othertocON) / surf
; unit : PmolO2/m2/an

;; ; horizontal integral
hadtocON_tot = moyenne(hadtocON, 'xy', /NAN, /INT)
zadtocON_tot = moyenne(zadtocON, 'xy', /NAN, /INT)
mldtocON_tot = moyenne(mldtocON, 'xy', /NAN, /INT)
zdftocON_tot = moyenne(zdftocON, 'xy', /NAN, /INT)
ldftocON_tot = moyenne(ldftocON, 'xy', /NAN, /INT)
;; zlftocON_tot = moyenne(zlftocON, 'xy', /NAN, /INT)
;; hlftocON_tot = moyenne(hlftocON, 'xy', /NAN, /INT)
gmwtocON_tot = moyenne(gmwtocON, 'xy', /NAN, /INT)

; totaux
kintocON = hadtocON + zadtocON + mldtocON + gmwtocON
diftocON = zdftocON + ldftocON
dyntocON = kintocON + diftocON

kintocON_tot = hadtocON_tot + zadtocON_tot + mldtocON_tot + gmwtocON_tot
diftocON_tot = zdftocON_tot + ldftocON_tot
dyntocON_tot = kintocON_tot + diftocON_tot

; ________________________________________
; SUB DOC OFFLINE

zfactbio   = -160./122.  * 1000. ; kmolC -> molO2  

subdu = restore_subdu(fdirOFF, fsubdocOFF, 1e11) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / maille / month

; selection de 12 mois
subh   = subh   (*, *, 12:23) 
subw   = subw   (*, *, 12:23)
submld = submld (*, *, 12:23) 
subzdf = subzdf (*, *, 12:23)
subgm  = subgm  (*, *, 12:23)
subldf = subldf (*, *, 12:23)

; TEMPORAL INTEGRAL
jpt = 12
time = findgen(jpt)
haddocOFF = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio / surf
zaddocOFF = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio / surf
mlddocOFF = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio / surf                   
zdfdocOFF = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio / surf                   
ldfdocOFF = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio / surf                   
gmwdocOFF  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio / surf                    
; unit : molO2 / m2 / ans

; to remove the overlap in the indian ocean and north pole
haddocOFF = orca2_2d_overlap_remove(haddocOFF)
zaddocOFF = orca2_2d_overlap_remove(zaddocOFF)
zdfdocOFF = orca2_2d_overlap_remove(zdfdocOFF)
mlddocOFF = orca2_2d_overlap_remove(mlddocOFF)
ldfdocOFF = orca2_2d_overlap_remove(ldfdocOFF)
gmwdocOFF = orca2_2d_overlap_remove(gmwdocOFF)
; unit : PmolO2/m2/year

; horizontal integral
haddocOFF_tot = moyenne(haddocOFF, 'xy', /NAN, /INT)
zaddocOFF_tot = moyenne(zaddocOFF, 'xy', /NAN, /INT)
mlddocOFF_tot = moyenne(mlddocOFF, 'xy', /NAN, /INT)
zdfdocOFF_tot = moyenne(zdfdocOFF, 'xy', /NAN, /INT)
gmwdocOFF_tot = moyenne(gmwdocOFF, 'xy', /NAN, /INT)
ldfdocOFF_tot = moyenne(ldfdocOFF, 'xy', /NAN, /INT)

; totaux
kindocOFF = haddocOFF + zaddocOFF + mlddocOFF + gmwdocOFF
difdocOFF = zdfdocOFF + ldfdocOFF
dyndocOFF = kindocOFF + difdocOFF

kindocOFF_tot = haddocOFF_tot + zaddocOFF_tot + mlddocOFF_tot + gmwdocOFF_tot
difdocOFF_tot = zdfdocOFF_tot + ldfdocOFF_tot
dyndocOFF_tot = kindocOFF_tot + difdocOFF_tot

; ________________________________________
; SUB DOC ONLINE

restorefile = fdirON + fsubdocON

;%%% Masques U et V (sur les bords uniquement)
uumask = fltarr(182, 149, 24)+1.
vvmask = fltarr(182, 149, 24)+1.

uumask(0, *, *) = 0.
uumask(181, *, *) = 0.
uumask(*, 148, *) = 0.
uumask(92:181, 147, *) = 0.
uumask(1, 147, *) = 0.

vvmask(0, *, *) = 0.
vvmask(181, *, *) = 0.
vvmask(*, 148, *) = 0.

ttmask =  fltarr(182, 149, 24)
for tt=0, 23 DO ttmask[*, *, tt] = tmask[*, *, 0]

;Advection by u,v,w
xad = read_ncdf('xad_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xad(where (xad GE 1.e19)) = 0.
yad = read_ncdf('yad_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
yad(where(yad GE 1.e19)) = 0.
xad = xad*uumask
yad = yad*vvmask
zad = read_ncdf('zad_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask

;advection totale= advection by u,v,w + adv GM
xeiv = read_ncdf('xei_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xeiv(where (xeiv GE 1.e19)) = 0.
yeiv = read_ncdf('yei_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
yeiv(where (yeiv GE 1.e19)) = 0.
xeiv = xeiv*uumask
yeiv = yeiv*vvmask
zeiv = read_ncdf('zei_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask
heiv=xeiv+yeiv
eiv=xeiv+yeiv+zeiv

;Diffusion laterale
;warning: signe - sur xlf et ylf pour corriger bug dans le code
xlf = read_ncdf('xlf_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
ylf = read_ncdf('ylf_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
zlf = read_ncdf('zlf_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)
xlf(where (xlf GE 1.e19)) = 0.
ylf(where (ylf GE 1.e19)) = 0.
lf = -xlf*uumask-ylf*vvmask+zlf*ttmask
;; hlf =  -xlf*uumask-ylf*vvmask
;; zlf =  zlf*ttmask

;entrainement de la MLD
mld = read_ncdf('mld_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion verticale
zdf = read_ncdf('zdf_sub_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;surface boundary condition: effet des precipitations/evaporation
sbc = read_ncdf('sbc_mld_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion dans bottom boundary layer
bbl = read_ncdf('bbl_mld_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement de la surface libre (compte en negatif car dans la subdu)
rof = read_ncdf('rof_mld_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement correction des valeurs negatives 
rad = read_ncdf('rad_mld_DOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;Facteur de conversion pour avoir des unites en molO2/an 
zfactbio   = -160./122. ; molC -> molO2  
zfact = 365. * 86400. * 1000. * zfactbio
xad = zfact * xad;(*, *, 0:11)
yad = zfact * yad;(*, *, 0:11)
zad = zfact * zad;(*, *, 0:11)
eiv = zfact * eiv;(*, *, 0:11)
lf  = zfact * lf ;(*, *, 0:11)
;; hlf = zfact * hlf;(*, *, 0:11)
;; zlf = zfact * zlf;(*, *, 0:11)
mld = zfact * mld;(*, *, 0:11)
zdf = zfact * zdf;(*, *, 0:11)
rof = zfact * rof;(*, *, 0:11)
sbc = zfact * sbc;(*, *, 0:11)
rad = zfact * rad;(*, *, 0:11)
bbl = zfact * bbl;(*, *, 0:11)

;quelques bilans:
haddocON = (xad + yad)
zaddocON = zad
mlddocON = mld
zdfdocON = zdf
ldfdocON = lf
;; zlfdocON = zlf
;; hlfdocON = hlf
gmwdocON = (eiv - xad - yad - zad)
otherdocON = rof + sbc + rad + bbl

;; ; temporal integral
jpt=12
time = findgen(jpt)
haddocON = grossemoyenne(haddocON, 't', /NAN)
zaddocON = grossemoyenne(zaddocON, 't', /NAN)
mlddocON = grossemoyenne(mlddocON, 't', /NAN)
zdfdocON = grossemoyenne(zdfdocON, 't', /NAN)
ldfdocON = grossemoyenne(ldfdocON, 't', /NAN)
;; zlfdocON = grossemoyenne(zlfdocON, 't', /NAN)
;; hlfdocON = grossemoyenne(hlfdocON, 't', /NAN)
gmwdocON = grossemoyenne(gmwdocON, 't', /NAN)
otherdocON = grossemoyenne(otherdocON, 't', /INT, /NAN)

haddocON = orca2_2d_overlap_remove(haddocON) / surf
zaddocON = orca2_2d_overlap_remove(zaddocON) / surf
mlddocON = orca2_2d_overlap_remove(mlddocON) / surf
zdfdocON = orca2_2d_overlap_remove(zdfdocON) / surf
ldfdocON = orca2_2d_overlap_remove(ldfdocON) / surf
;; zlfdocON = orca2_2d_overlap_remove(zlfdocON) / surf
;; hlfdocON = orca2_2d_overlap_remove(hlfdocON) / surf
gmwdocON = orca2_2d_overlap_remove(gmwdocON) / surf
otherdocON = orca2_2d_overlap_remove(otherdocON) / surf
; unit : PmolO2/m2/an

;; ; horizontal integral
haddocON_tot = moyenne(haddocON, 'xy', /NAN, /INT)
zaddocON_tot = moyenne(zaddocON, 'xy', /NAN, /INT)
mlddocON_tot = moyenne(mlddocON, 'xy', /NAN, /INT)
zdfdocON_tot = moyenne(zdfdocON, 'xy', /NAN, /INT)
ldfdocON_tot = moyenne(ldfdocON, 'xy', /NAN, /INT)
;; zlfdocON_tot = moyenne(zlfdocON, 'xy', /NAN, /INT)
;; hlfdocON_tot = moyenne(hlfdocON, 'xy', /NAN, /INT)
gmwdocON_tot = moyenne(gmwdocON, 'xy', /NAN, /INT)

; totaux
kindocON = haddocON + zaddocON + mlddocON + gmwdocON
difdocON = zdfdocON + ldfdocON
dyndocON = kindocON + difdocON

kindocON_tot = haddocON_tot + zaddocON_tot + mlddocON_tot + gmwdocON_tot
difdocON_tot = zdfdocON_tot + ldfdocON_tot
dyndocON_tot = kindocON_tot + difdocON_tot

; ________________________________________
; SUB NH4 OFFLINE

zfactnit   = - 40./122.  * 1000. ; kmolNH4 -> molO2  nitrification amonium

subdu = restore_subdu(fdirOFF, fsubnh4OFF, 1e11) ; varying mld, correction diffusion isoneutral

subh   = subdu.subh
subw   = subdu.subw  
subgm  = subdu.subgm 
subldf = subdu.subldf
subzdf = subdu.subzdf
submld = subdu.submld
; unit : kmolO2 / maille / month

; selection de 12 mois
subh   = subh   (*, *, 12:23) 
subw   = subw   (*, *, 12:23)
submld = submld (*, *, 12:23) 
subzdf = subzdf (*, *, 12:23)
subgm  = subgm  (*, *, 12:23)
subldf = subldf (*, *, 12:23)

; TEMPORAL INTEGRAL
jpt = 12
time = findgen(jpt)
hadnh4OFF = grossemoyenne(subh  , 't', /INT, /NAN) * zfactbio / surf
zadnh4OFF = grossemoyenne(subw  , 't', /INT, /NAN) * zfactbio / surf
mldnh4OFF = grossemoyenne(submld, 't', /INT, /NAN) * zfactbio / surf                   
zdfnh4OFF = grossemoyenne(subzdf, 't', /INT, /NAN) * zfactbio / surf                   
ldfnh4OFF = grossemoyenne(subldf, 't', /INT, /NAN) * zfactbio / surf                   
gmwnh4OFF  = grossemoyenne(subgm , 't', /INT, /NAN) * zfactbio / surf                    
; unit : molO2 / m2 / ans

; to remove the overlap in the indian ocean and north pole
hadnh4OFF = orca2_2d_overlap_remove(hadnh4OFF)
zadnh4OFF = orca2_2d_overlap_remove(zadnh4OFF)
zdfnh4OFF = orca2_2d_overlap_remove(zdfnh4OFF)
mldnh4OFF = orca2_2d_overlap_remove(mldnh4OFF)
ldfnh4OFF = orca2_2d_overlap_remove(ldfnh4OFF)
gmwnh4OFF = orca2_2d_overlap_remove(gmwnh4OFF)
; unit : PmolO2/m2/year

; horizontal integral
hadnh4OFF_tot = moyenne(hadnh4OFF, 'xy', /NAN, /INT)
zadnh4OFF_tot = moyenne(zadnh4OFF, 'xy', /NAN, /INT)
mldnh4OFF_tot = moyenne(mldnh4OFF, 'xy', /NAN, /INT)
zdfnh4OFF_tot = moyenne(zdfnh4OFF, 'xy', /NAN, /INT)
gmwnh4OFF_tot = moyenne(gmwnh4OFF, 'xy', /NAN, /INT)
ldfnh4OFF_tot = moyenne(ldfnh4OFF, 'xy', /NAN, /INT)

; totaux
kinnh4OFF = hadnh4OFF + zadnh4OFF + mldnh4OFF + gmwnh4OFF
difnh4OFF = zdfnh4OFF + ldfnh4OFF
dynnh4OFF = kinnh4OFF + difnh4OFF

kinnh4OFF_tot = hadnh4OFF_tot + zadnh4OFF_tot + mldnh4OFF_tot + gmwnh4OFF_tot
difnh4OFF_tot = zdfnh4OFF_tot + ldfnh4OFF_tot
dynnh4OFF_tot = kinnh4OFF_tot + difnh4OFF_tot

; ________________________________________
; SUB NH4 ONLINE

restorefile = fdirON + fsubnh4ON

;%%% Masques U et V (sur les bords uniquement)
uumask = fltarr(182, 149, 24)+1.
vvmask = fltarr(182, 149, 24)+1.

uumask(0, *, *) = 0.
uumask(181, *, *) = 0.
uumask(*, 148, *) = 0.
uumask(92:181, 147, *) = 0.
uumask(1, 147, *) = 0.

vvmask(0, *, *) = 0.
vvmask(181, *, *) = 0.
vvmask(*, 148, *) = 0.

ttmask =  fltarr(182, 149, 24)
for tt=0, 23 DO ttmask[*, *, tt] = tmask[*, *, 0]

;Advection by u,v,w
xad = read_ncdf('xad_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
xad(where (xad GE 1.e19)) = 0.
yad = read_ncdf('yad_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
yad(where(yad GE 1.e19)) = 0.
xad = xad*uumask
yad = yad*vvmask
zad = read_ncdf('zad_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask

;advection totale= advection by u,v,w + adv GM
xeiv = read_ncdf('xei_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
xeiv(where (xeiv GE 1.e19)) = 0.
yeiv = read_ncdf('yei_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
yeiv(where (yeiv GE 1.e19)) = 0.
xeiv = xeiv*uumask
yeiv = yeiv*vvmask
zeiv = read_ncdf('zei_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask
heiv=xeiv+yeiv
eiv=xeiv+yeiv+zeiv

;Diffusion laterale
;warning: signe - sur xlf et ylf pour corriger bug dans le code
xlf = read_ncdf('xlf_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
ylf = read_ncdf('ylf_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
zlf = read_ncdf('zlf_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)
xlf(where (xlf GE 1.e19)) = 0.
ylf(where (ylf GE 1.e19)) = 0.
lf = -xlf*uumask-ylf*vvmask+zlf*ttmask
;; hlf =  -xlf*uumask-ylf*vvmask
;; zlf =  zlf*ttmask

;entrainement de la MLD
mld = read_ncdf('mld_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion verticale
zdf = read_ncdf('zdf_sub_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;surface boundary condition: effet des precipitations/evaporation
sbc = read_ncdf('sbc_mld_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;diffusion dans bottom boundary layer
bbl = read_ncdf('bbl_mld_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement de la surface libre (compte en negatif car dans la subdu)
rof = read_ncdf('rof_mld_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 
;ajustement correction des valeurs negatives 
rad = read_ncdf('rad_mld_NH4', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;Facteur de conversion pour avoir des unites en molO2/an 
zfactnit   = -40./122. ; molNH4 -> molO2  
zfact = 365. * 86400. * 1000. * zfactnit
xad = zfact * xad;(*, *, 0:11)
yad = zfact * yad;(*, *, 0:11)
zad = zfact * zad;(*, *, 0:11)
eiv = zfact * eiv;(*, *, 0:11)
lf  = zfact * lf ;(*, *, 0:11)
;; hlf = zfact * hlf;(*, *, 0:11)
;; zlf = zfact * zlf;(*, *, 0:11)
mld = zfact * mld;(*, *, 0:11)
zdf = zfact * zdf;(*, *, 0:11)
rof = zfact * rof;(*, *, 0:11)
sbc = zfact * sbc;(*, *, 0:11)
rad = zfact * rad;(*, *, 0:11)
bbl = zfact * bbl;(*, *, 0:11)

;quelques bilans:
hadnh4ON = (xad + yad)
zadnh4ON = zad
mldnh4ON = mld
zdfnh4ON = zdf
ldfnh4ON = lf
;; zlfnh4ON = zlf
;; hlfnh4ON = hlf
gmwnh4ON = (eiv - xad - yad - zad)
othernh4ON = rof + sbc + rad + bbl

;; ; temporal integral
jpt=12
time = findgen(jpt)
hadnh4ON = grossemoyenne(hadnh4ON, 't', /NAN)
zadnh4ON = grossemoyenne(zadnh4ON, 't', /NAN)
mldnh4ON = grossemoyenne(mldnh4ON, 't', /NAN)
zdfnh4ON = grossemoyenne(zdfnh4ON, 't', /NAN)
ldfnh4ON = grossemoyenne(ldfnh4ON, 't', /NAN)
;; zlfnh4ON = grossemoyenne(zlfnh4ON, 't', /NAN)
;; hlfnh4ON = grossemoyenne(hlfnh4ON, 't', /NAN)
gmwnh4ON = grossemoyenne(gmwnh4ON, 't', /NAN)
othernh4ON = grossemoyenne(othernh4ON, 't', /INT, /NAN)

hadnh4ON = orca2_2d_overlap_remove(hadnh4ON) / surf
zadnh4ON = orca2_2d_overlap_remove(zadnh4ON) / surf
mldnh4ON = orca2_2d_overlap_remove(mldnh4ON) / surf
zdfnh4ON = orca2_2d_overlap_remove(zdfnh4ON) / surf
ldfnh4ON = orca2_2d_overlap_remove(ldfnh4ON) / surf
;; zlfnh4ON = orca2_2d_overlap_remove(zlfnh4ON) / surf
;; hlfnh4ON = orca2_2d_overlap_remove(hlfnh4ON) / surf
gmwnh4ON = orca2_2d_overlap_remove(gmwnh4ON) / surf
othernh4ON = orca2_2d_overlap_remove(othernh4ON) / surf
; unit : PmolO2/m2/an

;; ; horizontal integral
hadnh4ON_tot = moyenne(hadnh4ON, 'xy', /NAN, /INT)
zadnh4ON_tot = moyenne(zadnh4ON, 'xy', /NAN, /INT)
mldnh4ON_tot = moyenne(mldnh4ON, 'xy', /NAN, /INT)
zdfnh4ON_tot = moyenne(zdfnh4ON, 'xy', /NAN, /INT)
ldfnh4ON_tot = moyenne(ldfnh4ON, 'xy', /NAN, /INT)
;; zlfnh4ON_tot = moyenne(zlfnh4ON, 'xy', /NAN, /INT)
;; hlfnh4ON_tot = moyenne(hlfnh4ON, 'xy', /NAN, /INT)
gmwnh4ON_tot = moyenne(gmwnh4ON, 'xy', /NAN, /INT)

; totaux
kinnh4ON = hadnh4ON + zadnh4ON + mldnh4ON + gmwnh4ON
difnh4ON = zdfnh4ON + ldfnh4ON
dynnh4ON = kinnh4ON + difnh4ON

kinnh4ON_tot = hadnh4ON_tot + zadnh4ON_tot + mldnh4ON_tot + gmwnh4ON_tot
difnh4ON_tot = zdfnh4ON_tot + ldfnh4ON_tot
dynnh4ON_tot = kinnh4ON_tot + difnh4ON_tot

;__________________________________________________
; EXPORT OFFLINE

zfactexp   =  160./122.  * 1.e9  ; GmolC -> molO2

; restore export
restorefile = fdirOFF + fexpOFF
print, 'file restored : ', restorefile
restore, file = restorefile, /verbose

; suppression des valeurs aberrantes
export( where( abs(export) GT 1e15 ) ) = 0.
; selection de 12 mois
export   = export   (*, *, 12:23) 

; TEMPORAL INTEGRAL
jpt = 12
time = findgen(jpt)
expOFF = grossemoyenne(export, 't', /NAN, /INT) * zfactexp / surf
; unit : PmolO2 / m2 / year
; to remove the overlap in the indian ocean and north pole
expOFF = orca2_2d_overlap_remove(expOFF)

; horizontal integral
expOFF_tot = moyenne(expOFF, 'xy', /NAN, /INT)

; ________________________________________
; EXPORT ONLINE

restorefile = fdirON + fsubtocON

;%%% Masques U et V (sur les bords uniquement)
uumask = fltarr(182, 149, 24)+1.
vvmask = fltarr(182, 149, 24)+1.

uumask(0, *, *) = 0.
uumask(181, *, *) = 0.
uumask(*, 148, *) = 0.
uumask(92:181, 147, *) = 0.
uumask(1, 147, *) = 0.

vvmask(0, *, *) = 0.
vvmask(181, *, *) = 0.
vvmask(*, 148, *) = 0.

ttmask =  fltarr(182, 149, 24)
for tt=0, 23 DO ttmask[*, *, tt] = tmask[*, *, 0]

;sedim at mld base
sin = read_ncdf('sin_mld_TOC', 12, 23, /timestep, file = restorefile,  /nostruct)*ttmask 

;Facteur de conversion pour avoir des unites en molO2/an 
zfactexp   = 160./122. ; molC -> molO2  
zfact = 365. * 86400. * 1000. * zfactexp
sin = zfact * sin 

; temporal integral
jpt=12
time = findgen(jpt)
expON = grossemoyenne(sin, 't', /NAN)

expON = orca2_2d_overlap_remove(expON) / surf
; unit : PmolO2/m2/an

;; ; horizontal integral
expON_tot = moyenne(expON, 'xy', /NAN, /INT)

; ________________________________________
; TOTAUX 

dynomON = dynnh4ON + dyndocON + dyntocON
dynomON_tot = dynnh4ON_tot + dyndocON_tot + dyntocON_tot

dynomOFF = dynnh4OFF + dyndocOFF + dyntocOFF
dynomOFF_tot = dynnh4OFF_tot + dyndocOFF_tot + dyntocOFF_tot

respON = dynomON + expON
respON_tot = dynomON_tot + expON_tot

respOFF = dynomOFF + expOFF
respOFF_tot = dynomOFF_tot + expOFF_tot

totON_tot = respON_tot + dynON_tot
totOFF_tot = respOFF_tot + dynOFF_tot

OFF = [respOFF_tot, dynOFF_tot, kinOFF_tot, difOFF_tot, hadOFF_tot, zadOFF_tot, gmwOFF_tot, mldOFF_tot, zdfOFF_tot, ldfOFF_tot, dynomOFF_tot, expOFF_tot]
ON  = [respON_tot , dynON_tot , kinON_tot , difON_tot , hadON_tot , zadON_tot , gmwON_tot , mldON_tot , zdfON_tot , ldfON_tot , dynomON_tot , expON_tot ]

; ________________________________________
; DELTA

Dhad = hadON - hadOFF
Dzad = zadON - zadOFF
Dmld = mldON - mldOFF
Dzdf = zdfON - zdfOFF
Dldf = ldfON - ldfOFF
;; Dzlf = zlfON - zlfOFF
;; Dhlf = hlfON - hlfOFF
Dgmw = gmwON - gmwOFF

Dkin = kinON - kinOFF
Ddif = difON - difOFF
Ddyn = dynON - dynOFF

;==================================================
; PLOT
;==================================================

stop


w = WINDOW(/BUFFER, DIMENSIONS = [600, 200]) ; dimension = [width, height]

props = {NBARS:2, LINESTYLE:'-', THICK:0.5, OVERPLOT:1, BUFFER:1}

zfact = 1e-15
p0 = plot(0.*findgen(13), /CURRENT, /BUFFER, LINESTYLE = 'none')
b01 = barplot(zfact * ON   , index=0, fill_color='black'      , _EXTRA = props)
b02 = barplot(zfact * OFF  , index=1, fill_color='royal blue' , _EXTRA = props)

b01.NAME  = 'ONLINE'
b02.NAME  = 'OFFLINE'

ax               = p0.AXES
p0.TITLE         = 'O2 subduction'
p0.FONT_SIZE     = 8
p0.XTICKVALUES       = [0,1,2,3,4,5,6,7,8,9,10,11]
p0.XTICKNAME         = ['RES','DYN', 'KIN', 'DIF', 'HAD', 'ZAD', 'GMW', 'MLD', 'ZDF', 'LDF','DYNOM','SEDIM']
p0.XTEXT_ORIENTATION = 25
p0.XTHICK            = 0
p0.XTICKLEN          = 0
p0.YTITLE        = 'Deoxygenation       [PmolO2/year]       Oxygenation'
p0.YRANGE        = [-1.2, 1.2]
p0.YTICKINTERVAL = .2
p0.YMINOR        = 0
p0.YTICKLEN      = 1.0
p0.YTHICK        = 0.5
p0.YGRIDSTYLE    = 2
ax[2].HIDE  = 1
ax[3].MAJOR = 0
ax[3].MINOR = 0

leg = legend(TARGET = [b01, b02], /AUTO_TEXT_COLOR, LINESTYLE = '-', $
             TRANSPARENCY = 0, SHADOW = 0, SAMPLE_WIDTH = 0.07, $
             POSITION = [0.5, 0.3], FONT_SIZE = 6)

filename = '/data/daclod/FIG/' + namefigbase + '.eps'
b01.save, filename

stop

@init_orca2_shift

unit2 = 'Deoxygenation            molO2/m2            Oxygenation'
props2 = {CB_TITLE:unit2, LCT:42, NOCONTOUR:1, REALCONT:1, SUBTITLE:'', FORMAT:'(e11.2)', INV:1}

saveplot = namefigbase + '_dyn_kin_dif'
print, 'plot saved: ', saveplot
openps, file = saveplot, /landscape

plt, shift(dynON (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 1], title = 'dyn ON'   , _extra = props2, win = 1, /landscape
plt, shift(dynOFF(1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 2], title = 'dyn OFF'  , _extra = props2, /noerase
plt, shift(Ddyn  (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 3], title = 'Delta dyn', _extra = props2, /noerase

plt, shift(kinON (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 4], title = 'kin ON'   , _extra = props2, /noerase
plt, shift(kinOFF(1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 5], title = 'kin OFF'  , _extra = props2, /noerase
plt, shift(Dkin  (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 6], title = 'Delta kin', _extra = props2, /noerase

plt, shift(difON (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 7], title = 'dif ON'   , _extra = props2, /noerase
plt, shift(difOFF(1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 8], title = 'dif OFF'  , _extra = props2, /noerase
plt, shift(Ddif  (1:180, 0:147), [30, 0]), -20, 20, small=[3, 3, 9], title = 'Delta dif', _extra = props2, /noerase

closeps


saveplot = namefigbase + '_had_zad_gmw_mld_zdf_ldf_OFF'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait

plt, shift(hadOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 1], title = 'had OFF', _extra = props2, win = 1, /portrait
plt, shift(zadOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 2], title = 'zad OFF', _extra = props2, /noerase
plt, shift(gmwOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 3], title = 'gmw OFF', _extra = props2, /noerase
plt, shift(mldOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 4], title = 'mld OFF', _extra = props2, /noerase
plt, shift(zdfOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 5], title = 'zdf OFF', _extra = props2, /noerase
plt, shift(ldfOFF(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 6], title = 'ldf OFF', _extra = props2, /noerase

closeps

saveplot = namefigbase + '_had_zad_gmw_mld_zdf_ldf_ON'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait

plt, shift(hadON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 1], title = 'had ON', _extra = props2, win = 1, /portrait
plt, shift(zadON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 2], title = 'zad ON', _extra = props2, /noerase
plt, shift(gmwON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 3], title = 'gmw ON', _extra = props2, /noerase
plt, shift(mldON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 4], title = 'mld ON', _extra = props2, /noerase
plt, shift(zdfON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 5], title = 'zdf ON', _extra = props2, /noerase
plt, shift(ldfON(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 6], title = 'ldf ON', _extra = props2, /noerase

closeps

saveplot = namefigbase + '_had_zad_gmw_mld_zdf_ldf_DELTA'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait

plt, shift(Dhad(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 1], title = 'had, Delta=ON-OFF', _extra = props2, win = 1, /portrait
plt, shift(Dzad(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 2], title = 'zad, Delta=ON-OFF', _extra = props2, /noerase
plt, shift(Dgmw(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 3], title = 'gmw, Delta=ON-OFF', _extra = props2, /noerase
plt, shift(Dmld(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 4], title = 'mld, Delta=ON-OFF', _extra = props2, /noerase
plt, shift(Dzdf(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 5], title = 'zdf, Delta=ON-OFF', _extra = props2, /noerase
plt, shift(Dldf(1:180, 0:147), [30, 0]), -20, 20, small=[2, 3, 6], title = 'ldf, Delta=ON-OFF', _extra = props2, /noerase

closeps


saveplot = namefigbase + '_respiration'
print, 'plot saved: ', saveplot
openps, file = saveplot, /portrait

plt, shift(respON  (1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 1], title = 'resp ON'  , _extra = props2, win = 1, /portrait
plt, shift(respOFF (1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 2], title = 'resp OFF' , _extra = props2, /noerase
plt, shift(expON   (1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 3], title = 'sedi ON'  , _extra = props2, /noerase
plt, shift(expOFF  (1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 4], title = 'sedi OFF' , _extra = props2, /noerase
plt, shift(dynomON (1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 5], title = 'subOM ON' , _extra = props2, /noerase
plt, shift(dynomOFF(1:180, 0:147), [30, 0]), -10, 10, small=[2, 3, 6], title = 'subOM OFF', _extra = props2, /noerase

closeps


stop
END
