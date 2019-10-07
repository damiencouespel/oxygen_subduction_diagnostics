FUNCTION bilan_global_subdu_vs_time, vsubdu, vzfact
  ; return the global integral of subduction terms

  @common  
  @init_run_clim_dc
  surf = e1t*e2t
  
  subh   = vsubdu.subh
  subw   = vsubdu.subw  
  subgm  = vsubdu.subgm 
  subldf = vsubdu.subldf
  subzdf = vsubdu.subzdf
  submld = vsubdu.submld
  ; unit : kmolO2 / maille / month
  
  ; annual integral
  jpt = 12
  had = FLTARR(jpi, jpj, 110) ; lateral advection                               
  zad = FLTARR(jpi, jpj, 110) ; vertical advection                              
  zdf = FLTARR(jpi, jpj, 110) ; entrainment                                     
  ldf = FLTARR(jpi, jpj, 110) ; vertical diffusion                              
  gm  = FLTARR(jpi, jpj, 110) ; advection by bolus velocities (Gent McWilliams) 
  mld = FLTARR(jpi, jpj, 110) ; isopycnal diffusion                             
  FOR yy  = 0, 109 DO BEGIN
     had(*, *, yy) = grossemoyenne(subh(*, *, yy*12 : yy*12 + 11)  , 't', /INT, /NAN) * vzfact / surf
     zad(*, *, yy) = grossemoyenne(subw(*, *, yy*12 : yy*12 + 11)  , 't', /INT, /NAN) * vzfact / surf
     mld(*, *, yy) = grossemoyenne(submld(*, *, yy*12 : yy*12 + 11), 't', /INT, /NAN) * vzfact / surf                   
     zdf(*, *, yy) = grossemoyenne(subzdf(*, *, yy*12 : yy*12 + 11), 't', /INT, /NAN) * vzfact / surf                   
     gm(*, *, yy)  = grossemoyenne(subgm(*, *, yy*12 : yy*12 + 11) , 't', /INT, /NAN) * vzfact / surf                    
     ldf(*, *, yy) = grossemoyenne(subldf(*, *, yy*12 : yy*12 + 11), 't', /INT, /NAN) * vzfact / surf                   
  ENDFOR
  ; unit : kmolO2 * vzfact / m2 / y

  ; remove overlap on the north pole and east-west
  had = orca2_3d_overlap_remove(had)
  zad = orca2_3d_overlap_remove(zad)
  mld = orca2_3d_overlap_remove(mld)
  zdf = orca2_3d_overlap_remove(zdf)
  gm  = orca2_3d_overlap_remove(gm )
  ldf = orca2_3d_overlap_remove(ldf)

  ; horizontal integral
  jpt = 110
  had = grossemoyenne(had, 'xy', /NAN, /INT)
  zad = grossemoyenne(zad, 'xy', /NAN, /INT)
  mld = grossemoyenne(mld, 'xy', /NAN, /INT)
  zdf = grossemoyenne(zdf, 'xy', /NAN, /INT)
  gm  = grossemoyenne(gm , 'xy', /NAN, /INT)
  ldf = grossemoyenne(ldf, 'xy', /NAN, /INT)
  ; unit : kmolO2 * vzfact / y
  
  ; totaux
  edd = ldf + gm  ; eddies          
  zmx = zdf + mld ; vertical mixing 
  adv = had + zad ; advection       
  dyn = adv + zmx + edd
  
  dataout = {dyn:dyn, adv:adv, zmx:zmx, edd:edd, had:had, zad:zad, mld:mld, zdf:zdf, ldf:ldf, gm:gm}
  RETURN, dataout

END
