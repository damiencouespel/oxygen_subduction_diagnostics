FUNCTION restore_subdu,  vfdir, vfile, nanval, FILEZDF1 = filezdf1, FILEZDF2 = filezdf2
  ; Function to read subduction
  ; diagnostics and make them ready to
  ; plot
    
  ; restore subdu
  restorefile = vfdir + vfile
  print, 'file restored : ', restorefile
  restore, file = restorefile, /verbose
  ; to adjust kz to 1.2e-4
  subzdf = subzdf * 1.2
  
  ; restore filezdf1 in case subduction
  ; vertical diffusion is in a separated file
  IF KEYWORD_SET(filezdf1) THEN BEGIN
     restorefile = vfdir + filezdf1
     print, 'file restored : ', restorefile
     restore, file = restorefile, /verbose
     ; to adjust kz to 1.2e-4
     subzdf = subzdf * 1.2
  ENDIF

  ; restore subzdf2 to remove the
  ; vertical component of isoneutral
  ; mixing from subzdf and to add it in subldf
  IF KEYWORD_SET(filezdf2) THEN BEGIN
     subzdf_old = subzdf        ; save old subzdf
     restorefile = vfdir + filezdf2
     print, 'file restored : ', restorefile
     restore, file = restorefile, /verbose
     subldf = subldf + (subzdf_old - subzdf) ; add vertical component of iso-neutral mixing to subldf
  ENDIF

  ; remove outliers
  subw(   where( abs(subw)   gt nanval ) ) = !VALUES.F_NAN
  subh(   where( abs(subh)   gt nanval ) ) = !VALUES.F_NAN
  subgm(  where( abs(subgm)  gt nanval ) ) = !VALUES.F_NAN
  subldf( where( abs(subldf) gt nanval ) ) = !VALUES.F_NAN
  subzdf( where( abs(subzdf) gt nanval ) ) = !VALUES.F_NAN
  submld( where( abs(submld) gt nanval ) ) = !VALUES.F_NAN

  ; to remove the overlap in the indian ocean DC2
  ; and north pole
  subh   = orca2_3d_overlap_remove(subh  )
  subw   = orca2_3d_overlap_remove(subw  )
  subzdf = orca2_3d_overlap_remove(subzdf)
  submld = orca2_3d_overlap_remove(submld)
  subldf = orca2_3d_overlap_remove(subldf)
  subgm  = orca2_3d_overlap_remove(subgm )
  
  dataout = {subh:subh, subw:subw, submld:submld, subzdf:subzdf, subldf:subldf, subgm:subgm}
  RETURN, dataout

END
