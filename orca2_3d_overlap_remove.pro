FUNCTION orca2_3d_overlap_remove,  tr
  ; function to remove the cyclical part of the mesh grid
  ; specific to ORCA grid
  tr[0, *, *]      = !VALUES.F_NAN
  tr[181, *, *]    = !VALUES.F_NAN
  tr[*, 148, *]    = !VALUES.F_NAN
  tr[0:91, 147, *] = !VALUES.F_NAN
  RETURN, tr
END
