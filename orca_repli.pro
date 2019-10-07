FUNCTION orca_repli,  tr
; reconstruction of cyclical
; conditions on x, y (specific to ORCA
; grid)
@common
; cyclical condition in E-W:
tr(jpi-1, *, *) = tr(1, *, *)
;cyclical condtion at north pole:
FOR ji = 1, jpj/2 DO BEGIN 
tr(jpi-1-ji, jpj-2, *) = tr(ji+1, jpj-2, *)
tr(ji+1, jpj-1, *) = tr(jpi-1-ji, jpj-3, *)
ENDFOR

RETURN, tr
END
