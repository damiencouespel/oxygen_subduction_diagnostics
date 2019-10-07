; Wed 7 May 2014 : DC comment xnotice, widget_control
;+
;
; @file_comments
; read NetCDF meshmask file created by OPA
;
; @categories
; Grid
;
; @examples
;
;   IDL> ncdf_meshread [,' filename']
;
; @param filename {in}{optional}{default='meshmask.nc'}{type=scalar string}
;    Name of the meshmask file to read. If this name does not contain any "/"
;    and if iodirectory keyword is not specified, then the common variable
;    iodir will be used to define the mesh file path.
;
; @keyword GLAMBOUNDARY {default=those defined in the file}{type=2 elements vector}
;    Longitude boundaries that should be used to visualize the data.
;      lon2 > lon1
;      lon2 - lon1 le 360
;    By default, the common (cm_4mesh) variable key_shift will be automatically
;    defined according to GLAMBOUNDARY.
;
; @keyword CHECKDAT
; Suppressed. Use <pro>micromeshmask</pro> to create an appropriate meshmask.
;
; @keyword ONEARTH {default=1}{type=scalar: 0 or 1}
;    Force the manual definition of data localization on the earth or not
;       0) if the data are not on the earth
;       1) if the data are on earth (in that case we can for example use
;          the labels 'longitude', 'latitude' in plots).
;    The resulting value will be stored in the common (cm_4mesh) variable key_onearth
;    ONEARTH=0 forces PERIODIC=0, SHIFT=0 and is cancelling GLAMBOUNDARY
;
; @keyword GETDIMENSIONS {default=0}{type=scalar: 0 or 1}
;    Activate this keywords if you only want to know the dimension
;    of the domain stored in the mesh file. This dimension will be
;    defined in jpiglo, jpjglo, jpkglo (cm_4mesh common variables)
;
; @keyword PERIODIC {default=computed by using the first line of glamt}{type=scalar: 0 or 1}
;    Force the manual definition of the grid zonal periodicity.
;    The resulting value will be stored in the common (cm_4mesh) variable key_periodic
;    PERIODIC=0 forces SHIFT=0
;
; @keyword SHIFT {default=computed according to glamboundary}{type=scalar}
;    Force the manual definition of the zonal shift that must be apply to the data.
;    The resulting value will be stored in the common (cm_4mesh) variable key_shift
;    Note that if key_periodic=0 then in any case key_shift = 0.
;
; @keyword STRCALLING {type=scalar string}
;    the calling command used to call <pro>computegrid</pro> (this is used by <pro>xxx</pro>)
;
; @keyword STRIDE {default=[1, 1, 1]}{type=3 elements vector}
;    Specify the stride in x, y and z direction. The resulting
;    value will be stored in the common (cm_4mesh) variable key_stride
;
; @keyword _EXTRA
; Used to pass keywords to <pro>isafile</pro>
;
; @uses
; <pro>cm_4mesh</pro>
; <pro>cm_4data</pro>
; <pro>cm_4cal</pro>
;
; @restrictions
; ixminmesh, ixmaxmesh, iyminmesh, iymaxmesh, izminmesh, izmaxmesh must
; be defined before calling ncdf_meshread. If some of those values
; are equal to -1 they will be automatically defined
;
; @history
; Sebastien Masson (smasson\@lodyc.jussieu.fr)
;                      12/1999
; July 2004, Sebastien Masson: Several modifications (micromeshmask,
; clean partial steps, clean use of key_stride, automatic definition
; of key_shift, ...)
; Oct. 2004, Sebastien Masson: add PERIODIC and SHIFT
; Aug. 2005, Sebastien Masson: some cleaning + english
;
; @version
; $Id: ncdf_meshread.pro 458 2011-04-12 07:45:09Z smasson $
;
;-
PRO ncdf_meshread_dc1, filename, GLAMBOUNDARY=glamboundary, CHECKDAT=checkdat $
                  , ONEARTH=onearth, GETDIMENSIONS=getdimensions $
                  , PERIODIC=periodic, SHIFT=shift, STRIDE=stride $
                  , STRCALLING=strcalling, _EXTRA=ex
;
  compile_opt idl2, strictarrsubs
;
@cm_4mesh
@cm_4data
@cm_4cal
  IF NOT keyword_set(key_forgetold) THEN BEGIN
@updatenew
@updatekwd
  ENDIF
;
  tempsun = systime(1)          ; for key_performance
  IF keyword_set(CHECKDAT) THEN BEGIN
    ras = report([' The keyword CHECKDAT has been suppressed (it could create bugs).', $
    ' Remove it from the call of ncdf_meshread', $
    ' Please use smallmeshmask.pro or micromeshmask.pro to create a', $
    ' meshmask that has manageable size'])
    return
  ENDIF
;-------------------------------------------------------
; find meshfile name and open it!
;-------------------------------------------------------
; def of filename by default
  IF n_params() EQ 0 then filename = 'meshmask.nc'
  meshname = isafile(file = filename, iodirectory = iodir, _EXTRA = ex)
  meshname = meshname[0]
;
;  noticebase = xnotice('Reading file !C '+meshname+'!C ...') ;DC
; if the meshmask is on tape archive ... get it back
  IF !version.OS_FAMILY EQ 'unix' THEN spawn, '\file '+meshname+' > /dev/null'
  cdfid = ncdf_open(meshname)
  inq = ncdf_inquire(cdfid)
;------------------------------------------------------------
; dimensions
;------------------------------------------------------------
  ncdf_diminq, cdfid, 'x', name, jpiglo
  ncdf_diminq, cdfid, 'y', name, jpjglo
  listdims = strlowcase(ncdf_listdims(cdfid))
  IF (where(listdims EQ 'z'))[0] NE -1 THEN ncdf_diminq, cdfid, 'z', name, jpkglo ELSE BEGIN
    dimid = (where(strmid(listdims, 0, 5) EQ 'depth'))[0]
    IF dimid NE -1 THEN ncdf_diminq, cdfid, dimid, name, jpkglo ELSE BEGIN
      dummy = report('We could not find the vertical dimension..., its name must be z or start with depth')
      stop
    ENDELSE
  ENDELSE
;
  if keyword_set(getdimensions) then begin
;DC    widget_control, noticebase, bad_id = nothing, /destroy
    ncdf_close,  cdfid
    return
  endif
;-------------------------------------------------------
; check that all i[xyz]min[ax]mesh are well defined
;-------------------------------------------------------
  if n_elements(ixminmesh) EQ 0 THEN ixminmesh = 0
  if n_elements(ixmaxmesh) EQ 0 then ixmaxmesh = jpiglo-1
  if ixminmesh EQ -1 THEN ixminmesh = 0
  IF ixmaxmesh EQ -1 then ixmaxmesh = jpiglo-1
  if n_elements(iyminmesh) EQ 0 THEN iyminmesh = 0
  IF n_elements(iymaxmesh) EQ 0 then iymaxmesh = jpjglo-1
  if iyminmesh EQ -1 THEN iyminmesh = 0
  IF iymaxmesh EQ -1 then iymaxmesh = jpjglo-1
  if n_elements(izminmesh) EQ 0 THEN izminmesh = 0
  IF n_elements(izmaxmesh) EQ 0 then izmaxmesh = jpkglo-1
  if izminmesh EQ -1 THEN izminmesh = 0
  IF izmaxmesh EQ -1 then izmaxmesh = jpkglo-1
; definition of jpi,jpj,jpj
  jpi = long(ixmaxmesh-ixminmesh+1)
  jpj = long(iymaxmesh-iyminmesh+1)
  jpk = long(izmaxmesh-izminmesh+1)
;-------------------------------------------------------
; check onearth and its consequences
;-------------------------------------------------------
  IF n_elements(onearth) EQ 0 THEN key_onearth = 1 $
  ELSE key_onearth = keyword_set(onearth)
  IF NOT key_onearth THEN BEGIN
    periodic = 0
    shift = 0
  ENDIF
;-------------------------------------------------------
; automatic definition of key_periodic
;-------------------------------------------------------
  IF n_elements(periodic) EQ 0 THEN BEGIN
    IF jpi GT 1 THEN BEGIN
      varinq = ncdf_varinq(cdfid, 'glamt')
      CASE varinq.ndims OF
        2:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh], count = [jpi, 1]
        3:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh, 0], count = [jpi, 1, 1]
        4:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh, 0, 0], count = [jpi, 1, 1, 1]
      ENDCASE
      xaxis = (xaxis+720) MOD 360
      xaxis = xaxis[sort(xaxis)]
      key_periodic = (xaxis[jpi-1]+2*(xaxis[jpi-1]-xaxis[jpi-2])) $
                     GE (xaxis[0]+360)
    ENDIF ELSE key_periodic = 0
  ENDIF ELSE key_periodic = keyword_set(periodic)
;-------------------------------------------------------
; automatic definition of key_shift
;-------------------------------------------------------
  IF n_elements(shift) EQ 0 THEN BEGIN
    key_shift = long(testvar(var = key_shift))
;  key_shift will be defined according to the first line of glamt.
    if keyword_set(glamboundary) AND jpi GT 1 AND key_periodic EQ 1 $
    THEN BEGIN
      varinq = ncdf_varinq(cdfid, 'glamt')
      CASE varinq.ndims OF
        2:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh], count = [jpi, 1]
        3:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh, 0], count = [jpi, 1, 1]
        4:ncdf_varget, cdfid, 'glamt', xaxis $
                       , offset = [ixminmesh, iyminmesh, 0, 0], count = [jpi, 1, 1, 1]
      ENDCASE
; xaxis between glamboundary[0] and glamboundary[1]
      xaxis = xaxis MOD 360
      smaller = where(xaxis LT glamboundary[0])
      if smaller[0] NE -1 then xaxis[smaller] = xaxis[smaller]+360
      bigger = where(xaxis GE glamboundary[1])
      if bigger[0] NE -1 then xaxis[bigger] = xaxis[bigger]-360
;
      key_shift = (where(xaxis EQ min(xaxis)))[0]
      IF key_shift NE 0 THEN BEGIN
        key_shift = jpi-key_shift
        xaxis = shift(xaxis, key_shift)
      ENDIF
;
      IF array_equal(sort(xaxis), lindgen(jpi)) NE 1 THEN BEGIN
        ras = report(['the x axis (1st line of glamt) is not sorted in the increasing order after the automatic definition of key_shift', $
        'please use the keyword shift (and periodic) to suppress the automatic definition of key_shift (and key_periodic) and define by hand a more suitable value...'])
;DC        widget_control, noticebase, bad_id = nothing, /destroy
        return
      ENDIF
;
    ENDIF ELSE key_shift = 0
  ENDIF ELSE key_shift = long(shift)*(key_periodic EQ 1)
;-------------------------------------------------------
; check key_stride and related things
;-------------------------------------------------------
  if n_elements(stride) eq 3 then key_stride = stride
  if n_elements(key_stride) LE 2 then key_stride = [1, 1, 1]
  key_stride = 1l > long(key_stride)
  IF total(key_stride) NE 3  THEN BEGIN
    IF key_shift NE 0 THEN BEGIN
; for explanation, see header of read_ncdf_varget.pro
      jpiright = key_shift
      jpileft = jpi - key_shift - ( (key_stride[0]-1)-((key_shift-1) MOD key_stride[0]) )
      jpi = ((jpiright-1)/key_stride[0]+1) + ((jpileft-1)/key_stride[0]+1)
    ENDIF ELSE jpi = (jpi-1)/key_stride[0]+1
    jpj = (jpj-1)/key_stride[1]+1
    jpk = (jpk-1)/key_stride[2]+1
  ENDIF
;-------------------------------------------------------
; default definitions to be able to use read_ncdf_varget
;-------------------------------------------------------
; default definitions to be able to use read_ncdf_varget
  ixmindtasauve = testvar(var = ixmindta)
  iymindtasauve = testvar(var = iymindta)
  izmindtasauve = testvar(var = izmindta)
;
  ixmindta = 0l
  iymindta = 0l
  izmindta = 0l
;
  jpt = 1
  time = 1
  firsttps = 0
;
  firstx = 0
  lastx = jpi-1
  firsty = 0
  lasty = jpj-1
  firstz = 0
  lastz = jpk-1
  nx = jpi
  ny = jpj
  nz = 1
  izminmeshsauve = izminmesh
  izminmesh = 0
;
  key_yreverse = 0
  key_zreverse = 0
  key_gridtype = 'c'
;-------------------------------------------------------
; 2d arrays:
;-------------------------------------------------------
; list the 2d variables that must be read
  namevar = ['glamt', 'glamu', 'glamv', 'glamf' $
             , 'gphit', 'gphiu', 'gphiv', 'gphif' $
             , 'e1t', 'e1u', 'e1v', 'e1f' $
             , 'e2t', 'e2u', 'e2v', 'e2f']
; for the variables related to the partial steps
  allvarname =  ncdf_listvars(cdfid)
;
  key_partialstep = 0
  hdept = -1
  hdepw = -1
  IF (where(allvarname EQ 'hdept'))[0] NE -1 AND jpk EQ jpkglo THEN BEGIN
    key_partialstep = 1
    namevar = [namevar, 'hdept', 'hdepw']
  ENDIF
; is gdept a real 3D array?
  IF (where(allvarname EQ 'gdept' ))[0] NE -1 AND jpk EQ jpkglo THEN BEGIN
    varinq = ncdf_varinq(cdfid, 'gdept')
    if varinq.ndims GE 3 THEN BEGIN
      ncdf_diminq, cdfid, varinq.dim[0], name, iii
      ncdf_diminq, cdfid, varinq.dim[1], name, jjj
      ncdf_diminq, cdfid, varinq.dim[2], name, kkk
      IF iii EQ jpiglo AND jjj EQ jpjglo AND kkk EQ jpkglo THEN key_gdep_3d = 1
    ENDIF   
  ENDIF   
; for compatibility with old versions of meshmask/partial steps
  e3t_ps = -1
  e3w_ps = -1
  IF (where(allvarname EQ 'e3tp'  ))[0] NE -1 THEN namevar = [namevar, 'e3tp', 'e3wp']
  IF (where(allvarname EQ 'e3t_ps'))[0] NE -1 THEN namevar = [namevar, 'e3t_ps', 'e3w_ps' ]
; is e3t a real 3D array?
  IF (where(allvarname EQ 'e3t'   ))[0] NE -1 AND jpk EQ jpkglo THEN BEGIN
    varinq = ncdf_varinq(cdfid, 'e3t')
    if varinq.ndims GE 3 THEN BEGIN
      ncdf_diminq, cdfid, varinq.dim[0], name, iii
      ncdf_diminq, cdfid, varinq.dim[1], name, jjj
      ncdf_diminq, cdfid, varinq.dim[2], name, kkk
      IF iii EQ jpiglo AND jjj EQ jpjglo AND kkk EQ jpkglo THEN BEGIN
        key_e3_3d = 1
        key_partialstep = 1
      ENDIF
    ENDIF   
  ENDIF
;
; read all the 2d variables
;
  for i = 0, n_elements(namevar)-1 do begin
    varinq = ncdf_varinq(cdfid, namevar[i])
    name = varinq.name
@read_ncdf_varget
    CASE namevar[i] OF
      'glamt':glamt = float(temporary(res))
      'glamu':glamu = float(temporary(res))
      'glamv':glamv = float(temporary(res))
      'glamf':glamf = float(temporary(res))
      'gphit':gphit = float(temporary(res))
      'gphiu':gphiu = float(temporary(res))
      'gphiv':gphiv = float(temporary(res))
      'gphif':gphif = float(temporary(res))
      'e1t':e1t = float(temporary(res))
      'e1u':e1u = float(temporary(res))
      'e1v':e1v = float(temporary(res))
      'e1f':e1f = float(temporary(res))
      'e2t':e2t = float(temporary(res))
      'e2u':e2u = float(temporary(res))
      'e2v':e2v = float(temporary(res))
      'e2f':e2f = float(temporary(res))
      'hdept':hdept = float(temporary(res))
      'hdepw':hdepw = float(temporary(res))
      'e3t_ps':e3t_ps = temporary(res)
      'e3tp':e3t_ps = temporary(res)
      'e3w_ps':e3w_ps = temporary(res)
      'e3wp':e3w_ps = temporary(res)
    ENDCASE
  ENDFOR
; in the case of key_stride ne [1, 1, 1] redefine f points
; coordinates: they must be in the middle of 3 T points
  if key_stride[0] NE 1 OR key_stride[1] NE 1 then BEGIN
; we must recompute glamf and gphif...
    IF jpi GT 1 THEN BEGIN
      if (keyword_set(key_onearth) AND keyword_set(xnotsorted)) $
        OR (keyword_set(key_periodic) AND key_irregular) then BEGIN
        stepxf = (glamt + 720) MOD 360
        stepxf = shift(stepxf, -1, -1) - stepxf
        stepxf = [ [[stepxf]], [[stepxf + 360]], [[stepxf - 360]] ]
        stepxf = min(abs(stepxf), dimension = 3)
        IF NOT keyword_set(key_periodic) THEN $
          stepxf[jpi-1, *] = stepxf[jpi-2, *]
      ENDIF ELSE BEGIN
        stepxf = shift(glamt, -1, -1) - glamt
        IF keyword_set(key_periodic) THEN $
          stepxf[jpi-1, *] = 360 + stepxf[jpi-1, *] $
        ELSE stepxf[jpi-1, *] = stepxf[jpi-2, *]
      ENDELSE
      IF jpj GT 1 THEN BEGIN
        stepxf[*, jpj-1] = stepxf[*, jpj-2]
        stepxf[jpi-1, jpj-1] = stepxf[jpi-2, jpj-2]
      ENDIF
      glamf = glamt + 0.5 * stepxf
    ENDIF ELSE glamf = glamt + 0.5
    IF jpj GT 1 THEN BEGIN
; we must compute stepyf: y distance between T(i,j) T(i+1,j+1)
      stepyf = shift(gphit, -1, -1) - gphit
      stepyf[*, jpj-1] = stepyf[*, jpj-2]
      IF jpi GT 1 THEN BEGIN
        if NOT keyword_set(key_periodic) THEN $
          stepyf[jpi-1, *] = stepyf[jpi-2, *]
        stepyf[jpi-1, jpj-1] = stepyf[jpi-2, jpj-2]
      ENDIF
      gphif = gphit + 0.5 * stepyf
    ENDIF ELSE gphif = gphit + 0.5
  ENDIF
;-------------------------------------------------------
; 3d arrays:
;-------------------------------------------------------
  nz = jpk
  izminmesh = izminmeshsauve
;
  listdims = ncdf_listdims(cdfid)
  micromask = (where(listdims EQ 'y_m'))[0]
;
  varinq = ncdf_varinq(cdfid, 'tmask')
  name = varinq.name
  IF micromask NE -1 THEN BEGIN
; keep original values
    iyminmeshtrue = iyminmesh
    key_stridetrue = key_stride
    yyy1 = firsty*key_stridetrue[1]+iyminmeshtrue
    yyy2 = lasty*key_stridetrue[1]+iyminmeshtrue
; the mask is stored as the bit values of the byte array (along the y
; dimension, see micromeshmask.pro)...
; we must modify several parameters...
    iyminmesh = 0L
    firsty = yyy1/8
    lasty = yyy2/8
    ny = lasty-firsty+1
    key_stride = [key_stride[0], 1, key_stride[2]]
@read_ncdf_varget
    tmask = bytarr(jpi, jpj, jpk)
; now we must get back the mask
; loop on the level to save memory (the loop is short and, thus,
; should be fast enough)
    FOR k = 0, jpk-1 DO BEGIN
      zzz = transpose(res[*, *, k])
      zzz = reform(binary(zzz), 8*ny, nx, /over)
      zzz = transpose(temporary(zzz))
      zzz = zzz[*, yyy1 MOD 8: 8*ny - 8 + yyy2 MOD 8]
      IF key_stridetrue[1] NE 1 THEN BEGIN
;        IF float(strmid(!version.release,0,3)) LT 5.6 THEN BEGIN
        nnny = (size(zzz))[2]
        yind = key_stridetrue[1]*lindgen((nnny-1)/key_stridetrue[1]+1)
        tmask[*, *, k] = temporary(zzz[*, yind])
;        ENDIF ELSE tmask[*, *, k] = temporary(zzz[*, 0:*:key_stridetrue[1]])
      ENDIF ELSE tmask[*, *, k] = temporary(zzz)
    ENDFOR
  ENDIF ELSE BEGIN
@read_ncdf_varget
    tmask = byte(res)
  ENDELSE
; if e3t is true 3D array...
  IF keyword_set(key_e3_3d)  THEN BEGIN
    bottom = 0 > ( total(tmask, 3)-1 )
    bottom = lindgen(jpi, jpj) + jpi*jpj*temporary(bottom)
;
    varinq = ncdf_varinq(cdfid, 'e3t')
    name = varinq.name
@read_ncdf_varget
    e3t_ps = (temporary(res))[bottom] * tmask[*, *, 0]
;
    varinq = ncdf_varinq(cdfid, 'e3w')
    name = varinq.name
@read_ncdf_varget
    e3w_ps = (temporary(res))[bottom] * tmask[*, *, 0]
  ENDIF
; if gdep is true 3D array...
  IF keyword_set(key_gdep_3d) THEN BEGIN
    bottom = 0 > ( total(tmask, 3)-1 )
    bottom = lindgen(jpi, jpj) + jpi*jpj*temporary(bottom)
;
    varinq = ncdf_varinq(cdfid, 'gdept')
    name = varinq.name
@read_ncdf_varget
    hdept = (temporary(res))[bottom] * tmask[*, *, 0]
;
    bottom = jpi*jpj + temporary(bottom)
    varinq = ncdf_varinq(cdfid, 'gdepw')
    name = varinq.name
@read_ncdf_varget
    hdepw = (temporary(res))[bottom] * tmask[*, *, 0]
  ENDIF
; boundary conditions used to compute umask.
  IF ncdf_varid(cdfid, 'umask') NE -1 THEN BEGIN 
     varinq = ncdf_varinq(cdfid, 'umask')
     name = varinq.name
     nx = 1L
     firstx = jpi-1
     lastx = jpi-1
     IF micromask NE -1 THEN BEGIN
@read_ncdf_varget
        umaskred = reform(binary(res), 8*ny, jpk, /over)
        umaskred = umaskred[yyy1 MOD 8: 8*ny - 8 + yyy2 MOD 8, *]
        IF key_stridetrue[1] NE 1 THEN umaskred = temporary(umaskred[yind, *])
     ENDIF ELSE BEGIN
@read_ncdf_varget
        umaskred = reform(byte(res), /over)
     ENDELSE
  ENDIF ELSE umaskred = bytarr(jpj, jpk)
; boundary conditions used to compute fmask (1).
  IF ncdf_varid(cdfid, 'fmask') NE -1 THEN BEGIN 
     varinq = ncdf_varinq(cdfid, 'fmask')
     name = varinq.name
     IF micromask NE -1 THEN BEGIN
@read_ncdf_varget
        fmaskredy = reform(binary(res), 8*ny, jpk, /over)
        fmaskredy = fmaskredy[yyy1 MOD 8: 8*ny - 8 + yyy2 MOD 8, *]
        IF key_stridetrue[1] NE 1 THEN fmaskredy = temporary(fmaskredy[yind, *])
     ENDIF ELSE BEGIN
@read_ncdf_varget
        fmaskredy = reform(byte(res), /over)
        fmaskredy = temporary(fmaskredy) MOD 2
     ENDELSE
  ENDIF ELSE fmaskredy = bytarr(jpj, jpk)
; boundary conditions used to compute vmask
  IF ncdf_varid(cdfid, 'vmask') NE -1 THEN BEGIN 
     varinq = ncdf_varinq(cdfid, 'vmask')
     name = varinq.name
     nx = jpi
     firstx = 0L
     lastx = jpi-1L
     ny = 1L
     firsty = jpj-1
     lasty = jpj-1
     IF micromask NE -1 THEN BEGIN
        yyy1 = firsty*key_stridetrue[1]+iyminmeshtrue
        yyy2 = lasty*key_stridetrue[1]+iyminmeshtrue
        iyminmesh = 0L
        firsty = yyy1/8
        lasty = yyy2/8
        ny = lasty-firsty+1
@read_ncdf_varget
        IF jpk EQ 1 THEN res = reform(res, jpi, 1, jpk, /over)
        vmaskred = transpose(temporary(res), [1, 0, 2])
        vmaskred = reform(binary(vmaskred), 8*ny, nx, nz, /over)
        vmaskred = transpose(temporary(vmaskred), [1, 0, 2])
        vmaskred = reform(vmaskred[*, yyy1 MOD 8: 8*ny - 8 + yyy2 MOD 8, *])
     ENDIF ELSE BEGIN
@read_ncdf_varget
        vmaskred = reform(byte(res), /over)
     ENDELSE
  ENDIF ELSE vmaskred = bytarr(jpi, jpk)
; boundary conditions used to compute fmask (2).
  IF ncdf_varid(cdfid, 'fmask') NE -1 THEN BEGIN 
     varinq = ncdf_varinq(cdfid, 'fmask')
     name = varinq.name
     IF micromask NE -1 THEN BEGIN
@read_ncdf_varget
        IF jpk EQ 1 THEN res = reform(res, jpi, 1, jpk, /over)
        fmaskredx = transpose(temporary(res), [1, 0, 2])
        fmaskredx = reform(binary(fmaskredx), 8*ny, nx, nz, /over)
        fmaskredx = transpose(temporary(fmaskredx), [1, 0, 2])
        fmaskredx = reform(fmaskredx[*, yyy1 MOD 8: 8*ny - 8 + yyy2 MOD 8, *])
;
        iyminmesh = iyminmeshtrue
        key_stride = key_stridetrue
     ENDIF ELSE BEGIN
@read_ncdf_varget
        fmaskredx = reform(byte(res), /over)
        fmaskredx = fmaskredx MOD 2
     ENDELSE
  ENDIF ELSE fmaskredx = bytarr(jpi, jpk)
;-------------------------------------------------------
; 1d arrays
;-------------------------------------------------------
  IF (where(allvarname EQ 'e3t_0'))[0] NE -1 THEN fnamevar = ['e3t_0', 'e3w_0', 'gdept_0', 'gdepw_0'] $
  ELSE                                            fnamevar = ['e3t', 'e3w', 'gdept', 'gdepw'] 
  for i = 0, n_elements(fnamevar)-1 do begin
    varinq = ncdf_varinq(cdfid, fnamevar[i])
    CASE n_elements(varinq.dim) OF
      4:ncdf_varget,cdfid,fnamevar[i],tmp,offset = [0,0,izminmesh,0], count = [1,1,jpk,1], stride=[1,1,key_stride[2],1]
      2:ncdf_varget,cdfid,fnamevar[i],tmp,offset = [izminmesh,0], count = [jpk,1], stride=[key_stride[2], 1]
      1:ncdf_varget,cdfid,fnamevar[i],tmp,offset = [izminmesh], count = [jpk], stride=key_stride[2]
    ENDCASE
    if size(tmp, /n_dimension) gt 0 then tmp = reform(tmp, /over)
    CASE fnamevar[i] OF
      'gdept':gdept = float(tmp)
      'gdept_0':gdept = float(tmp)
      'gdepw':gdepw = float(tmp)
      'gdepw_0':gdepw = float(tmp)
      'e3t':e3t = tmp
      'e3t_0':e3t = tmp
      'e3w':e3w = tmp
      'e3w_0':e3w = tmp
    ENDCASE
  ENDFOR
;-------------------------------------------------------
  ncdf_close,  cdfid
;-------------------------------------------------------
; Apply Glamboundary
;-------------------------------------------------------
  if keyword_set(glamboundary) AND key_onearth then BEGIN
    if glamboundary[0] NE glamboundary[1] then BEGIN
      glamt = temporary(glamt) MOD 360
      smaller = where(glamt LT glamboundary[0])
      if smaller[0] NE -1 then glamt[smaller] = glamt[smaller]+360
      bigger = where(glamt GE glamboundary[1])
      if bigger[0] NE -1 then glamt[bigger] = glamt[bigger]-360
      glamu = temporary(glamu) MOD 360
      smaller = where(glamu LT glamboundary[0])
      if smaller[0] NE -1 then glamu[smaller] = glamu[smaller]+360
      bigger = where(glamu GE glamboundary[1])
      if bigger[0] NE -1 then glamu[bigger] = glamu[bigger]-360
      glamv = temporary(glamv) MOD 360
      smaller = where(glamv LT glamboundary[0])
      if smaller[0] NE -1 then glamv[smaller] = glamv[smaller]+360
      bigger = where(glamv GE glamboundary[1])
      if bigger[0] NE -1 then glamv[bigger] = glamv[bigger]-360
      glamf = temporary(glamf) MOD 360
      smaller = where(glamf LT glamboundary[0])
      if smaller[0] NE -1 then glamf[smaller] = glamf[smaller]+360
      bigger = where(glamf GE glamboundary[1])
      if bigger[0] NE -1 then glamf[bigger] = glamf[bigger]-360
      toosmall = where(glamu EQ glamboundary[0])
      IF toosmall[0] NE -1 THEN glamu[toosmall] = glamu[toosmall] + 360
      toosmall = where(glamf EQ glamboundary[0])
      IF toosmall[0] NE -1 THEN glamf[toosmall] = glamf[toosmall] + 360
    endif
  endif
;-------------------------------------------------------
; make sure we do have 2d arrays when jpj eq 1
;-------------------------------------------------------
  IF jpj EQ 1 THEN BEGIN
    glamt = reform(glamt, jpi, jpj, /over)
    gphit = reform(gphit, jpi, jpj, /over)
    e1t = reform(e1t, jpi, jpj, /over)
    e2t = reform(e2t, jpi, jpj, /over)
    glamu = reform(glamu, jpi, jpj, /over)
    gphiu = reform(gphiu, jpi, jpj, /over)
    e1u = reform(e1u, jpi, jpj, /over)
    e2u = reform(e2u, jpi, jpj, /over)
    glamv = reform(glamv, jpi, jpj, /over)
    gphiv = reform(gphiv, jpi, jpj, /over)
    e1v = reform(e1v, jpi, jpj, /over)
    e2v = reform(e2v, jpi, jpj, /over)
    glamf = reform(glamf, jpi, jpj, /over)
    gphif = reform(gphif, jpi, jpj, /over)
    e1f = reform(e1f, jpi, jpj, /over)
    e2f = reform(e2f, jpi, jpj, /over)
    IF keyword_set(key_partialstep) THEN BEGIN
      hdept = reform(hdept, jpi, jpj, /over)
      hdepw = reform(hdepw, jpi, jpj, /over)
      e3t_ps = reform(e3t_ps, jpi, jpj, /over)
      e3w_ps = reform(e3w_ps, jpi, jpj, /over)
    ENDIF
  ENDIF
;-------------------------------------------------------
  ixmindta = ixmindtasauve
  iymindta = iymindtasauve
  izmindta = izmindtasauve
;-------------------------------------------------------
;  widget_control, noticebase, bad_id = nothing, /destroy ;DC
;
;====================================================
; grid parameters used by xxx
;====================================================
;
  IF NOT keyword_set(strcalling) THEN BEGIN
    IF n_elements(ccmeshparameters) EQ 0 THEN strcalling = 'ncdf_meshread' $
    ELSE strcalling = ccmeshparameters.filename
  ENDIF
  IF n_elements(glamt) GE 2 THEN BEGIN
    glaminfo = moment(glamt)
    IF finite(glaminfo[2]) EQ 0 THEN glaminfo = glaminfo[0:1]
    gphiinfo = moment(gphit)
    IF finite(gphiinfo[2]) EQ 0 THEN gphiinfo = gphiinfo[0:1]
  ENDIF ELSE BEGIN
    glaminfo = glamt
    gphiinfo = gphit
  ENDELSE
  romszinfos = {h:-1, zeta:-1, theta_s:-1, theta_b:-1, hc:-1}
  ccmeshparameters = {filename:strcalling  $
          , glaminfo:float(string(glaminfo, format = '(E11.4)')) $
          , gphiinfo:float(string(gphiinfo, format = '(E11.4)')) $
          , jpiglo:jpiglo, jpjglo:jpjglo, jpkglo:jpkglo $
          , jpi:jpi, jpj:jpj, jpk:jpk $
          , ixminmesh:ixminmesh, ixmaxmesh:ixmaxmesh $
          , iyminmesh:iyminmesh, iymaxmesh:iymaxmesh $
          , izminmesh:izminmesh, izmaxmesh:izmaxmesh $
          , key_shift:key_shift, key_periodic:key_periodic $
          , key_stride:key_stride, key_gridtype:key_gridtype $
          , key_yreverse:key_yreverse, key_zreverse:key_zreverse $
          , key_partialstep:key_partialstep, key_onearth:key_onearth}
;
  if keyword_set(key_performance) THEN $
    print, 'time ncdf_meshread', systime(1)-tempsun

;-------------------------------------------------------
   @updateold
;-------------------------------------------------------
   return
 end
