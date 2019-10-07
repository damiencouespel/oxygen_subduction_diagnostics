; Thu 24 april 2014 : DC
;Thu 21 Mar 2014 : DC

@cm_4mesh
@cm_4data
;----------------------------------------------------------
; boundaries of the read grid regarding to the original grid
;----------------------------------------------------------
ixminmesh  =-1
ixmaxmesh  =-1
;
iyminmesh  =-1
iymaxmesh  =-1
;
izminmesh  =-1
izmaxmesh  =-1
;------------------------------------------------------
; read the grid
;------------------------------------------------------
key_stride = [1, 1, 1]
iodir = isadirectory('/data/daclod/DATA_CMIP5/', title = 'Select the default IO directory') ;DC
;ncdf_meshread_orca2, 'mesh_mask.nc', glamboundary = [20, 380]
;ncdf_meshread, 'mesh_mask.nc', glamboundary = [20, 380]
;ncdf_meshread, '/data/daclod/mesh_mask.nc'
ncdf_meshread_dc1, '/data/daclod/mesh_mask.nc'

;-------------------------------------------------------------
domdef
triangles_list = triangule()
;----------------------------------------------------------
; boundaries of the data regarding to the original grid
;----------------------------------------------------------
jpidta = jpiglo
jpjdta = jpjglo
jpkdta = jpkglo
ixmindta = 0
ixmaxdta = jpidta-1
iymindta = 0
iymaxdta = jpjdta-1
izmindta = 0
izmaxdta = jpkdta-1
;----------------------------------------------------------
@updateold
;ONLINE_HELP,book="/Users/marina/WORK/SOURCES/SAXO_DIR/SRC/Documentation/idldoc_assistant_output/idldoc-lib.adp"


