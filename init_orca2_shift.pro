@cm_4mesh
@cm_4cal
@cm_4data
;---------------------------------------------------------
; grid boundaries
;----------------------------------------------------------
ixminmesh  =1
ixmaxmesh  =180
;
iyminmesh  =0
iymaxmesh  =147
;
izminmesh  =-1
izmaxmesh  =-1
;
jpt = 1
time = 0
;------------------------------------------------------
; read the grid
;------------------------------------------------------
iodir = '/data/daclod/DATA_CMIP5/'
key_stride = [1, 1, 1]
;; ncdf_meshread, '/data/daclod/mesh_mask.nc', glamboundary = [70, 430]
;; ncdf_meshread_dc1, '/data/daclod/mesh_mask.nc', glamboundary = [70, 430]
ncdf_meshread_dc1, '/data/daclod/mesh_mask.nc', glamboundary = [20, 380]
;ncdf_meshread, 'meshmask.orca.2d.nc', glamboundary = [-180, 180]
;ncdf_meshread, 'meshmask_2.nc', glamboundary = [20, 380]
;
;-------------------------------------------------------------
domdef
;
triangles_list = triangule()
;----------------------------------------------------------
; data boundaries
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
