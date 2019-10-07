;;;  Name : rhoeos
;;;  ----
;;;  Purpose :
;;;  --------
;;;	compute the in situ density  and the potential
;;;	volumic mass (Kg/m3) from now potential temperature and salinity
;;;	fields  using an equation of state defined
;;;	through the input parameter neos.
;;
;;   Method :
;;   -------
;;	neos = 0 : Jackett and McDougall (1994) equation of state.
;;         the now in situ density is computed directly as a function of
;;         potential temperature relative to the surface (the opa tn
;;         variable), salt and pressure (assuming no pressure variation
;;         along geopotential surfaces, i.e. the pressure p in decibars
;;         is approximated by the depth in meters.
;;              rho(t,s,p) 
;;              rhop(t,s)  = rho(t,s,0)
;;         with pressure                      p        decibars
;;              potential temperature         t        deg celsius
;;              salinity                      s        psu
;;              reference volumic mass        rau0     kg/m**3
;;              in situ volumic mass          rho      kg/m**3
;;
;;         Check value: rho = 1059.8204 kg/m**3 for p=10000 dbar,
;;          t = 40 deg celcius, s=40 psu
;;
;;      neos = 1 : linear equation of state function of temperature only
;;              rho(t) =  rau0*(0.028 - ralpha * tn)+rau0
;;              rhop(t,s)  = rho(t,s)
;;
;;      neos = 2 : linear equation of state function of temperature and
;;		   salinity
;;		rho(t,s) = (rbeta * sn - ralpha * tn - 1.)*rau0+rau0
;;              rhop(t,s)  = rho(t,s)
;;
;;      Note that no boundary condition problem occurs in this routine
;;      as tn and sn are defined over the whole domain.
;;
;;   Input : potential temperature tn, salinity sn
;;   ------
;;
;;   Output : potential density rho (default), or insitu  density
;;   -------  rhop if type is set to 'situ'
;;   
;;
;;   Referen;es :
;;   -----------
;;      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
;;
;;   Modifications :
;;   --------------
;;      original :  89-03 (o. Marti)
;;      additions : 94-08 (G. Madec)
;;                : 96-01 (G. Madec) statement function for e3
;;		  : 97-07 (G. Madec) introduction of neos, OPA8.1
;;		  : 97-07 (G. Madec) density instead of volumic mass
;;                : 98-12 (M. Levy) idl version
;; attention: tt et ss doivent etre maskés, rau doit etre maske en sortie
;;------------------------------------------------------------------------
;                     
;----------------------------------------------------------------------
FUNCTION rhoeos, tn, sn, rau0
@common


; 1) Jackett and McDougall (1994) formulation
;---------------------------------------------
  ralpha = 2.e-4
  rbeta = 7.7e-4
  rhop = rau0* (rbeta*sn - ralpha*tn)
  rho = rhop
  sigma = 1

        zrau0r = 1 / rau0;
;   ... depth
            zh = 1000.* sigma
;   ... square root salinity
            zsr= sqrt( abs( sn ) )
;   ... compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332e-9*tn-1.120083e-6 )*tn+1.001685e-4)*tn $  
                        -9.095290e-3 )*tn+6.793952e-2 )*tn+999.842594
;   ... seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875e-9*tn-8.2467e-7 ) *tn+7.6438e-5 ) *tn $  
                      -4.0899e-3 ) *tn+0.824493
            zr3= ( -1.6546e-6*tn+1.0227e-4 ) *tn-5.72466e-3
            zr4= 4.8314e-4
;   ... potential volumic mass (reference to the surface)
            rhop= (( zr4*sn + zr3*zsr + zr2 ) *sn + zr1)* tmask
;   ... add the compression terms
            ze = ( -3.508914e-8*tn-1.248266e-8 ) *tn-2.595994e-6
            zbw= (  1.296821e-6*tn-5.782165e-9 ) *tn+1.045941e-4
            zb = zbw + ze * sn
            zd = -2.042967e-2
            zc =   (-7.267926e-5*tn+2.598241e-3 ) *tn+0.1571896
            zaw= ( ( 5.939910e-6*tn+2.512549e-3 ) *tn-0.1028859 ) *tn  -4.721788
            za = ( zd*zsr + zc ) *sn + zaw
            zb1=   (-0.1909078*tn+7.390729 ) *tn-55.87545
            za1= ( ( 2.326469e-3*tn+1.553190)*tn-65.00517 ) *tn+1044.077
            zkw= ( ( (-1.361629e-4*tn-1.852732e-2 ) *tn-30.41638 ) *tn $  
                      +2098.925 ) *tn+190925.6
            zk0= ( zb1*zsr + za1 )*sn + zkw
;
;   ... masked in situ density
;            rho =  rhop  / (1.0 - zh/(zk0- zh *(za- zh*zb))) 
;
;   ... masked in situ density anomaly            
            prd = (rhop / (1.0- zh / ( zk0 - zh * ( za-zh*zb)))- rau0) * zrau0r * tmask           

return, prd
      

END

