;  Name: bvf (adapted from bn2)
;  -----
;                     
;  Purpose :
;  --------
;     Compute the local Brunt-Vaisala frequency 
;	
;  Method :
;  -------
;	neos = 0  : UNESCO sea water properties
;	   The brunt-vaisala frequency is computed using the polynomial
;	polynomial expression of McDougall (1987):
;		N^2 = g * beta * ( alpha/beta*dk[ tn ] - dk[ sn ] )/e3w
;
;	neos = 1  : linear equation of state (temperature only)
;		N^2 = g * ralpha * dk[ tn ]/e3w
;
;	neos = 2  : linear equation of state (temperature & salinity)
;		N^2 = g * (ralpha * dk[ tn ] - rbeta * dk[ sn ] ) / e3w
;
;	The use of potential density to compute N^2 introduces e r r o r
;      in the sign of N^2 at great depths. We recommand the use of 
;      neos = 0, except for academical studies.
;
;      N.B. N^2 is set to zero at the first level (JK=1) in inidtr
;      and is never used at this level.
;
;   Input : potential temperature tn, salinity sn
;   ------
;
;   Output :  brunt-vaisala frequency
;   -------
;	       
;   References :
;   -----------
;	McDougall, T. J., J. Phys. Oceanogr., 17, 1950-1964, 1987.
;
;   Modifications :
;   --------------
;      Original  : 94-07 (G. Madec, M. Imbard)
;      Additions : 97-07 (G. Madec) introduction of statement functions
;      Idl version : 98-12 (M. Levy)
;----------------------------------------------------------------------
FUNCTION bvf, tt,  ss
@common



   g = 9.81

   zgde3w = g/e3w

;
;  Interior points ( 2=< jk =< jpkm1 )
;  ----------------
;

            tn = tt*tmask
            sn = ss*tmask
; ... UNESCO seawater equation of state
;   ... temperature, salinity anomalie 
         zzt = 0.5*( tn + shift( tn, 0, 0,  1 ))
         zzs = 0.5*( sn + shift( sn, 0, 0,  1) ) - 35.0
         zzp =  gdepw
;   ... ratio alpha/beta and beta
         zalbet = fsalbt( zzt, zzs, zzp )
         zbeta  = fsbeta( zzt, zzs, zzp )
;   ... N^2 = g/e3w * beta * ( alpha/beta * dk[tn] - dk[sn])
         z = zgde3w * zbeta * tmask $ 
          * ( zalbet * ( shift(tn, 0, 0, 1) - tn ) $ 
              - ( shift(sn, 0, 0, 1) - sn ) )

;
;
;  first and last levels
; ------------------------
; bn^2=0. at first vertical point and jk=jpk 
        z(*, *, 0) = 0.
        z(*, *, jpk-1) = 0.

   return,  z

END

