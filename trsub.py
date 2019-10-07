def trsubwadv(wmld, trmld) : 
    """
    Subduction by vertical advection
    """
    surf = e1t * e2t
    subw = -wmld * trmld * surf
    return subw
#

def trsubhadv(mld, umld, vmld, trmldu, trmldv, umask_i,  vmask_i,  nmlnu,  nmlnv) : 
    """
    Subduction by horizontal advection across tilted MLD
    suadv at point mldu
    svadv at point mldv
    """

    suadv = fltarr(jpj, jpi)
    svadv = fltarr(jpj, jpi)

    for ii in range(jpi-1) :
        for jj in range(jpj) : 
            suadv[jj, ii] = -umld[jj, ii] * trmldu[jj, ii]* e2u[jj, ii] * ( mld[jj, ii+1] - mld[jj, ii] ) * umask_i[jj, ii, nmlnu[jj, ii]]
        #
    #

    for ii in range(jpi) :
        for jj in range(jpj-1) : 
            svadv[jj, ii] = -vmld[jj, ii] * trmldv[jj, ii] * e1v[jj, ii] * ( mld[jj+1, ii] - mld[jj, ii]) * vmask_i[jj, ii, nmlnv[jj, ii]]
        #
    #

    sadv = suadv + svadv

    return sadv
#

def trsubmld(mld, mldm1, mldp1, tr, trm1) : 
    """
    Subduction by entrainment/detrainment due to variation of the MLD
    """

    # computing d(MLD)/dt and center it in time
    # la subduction est égale à -( trcmld2 - trcmld1 ) = -tr(n) * ( mld(n+1) - mld(n) )
    trcmld1  = fltarr[jpj, jpi]       # =mld(n)*tr(n)
    trcmld2  = fltarr[jpj, jpi]       # =mld(n+1)*tr(n)
    mld1     = (mld + mldm1) * .5     # 1er du pas de temps en cours (mld en time - 1/2)
    mld2     = (mld + mldp1) * .5     # 1er du pas de temps suivant  (mld en time + 1/2)
    trcentre = (tr + trm1) * .5       # tr en time - 1/2

    # 1-2-1contenu de la ml: prise en compte de tous les niveaux
    # supérieurs+prorata du niveau dans lequel se trouve la ml
    kmldinf1 = fltarr[jpj, jpi]
    kmldinf2 = fltarr[jpj, jpi]
    for ii in range(jpi) :
        for jj in range(0, jpj) :
            aa = np.where(gdepw <= mld1[jj, ii])
            kmldinf1[jj, ii] = max(a)
            for  kk in range(kmldinf1[jj, ii]-1) :
                trcmld1[jj, ii] = trcmld1[jj, ii] + e1t[jj, ii] * e2t[jj, ii] * e3t * trcentre[kk, jj, ii] * tmask[kk, jj, ii]
            #
            b = where (gdepw(i, j, *) LE mld2[jj, ii])
            kmldinf2[jj, ii] = max(b)
            for  kk in range(kmldinf2[jj, ii]-1) :
                trcmld2[jj, ii] = trcmld2[jj, ii] + e1t[jj, ii] * e2t[jj, ii] * e3t * trcentre[kk, jj, ii] * tmask[kk, jj, ii]
            #
            trcmld1[jj, ii] = (   trcmld1[jj, ii] + e1t[jj, ii] * e2t[jj, ii] * ( mld1[jj, ii] - gdepw[kmldinf1[jj, ii]] ) 
                                  * trcentre[kmldinf1[jj, ii], jj, ii] * tmask[kmldinf1[jj, ii], jj, ii]   )
            trcmld2[jj, ii] = (   trcmld2[jj, ii] + e1t[jj, ii] * e2t[jj, ii] * ( mld2[jj, ii] - gdepw[kmldinf2[jj, ii]] )
                                  * trcentre[kmldinf2[jj, ii], jj, ii] * tmask[kmldinf2[jj, ii], jj, ii]   )        
        #
    #

    zz = trcmld1 - trcmld2
    return, zz

  



#
