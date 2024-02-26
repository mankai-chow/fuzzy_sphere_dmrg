function ITensors.space( :: SiteType"Fermion" ; o1 :: Int = 0,
    conserve_sz = true, conserve_nf = true, conserve_lz = true,
    qnname_sz = "Sz", qnname_nf = "Nf", qnname_lz = "Lz")
    ## Note that conserve nf and lz is always set on. Lz is TWICE the actual angular momentum
    if (!conserve_nf) 
        return 2 
    end
    m1 = div(o1 - 1, 4) + 1
    f1 = mod(o1 - 1, 4)
    sz1 = iseven(f1) ? 1 : -1
    if (conserve_sz && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0), (qnname_sz,   0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1), (qnname_sz, sz1)) => 1
        ]
    elseif ((!conserve_sz) && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1)) => 1
        ]        
    elseif (conserve_sz && (!conserve_lz))
        return [
            QN((qnname_nf, 0, -1), (qnname_sz,   0)) => 1
            QN((qnname_nf, 1, -1), (qnname_sz, sz1)) => 1
        ]
    else
        return [
            QN((qnname_nf, 0, -1)) => 1
            QN((qnname_nf, 1, -1)) => 1
        ]
    end
end