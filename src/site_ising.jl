function ITensors.space( :: SiteType"Fermion" ; o1 :: Int = 0,
    conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_z2 = true,
    qnname_sz = "Sz", qnname_nf = "Nf", qnname_lz = "Lz", qnname_z2 = "Z2")
    ## Note that conserve nf and lz is always set on. Lz is TWICE the actual angular momentum
    if (conserve_sz || (!conserve_nf)) 
        return 2 
    end
    m1 = div(o1 - 1, 2) + 1
    f1 = mod(o1 - 1, 2)
    if (conserve_z2 && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0), (qnname_z2,  0, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1), (qnname_z2, f1, 2)) => 1
        ]
    elseif ((!conserve_z2) && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1)) => 1
        ]        
    elseif (conserve_z2 && (!conserve_lz))
        return [
            QN((qnname_nf, 0, -1), (qnname_z2,  0, 2)) => 1
            QN((qnname_nf, 1, -1), (qnname_z2, f1, 2)) => 1
        ]
    else
        return [
            QN((qnname_nf, 0, -1)) => 1
            QN((qnname_nf, 1, -1)) => 1
        ]
    end
end