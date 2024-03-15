function ITensors.space( :: SiteType"Fermion" ; o1 :: Int = 0,
    conserve_sz = true, conserve_nf = true, conserve_lz = true,
    qnname_sz = "Sz", qnname_nf = "Nf", qnname_lz = "Lz", fac_sz = [ 10 ^ f for f :: Int = 0 : nf / 2 - 1])
    m1 = div(o1 - 1, nf) + 1
    f1 = mod(o1 - 1, nf)
    sz1 :: Int = (f1 < nf / 2 ? 1 : -1) * fac_sz[Int(1 + mod(f1, nf / 2))]
    if (conserve_sz && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0), (qnname_sz,      0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1), (qnname_sz, sz1[1])) => 1
        ]
    elseif ((!conserve_sz) && conserve_lz)
        return [
            QN((qnname_nf, 0, -1), (qnname_lz,  0)) => 1
            QN((qnname_nf, 1, -1), (qnname_lz, m1)) => 1
        ]        
    elseif (conserve_sz && (!conserve_lz))
        return [
            QN((qnname_nf, 0, -1), (qnname_sz,      0)) => 1
            QN((qnname_nf, 1, -1), (qnname_sz, sz1[1])) => 1
        ]
    else
        return [
            QN((qnname_nf, 0, -1)) => 1
            QN((qnname_nf, 1, -1)) => 1
        ]
    end
end