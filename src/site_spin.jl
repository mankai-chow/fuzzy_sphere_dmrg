function ITensors.space( :: SiteType"Spin" ; 
    conserve_sz = true, conserve_nf = true, conserve_lz = true,
    qnname_sz = "Sz", qnname_nf = "Nf", qnname_lz = "Lz")
    if (conserve_sz)
        return [ QN((qnname_nf, 0, -1), (qnname_lz, 0), (qnname_sz, ns + 1 - m * 2)) => 1 for m = 1 : ns ]
    else 
        return ns
    end
end

ITensors.state(::StateName"Up", ::SiteType"Spin") = [ i == 1 ? 1. : 0. for i = 1 : ns]
ITensors.state(::StateName"Dn", ::SiteType"Spin") = [ i == ns ? 1. : 0. for i = 1 : ns]

ITensors.op(::OpName"I",::SiteType"Spin") = [ m1 == m2 ? 1. : 0 for m1 = 1 : ns, m2 = 1 : ns ]
ITensors.op(::OpName"N",::SiteType"Spin") = [ m1 == m2 ? .5 * (ns + 1) - m1 : 0 for m1 = 1 : ns, m2 = 1 : ns ]
ITensors.op(::OpName"Cdag",::SiteType"Spin") = [ m1 == m2 - 1 ? sqrt((ns - m1) * m1) : 0 for m1 = 1 : ns, m2 = 1 : ns ]
ITensors.op(::OpName"C",::SiteType"Spin") = [ m2 == m1 - 1 ? sqrt((ns - m2) * m2) : 0 for m1 = 1 : ns, m2 = 1 : ns ]
# I define these names so that they are the same as fermions, need to revisit later