function generate_hmt_ops_fuzzy_ising_bd_ex(nm :: Int, nm0 :: Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16)
    os = OpSum()
    int_el = calc_int_ps(nm0, ps_pot_u)
    no = nm * nf
    fld_z = zeros(nm)
    chem_pot = zeros(no)
    for m1 = 1 : nm
        f1 = 0
        o1 = (m1 - 1) * nf + f1 + 1
        for m2 = 1 : nm 
            f2 = 1
            o2 = (m2 - 1) * nf + f2 + 1
            for m3 = 1 : nm
                f3 = 1
                o3 = (m3 - 1) * nf + f3 + 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                f4 = 0
                o4 = (m4 - 1) * nf + f4 + 1
                
                V1 = int_el[m1, m2, m3]
                os += V1 * 2, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
                if (o2 == o3)
                    chem_pot[o1] += -V1
                    chem_pot[o2] += -V1 
                end
            end
        end
        for m2 = nm + 1 : nm0 
            V1 = int_el[m1, m2, m2]
            fld_z[m1] += V1
        end
    end
    for m1 = 1 : nm 
        os += -fld_h, "Cdag", m1 * 2 - 1, "C", m1 * 2
        os += -fld_h, "Cdag", m1 * 2, "C", m1 * 2 - 1
        os +=  fld_z[m1] + chem_pot[m1 * 2 - 1], "N", m1 * 2 - 1
        os += -fld_z[m1] + chem_pot[m1 * 2    ], "N", m1 * 2
    end
    return os
end

function generate_hmt_mpo_fuzzy_ising_bd_ex(nm :: Int, nm0 ::Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16, old = false)
    if (isfile("hmt_n$nm.h5"))
        f = h5open("hmt_n$nm.h5","r")
        hmt = read(f, "hmt", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        sites = [siteind("Fermion", o1 = o, conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_z2 = false) for o :: Int in 1 : no]
        @time "GENERATE HAMILTONIAN OPERATOR STRING" os = generate_hmt_ops_fuzzy_ising_bd_ex(nm, nm0, nf)
        
        if (old)
            @time "GENERATE HAMILTONIAN MPO" hmt = MPO(os, sites)
        else
            operatorNames = [ "I", "C", "Cdag", "N" ]
            opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
            @time "GENERATE HAMILTONIAN MPO" hmt = MPO_new(os, sites; basisOpCacheVec = opCacheVec)
        end

        f = h5open("hmt_n$nm.h5","cw")
        write(f, "hmt", hmt)
        write(f, "sites", sites)
        close(f)
        println("FINISHED GENERATING HAMILTONIAN, BOND DIM $(maxlinkdim(hmt))")
    end 
    return hmt, sites
end