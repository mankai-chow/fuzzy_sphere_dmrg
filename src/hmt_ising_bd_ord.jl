function generate_hmt_ops_fuzzy_ising_bd_ord(nm :: Int, nm0 :: Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16)
    os = OpSum()
    int_el = calc_int_ps(nm0, ps_pot_u)
    no = nm * nf
    chem_pot = zeros(no)
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf)
        for o2 = o1 + 1 : no # consider only o1 < o2, o3 > o4
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf)
            for o3 = 1 : no 
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf)
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                f4 = mod(f1 + f2 + f3, 2)
                o4 = (m4 - 1) * nf + f4 + 1
                if (o3 <= o4) continue end
                
                V1 = int_el[m1, m2, m3] * (f1 == f4 ? 1 : -1) -int_el[m2, m1, m3] * (f2 == f4 ? 1 : -1)
                os += V1, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
                if (o2 == o3)
                    chem_pot[o1] += -V1 / 2
                    chem_pot[o2] += -V1 / 2
                end
            end
        end
        chem_pot[o1] += mod(o1, 2) == 0 ? fld_h : -fld_h
    end
    for o1 = 1 : no 
        os += chem_pot[o1], "N", o1
    end
    return os
end

function generate_hmt_mpo_fuzzy_ising_bd_ord(nm :: Int, nm0 ::Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16, old = false)
    if (isfile("hmt_n$nm.h5"))
        f = h5open("hmt_n$nm.h5","r")
        hmt = read(f, "hmt", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        sites = [siteind("Fermion", o1 = o, conserve_sz = false, conserve_nf = true, conserve_lz = true, conserve_z2 = true) for o :: Int in 1 : no]
        @time "GENERATE HAMILTONIAN OPERATOR STRING" os = generate_hmt_ops_fuzzy_ising_bd_ord(nm, nm0, nf)
        
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