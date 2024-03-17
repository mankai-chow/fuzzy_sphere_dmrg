function generate_hmt_ops_fuzzy_spn(nm :: Int, nf :: Int ; 
    ps_pot_u = [1.0], ps_pot_v = [.9])
    os = OpSum()
    int_el_u = calc_int_ps(nm, ps_pot_u)
    int_el_v = calc_int_ps(nm, ps_pot_v)
    no = nm * nf
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf)
        for o2 = 1 : no # consider only o1 < o2, o3 > o4
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf)
            if (f1 < f2) continue end # f1 >= f2
            if (f1 == f2 && m1 <= m2) continue end 
            for o3 = 1 : no 
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf)
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = f3 : nf - 1  # f3 <= f4
                    if (f3 == f4 && m3 >= m4) continue end
                    o4 = (m4 - 1) * nf + f4 + 1
                    if (!(f1 == f4 && f2 == f3)) 
                        if (!(f1 & ~f2 == 1 && f4 & ~f3 == 1)) continue end
                    end
                    val = 0
                    if (f1 == f2)
                        val = val + (int_el_u[m1, m2, m3] - int_el_u[m2, m1, m3]) * 2.
                    else
                        if (f1 == f4 && f2 == f3) 
                            val = val + (int_el_u[m1, m2, m3]) * 2. 
                        end
                        if (f1 & ~f2 == 1 && f4 & ~f3 == 1) 
                            val = val - int_el_v[m1, m2, m3] 
                        end
                    end
                    os += val, "Cdag", o1, "Cdag", o2, "C", o3, "C", o4
                end 
            end
        end
    end
    return os
end

function generate_hmt_mpo_fuzzy_spn(id, nm :: Int, nf :: Int ; 
    ps_pot_u = [1.0], ps_pot_v = [.65], old = false)
    if (isfile("hmt_$(id)_n$nm.h5"))
        f = h5open("hmt_$(id)_n$nm.h5","r")
        hmt = read(f, "hmt", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        no = nm * nf
        sites = [siteind("Fermion", o1 = o, conserve_sz = true, conserve_nf = true, conserve_lz = true) for o :: Int in 1 : no]
        @time "GENERATE HAMILTONIAN OPERATOR STRING" os = generate_hmt_ops_fuzzy_spn(nm, nf ; ps_pot_u, ps_pot_v)
        
        if (old)
            @time "GENERATE HAMILTONIAN MPO" hmt = MPO(os, sites)
        else
            operatorNames = [ "I", "C", "Cdag", "N" ]
            opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
            @time "GENERATE HAMILTONIAN MPO" hmt = MPO_new(os, sites; basisOpCacheVec = opCacheVec)
        end

        f = h5open("hmt_$(id)_n$nm.h5","cw")
        write(f, "hmt", hmt)
        write(f, "sites", sites)
        close(f)
        println("FINISHED GENERATING HAMILTONIAN, BOND DIM $(maxlinkdim(hmt))")
    end 
    return hmt, sites
end

function generate_s2_mpo_fuzzy_spn(nm :: Int, nf :: Int, ne :: Int, sites ; old = false)
    if (isfile("s2_n$nm.h5"))
        f = h5open("s2_n$nm.h5","r")
        s2 = read(f, "s2", MPO)
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        no = nm * nf
        ps_pot_u :: Vector{Float64} = [iseven(l) ? -.5 : 0 for l = 0 : nm - 1]
        ps_pot_v :: Vector{Float64} = [iseven(l) ? 1.  : 0 for l = 0 : nm - 1]
        @time "GENERATE HAMILTONIAN OPERATOR STRING" os = generate_hmt_ops_fuzzy_spn(nm, nf ; ps_pot_u, ps_pot_v)
        os += .25 * ne * (nf + ne), "I", 1
        
        if (old)
            @time "GENERATE HAMILTONIAN MPO" s2 = MPO(os, sites)
        else
            operatorNames = [ "I", "C", "Cdag", "N" ]
            opCacheVec = [ [OpInfo(ITensors.Op(name, n), sites[n]) for name in operatorNames] for n in eachindex(sites)  ]
            @time "GENERATE HAMILTONIAN MPO" s2 = MPO_new(os, sites; basisOpCacheVec = opCacheVec)
        end

        f = h5open("s2_n$nm.h5","cw")
        write(f, "s2", s2)
        write(f, "sites", sites)
        close(f)
        println("FINISHED GENERATING CASIMIR, BOND DIM $(maxlinkdim(s2))")
    end 
    return s2
end