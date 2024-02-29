using SpecialFunctions

function mono_harm_ss(nm :: Int, m1 :: Int, theta, phi) 
    # nm = 2s + 1 ; m1 = 1 -- nm ; phi = 0 
    logfac = logfactorial(nm) - logfactorial(m1 - 1) - logfactorial(nm - m1)
    return exp(logfac * .5) / sqrt(4 * pi) * cos(theta / 2) ^ (m1 - 1) * sin(theta / 2) ^ (nm - m1) * exp((m1 - 1) * im * phi)
end 

function generate_hmt_ops_fuzzy_ising_def_cusp(nm :: Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16, fld_def = [ 1E+4, 1E+4 ], angle_def = [ 0 0 ; pi 0 ])
    os = OpSum()
    int_el = calc_int_ps(nm, ps_pot_u)
    no = nm * nf
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
            end
        end
    end
    for m1 = 1 : nm 
        os += -fld_h, "Cdag", m1 * 2 - 1, "C", m1 * 2
        os += -fld_h, "Cdag", m1 * 2, "C", m1 * 2 - 1
        for m2 = 1 : nm 
            fld1 = sum([ fld_def[i] * conj(mono_harm_ss(nm, m1, angle_def[i, 1], angle_def[i, 2] )) * mono_harm_ss(nm, m2, angle_def[i, 1], angle_def[i, 2] ) for i = 1 : length(fld_def) ])
            if (abs(fld1) < 1E-8) continue end
            os += -fld1, "Cdag", m1 * 2,     "C", m2 * 2
            os +=  fld1, "Cdag", m1 * 2 - 1, "C", m2 * 2 - 1
        end 
    end
    return os
end

function generate_hmt_mpo_fuzzy_ising_def_cusp(id, nm :: Int, nf :: Int ; 
    ps_pot_u = [4.75, 1.0], fld_h = 3.16, old = false, fld_def = [ 1E+3, 1E+3 ], angle_def = [ 0 0 ; pi 0 ])
    if (isfile("hmt_$(id)_n$nm.h5"))
        f = h5open("hmt_$(id)_n$nm.h5","r")
        hmt = read(f, "hmt", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        sites = [siteind("Fermion", o1 = o, conserve_sz = false, conserve_nf = true, conserve_lz = false, conserve_z2 = false) for o :: Int in 1 : no]
        @time "GENERATE HAMILTONIAN OPERATOR STRING" os = generate_hmt_ops_fuzzy_ising_def_cusp(nm, nf ; ps_pot_u, fld_h, fld_def, angle_def)
        
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