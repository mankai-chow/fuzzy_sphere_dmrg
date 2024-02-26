include("hmt_nn.jl")

function generate_hmt_ops_fuzzy_o3_kondo(nm, nf, int_def, os0)
    os = os0
    no = nm * nf 
    os +=  .5 * int_def, "N", 1, "N", 2 # Sz Sz
    os += -.5 * int_def, "N", 1, "N", 3
    os +=  .5 * int_def, "N", 1, "N", 4
    os += -.5 * int_def, "N", 1, "N", 5
    os +=  .5 * int_def, "Cdag", 1, "Cdag", 3, "C", 2  # S+ S- / 2
    #os +=  .5 * int_def, "Cdag", 1, "Cdag", 5, "C", 4
    #os +=  .5 * int_def, "C", 1, "Cdag", 2, "C", 3  # S- S+ / 2
    #os +=  .5 * int_def, "C", 1, "Cdag", 4, "C", 5
    os +=  .5 * int_def, "N", no + 2, "N", no - 2 # Sz Sz
    os += -.5 * int_def, "N", no + 2, "N", no - 1
    os +=  .5 * int_def, "N", no + 2, "N", no
    os += -.5 * int_def, "N", no + 2, "N", no + 1
    #os +=  .5 * int_def, "Cdag", no + 2, "Cdag", no - 1, "C", no - 2  # S+ S- / 2
    #os +=  .5 * int_def, "Cdag", no + 2, "Cdag", no + 1, "C", no
    #os +=  .5 * int_def, "C", no + 2, "Cdag", no - 2, "C", no - 1  # S- S+ / 2
    #os +=  .5 * int_def, "C", no + 2, "Cdag", no, "C", no + 1
    return os
end

function generate_hmt_mpo_fuzzy_o3_kondo(nm :: Int, nf :: Int ; 
        U0 = 0.0795775, U1 = 0.0218838, U2 = 0.00377993, fld_h = -0.2248, int_def = 1., old = false)
    if (isfile("hmt_n$nm.h5"))
        f = h5open("hmt_n$nm.h5","r")
        hmt = read(f, "hmt", MPO)
        sites = read(f, "sites", Vector{<:Index})
        close(f)
        println("FINISHED READING HAMILTONIAN")
    else
        no = nm * nf
        sites = [ siteind("Spin", conserve_sz = true) ; [ siteind("Fermion", o1 = o, conserve_sz = true, conserve_nf = true, conserve_lz = true) for o :: Int in 1 : no ] ; siteind("Spin", conserve_sz = true)]
        @show sites
        mat_v = [[ 1  0  0  0 ;  0  1  0  0 ;  0  0  1  0 ;  0  0  0  1 ],
                 [ 0  1  0  0 ;  1  0  0  0 ;  0  0  0  1 ;  0  0  1  0 ],
                 [ 0  1  0  0 ; -1  0  0  0 ;  0  0  0  1 ;  0  0 -1  0 ],
                 [ 1  0  0  0 ;  0 -1  0  0 ;  0  0  1  0 ;  0  0  0 -1 ],
                 [ 0  1  0  0 ;  1  0  0  0 ;  0  0  0 -1 ;  0  0 -1  0 ],
                 [ 0  1  0  0 ; -1  0  0  0 ;  0  0  0 -1 ;  0  0  1  0 ],
                 [ 1  0  0  0 ;  0 -1  0  0 ;  0  0 -1  0 ;  0  0  0  1 ]]
        ps_pot = [[U0], [-U1, U2], [-U1, U2], [-U1, U2], [-U1, -U2], [-U1, -U2], [-U1, -U2]]
        fac = [1., 1, -1, 1, 1, -1, 1]
        ps_pot = [ ps_pot[i] * fac[i] for i = 1 : length(ps_pot) ]
        mat_h = [ 0. 0  1  0 ;  0  0  0  1 ;  1  0  0  0 ;  0  1  0  0]
        mat_h *= fld_h

        @time "GENERATE HAMILTONIAN OPERATOR STRING BULK" os = generate_hmt_ops_fuzzy_nn(nm, nf, mat_v, ps_pot, mat_h ; m0 = 1)
        @time "GENERATE HAMILTONIAN OPERATOR STRING DEF"  os = generate_hmt_ops_fuzzy_o3_kondo(nm, nf, int_def, os)
        
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