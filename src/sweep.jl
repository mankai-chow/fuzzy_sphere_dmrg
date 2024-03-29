function sweep_one(id, hmt, st0, nm, nb, dim1 ; 
        cutoff = [1E-9], maxdim = [dim1], nsweeps = 10, noise = [1E-6,1E-7,0], proj = [], e_tol = 1E-6)
    if (isfile("st$(id)_n$(nm).h5"))
        f = h5open("st$(id)_n$(nm).h5","r")
        if (haskey(f, "st_d$(dim1)"))
            st1 = read(f, "st_d$(dim1)", MPS)
            E1 = read(f, "E_d$(dim1)")
            close(f)
            println("FINISHED READING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
            return E1, st1
        end
        close(f)
    end

    if (isempty(proj))
        E1, st1 = dmrg(hmt, st0 ; nsweeps, maxdim, cutoff, noise, observer = EntanglementObserver(nb, e_tol), outputlevel = 1)  # ground state
    else
        fs = [ h5open("st$(fi)_n$(nm).h5", "r") for fi in proj ]
        grs = [ "st_d$(dim1)" for fi in proj ]
        for i = 1 : length(proj)
            if (haskey(fs[i], grs[i])) continue end
            grs[i] = "st_fin"
        end
        sts = [ read(fs[i], grs[i], MPS) for i = 1 : length(proj) ]
        for fi in fs
            close(fi) 
        end 
        # strategy for reading excited states : 
        # first try to read the same bond dimension
        # if not exist, read the final
        E1, st1 = dmrg(hmt, sts, st0 ; nsweeps, maxdim, cutoff, noise, observer = EntanglementObserver(nb, e_tol), outputlevel = 1, weight = nm)  # ground state
    end

    f = h5open("st$(id)_n$(nm).h5","cw")
    write(f, "st_d$(dim1)", st1)
    write(f, "E_d$(dim1)", E1)
    close(f)
    println("FINISHED GENERATING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
    return E1, st1
end

function sweep_full(id, hmt, st00, nm, nb ; 
        dim_list = [1000,2000,3000,4000,5000,6000], proj = [], e_tol = 1E-6, e_tol1 = 1E-7, 
        maxdim0 = [10,20,50,100,200,500], noise0 = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], noise = [1E-6,2E-7,5E-8,1E-8,0], nsweeps = 10)
    if (isfile("st$(id)_n$(nm).h5"))
        f = h5open("st$(id)_n$(nm).h5","r")
        if (haskey(f, "st_fin"))
            st1 = read(f, "st_fin", MPS)
            E1 = read(f, "E_fin")
            close(f)
            println("FINISHED READING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
            return E1, st1
        end
        close(f)
    end
    
    E0, st0 = sweep_one(id, hmt, st00, nm, nb, 0 ; maxdim = maxdim0, nsweeps = length(maxdim0), noise = noise0, e_tol = 0., proj)
    for i = 1 : length(dim_list)
        E1, st1 = sweep_one(id, hmt, st0, nm, nb, dim_list[i] ; proj, e_tol, noise, nsweeps)
        if (abs(E1 - E0) < e_tol1 || maxlinkdim(st1) < .9 * dim_list[i]) break end
        E0 = E1 
        st0 = st1
    end 

    f = h5open("st$(id)_n$(nm).h5","cw")
    write(f, "st_fin", st1)
    write(f, "E_fin", E1)
    close(f)

    println("FINISHED GENERATING STATE $id, FINAL, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
    return E1, st1
end
    