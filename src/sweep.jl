function sweep_one(id, hmt, st0, nm, nb, dim1 ; cutoff = [1E-9], maxdim = [dim1], nsweeps = 10, noise = [1E-6,1E-7,0], proj = [], e_tol = 1E-6)
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
        E1, st1 = dmrg(hmt, proj, st0 ; nsweeps, maxdim, cutoff, noise, observer = EntanglementObserver(nb, e_tol), outputlevel = 1, weight = nm)  # ground state
    end

    f = h5open("st$(id)_n$(nm).h5","cw")
    write(f, "st_d$(dim1)", st1)
    write(f, "E_d$(dim1)", E1)
    close(f)
    println("FINISHED GENERATING STATE $id, BOND DIM $(maxlinkdim(st1)), ENERGY $E1")
    return E1, st1
end

function sweep_full(id, hmt, st00, nm, nb ; dim_list = [1000, 1500, 2000, 2500, 3000, 3500, 4000, 4400, 4800, 5200, 5600, 6000], proj = [], e_tol = 1E-6, maxdim0 = [10,20,50,100,200,500], noise0 = [1E-4,3E-5,1E-5,3E-6,1E-6,3E-7], noise = [1E-6,1E-7,0], nsweeps = 10)
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
        if (abs(E1 - E0) < 1E-7 || maxlinkdim(st1) < .9 * dim_list[i]) break end
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
    