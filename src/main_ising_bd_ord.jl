using LinearAlgebra
using WignerSymbols
using ITensors
using ITensors.HDF5
using ITensorMPOConstruction

BLAS.set_num_threads(1);
NDTensors.Strided.disable_threads();
ITensors.enable_threaded_blocksparse();

include("threej.jl")
include("converge.jl")
include("sweep.jl")
include("site_ising.jl")
include("hmt_ising_bd_ord.jl")
include("init_st_ising.jl")

nm = parse(Int, ARGS[1])
nm0 = 2 * nm
nf = 2
no = nm * nf
nb = div(nm, 2) * nf

hmt, sites = generate_hmt_mpo_fuzzy_ising_bd_ord(nm, nm0, nf)

st0 = generate_init_st_ising(no, sites)
st1 = generate_init_st_ising(no, sites ; z2 = -1)

Eg, stg = sweep_full("g", hmt, st0, nm, nb)
Ed, std = sweep_full("d", hmt, st0, nm, nb ; proj = ["g"])
Eo, sto = sweep_full("o", hmt, st1, nm, nb)

if (!isfile("overlap_n$nm.dat"))
    inner_product = zeros(3, nm)
    for m = 1 : nm
        osz = OpSum()
        osz +=     "N", 2 * m - 1
        osz += -1, "N", 2 * m
        nmz = MPO(osz, sites)
        inner_product[1, m] = inner(stg', nmz, stg)
        inner_product[2, m] = inner(std', nmz, stg)

        osx = OpSum()
        osx += "Cdag", 2 * m - 1, "C", 2 * m 
        osx += "Cdag", 2 * m,     "C", 2 * m - 1
        nmx = MPO(osx, sites)
        inner_product[3, m] = inner(sto', nmx, stg)
    end 
    f = open("overlap_n$nm.dat", "a")
    for m = 1 : nm
        write(f, string(nm), "\t", string(m), "\t", 
            string(inner_product[1, m]), "\t", 
            string(inner_product[2, m]), "\t", 
            string(inner_product[3, m]), "\n")
    end
    close(f)
end

# Eo1, sto1 = sweep_full("o1", hmt, generate_init_st_ising(no, sites ; z2 = -1, lz = 1), nm, nb)
