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
include("hmt_ising_bd_ex.jl")
include("init_st_ising.jl")

nm = parse(Int, ARGS[1])
nm0 = 2 * nm
nf = 2
no = nm * nf
nb = div(nm, 2) * nf

hmt, sites = generate_hmt_mpo_fuzzy_ising_bd_ex(nm, nm0, nf)

st0 = generate_init_st_ising(no, sites)

Eg, stg = sweep_full("g", hmt, st0, nm, nb)
Ed, std = sweep_full("d", hmt, st0, nm, nb ; proj = [stg])
Ed1, std1 = sweep_full("d1", hmt, generate_init_st_ising(no, sites ; lz = 1), nm, nb)
