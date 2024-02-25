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
include("hmt_ising.jl")
include("init_st_ising.jl")

nm = parse(Int, ARGS[1])
nf = 2
no = nm * nf
nb = div(nm, 2) * nf

hmt, sites = generate_hmt_mpo_fuzzy_ising(nm, nf)

st0 = generate_init_st_ising(no, sites)
st1 = generate_init_st_ising(no, sites ; z2 = -1)

Eg, stg = sweep_full("g", hmt, st0, nm, nb)
Ee, ste = sweep_full("e", hmt, st0, nm, nb ; proj = [stg])
Es, sts = sweep_full("s", hmt, st1, nm, nb)
