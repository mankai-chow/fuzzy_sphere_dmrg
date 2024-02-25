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
Ed, std = sweep_full("d", hmt, st0, nm, nb ; proj = [stg])
Eo, sto = sweep_full("o", hmt, st1, nm, nb)
Eo1, sto1 = sweep_full("o1", hmt, generate_init_st_ising(no, sites ; z2 = -1, lz = 1), nm, nb)
