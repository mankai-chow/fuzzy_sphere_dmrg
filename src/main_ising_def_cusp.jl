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
include("hmt_ising_def_cusp.jl")
include("init_st_ising.jl")

nm = parse(Int, ARGS[1])
nf = 2
no = nm * nf
nb = div(nm, 2) * nf
theta = parse(Float64, ARGS[2])

hmt, sites = generate_hmt_mpo_fuzzy_ising_def_cusp("th$theta", nm, nf ; angle_def = [ 0 0 ; theta*pi 0 ])
st0 = generate_init_st_ising(no, sites)
Eg, stg = sweep_full("g_th$theta", hmt, st0, nm, nb)