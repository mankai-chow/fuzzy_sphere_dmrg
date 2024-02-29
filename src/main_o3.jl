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
include("site_o3.jl")
include("hmt_o3.jl")
include("init_st_o3.jl")

nm = parse(Int, ARGS[1])
nf = 4
no = nm * nf
nb = div(nm, 2) * nf

hmt, sites = generate_hmt_mpo_fuzzy_o3(nm, nf)

st0 = generate_init_st_o3(no, sites)

E0, st0 = sweep_full("0", hmt, st0, nm, nb)
E1, st1 = sweep_full("1", hmt, st0, nm, nb ; proj = ["0"])
E2, st2 = sweep_full("2", hmt, st0, nm, nb ; proj = ["0", "1"])
