using LinearAlgebra
using WignerSymbols
using ITensors
using ITensors.HDF5
using ITensorMPOConstruction

BLAS.set_num_threads(1);
NDTensors.Strided.disable_threads();
ITensors.enable_threaded_blocksparse();

nm = parse(Int, ARGS[1])
ns = parse(Int, ARGS[2])
nf = 4
no = nm * nf
nb = div(nm, 2) * nf + 1

include("threej.jl")
include("converge.jl")
include("sweep.jl")
include("site_o3.jl")
include("site_spin.jl")
include("hmt_o3_kondo.jl")
include("init_st_o3_kondo.jl")


hmt, sites = generate_hmt_mpo_fuzzy_o3_kondo(nm, nf)

st0 = generate_init_st_o3_kondo(no, sites)
E0, st0 = sweep_full("0", hmt, st0, nm, nb)

# sites = [ siteind("Spin", conserve_sz = true) ; [ siteind("Fermion", o1 = o, conserve_sz = true, conserve_nf = true, conserve_lz = true) for o :: Int in 1 : no ] ; siteind("Spin", conserve_sz = true)]
# @show  eachindex(sites)
# conf = [ "Up" ; [ (mod(o, 4) == 1 || mod(o, 4) == 2 ? "Occ" : "Emp") for o = 1 : no ] ; "Dn" ]

# st0 = MPS(sites, conf)