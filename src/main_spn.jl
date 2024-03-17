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
include("site_spn.jl")
include("hmt_spn.jl")

nm = parse(Int, ARGS[1])
nf = parse(Int, ARGS[2])
nb = div(nm, 2) * nf
ps_pot_u = [1.]
ps_pot_v = [parse(Float64, ARGS[3])]

@show nm 
@show nf 

hmt, sites = generate_hmt_mpo_fuzzy_spn("sp$(div(nf,2))_u0$(ps_pot_u[1])_v0$(ps_pot_v[1])", nm, nf ; ps_pot_u, ps_pot_v)
s2 = generate_s2_mpo_fuzzy_spn(nm, nf, nm * 2, sites)

conf = [(mod(o, nf) < 2 ? "Occ" : "Emp") for o = 0 : nm * nf - 1]
st0 = randomMPS(sites, conf ; linkdims = 10) 

Eg, stg = sweep_full("g", hmt, st0, nm, nb ; dim_list = [1000, 1500, 2000])
s2ip = inner(stg', s2, stg)
@show s2ip
Ef0, stf0 = sweep_full("f0", hmt, st0, nm, nb ; proj = [ "g" ], dim_list = [1000, 1500, 2000])
s2ip = inner(stf0', s2, stf0)
@show s2ip

conf[1] = "Emp"
conf[3] = "Occ"
st1 = randomMPS(sites, conf ; linkdims = 10) 
Ef1, stf1 = sweep_full("f1", hmt, st1, nm, nb ; dim_list = [1000, 1500, 2000])
s2ip = inner(stf1', s2, stf1)
@show s2ip

conf[nf + 1] = "Emp"
conf[nf + 3] = "Occ"
st2 = randomMPS(sites, conf ; linkdims = 10) 
ET, stT = sweep_full("T", hmt, st2, nm, nb ; dim_list = [1000, 2000])
s2ip = inner(stT', s2, stT)
@show s2ip