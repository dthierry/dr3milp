
include("../mod/create_dummy_data.jl")
include("../mod/sets_v10_4.jl")
include("../mod/bm_v1213.jl")

s = sets(p.t_horizon, p.n_location, p.n_rtft, p.n_new, p.n_fu)

m = createBlockMod(s.L, p, s)



