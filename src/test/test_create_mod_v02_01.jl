
include("../mod/create_dummy_data_v02_01.jl")
include("../mod/sets_v01_30.jl")
include("../mod/bm_v02_02.jl")

s = sets(p.n_periods, p.n_years, p.n_location, p.n_rtft, p.n_new, p.n_fu)

m = createBlockMod(s.P, s.L, p, s)
D, vv, vp = complicatingMatrix(m, p, s)


