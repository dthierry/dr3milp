


include("../ins/create_data_v02_12.jl")
include("../mod/sets_v01_30.jl")
include("../mod/bm_v02_20.jl")

s = sets(p.n_periods, p.n_years, p.n_location, p.n_rtft, p.n_new, p.n_fu)

m = createBlockMod(s.P, s.L, p, s)

