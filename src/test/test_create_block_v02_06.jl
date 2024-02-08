
include("../mod/create_dummy_data_v02_01.jl")
include("../mod/sets_v01_30.jl")
include("../mod/bm_v02_02.jl")

s = sets(p.n_periods, p.n_years, p.n_location, p.n_rtft, p.n_new, p.n_fu)

# create full model
@printf "Generating full model\n"
@time m = createBlockMod(s.P, s.L, p, s)
# D, var_vect, var_ptr, nvars_block, rhs
@printf "generating matrices"
@time D, vv, vp, nvblk, rhs = complicatingMatrix(m, p, s)

period = rand(s.P)
location = rand(s.L)

@printf "period=\t%i\t" period
@printf "location=\t%i\n" location

@printf "Generating block model\n"
@time mb = createBlockMod(period, location, p, s)
# generate the variable block vector
@time vvb, vpb = genBlockVarVec(mb, p, s, period, location)

# the respective index is (p.n_location + 1)*period + location
# we need a mechanism to compute this?
id0 = (p.n_location+1)*period + location
id1 = id0 * nvblk
id2 = id1 + nvblk
#
D_ij = D[:, id1+1:id2]
pi_ = ones(size(D_ij)[1])
@time attachBlockObjective(mb, p, s, D_ij, vvb, pi_, period, location)
#
@printf "Done\n"
