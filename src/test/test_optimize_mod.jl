
using HiGHS

include("test_create_mod.jl")

# /Users/dthierry/Projects/dr3milp/src/test/test_create_mod.jl
#
@info "HiGHS attached"

reattachBlockMod!(m, s.L, p, s)

# unset_binary.(m[:y_o][:, :])
# unset_binary.(m[:y_e][:, :])
# 
# unset_binary.(m[:y_r][:, :, :])
# unset_binary.(m[:y_n][:, :, :])
#  
# unset_binary.(m[:e_yps][:, :])
# unset_binary.(m[:r_yps][:, :])
# unset_binary.(m[:n_yps][:, :])


set_optimizer(m, HiGHS.Optimizer)
unset_silent(m)
optimize!(m)


