
using HiGHS

include("test_create_mod_v02_01.jl")



attachPeriodBlock(m, p, s)
attachLocationBlock(m, p, s)



set_optimizer(m, HiGHS.Optimizer)
unset_silent(m)
optimize!(m)


