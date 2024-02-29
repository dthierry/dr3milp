import MathOptInterface as MOI
using Gurobi

include("./create_data_v02_26.jl")
# create_dummy_data_v02_26.jl
include("../mod/sets_v01_30.jl")
include("../mod/bm_v02_26.jl")
include("../mod/compMat_v02_28.jl")

s = sets(p.n_periods, p.n_years, p.n_location, p.n_rtft, p.n_new, p.n_fu)

m = createBlockMod(s.P, s.L, p, s)

mv =  create_vector_mods(p, s)

D, _, _, nvb, rhs, consense = complicatingMatrix(m, p, s)

pi_ = zeros(size(D)[1])

blockObjAttach!(mv, p, s, D, pi_)

for m in mv
    set_optimizer(m, Gurobi.Optimizer)
    optimize!(m)
end

xk = primalValVec(mv, p, s)

#  
rmp = Model()

ncomp_con_i = size(consense[consense.>0])[1]
ncomp_con_e = size(consense[consense.==0])[1]
ncomp_con = ncomp_con_i + ncomp_con_e

K = length(mv)

# lambda is K by ncolumns
ncolumns = 1

lambda = Matrix{VariableRef}(undef, K, ncolumns)

for col in 1:ncolumns
    for k in 1:K
        lambda[k, col] = @variable(rmp, base_name="lm_[$k,$col]", 
                                   lower_bound=0)
    end
end
rmp[:lambda] = lambda


# each column has ncomp_con_x rows
DixG = Array{Float64, 3}(undef, ncomp_con_i, K, ncolumns)
DexG = Array{Float64, 3}(undef, ncomp_con_e, K, ncolumns)

rhs_i = rhs[consense.>0]
rhs_e = rhs[consense.==0]

for col in 1:ncolumns
    k = 1
    for i_ in s.P
        for l_ in s.L
            c0 = nvb * (k-1) + 1
            c1 = nvb * k
            di_ = D[consense.>0, c0:c1]
            de_ = D[consense.==0, c0:c1]
            x = xk[1:nvb, k, col]  # this guy is a matrix per column
            # calculate column (vector)
            DixG[:, k, col] = di_*x  # ncomp_con_i by K
            DexG[:, k, col] = de_*x  # ncomp_con_e by K
            k += 1
        end
    end
end
 
# DixG
@constraint(rmp, comp_ineq, 
            sum(DixG[:, :, c]*lambda[:,c] for c in 1:ncolumns) - rhs_i
            in MOI.Nonnegatives(ncomp_con_i))
# DexG
@constraint(rmp, comp_eq, 
            sum(DexG[:, :, c]*lambda[:,c] for c in 1:ncolumns) - rhs_e
            in MOI.Zeros(ncomp_con_e))

