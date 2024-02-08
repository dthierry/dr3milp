
using HiGHS

include("test_create_block_v02_06.jl")


pi_ = 1e-03.*ones(size(D_ij)[1])

for i in s.P
    for l in s.L
        @printf "period=\t%i\t" i
        @printf "location=\t%i\n" l
        @time mb = createBlockMod(i, l, p, s)
        # generate the variable block vector
        @time vvb, vpb = genBlockVarVec(mb, p, s, i, l)
        id0 = (p.n_location+1)*i + l
        id1 = id0 * nvblk
        id2 = id1 + nvblk
        #
        D_ij = D[:, id1+1:id2]
        @time attachBlockObjective(mb, p, s, D_ij, vvb, pi_, i, l)
        set_optimizer(mb, HiGHS.Optimizer)
        unset_silent(mb)
        optimize!(mb)
    end
end




