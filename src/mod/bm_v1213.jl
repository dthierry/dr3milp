 # 7n8
using JuMP
using Printf


include("sets_v10_4.jl")
include("params_v10_3.jl")
#include("info_i_v10_3.jl")

# 80 
# 80 ###########################################################################
# defining the model object
#
function createBlockMod(index_l, p::parms, s::sets)
    # 76 
    # 76 #######################################################################
    @info "Generating the damn block"
    #m = Model(Cbc.Optimizer)
    m = Model()
    set_silent(m)
    if index_l isa Int
        L = [index_l]
    else
        @printf "Set passed as a collection\n"
        L = index_l
    end

    T = s.T
    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu
    
    # 76 
    # 76 #######################################################################
    # side-effects from the existing plant:
    # capacity factor
    # heating requirement
    # electricity requirement
    # intrinsic emissions
    # exstrinsic emissions
    # scope 1, emitted
    # scope 1, captured
    # scope 1, stored
    # o&m
    #
    # states: 
    # loan
    # triggered:
    # annuity
    # additional capital (ladd)
    ##
    # variables
    ##
    # tier 0: online status
    @variable(m, y_o[t=T, l=L], Bin)  # online
    # tier 1: retrofit variable
    @variable(m, y_r[t=T, l=L, k=Kr], Bin)
    # tier 2: retirement
    # there is no retirement binary variable per se.  
    

    # d082923
    # expansion
    @variable(m, y_e[t=T, l=L], Bin)  # 1 if yes expand else no
    @variable(m, e_c_d_[t=T, l=L, (0,1)] >= 0)
    @variable(m, e_c[t=T, l=L])

    @variable(m, e_l[t=T, l=L])
    @variable(m, e_l_d_[t=T, l=L, k=(0,1)]) # disaggregated

    @variable(m, 0 <= x[l=L] <= 10) # this should be int
    @variable(m, x_d_[t=T, l=L, (0, 1)] >= 0)
    
    # d082223
    # capacity
    @variable(m, r_cp[t=T, l=L])
    @variable(m, r_cp_d_[t=T, l=L, k=Kr])  # retrofit capacity disagg
    @variable(m, cpb[t=T, l=L] >= 0.0)  # base capacity
    @variable(m, r_cpb_d_[t=T, l=L, k=Kr] >= 0.0)  # base capacity disagg
    
    # d082923
    # heating requirement
    @variable(m, r_eh[t=T, l=L]) # 0
    @variable(m, r_eh_d_[t=T, l=L, k=Kr]) # 1
    # fuel requirement 
    @variable(m, r_ehf[t=T, l=L, f=Fu]) # 2
    @variable(m, r_ehf_d_[t=T, l=L, k=Kr, f=Fu]) # 3
    # electricity requirement
    @variable(m, r_u[t=T, l=L]) # 4
    @variable(m, r_u_d_[t=T, l=L, k=Kr]) # 5
    # process (intrinsic) emissions
    @variable(m, r_cp_e[t=T, l=L]) # 6
    @variable(m, r_cp_e_d_[t=T, l=L, k=Kr]) # 7
    # process (total) disaggregated emissions, e.g. process + fuel, etc.
    @variable(m, r_ep0[t=T, l=L]) # 8
    @variable(m, r_ep0_d_[t=T, l=L, k=Kr]) # 9
    # scope 1, emitted
    @variable(m, r_ep1ge[t=T, l=L]) # 10
    @variable(m, r_ep1ge_d_[t=T, l=L, k=Kr]) # 11
    # scope 1, captured
    @variable(m, r_ep1gce[t=T, l=L]) # 12
    @variable(m, r_ep1gce_d_[t=T, l=L, k=Kr]) # 13
    # scope 1,  stored
    @variable(m, r_ep1gcs[t=T, l=L]) # 14
    @variable(m, r_ep1gcs_d_[t=T, l=L, k=Kr]) # 15
    # operating and maintenance
    @variable(m, r_conm_d_[t=T, l=L, k=Kr])  # 16
    @variable(m, r_conm[t=T, l=L])  # 17

    # d083023
    # loan state
    @variable(m, r_loan[t=T, l=L])
    @variable(m, r_loan_p[t=T, l=L] >= 0)
    @variable(m, r_loan_n[t=T, l=L] >= 0)
    # payment
    @variable(m, r_pay[t=T, l=L])
    @variable(m, r_pay_1[t=T, l=L] >= 0)
    # add cost 
    @variable(m, r_ladd[t=T, l=L])
    @variable(m, r_ladd_d_[t=T, l=L, (0,1)] >= 0)

    @variable(m, r_yps[t=T, l=L], Bin)  # paid or not
    #
    @variable(m, r_l[t=T, l=L]) # capital (loan)
    @variable(m, r_l_md_[t=T, l=L, k=Kr]) # 
    @variable(m, r_l_pd_[t=T, l=L, (0,1)])  # payd or not paid
    #
    @variable(m, r_ann[t=T, l=L])
    @variable(m, r_ann_0[t=T, l=L] >= 0)
    @variable(m, r_ann_1[t=T, l=L] >= 0)
    @variable(m, r_ann_md_[t=T, l=L, k=Kr])
    
    
    # expansion capacity
    @variable(m, e_ladd[t=T, l=L])
    @variable(m, e_ladd_d_[t=T, l=L, (0,)] >= 0)
    @variable(m, e_l_pd_[t=T, l=L, (0,1)])
    @variable(m, e_yps[t=T, l=L], Bin) # 1 if paid, 0 otw

    @variable(m, e_loan[t=T, l=L])
    @variable(m, e_loan_p[t=T, l=L] >= 0)
    @variable(m, e_loan_n[t=T, l=L] >= 0)

    @variable(m, e_ann[t=T, l=L])

    @variable(m, e_ann_d_[t=T, l=L, (0,)])
    @variable(m, e_ann_0[t=T, l=L]>=0)
    @variable(m, e_ann_1[t=T, l=L]>=0)

    @variable(m, e_pay[t=T, l=L])
    @variable(m, e_pay_1[t=T, l=L])



    # d090423
    # retirement cost
    @variable(m, t_loan_d_[t=T, l=L, (0,1)] >= 0)  # total loan
    @variable(m, t_ret_cost[t=T, l=L])
    @variable(m, t_ret_cost_d_[t=T, l=L, k=(0,)] >= 0) # only 0th needed

    # 76 
    # 76 #######################################################################
    @variable(m, y_n[t=T, l=L, k=Kn], Bin) # this means plant kind k
    @variable(m, n_c0[l=L])  # the new capacity
    @variable(m, n_c0_d_[t=T, l=L, k=Kn] >= 0) # disaggregated

    @variable(m, n_cp[t=T, l=L])
    @variable(m, n_cp_d_[t=T, l=L, k=Kn] >= 0) # disaggregated variable

    #
    @variable(m, n_l_d_[t=T, l=L, k=Kn])
    @variable(m, n_l[t=T, l=L])

    @variable(m, n_ann_d_[t=T, l=L, k=Kn] >= 0)
    @variable(m, n_ann[t=T, l=L])
    @variable(m, n_ann_0[t=T, l=L]>=0)
    @variable(m, n_ann_1[t=T, l=L]>=0)

    @variable(m, n_ladd_d_[t=T, l=L, (0,1)] >= 0)
    @variable(m, n_ladd[t=T, l=L])

    @variable(m, n_l_pd_[t=T, l=L, (0,1)])

    @variable(m, n_loan[t=T, l=L])
    @variable(m, n_loan_p[t=T, l=L] >= 0)
    @variable(m, n_loan_n[t=T, l=L] >= 0)
    @variable(m, n_yps[t=T, l=L], Bin)

    @variable(m, n_pay[t=T, l=L])
    @variable(m, n_pay_1[t=T, l=L] >= 0)
    
    # 76 
    # 76 #######################################################################
    
    # heating requirement
    @variable(m, n_eh[t=T, l=L]) # 0
    @variable(m, n_eh_d_[t=T, l=L, k=Kn]) # 1
    # fuel requirement 
    @variable(m, n_ehf[t=T, l=L, f=Fu]) # 2
    @variable(m, n_ehf_d_[t=T, l=L, k=Kn, f=Fu]) # 3
    # electricity requirement
    @variable(m, n_u[t=T, l=L]) # 4
    @variable(m, n_u_d_[t=T, l=L, k=Kn]) # 5
    # process (intrinsic) emissions
    @variable(m, n_cp_e[t=T, l=L]) # 6
    @variable(m, n_cp_e_d_[t=T, l=L, k=Kn]) # 7
    # process (extrinsic) disaggregated emissions, e.g. scope 0, etc.
    @variable(m, n_ep0[t=T, l=L]) # 8
    @variable(m, n_ep0_d_[t=T, l=L, k=Kn]) # 9
    # scope 1, emitted
    @variable(m, n_ep1ge[t=T, l=L]) # 10
    @variable(m, n_ep1ge_d_[t=T, l=L, k=Kn]) # 11
    # scope 1, captured
    @variable(m, n_ep1gce[t=T, l=L]) # 12
    @variable(m, n_ep1gce_d_[t=T, l=L, k=Kn]) # 13
    # scope 1,  stored
    @variable(m, n_ep1gcs[t=T, l=L]) # 14
    @variable(m, n_ep1gcs_d_[t=T, l=L, k=Kn]) # 15
    # operating and maintenance
    @variable(m, n_conm_d_[t=T, l=L, k=Kn])  # 16
    @variable(m, n_conm[t=T, l=L])  # 17

    # 76 
    # 76 #######################################################################
    ##
    # k = 0 means on
    @variable(m, o_pay[t=T, l=L])  # we put this in a disjunction so it can be
    @variable(m, o_pay_d_[t=T, l=L, (1,)] >= 0.0)
    # k = 0 means on
    @variable(m, o_conm[t=T, l=L])
    @variable(m, o_conm_d_[t=T, l=L, (1,)] >= 0.0)
    @variable(m, o_tconm_d_[t=T, l=L, (0,1)])
    #
    @variable(m, o_cp[t=T, l=L])
    @variable(m, o_cp_d_[t=T, l=L, (1,)] >= 0) # disaggregated variable

    @variable(m, o_tcp_d_[t=T, l=L, (0,1)] >= 0)

    @variable(m, o_tpay_d_[t=T, l=L, (0,1)] >= 0)
    #
    @variable(m, o_cp_e_d_[t=T, l=L, (1,)] >= 0)
    @variable(m, o_tcp_e_d_[t=T, l=L, (0, 1)] >= 0)

    @variable(m, o_cp_e[t=T, l=L])
    @variable(m, o_tcp_e[t=T, l=L])
    #
    @variable(m, o_ep0_d_[t=T, l=L, (1,)] >= 0)
    @variable(m, o_ep0[t=T, l=L])

    @variable(m, o_tep0_d_[t=T, l=L, (0,1)] >= 0)
    #
    @variable(m, o_ep1ge_d_[t=T, l=L, (1,)] >= 0)
    @variable(m, o_ep1ge[t=T, l=L])

    @variable(m, o_tep1ge_d_[t=T, l=L, (0,1)] >= 0)
    #
    @variable(m, o_ep1gcs[t=T, l=L]) # 14
    @variable(m, o_ep1gcs_d_[t=T, l=L, k=Kr]) # 15

    @variable(m, o_tep1gcs_d_[t=T, l=L, (0,1)] >= 0)
    #


    # 76 
    # 76 #######################################################################
    ##
    # tier 0 logic
    @constraint(m, o_logic_1[t=T, l=L; t<p.t_horizon],
                y_o[t+1, l] <= y_o[t, l]) # this can only go offline

    @constraint(m, o_logic_init[l=L],
                y_o[0, l] == 1)  # plants must start online 


    # tier 1
    #@constraint(m, logic_tier01_0_e[t=T, l=L],  # only one option
    #            sum(y_r[t, l, k] for k in Kr) >= y_o[t, l])

    #@constraint(m, logic_tier01_1_e[t=T, l=L, k=Kr],
    #            y_o[t, l] >= y_r[t, l, k]
    #           )

    #
    @constraint(m, r_logic_init_0[l=L], # all non 0 modes
                y_r[0, l, 0] == 1
               )
    #
    @constraint(m, r_logic_init_1[l=L, k=Kr; k > 0], # all non 0 modes
                y_r[0, l, k] == 0
               )
    #
    @constraint(m, r_budget_s[t=T, l=L; t>0], # only one mode
                sum(y_r[t, l, k] for k in Kr if k > 0) <= 1
               )
    #
    @constraint(m, logic_tier_1_0m_e_[t=T, l=L, k=Kr; k>0], # either 0 or >0 
                1 >= y_r[t, l, k] + y_r[t, l, 0]
               )
    #
    #@constraint(m, r_budget[t=T, l=L, k=Kr; k>0 && t<p.t_horizon],
    #            y_r[t+1, l, k] + (1 - y_o[t+1, l]) >= y_r[t, l, k]
    #           )
    #
    @constraint(m, r_budget[t=T, l=L, k=Kr; k>0 && t<p.t_horizon],
                y_r[t+1, l, k]  >= y_r[t, l, k]
               )
    # continuity


    # 76 
    # 76 #######################################################################
    ##
    # d082923
    # -> expansion
    # p.e_C[l], the capacity per unit of allocation
    @constraint(m, exp_d0_e_[t=T, l=L], # only 0 counts
                e_c_d_[t, l, 0] == p.e_C[l+1] * x_d_[t, l, 0]
               )
    # e_c_bM
    @constraint(m, exp_ncw_m_i0_[t=T, l=L],
                e_c_d_[t, l, 0] <= p.e_c_bM * y_e[t, l]
               )
    @constraint(m, exp_ncw_s_[t=T, l=L],
                e_c[t, l] == e_c_d_[t, l, 0]
               )
    # x_bM
    @constraint(m, exp_x_m_i0_[t=T, l=L],
                x_d_[t, l, 0] <= p.x_bM * y_e[t, l]
               )
    # x_bM (relax)
    @constraint(m, exp_x_m_i1_[t=T, l=L],
                x_d_[t, l, 1] <= p.x_bM * (1 - y_e[t, l])
               )
    @constraint(m, exp_x_s_[t=T, l=L],
                x[l] == x_d_[t, l, 0] + x_d_[t, l, 1]
               )

    # -> expansion cost (loan)
    # p.e_loanFact
    @constraint(m, e_l_d0_e_[t=T, l=L], # only zero counts
                e_l_d_[t, l, 0] == p.e_loanFact[l+1] * x_d_[t, l, 0]
               )
    #@constraint(m, e_l_d1_e_[t=T, l=L], # only zero counts
    #            e_l_d_[t, l, 1] == 0.0 
    #           )
    
    # e_l_bM
    @constraint(m, e_l_m_i0_[t=T, l=L],
                e_l_d_[t, l, 0] <= p.e_l_bM * y_e[t, l]
               )
    # e_l_bM
    #@constraint(m, e_l_m_i1_[t=T, l=L],
    #            e_l_d_[t, l, 1] <= e_l_bM * (1 - y_e[t, l])
    #           )
    @constraint(m, e_l_s_[t=T, l=L],
                #e_l[t, l] == sum(e_l_d_[t, l, k] for k in (0,1))
                e_l[t, l] == e_l_d_[t, l, 0]
               )

    # -> expansion cost annuity (annual payment)
    # p.e_Ann
    @constraint(m, e_ann_d0_e_[t=T, l=L], # only zero counts
                e_ann_d_[t, l, 0] == p.e_Ann[l+1] * x_d_[t, l, 0]
               )
    # e_ann_bM
    @constraint(m, e_ann_m_i0_[t=T, l=L],
                e_ann_d_[t, l, 0] <= p.e_ann_bM * y_e[t, l]
               )
    @constraint(m, e_ann_s_[t=T, l=L],
                e_ann[t, l] == e_ann_d_[t, l, 0]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion logic
    @constraint(m, e_logic_init[l=L], 
                y_e[0, l] == 0  # start not-expanded
               )
    @constraint(m, e_logic_1[t=T, l=L; t<p.t_horizon],
                y_e[t+1, l] >= y_e[t, l]  # only expand in the future
               )

    # 76 
    # 76 #######################################################################
    ##
    # capacity expansion loans (they should be agnostic to retrofit or rf)
    # components: e_ladd, e_ann, e_pslack, e_ploan
    # -> e_add
    # exp add loan
    @constraint(m, e_ladd_d0_e_[t=T, l=L],
                e_ladd_d_[t, l, 0] == e_l_pd_[t, l, 0]
               )
    @constraint(m, e_ladd_m_i0_[t=T, l=L; t>0], # y_e goes from 0 to 1
                e_ladd_d_[t, l, 0] <= 
                p.e_ladd_bM * (y_e[t, l] - y_e[t-1, l])
               )  # 
    @constraint(m, e_l_m_i1_[t=T, l=L; t>0],  # we need to set this to 0
                e_l_pd_[t, l, 1] <= 
                p.e_l_bM * (1-y_e[t, l]+y_e[t-1, l])
               )  #  the 0th component is implied by the e_ladd_m constr
    @constraint(m, e_ladd_s_e[t=T, l=L],
                e_ladd[t, l] == e_ladd_d_[t, l, 0]
               )
    @constraint(m, e_l_s_e_[t=T, l=L],
                e_l[t, l] == 
                #sum(e_l_pd_[t, l, k] for k in (0, 1))
                e_l_pd_[t, l, 0] + e_l_pd_[t, l, 1]
               )

    # load
    @constraint(m, e_loan_s_e_[t=T, l=L],
                e_loan[t, l] == e_loan_p[t, l] - e_loan_n[t, l]
               )
    @constraint(m, e_loan_p_m0_i_[t=T, l=L],
                e_loan_p[t, l] <= p.e_loan_bM * (1 - e_yps[t, l])
               )
    @constraint(m, e_loan_n_m0_i_[t=T, l=L],
                e_loan_n[t, l] <= p.e_loan_bM * e_yps[t, l]
               )
    # pay
    @constraint(m, e_pay_s_e_[t=T, l=L],
                e_pay[t, l] == e_pay_1[t, l]
               )
    @constraint(m, e_pay_n_m0_i_[t=T, l=L],
                e_pay_1[t, l] <= p.e_pay_bM * (1 - e_yps[t, l])
               )
    # annuity
    @constraint(m, e_ann_s_e_[t=T, l=L],
                e_ann[t, l] == e_ann_0[t, l] + e_ann_1[t, l]
               )
    @constraint(m, e_ann_0_m_i_[t=T, l=L],
                e_ann_0[t, l] <= p.e_ann_bM * e_yps[t, l]
               )
    @constraint(m, e_ann_1_m_i_[t=T, l=L],
                e_ann_1[t, l] <= p.e_ann_bM * (1 - e_yps[t, l])
               )
    @constraint(m, e_ann_dpay0_e_[t=T, l=L],
                e_ann_1[t, l] == e_pay_1[t, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    @constraint(m, e_loan_bal_e_[t=T, l=L; t<last(T)],
                e_loan[t+1, l] == e_loan[t, l]
                - e_pay[t, l]
                + e_ladd[t, l]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    # d082923
    # base capacity
    @constraint(m, cpb_e_[t=T, l=L],
                cpb[t, l] == p.c0[l+1] + e_c[t, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    # retrofit 
    # -> capacity.
    # p.Kr (capacity factor, e.g. prod in retrofit r / base prod)
    # viz. cap_mod = factor * cap_base
    @constraint(m, r_cp_d_e_[t=T, l=L, k=Kr],
                r_cp_d_[t, l, k] == p.r_c_C[l+1, k+1] * r_cpb_d_[t, l, k] 
                + p.r_rhs_C[l+1, k+1] * y_r[t, l, k]
               ) 
    # r_cp_bM, retrofit capacity
    @constraint(m, r_cp_d_m_e_[t=T, l=L, k=Kr],
                r_cp_d_[t, l, k] <= p.r_cp_bM * y_r[t, l, k]
               )
    # cpb_bM, base capacity
    @constraint(m, r_cpb_d_m_i_[t=T, l=L, k=Kr],
                r_cpb_d_[t, l, k] <= p.r_cpb_bM * y_r[t, l, k]
               )
    @constraint(m, r_cp_s_e_[t=T, l=L],
                r_cp[t, l] == sum(r_cp_d_[t, l, k] for k in Kr)
               )
    @constraint(m, r_cpb_s_e_[t=T, l=L],
                cpb[t, l] == sum(r_cpb_d_[t, l, k] for k in Kr)
               )
    # -> heating requirement.
    # p.Hm (heating factor, i.e. heat / product)
    @constraint(m, r_eh_d_e_[t=T, l=L, k=Kr],
                r_eh_d_[t, l, k] == p.r_c_H[l+1, k+1] * r_cp_d_[t, l, k] 
                + p.r_rhs_H[l+1, k+1] * y_r[t, l, k]
               )
    @constraint(m, r_eh_d_m_i_[t=T, l=L, k=Kr],
                r_eh_d_[t, l, k] <= p.r_eh_bM * y_r[t, l, k]
               )
    @constraint(m, r_eh_s_e_[t=T, l=L],
                r_eh[t, l] == sum(r_eh_d_[t, l, k] for k in Kr)
               )
    # d082923
    # -> fuel required for heat.
    # p.r_c_F in [0,1], (i.e. a fraction, heat by fuel /tot heat)
    @constraint(m, r_ehf_d_e_[t=T, l=L, k=Kr, f=Fu],
                r_ehf_d_[t, l, k, f] ==
                p.r_c_F[l+1, k+1, f+1] * r_eh_d_[t, l, k] 
                + p.r_rhs_F[l+1, k+1, f+1] * y_r[t, l, k]
               )
    # r_ehf_bM
    @constraint(m, r_ehf_m_i_[t=T, l=L, k=Kr, f=Fu],
                r_ehf_d_[t, l, k, f] <= p.r_ehf_bM * y_r[t, l, k]
               )
    @constraint(m, r_ehf_s_e_[t=T, l=L, f=Fu],
                r_ehf[t, l, f] == sum(r_ehf_d_[t, l, k, f] for k in Kr)
               )

    # -> electricity requirement
    # r_u and m_ud_, p.Um & p.UmRhs
    @constraint(m, r_u_d_e_[t=T, l=L, k=Kr],
                r_u_d_[t, l, k] == p.r_c_U[l+1, k+1] * r_cp_d_[t, l, k] 
                + p.r_rhs_U[l+1, k+1] * y_r[t, l, k]
               )
    # r_u_bM (big-M)
    @constraint(m, r_u_i_[t=T, l=L, k=Kr],
                r_u_d_[t, l, k] <= p.r_u_bM * y_r[t, l, k]
               )
    @constraint(m, r_u_s_e_[t=T, l=L],
                r_u[t, l] == sum(r_u_d_[t, l, k] for k in Kr)
               )

    # -> process (intrinsic) emissions
    # r_cp_e, r_cp_e_d_. p.Cp & p.CpRhs
    @constraint(m, r_cp_e_d_e_[t=T, l=L, k=Kr],
                r_cp_e_d_[t, l, k] == p.r_c_cp_e[l+1, k+1] * r_cp_d_[t, l, k]
                + p.r_rhs_cp_e[l+1, k+1] * y_r[t, l, k]
               )
    # r_cp_e_bM
    @constraint(m, r_cp_e_m_i_[t=T, l=L, k=Kr],
                r_cp_e_d_[t, l, k] <= p.r_cp_e_bM * y_r[t, l, k]
               )
    @constraint(m, r_cp_e_s_e_[t=T, l=L],
                r_cp_e[t, l] == sum(r_cp_e_d_[t, l, k] for k in Kr))

    # -> -> process (disaggregated) emissions

    # -> scope 0 emission
    # Fef (fuel emission factor)
    @constraint(m, r_ep0_d_e_[t=T, l=L, k=Kr],
                r_ep0_d_[t, l, k] == 
                sum(p.r_Fef[l+1, k+1, f+1] * r_ehf_d_[t, l, k, f] for f in Fu) +
                r_cp_e_d_[t, l, k]
               )
    # r_ep0_bM
    @constraint(m, r_ep0_m_i_[t=T, l=L, k=Kr],
                r_ep0_d_[t, l, k] <= p.r_ep0_bM * y_r[t, l, k]
               )
    @constraint(m, r_ep0_s_e_[t=T, l=L],
                r_ep0[t, l] == sum(r_ep0_d_[t, l, k] for k in Kr)
               )

    # -> scope 1 emitted
    # p.r_chi
    @constraint(m, r_ep1ge_d_e_[t=T, l=L, k=Kr],
                r_ep1ge_d_[t, l, k] == (1.0 - p.r_chi[l+1, k+1]) 
                * r_ep0_d_[t, l, k]
               )
    # r_ep1ge_bM
    @constraint(m, r_ep1ge_m_i_[t=T, l=L, k=Kr],
                r_ep1ge_d_[t, l, k] <= p.r_ep1ge_bM * y_r[t, l, k]
               )
    @constraint(m, r_ep1ge_s_e_[t=T, l=L],
                r_ep1ge[t, l] == sum(r_ep1ge_d_[t, l, k] for k in Kr)
               )

    # -> scope 1 captured
    # p.r_sigma
    @constraint(m, r_ep1gce_d_e_[t=T, l=L, k=Kr],
                r_ep1gce_d_[t, l, k] == 
                p.r_chi[l+1, k+1] * (1 - p.r_sigma[l+1, k+1]) 
                * r_ep0_d_[t, l, k]
               )
    @constraint(m, r_ep1gce_m_i_[t=T, l=L, k=Kr],
                r_ep1gce_d_[t, l, k] <= p.r_ep1gce_bM * y_r[t, l, k]
               )
    @constraint(m, r_ep1gce_s_e_[t=T, l=L],
                r_ep1gce[t, l] == sum(r_ep1gce_d_[t, l, k] for k in Kr)
               )

    # -> scope 1 stored
    # p.sigma ?
    @constraint(m, r_ep1gcs_d_e_[t=T, l=L, k=Kr],
                r_ep1gcs_d_[t, l, k] ==
                p.r_chi[l+1, k+1] * p.r_sigma[l+1, k+1] * r_ep0_d_[t, l, k]
               )
    # ep1gcsm_bM
    @constraint(m, r_ep1gcs_m_i_[t=T, l=L, k=Kr],
                r_ep1gcs_d_[t, l, k] <= p.r_ep1gcs_bM * y_r[t, l, k]
               )

    @constraint(m, r_ep1gcs_s_e_[t=T, l=L],
                r_ep1gcs[t, l] == sum(r_ep1gcs_d_[t, l, k] for k in Kr)
               )

    
    # params include carbon intensity, chi and sigma.
    # chi : carbon capture
    # sigma : carbon storage
    
    # r_cp_e and r_cp_e_d_
    
    # ep_0 = (sum(hj*chj)+cp) * prod
    # ep_1ge = ep_0 * (1-chi)
    # ep_1gce = ep_0 * chi * (1-sigma)
    # ep_1gcs # stored ?
    # ep_1nge # not generated emmited
    # ep_2 = (sum(uj*cuj))

    # 76 
    # 76 #######################################################################
    ##
    # -> operating and maintenance
    # p.m_Onm
    @constraint(m, r_conm_d_e_[t=T, l=L, k=Kr],
                r_conm_d_[t, l, k] == p.r_c_Onm[l+1, k+1] * r_cp_d_[t, l, k]
                + p.r_rhs_Onm[l+1, k+1] * y_r[t, l, k]
               )
    @constraint(m, r_conm_m_i_[t=T, l=L, k=Kr],
                r_conm_d_[t, l, k] <= p.r_conm_bM * y_r[t, l, k]
               )
    @constraint(m, r_conm_s_e_[t=T, l=L],
                r_conm[t, l] == sum(r_conm_d_[t, l, k] for k in Kr)
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> retrofit disagg LOAN added amount (ladd), p.r_loanFact
    @constraint(m, r_l_md_e_[t=T, l=L, k=Kr], # retrofit disagg
                r_l_md_[t, l, k] == p.r_loanFact[t+1, l+1, k+1]*r_cp_d_[t, l, k]
               ) # uses c0+expc
    @constraint(m, r_l_mm_i_[t=T, l=L, k=Kr], # big-M
                r_l_md_[t, l, k] <= p.r_l_bM * y_r[t, l, k]
               )
    @constraint(m, r_l_ms_e_[t=T, l=L],
                r_l[t, l] == sum(r_l_md_[t, l, k] for k in Kr)
               )
    # 76 
    # 76 #######################################################################
    ##
    # -> retrofit disagg payment (associated with the payment)
    @constraint(m, r_ann_md_e_[t=T, l=L, k=Kr],
                r_ann_md_[t, l, k] == p.r_annf[t+1, l+1, k+1] * r_cp_d_[t, l, k]
               ) # uses c0+expc
    @constraint(m, r_ann_mm_i_[t=T, l=L, k=Kr],
                r_ann_md_[t, l, k] <= p.r_ann_bM * y_r[t, l, k]
               )
    @constraint(m, r_ann_s_e_[t=T, l=L],
                r_ann[t, l] == sum(r_ann_md_[t, l, k] for k in Kr)
               )

    # 76 
    # 76 #######################################################################
    ##
    # --> padd switch (loan-switch), associated with the loan
    @constraint(m, r_ladd_d_e_[t=T, l=L], 
                r_ladd_d_[t, l, 0] == r_l_pd_[t, l, 0]
               )
    # p_add_bM
    @constraint(m, r_ladd_m_i0_[t=T, l=L; t>0],
                r_ladd_d_[t, l, 0] <= # goes from 1 to 0 only
                p.r_ladd_bM * (y_r[t-1,l,0] - y_r[t,l,0])
               )
    # (zero otw.)
    #@constraint(m, m_ladd_m_i1_[t=T, l=L],
    #            r_ladd_d_[t, l, 1] <=
    #            maximum(r_ladd_bM[t, l, k])*(1-y_r[t-1,l,0]+y_r[t,l,0])
    #           )
    #
    @constraint(m, r_ladd_s_e_[t=T, l=L],
                r_ladd[t, l] == 
                #sum(r_ladd_d_[t, l, k] for k in (0, 1))
                r_ladd_d_[t, l, 0]
               )
    # padd switch (loan-retrofit)
    @constraint(m, r_l_m_i0_[t=T, l=L; t>0],
                r_l_pd_[t, l, 0] <=
                p.r_l_bM * (y_r[t-1,l,0] - y_r[t,l,0])
               )
    #
    @constraint(m, r_l_m_i1_[t=T, l=L; t>0],
                r_l_pd_[t, l, 1] <=
                p.r_l_bM * (1.0 - y_r[t-1,l,0] + y_r[t,l,0])
               )
    # connection to the retrofit 
    @constraint(m, m_l_s_e_[t=T, l=L], # slack (1) is necessary here
                r_l[t, l] == sum(r_l_pd_[t, l, k] for k in (0, 1))
               )
    
    # 76 
    # 76 #######################################################################
    ##
    # -> loan disaggregation, r_yps = 1 if paid, 0 otw
    @constraint(m, r_loan_ps_s_e_[t=T, l=L],
                r_loan[t, l] == r_loan_p[t, l] - r_loan_n[t, l]
               )
    @constraint(m, r_loan_p_m0_i_[t=T, l=L],
                r_loan_p[t, l] <= p.r_loan_bM * (1 - r_yps[t, l])
               )
    @constraint(m, r_loan_n_m0_i_[t=T, l=L],
                r_loan_n[t, l] <= p.r_loan_bM * r_yps[t, l]
               )
    # payment
    @constraint(m, r_pay_s_e_[t=T, l=L],
                r_pay[t, l] == r_pay_1[t, l]
               )
    @constraint(m, r_pay_1_m0_i_[t=T, l=L],
                r_pay_1[t, l] <= p.r_pay_bM * (1 - r_yps[t, l])
               )
    # annuity
    @constraint(m, r_ann_01_s_e_[t=T, l=L],
                r_ann[t, l] == r_ann_0[t, l] + r_ann_1[t, l]
               )
    @constraint(m, r_ann_0_m_i_[t=T, l=L],
                r_ann_0[t, l] <= p.r_ann_bM * r_yps[t, l]
               )
    @constraint(m, r_ann_1_m_i_[t=T, l=L],
                r_ann_1[t, l] <= p.r_ann_bM * (1 - r_yps[t, l])
               )
    @constraint(m, r_ann_d0_e_[t=T, l=L],
                r_ann_1[t, l] == r_pay_1[t, l]
               )

    # if the plant becomes retired r_pay_p might still be positive
    # 76 
    # 76 #######################################################################
    ##
    # -> loan balance
    @constraint(m, r_loan_bal_e_[t=T, l=L;t<last(T)],
                r_loan[t+1, l] == r_loan[t, l] 
                - r_pay[t, l] 
                + r_ladd[t, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    # -> retirement
    # retirement is just r_loan, we just need a way to activate it. 
    # a) have a switch using y_o going from 0 to 1
    # b) take the snapshot of the current value of r_loan and use it as the cost
    #
    @constraint(m, t_ret_c_d_e_[t=T, l=L],  # only enforceable at the switch
                t_ret_cost_d_[t, l, 0] == t_loan_d_[t, l, 0]
               )
    # t_ret_c_bM
    @constraint(m, t_ret_c_bm0_i_[t=T, l=L; t>0],
                t_ret_cost_d_[t, l, 0] <= 
                p.t_ret_c_bM * (y_o[t-1, l] - y_o[t, l])
               )
    # m_loan_d_bM
    @constraint(m, r_loan_d_bm0_i_[t=T, l=L; t>0],
                t_loan_d_[t, l, 0] <=  # retired
                p.t_loan_bM * (y_o[t-1, l] - y_o[t, l])
               )
    # m_loan_d_bM
    @constraint(m, r_loan_d_bm1_i_[t=T, l=L; t>0], # not retired
                t_loan_d_[t, l, 1] <= p.t_loan_bM * 
                (1 + y_o[t, l] - y_o[t-1, l])
               )
    @constraint(m, t_ret_c_s_e_[t=T, l=L],
                #t_ret_c[t, l] == sum(t_ret_cost_d_[t, l, k] for k in (0, 1))
                t_ret_cost[t, l] == t_ret_cost_d_[t, l, 0]  # only one needed
               )
               
    @constraint(m, r_loan_s_e_[t=T, l=L],  # total loan
                r_loan_p[t, l] + e_loan_p[t, l]
                == t_loan_d_[t, l, 0] + t_loan_d_[t, l, 1]
               )
    # -> total payment (retrofit/existing + )
    @constraint(m, o_pay_s_e_[t=T, l=L],
                o_pay[t, l] == o_pay_d_[t, l, 1]
               )
    @constraint(m, o_pay_d1_e_[t=T, l=L], # online
                # there is pay for retrof but not for expansion
                o_pay_d_[t, l, 1] == o_tpay_d_[t, l, 1]
               )
    @constraint(m, o_pay_m_1_i_[t=T, l=L], # on
                o_pay_d_[t, l, 1] <= p.o_pay_bM * y_o[t, l]
               )
    # this one includes both expansion and retrof
    @constraint(m, o_tpay_s_e_[t=T, l=L],
                r_pay[t, l] + e_pay[t, l] == 
                o_tpay_d_[t, l, 0] + o_tpay_d_[t, l, 1]
               )
    #
    @constraint(m, o_tpay_m1_i_[t=T, l=L], # on
                o_tpay_d_[t, l, 1] <= p.o_pay_bM * y_o[t, l]
               )
    #
    @constraint(m, o_tpay_m0_i_[t=T, l=L], # off
                o_tpay_d_[t, l, 0] <= p.o_pay_bM * (1 - y_o[t, l])
               )
    # -> o&m
    @constraint(m, o_conm_s_e_[t=T, l=L],
                o_conm[t, l] == o_conm_d_[t, l, 1]
               )
    @constraint(m, o_conm_d1_e_[t=T, l=L],
                o_conm_d_[t, l, 1] == o_tconm_d_[t, l, 1]
               )

    @constraint(m, o_conm_m_1_i_[t=T, l=L], # on
                o_conm_d_[t, l, 1] <= p.o_conm_bM * y_o[t, l]
               )
    #
    @constraint(m, o_tconm_m1_i_[t=T, l=L],
                o_tconm_d_[t, l, 1] <= p.o_conm_bM * y_o[t, l]
               )
    @constraint(m, o_tconm_m0_i_[t=T, l=L],
                o_tconm_d_[t, l, 0] <= p.o_conm_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tconm_s_e_[t=T, l=L],
                r_conm[t, l] == o_tconm_d_[t, l, 0] + o_tconm_d_[t, l, 1]
               )


    # 76 
    # 76 #######################################################################
    ##
    # -> new plant capacity
    #@constraint(m, n_c0_e_[l=L],
    #            n_c0[l] >= p.c0[l+1]
    #           )
    # -> new plant capacity disaggregation
    @constraint(m, n_cp0_d_e_[t=T, l=L],
                n_cp_d_[t, l, 0] == 0
               ) # how do we make this 0 at k=0?
    @constraint(m, n_cp_d_e_[t=T, l=L, k=Kn; k>0],
                n_cp_d_[t, l, k] == n_c0_d_[t, l, k]
               )
    #@constraint(m, n_c0_d0_i_[t=T, l=L, k=Kn; k>0],
    #            n_c0_d_[t, l, k] >= p.c0[l+1] * y_n[t, l, k]
    #           )

    @constraint(m, n_cp_m_i_[t=T, l=L, k=Kn; k>0], # skip the 0-th
                n_cp_d_[t, l, k] <= p.n_cp_bM * y_n[t, l, k]
               )

    @constraint(m, n_cap_add_d_m_i_[t=T, l=L, k=Kn],
                n_c0_d_[t, l, k] <= p.n_c0_bM * y_n[t, l, k]
               )
    #
    @constraint(m, n_cp_s_e_[t=T, l=L],
                n_cp[t, l] == sum(n_cp_d_[t, l, k] for k in Kn)
               )

    @constraint(m, n_c0_s_e_[t=T, l=L],
                n_c0[l] == sum(n_c0_d_[t, l, k] for k in Kn)
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> base loan (proportional to the capacity)
    @constraint(m, n_l_d0_e_[t=T, l=L, k=Kn],
                n_l_d_[t, l, k] == p.n_loanFact[l+1, k+1] * n_cp_d_[t, l, k]
               ) # this should be 0 at 0th  
    # perhaps this should be n_c0_d_
    @constraint(m, n_l_m_i0_[t=T, l=L, k=Kn],
                n_l_d_[t, l, k] <= p.n_l_bM * y_n[t, l, k]
               )
    @constraint(m, n_l_s_[t=T, l=L],
                n_l[t, l] == sum(n_l_d_[t, l, k] for k in Kn)
               )
    # -> annuity (how much we pay)
    @constraint(m, n_ann_d0_e_[t=T, l=L, k=Kn],
                n_ann_d_[t, l, k] == p.n_Ann[l+1, k+1] * n_cp_d_[t, l, k]
               )
    @constraint(m, n_ann_m_i0_[t=T, l=L, k=Kn],
                n_ann_d_[t, l, k] <= p.n_ann_bM * y_n[t, l, k]
               )
    @constraint(m, n_ann_s_[t=T, l=L],
                n_ann[t, l] == sum(n_ann_d_[t, l, k] for k in Kn)
               )
    
    # 76 
    # 76 #######################################################################
    ##
    # -> 
    # link the n_cost_ to the ladd in a single time period.
    @constraint(m, n_ladd_m_i0_[t=T, l=L; t>0], 
                n_ladd_d_[t, l, 0] <= 
                p.n_ladd_bM * (y_n[t-1, l, 0] - y_n[t, l, 0])
               ) # y_n goes from 1 to 0
    @constraint(m, n_ladd_m_i1_[t=T, l=L;t>0],
                n_l_pd_[t, l, 1] <= 
                p.n_l_bM * (1 + y_n[t, l, 0] - y_n[t-1, l, 0])
               )
    #
    @constraint(m, n_ladd_d0_e_[t=T, l=L],
                n_ladd_d_[t, l, 0] == n_l_pd_[t, l, 0]
               )
    @constraint(m, n_ladd_s_e_[t=T, l=L],
                n_ladd[t, l] == n_ladd_d_[t, l, 0]
               )
    @constraint(m, n_l_s_e_[t=T, l=L],
                n_l[t, l] ==
                n_l_pd_[t, l, 0] + n_l_pd_[t, l, 1]
               )

    # loan
    @constraint(m, n_loan_s_e_[t=T, l=L],
                n_loan[t, l] == n_loan_p[t, l] - n_loan_n[t, l]
               )
    @constraint(m, n_loan_p_m0_i_[t=T, l=L],
                n_loan_p[t, l] <= p.n_loan_bM * (1 - n_yps[t, l])
               )
    @constraint(m, n_loan_n_m0_i_[t=T, l=L],
                n_loan_n[t, l] <= p.n_loan_bM * n_yps[t, l]
               )
    # pay
    @constraint(m, n_pay_s_e_[t=T, l=L],
                n_pay[t, l] == n_pay_1[t, l]
               )
    # note: this mechanism seems to have a lag of 1 period.
    @constraint(m, n_pay_n_m0_i_[t=T, l=L],
                n_pay_1[t, l] <= p.n_pay_bM * (1 - n_yps[t, l])
               )
    # annuity
    @constraint(m, n_ann_s_e_[t=T, l=L],
                n_ann[t, l] == n_ann_0[t, l] + n_ann_1[t, l]
               )
    @constraint(m, n_ann_0_m_i_[t=T, l=L],
                n_ann_0[t, l] <= p.n_ann_bM * n_yps[t, l]
               )
    @constraint(m, n_ann_1_m_i_[t=T, l=L],
                n_ann_1[t, l] <= p.n_ann_bM * (1 - n_yps[t, l])
               )
    @constraint(m, n_pay_d0_e_[t=T, l=L],
                n_pay_1[t, l] == n_ann_1[t, l]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    @constraint(m, n_loan_bal_e_[t=T, l=L; t<last(T)],
                n_loan[t+1, l] == n_loan[t, l]
                - n_pay[t, l]
                + n_ladd[t, l]
               )
    # 76 
    # 76 #######################################################################
    ##

    # -> heating requirement.
    # p.Hm (heating factor, i.e. heat / product)
    @constraint(m, n_eh_d_e_[t=T, l=L, k=Kn],
                n_eh_d_[t, l, k] == p.n_c_H[l+1, k+1] * n_cp_d_[t, l, k] 
                + p.n_rhs_h[l+1, k+1] * y_n[t, l, k]
               )
    @constraint(m, n_eh_d_m_i_[t=T, l=L, k=Kn],
                n_eh_d_[t, l, k] <= p.n_eh_bM * y_n[t, l, k]
               )
    @constraint(m, n_eh_s_e_[t=T, l=L],
                n_eh[t, l] == sum(n_eh_d_[t, l, k] for k in Kn)
               )
    # -> fuel required for heat.
    # p.n_c_F in [0,1], (i.e. a fraction, heat by fuel /tot heat)
    @constraint(m, n_ehf_d_e_[t=T, l=L, k=Kn, f=Fu],
                n_ehf_d_[t, l, k, f] ==
                p.n_c_F[l+1, k+1, f+1] * n_eh_d_[t, l, k] 
                + p.n_rhs_F[l+1, k+1, f+1] * y_n[t, l, k]
               )
    # n_ehf_bM
    @constraint(m, n_ehf_m_i_[t=T, l=L, k=Kn, f=Fu],
                n_ehf_d_[t, l, k, f] <= p.n_ehf_bM * y_n[t, l, k]
               )
    @constraint(m, n_ehf_s_e_[t=T, l=L, f=Fu],
                n_ehf[t, l, f] == sum(n_ehf_d_[t, l, k, f] for k in Kn)
               )

    # -> electricity requirement
    # n_u and m_ud_, p.Um & p.UmRhs
    @constraint(m, n_u_d_e_[t=T, l=L, k=Kn],
                n_u_d_[t, l, k] == p.n_c_U[l+1, k+1] * n_cp_d_[t, l, k] 
                + p.n_rhs_U[l+1, k+1] * y_n[t, l, k]
               )
    # n_u_bM (big-M)
    @constraint(m, n_u_i_[t=T, l=L, k=Kn],
                n_u_d_[t, l, k] <= p.n_u_bM * y_n[t, l, k]
               )
    @constraint(m, n_u_s_e_[t=T, l=L],
                n_u[t, l] == sum(n_u_d_[t, l, k] for k in Kn)
               )
    
    # -> process (intrinsic) emissions
    # n_cp_e, n_cp_e_d_. p.Cp & p.CpRhs
    @constraint(m, n_cp_e_d_e_[t=T, l=L, k=Kn],
                n_cp_e_d_[t, l, k] == p.n_c_cp_e[l+1, k+1] * n_cp_d_[t, l, k]
                + p.n_rhs_cp_e[l+1, k+1] * y_n[t, l, k]
               )
    # n_cp_e_bM
    @constraint(m, n_cp_e_m_i_[t=T, l=L, k=Kn],
                n_cp_e_d_[t, l, k] <= p.n_cp_e_bM * y_n[t, l, k]
               )
    @constraint(m, n_cp_e_s_e_[t=T, l=L],
                n_cp_e[t, l] == sum(n_cp_e_d_[t, l, k] for k in Kn)
               )

    # -> -> process (disaggregated) emissions

    # -> scope 0 emission
    # Fef (fuel emission factor)
    @constraint(m, n_ep0_d_e_[t=T, l=L, k=Kn],
                n_ep0_d_[t, l, k] == 
                sum(p.n_Fef[l+1, k+1, f+1] * n_ehf_d_[t, l, k, f] for f in Fu) +
                n_cp_e_d_[t, l, k]
               )
    # n_ep0_bM
    @constraint(m, n_ep0_m_i_[t=T, l=L, k=Kn],
                n_ep0_d_[t, l, k] <= p.n_ep0_bM * y_n[t, l, k]
               )
    @constraint(m, n_ep0_s_e_[t=T, l=L],
                n_ep0[t, l] == sum(n_ep0_d_[t, l, k] for k in Kn)
               )
    # -> scope 1 emitted
    # p.n_chi
    @constraint(m, n_ep1ge_d_e_[t=T, l=L, k=Kn],
                n_ep1ge_d_[t, l, k] == (1.0 - p.n_chi[l+1, k+1]) 
                * n_ep0_d_[t, l, k]
               )
    # n_ep1ge_bM
    @constraint(m, n_ep1ge_m_i_[t=T, l=L, k=Kn],
                n_ep1ge_d_[t, l, k] <= p.n_ep1ge_bM * y_n[t, l, k]
               )
    @constraint(m, n_ep1ge_s_e_[t=T, l=L],
                n_ep1ge[t, l] == sum(n_ep1ge_d_[t, l, k] for k in Kn)
               )

    # -> scope 1 captured
    # p.n_sigma
    @constraint(m, n_ep1gce_d_e_[t=T, l=L, k=Kn],
                n_ep1gce_d_[t, l, k] == 
                p.n_chi[l+1, k+1] * (1 - p.n_sigma[l+1, k+1]) 
                * n_ep0_d_[t, l, k]
               )
    @constraint(m, n_ep1gce_m_i_[t=T, l=L, k=Kn],
                n_ep1gce_d_[t, l, k] <= p.n_ep1gce_bM * y_n[t, l, k]
               )
    @constraint(m, n_ep1gce_s_e_[t=T, l=L],
                n_ep1gce[t, l] == sum(n_ep1gce_d_[t, l, k] for k in Kn)
               )

    # -> scope 1 stored
    # p.sigma ?
    @constraint(m, n_ep1gcs_d_e_[t=T, l=L, k=Kn],
                n_ep1gcs_d_[t, l, k] ==
                p.n_chi[l+1, k+1] * p.n_sigma[l+1, k+1] * n_ep0_d_[t, l, k]
               )
    # ep1gcsm_bM
    @constraint(m, n_ep1gcs_m_i_[t=T, l=L, k=Kn],
                n_ep1gcs_d_[t, l, k] <= p.n_ep1gcs_bM * y_n[t, l, k]
               )

    @constraint(m, n_ep1gcs_s_e_[t=T, l=L],
                n_ep1gcs[t, l] == sum(n_ep1gcs_d_[t, l, k] for k in Kn)
               )
    
    # params include carbon intensity, chi and sigma.
    # chi : carbon capture
    # sigma : carbon storage
    
    # n_cp_e and n_cp_e_d_
    
    # ep_0 = (sum(hj*chj)+cp) * prod
    # ep_1ge = ep_0 * (1-chi)
    # ep_1gce = ep_0 * chi * (1-sigma)
    # ep_1gcs # stored ?
    # ep_1nge # not generated emmited
    # ep_2 = (sum(uj*cuj))

    # 76 
    # 76 #######################################################################
    ##
    # -> operating and maintenance
    # p.m_Onm
    @constraint(m, n_conm_d_e_[t=T, l=L, k=Kn],
                n_conm_d_[t, l, k] == p.n_c_Onm[l+1, k+1] * n_cp_d_[t, l, k]
                + p.n_rhs_Onm[l+1, k+1] * y_n[t, l, k]
               )
    @constraint(m, n_conm_m_i_[t=T, l=L, k=Kn],
                n_conm_d_[t, l, k] <= p.n_conm_bM * y_n[t, l, k]
               )
    @constraint(m, n_conm_s_e_[t=T, l=L],
                n_conm[t, l] == sum(n_conm_d_[t, l, k] for k in Kn)
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> new alloc logic
    @constraint(m, n_logic0_init[l=L], 
                y_n[0, l, 0] == 1  # start with no new
               )
    @constraint(m, n_logic1_init[l=L, k=Kn; k>0], 
                y_n[0, l, k] == 0  # start with no new facility
               )
    @constraint(m, n_logic_0[t=T, l=L; t<p.t_horizon],
                y_n[t+1, l, 0] <= y_n[t, l, 0]  # only build in the future
               )
    @constraint(m, n_logic_1[t=T, l=L, k=Kn; t<p.t_horizon && k>0],
                y_n[t+1, l, k] >= y_n[t, l, k]  # only build in the future
               )
    @constraint(m, n_log_s_e[t=T, l=L], 
                sum(y_n[t, l, k] for k in Kn) == 1
               )
    @constraint(m, n_logic_2[t=T, l=L, k=Kn; k>0 && t>0],
                y_n[t, l, k] + y_o[t, l] <= 1
               )
    
    # 76
    # 76 #######################################################################
    ##
    # -> initial loans
    @constraint(m, r_loan_initial_cond[l=L],
                r_loan[0, l] == p.r_loan0[l+1]
               )
    @constraint(m, e_loan_init_cond[l=L],
                e_loan[0, l] == 0.0
               )
    @constraint(m, n_loan_init_cond[l=L],
                n_loan[0, l] == 0.0
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> output capacity
    @constraint(m, o_cp_d_e[t=T, l=L],
                o_cp_d_[t, l, 1] == o_tcp_d_[t, l, 1]
               )
    @constraint(m, o_cp_m1_i_[t=T, l=L], # on
                o_cp_d_[t, l, 1] <= p.o_cp_bM * y_o[t, l]
               )
    #
    @constraint(m, o_cp_s_e_[t=T, l=L],
                o_cp[t, l] == o_cp_d_[t, l, 1]
               )
    #
    @constraint(m, o_tcp_d_m1_i_[t=T, l=L], # on
                o_tcp_d_[t, l, 1] <= p.o_cp_bM * y_o[t, l]
               )
    @constraint(m, o_tcp_d_m0_i_[t=T, l=L], # off
                o_tcp_d_[t, l, 0] <= p.o_cp_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tcp_d_s_e_[t=T, l=L],
                r_cp[t, l] == o_tcp_d_[t, l, 0] + o_tcp_d_[t, l, 1]
               ) # total cap
    # 76 
    # 76 #######################################################################
    ##
    # -> output emissions
    @constraint(m, o_cp_e_d_e_[t=T, l=L],
                o_cp_e_d_[t, l, 1] == o_tcp_e_d_[t, l, 1]
               )
    @constraint(m, o_cp_e_m_i_[t=T, l=L], # on
                o_cp_e_d_[t, l, 1] <= p.o_cp_e_bM * y_o[t, l]
               )
    @constraint(m, o_cp_e_s_[t=T, l=L],
                o_cp_e[t, l] == o_cp_e_d_[t, l, 1]
               )
    @constraint(m, o_tcp_e_d_m1_i_[t=T, l=L], # on
                o_tcp_e_d_[t, l, 1] <= p.o_cp_e_bM * y_o[t, l]
               )
    @constraint(m, o_tcp_e_d_m0_i_[t=T, l=L], # off
                o_tcp_e_d_[t, l, 0] <= p.o_cp_e_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tcp_e_s_e_[t=T, l=L],
                r_cp_e[t, l] == o_tcp_e_d_[t, l, 0] + o_tcp_e_d_[t, l, 1]
               )
    # ->  
    @constraint(m, o_ep0_d_e_[t=T, l=L],
                o_ep0_d_[t, l, 1] == o_tep0_d_[t, l, 1]
               )
    @constraint(m, o_ep0_m_i_[t=T, l=L], # on
                o_ep0_d_[t, l, 1] <= p.o_ep0_bM * y_o[t, l]
               )
    @constraint(m, o_ep0_s_[t=T, l=L],
                o_ep0[t, l] == o_ep0_d_[t, l, 1]
               )
    @constraint(m, o_tep0_d_m1_i_[t=T, l=L], # on
                o_tep0_d_[t, l, 1] <= p.o_ep0_bM * y_o[t, l]
               )
    @constraint(m, o_tep0_d_m0_i_[t=T, l=L], # off
                o_tep0_d_[t, l, 0] <= p.o_ep0_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tep0_s_e_[t=T, l=L],
                r_ep0[t, l] == o_tep0_d_[t, l, 0] + o_tep0_d_[t, l, 1]
               )
    # ->  
    @constraint(m, o_ep1ge_d_e_[t=T, l=L],
                o_ep1ge_d_[t, l, 1] == o_tep1ge_d_[t, l, 1]
               )
    @constraint(m, o_ep1ge_m_i_[t=T, l=L], # on
                o_ep1ge_d_[t, l, 1] <= p.o_ep1ge_bM * y_o[t, l]
               )
    @constraint(m, o_ep1ge_s_[t=T, l=L],
                o_ep1ge[t, l] == o_ep1ge_d_[t, l, 1]
               )
    @constraint(m, o_tep1ge_d_m1_i_[t=T, l=L], # on
                o_tep1ge_d_[t, l, 1] <= p.o_ep1ge_bM * y_o[t, l]
               )
    @constraint(m, o_tep1ge_d_m0_i_[t=T, l=L], # off
                o_tep1ge_d_[t, l, 0] <= p.o_ep1ge_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tep1ge_s_e_[t=T, l=L],
                r_ep1ge[t, l] == o_tep1ge_d_[t, l, 0] + o_tep1ge_d_[t, l, 1]
               )
    # ->
    @constraint(m, o_ep1gcs_d_e_[t=T, l=L],
                o_ep1gcs_d_[t, l, 1] == o_tep1gcs_d_[t, l, 1]
               )
    @constraint(m, o_ep1gcs_m_i_[t=T, l=L], # on
                o_ep1gcs_d_[t, l, 1] <= p.o_ep1gcs_bM * y_o[t, l]
               )
    @constraint(m, o_ep1gcs_s_[t=T, l=L],
                o_ep1gcs[t, l] == o_ep1gcs_d_[t, l, 1]
               )
    @constraint(m, o_tep1gcs_d_m1_i_[t=T, l=L], # on
                o_tep1gcs_d_[t, l, 1] <= p.o_ep1gcs_bM * y_o[t, l]
               )
    @constraint(m, o_tepgcs_d_m0_i_[t=T, l=L], # off
                o_tep1gcs_d_[t, l, 0] <= p.o_ep1gcs_bM * (1 - y_o[t, l])
               )
    @constraint(m, o_tep1gcs_s_e_[t=T, l=L],
                r_ep1gcs[t, l] == o_tep1gcs_d_[t, l, 0] 
                + o_tep1gcs_d_[t, l, 1]
               )
    # ->
    @variable(m, 
              o_loan_last[l=L, (0,1)] >= 0e0
             )
    @constraint(m, o_last_loan_s_e_[l=L], 
                r_loan_p[last(T), l] + e_loan_p[last(T), l] == 
                o_loan_last[l, 0] + o_loan_last[l, 1]
               )
    @constraint(m, o_last_loan_d_m1_i_[l=L],
                o_loan_last[l, 1] <= 1e5 * y_o[last(T), l]
               )
    @constraint(m, o_last_loan_d_m0_i_[l=L],
                o_loan_last[l, 0] <= 1e5 * (1 - y_o[last(T), l])
               )


    # 76 
    # 76 #######################################################################
    ##
    # objective
    @objective(m, Min,
               sum(
               # loan
               sum(p.discount[t+1] * o_pay[t, l] for l in L)
               for t in T)
               # o&m
               + sum(sum(p.discount[t+1] * o_conm[t, l] for l in L)
                     for t in T)
               # loan new
               + sum(sum(p.discount[t+1] * n_pay[t, l] for l in L)
                     for t in T)
               # o&m new
               + sum(sum(p.discount[t+1] * n_conm[t, l] for l in L)
                     for t in T)
               # retirement
               + sum(sum(p.discount[t+1] * t_ret_cost[t, l] for l in L)
                     for t in T)
               #+ sum(p.discount[last(T)+1] * r_loan_p[last(T), l] for l in L)
               #+ sum(p.discount[last(T)+1] * e_loan_p[last(T), l] for l in L)
               #+ sum(p.discount[last(T)+1] * n_loan_p[last(T), l] for l in L)
               + sum(p.discount[last(T)+1] * o_loan_last[l, 0] for l in L)
               # if you retire but still have unpayed loan it is gonna be
               # reflected here :()
              )
    
    # pricing-problem
    return m
end


function reattachBlockMod!(m::Model, index_l, p::parms,
        s::sets)
    L = index_l 
    
    T = s.T
    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu

    o_cp = m[:o_cp]
    n_cp = m[:n_cp]
    # aggregated demand.
    @constraint(m, cp_con[t=2:last(T)],
                sum(o_cp[t, l] for l in L) + sum(n_cp[t, l] for l in L) 
                >= p.demand[t+1]
               )
    r_ep0 = m[:r_ep0]
    n_ep0 = m[:n_ep0]
    # emissions have to be aggregated.
    #@variable(m, t_ep0[t=T])
    #@constraint(m, em_con[t=T],
    #            sum(r_ep0[t, l] for l in L) 
    #            + sum(n_ep0[t, l] for l in L)
    #            == t_ep0[t])

end

