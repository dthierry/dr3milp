 # 7n8
using JuMP
using Printf


include("sets_v01_30.jl")
include("params_v02_10.jl")

# 80 
# 80 ###########################################################################
# defining the model object
"""
    createBlockMod(index_p, index_l, p::parms, s::sets])

Create a block of constraints for period `index_p` and location `index_l`.

If `index_p` or `index_l` are collections, then this would create the range
of constraints.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function createBlockMod(index_p, index_l, p::parms, s::sets)
    # 76 
    # 76 #######################################################################
    @info "Generating the block"
    #m = Model(Cbc.Optimizer)
    m = Model()
    #set_silent(m)
    if index_p isa Int
        P = [index_p]
    else
        @printf "Set `period` passed as a collection\n"
        P = index_p
    end
    if index_l isa Int
        L = [index_l]
    else
        @printf "Set `location` passed as a collection\n"
        L = index_l
    end

    #P = s.P
    Y = s.Y

    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu

    n_periods = p.n_periods
    n_years = p.n_years
    
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
    @variable(m, y_o[i=P, j=Y, l=L], Bin)  # online
    # tier 1: retrofit variable
    @variable(m, y_r[i=P, j=Y, l=L, k=Kr], Bin)
    # tier 2: retirement
    # there is no retirement binary variable per se.  
    

    # d082923
    # expansion
    @variable(m, y_e[i=P, j=Y, l=L], Bin)  # 1 if yes expand else no
    @variable(m, e_c_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    @variable(m, e_c[i=P, j=Y, l=L])

    @variable(m, e_l[i=P, j=Y, l=L])
    @variable(m, e_l_d_[i=P, j=Y, l=L, k=(0,1)]) # disaggregated

    @variable(m, 0 <= x[l=L] <= 10) # this should be int
    @variable(m, x_d_[i=P, j=Y, l=L, (0, 1)] >= 0)
    
    # d082223
    # capacity
    @variable(m, r_cp[i=P, j=Y, l=L])
    @variable(m, r_cp_d_[i=P, j=Y, l=L, k=Kr])  # retrofit capacity disagg
    @variable(m, cpb[i=P, j=Y, l=L] >= 0.0)  # base capacity
    @variable(m, r_cpb_d_[i=P, j=Y, l=L, k=Kr] >= 0.0)  # base capacity disagg
    
    # d082923
    # heating requirement
    @variable(m, r_eh[i=P, j=Y, l=L]) # 0
    @variable(m, r_eh_d_[i=P, j=Y, l=L, k=Kr]) # 1
    # fuel requirement 
    @variable(m, r_ehf[i=P, j=Y, l=L, f=Fu]) # 2
    @variable(m, r_ehf_d_[i=P, j=Y, l=L, k=Kr, f=Fu]) # 3
    # electricity requirement (on-site)
    @variable(m, r_u[i=P, j=Y, l=L]) # 4
    @variable(m, r_u_d_[i=P, j=Y, l=L, k=Kr]) # 5
    
    # electricity emissions (off-site)
    @variable(m, r_u_onsite[i=P, j=Y, l=L]) #
    @variable(m, r_u_onsite_d_[i=P, j=Y, l=L, k=Kr]) #
    #
    # process (intrinsic) emissions
    @variable(m, r_cp_e[i=P, j=Y, l=L]) # 
    @variable(m, r_cp_e_d_[i=P, j=Y, l=L, k=Kr]) # 
    # process (total) disaggregated emissions, e.g. process + fuel, etc.
    @variable(m, r_ep0[i=P, j=Y, l=L]) # 
    @variable(m, r_ep0_d_[i=P, j=Y, l=L, k=Kr]) # 
    # scope 1, emitted
    @variable(m, r_ep1ge[i=P, j=Y, l=L]) # 
    @variable(m, r_ep1ge_d_[i=P, j=Y, l=L, k=Kr]) # 
    # scope 1, captured
    @variable(m, r_ep1gce[i=P, j=Y, l=L]) # 
    @variable(m, r_ep1gce_d_[i=P, j=Y, l=L, k=Kr]) # 
    # scope 1,  stored
    @variable(m, r_ep1gcs[i=P, j=Y, l=L]) # 
    @variable(m, r_ep1gcs_d_[i=P, j=Y, l=L, k=Kr]) # 
    # operating and maintenance
    @variable(m, r_conm_d_[i=P, j=Y, l=L, k=Kr])  # 
    @variable(m, r_conm[i=P, j=Y, l=L])  # 

    # d083023
    # loan state
    @variable(m, r_loan[i=P, j=Y, l=L])
    @variable(m, r_loan_p[i=P, j=Y, l=L] >= 0)
    @variable(m, r_loan_n[i=P, j=Y, l=L] >= 0)
    # payment
    @variable(m, r_pay[i=P, j=Y, l=L])
    @variable(m, r_pay_1[i=P, j=Y, l=L] >= 0)
    # add cost 
    @variable(m, r_ladd[i=P, j=Y, l=L])
    @variable(m, r_ladd_d_[i=P, j=Y, l=L, (0,)] >= 0)

    @variable(m, r_yps[i=P, j=Y, l=L], Bin)  # paid or not
    #
    @variable(m, r_l[i=P, j=Y, l=L]) # capital (loan)
    @variable(m, r_l_md_[i=P, j=Y, l=L, k=Kr]) # 
    @variable(m, r_l_pd_[i=P, j=Y, l=L, (0,1)])  # payd or not paid
    #
    @variable(m, r_ann[i=P, j=Y, l=L])
    @variable(m, r_ann_0[i=P, j=Y, l=L] >= 0)
    @variable(m, r_ann_1[i=P, j=Y, l=L] >= 0)
    @variable(m, r_ann_md_[i=P, j=Y, l=L, k=Kr])
    
    
    # expansion capacity
    @variable(m, e_ladd[i=P, j=Y, l=L])
    @variable(m, e_ladd_d_[i=P, j=Y, l=L, (0,)] >= 0)
    @variable(m, e_l_pd_[i=P, j=Y, l=L, (0,1)])
    @variable(m, e_yps[i=P, j=Y, l=L], Bin) # 1 if paid, 0 otw

    @variable(m, e_loan[i=P, j=Y, l=L])
    @variable(m, e_loan_p[i=P, j=Y, l=L] >= 0)
    @variable(m, e_loan_n[i=P, j=Y, l=L] >= 0)

    @variable(m, e_ann[i=P, j=Y, l=L])

    @variable(m, e_ann_d_[i=P, j=Y, l=L, (0,)])
    @variable(m, e_ann_0[i=P, j=Y, l=L]>=0)
    @variable(m, e_ann_1[i=P, j=Y, l=L]>=0)

    @variable(m, e_pay[i=P, j=Y, l=L])
    @variable(m, e_pay_1[i=P, j=Y, l=L])



    # d090423
    # retirement cost
    @variable(m, t_loan_d_[i=P, j=Y, l=L, (0,1)] >= 0)  # total loan
    @variable(m, t_ret_cost[i=P, j=Y, l=L])
    @variable(m, t_ret_cost_d_[i=P, j=Y, l=L, k=(0,)] >= 0) # only 0th needed

    # 76 
    # 76 #######################################################################
    @variable(m, y_n[i=P, j=Y, l=L, k=Kn], Bin) # this means plant kind k
    @variable(m, n_c0[l=L])  # the new capacity
    @variable(m, n_c0_d_[i=P, j=Y, l=L, k=Kn] >= 0) # disaggregated

    @variable(m, n_cp[i=P, j=Y, l=L])
    @variable(m, n_cp_d_[i=P, j=Y, l=L, k=Kn] >= 0) # disaggregated variable

    #
    @variable(m, n_l_d_[i=P, j=Y, l=L, k=Kn])
    @variable(m, n_l[i=P, j=Y, l=L])

    @variable(m, n_ann_d_[i=P, j=Y, l=L, k=Kn] >= 0)
    @variable(m, n_ann[i=P, j=Y, l=L])
    @variable(m, n_ann_0[i=P, j=Y, l=L]>=0)
    @variable(m, n_ann_1[i=P, j=Y, l=L]>=0)

    @variable(m, n_ladd_d_[i=P, j=Y, l=L, (0,)] >= 0)
    @variable(m, n_ladd[i=P, j=Y, l=L])

    @variable(m, n_l_pd_[i=P, j=Y, l=L, (0,1)] >= 0)

    @variable(m, n_loan[i=P, j=Y, l=L])
    @variable(m, n_loan_p[i=P, j=Y, l=L] >= 0)
    @variable(m, n_loan_n[i=P, j=Y, l=L] >= 0)
    @variable(m, n_yps[i=P, j=Y, l=L], Bin)

    @variable(m, n_pay[i=P, j=Y, l=L])
    @variable(m, n_pay_1[i=P, j=Y, l=L] >= 0)
    
    # 76 
    # 76 #######################################################################
    
    # heating requirement
    @variable(m, n_eh[i=P, j=Y, l=L]) # 0
    @variable(m, n_eh_d_[i=P, j=Y, l=L, k=Kn]) # 1
    # fuel requirement 
    @variable(m, n_ehf[i=P, j=Y, l=L, f=Fu]) # 2
    @variable(m, n_ehf_d_[i=P, j=Y, l=L, k=Kn, f=Fu]) # 3
    # electricity requirement
    @variable(m, n_u[i=P, j=Y, l=L]) # 4
    @variable(m, n_u_d_[i=P, j=Y, l=L, k=Kn]) # 5

    # electricity emissions (off-site)
    @variable(m, n_u_onsite[i=P, j=Y, l=L]) #
    @variable(m, n_u_onsite_d_[i=P, j=Y, l=L, k=Kn]) #

    # process (intrinsic) emissions
    @variable(m, n_cp_e[i=P, j=Y, l=L]) # 6
    @variable(m, n_cp_e_d_[i=P, j=Y, l=L, k=Kn]) # 7

    # process (extrinsic) disaggregated emissions, e.g. scope 0, etc.
    @variable(m, n_ep0[i=P, j=Y, l=L]) # 8
    @variable(m, n_ep0_d_[i=P, j=Y, l=L, k=Kn]) # 9

    # scope 1, emitted
    @variable(m, n_ep1ge[i=P, j=Y, l=L]) # 10
    @variable(m, n_ep1ge_d_[i=P, j=Y, l=L, k=Kn]) # 11
    # scope 1, captured
    @variable(m, n_ep1gce[i=P, j=Y, l=L]) # 12
    @variable(m, n_ep1gce_d_[i=P, j=Y, l=L, k=Kn]) # 13
    # scope 1,  stored
    @variable(m, n_ep1gcs[i=P, j=Y, l=L]) # 14
    @variable(m, n_ep1gcs_d_[i=P, j=Y, l=L, k=Kn]) # 15
    # operating and maintenance
    @variable(m, n_conm_d_[i=P, j=Y, l=L, k=Kn] >= 0.0)  # 16
    @variable(m, n_conm[i=P, j=Y, l=L])  # 17

    # 76 
    # 76 #######################################################################
    ##
    # k = 0 means on
    @variable(m, o_pay[i=P, j=Y, l=L])  # we put this in a disjunction
    @variable(m, o_pay_d_[i=P, j=Y, l=L, (1,)] >= 0.0)
    # k = 0 means on
    @variable(m, o_conm[i=P, j=Y, l=L])
    @variable(m, o_conm_d_[i=P, j=Y, l=L, (1,)] >= 0.0)
    @variable(m, o_tconm_d_[i=P, j=Y, l=L, (0,1)])
    #
    @variable(m, o_cp[i=P, j=Y, l=L])
    @variable(m, o_cp_d_[i=P, j=Y, l=L, (1,)] >= 0) # disaggregated variable

    @variable(m, o_tcp_d_[i=P, j=Y, l=L, (0,1)] >= 0)

    @variable(m, o_tpay_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    #
    @variable(m, o_u_d_[i=P, j=Y, l=L, (1,)] >= 0)
    @variable(m, o_u[i=P, j=Y, l=L] >= 0)

    @variable(m, o_tu_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    #
    @variable(m, o_ep0_d_[i=P, j=Y, l=L, (1,)] >= 0)
    @variable(m, o_ep0[i=P, j=Y, l=L])

    @variable(m, o_tep0_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    #
    @variable(m, o_ep1ge_d_[i=P, j=Y, l=L, (1,)] >= 0)
    @variable(m, o_ep1ge[i=P, j=Y, l=L])

    @variable(m, o_tep1ge_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    #
    @variable(m, o_ep1gcs[i=P, j=Y, l=L]) # 14
    @variable(m, o_ep1gcs_d_[i=P, j=Y, l=L, k=Kr]) # 15

    @variable(m, o_tep1gcs_d_[i=P, j=Y, l=L, (0,1)] >= 0)
    #
    # -> had to put a period here
    @variable(m, 
              o_loan_last[i=P, l=L, (0,1)] >= 0e0
             )


    # 76 
    # 76 #######################################################################
    ##
    # tier 0 logic
    @constraint(m, o_logic_1_y_i_[i=P, j=Y, l=L; j<n_years],
                y_o[i, j+1, l] <= y_o[i, j, l]) # this can only go offline

    #if index_p isa Int && index_p == 0
        @constraint(m, o_logic_init_0_[l=L],
                    y_o[0, 0, l] == 1)  # plants must start online 
    #end


    # tier 1
    #@constraint(m, logic_tier01_0_e[i=P, j=Y, l=L],  # only one option
    #            sum(y_r[i, j, l, k] for k in Kr) >= y_o[i, j, l])

    #@constraint(m, logic_tier01_1_e[i=P, j=Y, l=L, k=Kr],
    #            y_o[i, j, l] >= y_r[i, j, l, k]
    #           )

    #
    #if index_p isa Int && index_p == 0
        @constraint(m, r_logic_init_0[l=L], # all non 0 rtrf
                    y_r[0, 0, l, 0] == 1
                   )
    #end
    #
    #if index_p isa Int && index_p == 0
        @constraint(m, r_logic_init_1[l=L, k=Kr; k > 0], # all non 0 rtrf
                    y_r[0, 0, l, k] == 0
               )
    #end
    #
    @constraint(m, r_logic_budget_s[i=P, j=Y, l=L;
                                    (i>0 && j>0)], # only one mode
                sum(y_r[i, j, l, k] for k in Kr if k > 0) <= 1
               )
    #
    @constraint(m, logic_tier_1_0m_e_[i=P, j=Y, l=L, k=Kr; 
                                      k>0], # either 0 or >0 
                1 >= y_r[i, j, l, k] + y_r[i, j, l, 0]
               )
    #
    @constraint(m, r_logic_budget_y_i_[i=P, j=Y, l=L, k=Kr; 
                                       k>0 && j<n_years],
                y_r[i, j+1, l, k] >= y_r[i, j, l, k]
               )

    @constraint(m, r_logic_onoff_1_y_i_[i=P, j=Y, l=L, k=Kr;
                                        j<n_years],
                y_o[i, j, l] + 1 - y_r[i, j, l, k] 
                + y_r[i, j+1, l, k] >= 1
               )

    @constraint(m, r_logic_onoff_2_y_i_[i=P, j=Y, l=L, k=Kr;
                                        j<n_years],
                y_o[i, j, l] + 1 - y_r[i, j+1, l, k] 
                + y_r[i, j, l, k] >= 1
               )

    # continuity


    # 76 
    # 76 #######################################################################
    ##
    # d082923
    # -> expansion
    # p.e_C[l], the capacity per unit of allocation
    @constraint(m, exp_d0_e_[i=P, j=Y, l=L], # only 0 counts
                e_c_d_[i, j, l, 0] == p.e_C[l+1] * x_d_[i, j, l, 0]
               )
    # e_c_bM
    @constraint(m, exp_ncw_m_i0_[i=P, j=Y, l=L],
                e_c_d_[i, j, l, 0] <= p.e_c_bM * y_e[i, j, l]
               )
    @constraint(m, exp_ncw_s_[i=P, j=Y, l=L],
                e_c[i, j, l] == e_c_d_[i, j, l, 0]
               )
    # x_bM
    @constraint(m, exp_x_m_i0_[i=P, j=Y, l=L],
                x_d_[i, j, l, 0] <= p.x_bM * y_e[i, j, l]
               )
    # x_bM (relax)
    @constraint(m, exp_x_m_i1_[i=P, j=Y, l=L],
                x_d_[i, j, l, 1] <= p.x_bM * (1 - y_e[i, j, l])
               )
    @constraint(m, exp_x_s_[i=P, j=Y, l=L],
                x[l] == x_d_[i, j, l, 0] + x_d_[i, j, l, 1]
               )

    # -> expansion cost (loan)
    # p.e_loanFact
    @constraint(m, e_l_d0_e_[i=P, j=Y, l=L], # only zero counts
                e_l_d_[i, j, l, 0] == p.e_loanFact[l+1] * x_d_[i, j, l, 0]
               )
    #@constraint(m, e_l_d1_e_[i=P, j=Y, l=L], # only zero counts
    #            e_l_d_[i, j, l, 1] == 0.0 
    #           )
    
    # e_l_bM
    @constraint(m, e_l_m_i0_[i=P, j=Y, l=L],
                e_l_d_[i, j, l, 0] <= p.e_l_bM * y_e[i, j, l]
               )
    # e_l_bM
    #@constraint(m, e_l_m_i1_[i=P, j=Y, l=L],
    #            e_l_d_[i, j, l, 1] <= e_l_bM * (1 - y_e[i, j, l])
    #           )
    @constraint(m, e_l_s_[i=P, j=Y, l=L],
                #e_l[i, j, l] == sum(e_l_d_[i, j, l, k] for k in (0,1))
                e_l[i, j, l] == e_l_d_[i, j, l, 0]
               )

    # -> expansion cost annuity (annual payment)
    # p.e_Ann
    @constraint(m, e_ann_d0_e_[i=P, j=Y, l=L], # only zero counts
                e_ann_d_[i, j, l, 0] == p.e_Ann[l+1] * x_d_[i, j, l, 0]
               )
    # e_ann_bM
    @constraint(m, e_ann_m_i0_[i=P, j=Y, l=L],
                e_ann_d_[i, j, l, 0] <= p.e_ann_bM * y_e[i, j, l]
               )
    @constraint(m, e_ann_s_[i=P, j=Y, l=L],
                e_ann[i, j, l] == e_ann_d_[i, j, l, 0]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion logic
    #if index_p isa Int && index_p == 0
        @constraint(m, e_logic_init[l=L], 
                    y_e[0, 0, l] == 0  # start not-expanded
                   )
    #end
    @constraint(m, e_logic_1_y_i_[i=P, j=Y, l=L; j<n_years],
                y_e[i, j+1, l] >= y_e[i, j, l]  # only expand in the future
               )

    # 76 
    # 76 #######################################################################
    ##
    # capacity expansion loans (they should be agnostic to retrofit or rf)
    # components: e_ladd, e_ann, e_pslack, e_ploan
    # -> e_add
    # exp add loan
    @constraint(m, e_ladd_d0_e_[i=P, j=Y, l=L],
                e_ladd_d_[i, j, l, 0] == e_l_pd_[i, j, l, 0]
               )
    @constraint(m, e_ladd_m_0_y_i_[i=P, j=Y, l=L; j<n_years], 
                # y_e goes from 0 to 1
                e_ladd_d_[i, j, l, 0] <= 
                p.e_ladd_bM * (y_e[i, j+1, l] - y_e[i, j, l])
               )
    # 
    @constraint(m, e_l_m_1_y_i_[i=P, j=Y, l=L; j<n_years], 
                # need to set this to 0
                e_l_pd_[i, j, l, 1] <= 
                p.e_l_bM * (1 - y_e[i, j+1, l] + y_e[i, j, l])
               )  #  the 0th component is implied by the e_ladd_m constr
    # 

    @constraint(m, e_ladd_s_e[i=P, j=Y, l=L],
                e_ladd[i, j, l] == e_ladd_d_[i, j, l, 0]
               )

    @constraint(m, e_l_s_e_[i=P, j=Y, l=L; j<n_years],
                e_l[i, j+1, l] == 
                e_l_pd_[i, j, l, 0] + e_l_pd_[i, j, l, 1]
               )

    # load
    @constraint(m, e_loan_s_e_[i=P, j=Y, l=L],
                e_loan[i, j, l] == e_loan_p[i, j, l] - e_loan_n[i, j, l]
               )
    @constraint(m, e_loan_p_m0_i_[i=P, j=Y, l=L],
                e_loan_p[i, j, l] <= p.e_loan_bM * (1 - e_yps[i, j, l])
               )
    @constraint(m, e_loan_n_m0_i_[i=P, j=Y, l=L],
                e_loan_n[i, j, l] <= p.e_loan_bM * e_yps[i, j, l]
               )
    # pay
    @constraint(m, e_pay_s_e_[i=P, j=Y, l=L],
                e_pay[i, j, l] == e_pay_1[i, j, l]
               )
    @constraint(m, e_pay_n_m0_i_[i=P, j=Y, l=L],
                e_pay_1[i, j, l] <= p.e_pay_bM * (1 - e_yps[i, j, l])
               )
    # annuity
    @constraint(m, e_ann_s_e_[i=P, j=Y, l=L],
                e_ann[i, j, l] == e_ann_0[i, j, l] + e_ann_1[i, j, l]
               )
    @constraint(m, e_ann_0_m_i_[i=P, j=Y, l=L],
                e_ann_0[i, j, l] <= p.e_ann_bM * e_yps[i, j, l]
               )
    @constraint(m, e_ann_1_m_i_[i=P, j=Y, l=L],
                e_ann_1[i, j, l] <= p.e_ann_bM * (1 - e_yps[i, j, l])
               )
    @constraint(m, e_ann_dpay0_e_[i=P, j=Y, l=L],
                e_ann_1[i, j, l] == e_pay_1[i, j, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    @constraint(m, e_loan_bal_y_e_[i=P, j=Y, l=L; j<n_years],
                e_loan[i, j+1, l] == e_loan[i, j, l]
                - e_pay[i, j, l]
                + e_ladd[i, j, l]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    # d082923
    # base capacity
    @constraint(m, cpb_e_[i=P, j=Y, l=L],
                cpb[i, j, l] == p.c0[l+1] + e_c[i, j, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    # retrofit 
    # -> capacity.
    # p.Kr (capacity factor, e.g. prod in retrofit r / base prod)
    # viz. cap_mod = factor * cap_base
    @constraint(m, r_cp_d_e_[i=P, j=Y, l=L, k=Kr],
                r_cp_d_[i, j, l, k] == p.r_c_C[l+1, k+1] * r_cpb_d_[i, j, l, k] 
                + p.r_rhs_C[l+1, k+1] * y_r[i, j, l, k]
               )
    # r_cp_bM, retrofit capacity
    @constraint(m, r_cp_d_m_e_[i=P, j=Y, l=L, k=Kr],
                r_cp_d_[i, j, l, k] <= p.r_cp_bM * y_r[i, j, l, k]
               )
    # cpb_bM, base capacity
    @constraint(m, r_cpb_d_m_i_[i=P, j=Y, l=L, k=Kr],
                r_cpb_d_[i, j, l, k] <= p.r_cpb_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_cp_s_e_[i=P, j=Y, l=L],
                r_cp[i, j, l] == sum(r_cp_d_[i, j, l, k] for k in Kr)
               )
    @constraint(m, r_cpb_s_e_[i=P, j=Y, l=L],
                cpb[i, j, l] == sum(r_cpb_d_[i, j, l, k] for k in Kr)
               )
    # -> heating requirement.
    # p.Hm (heating factor, i.e. heat / product)
    @constraint(m, r_eh_d_e_[i=P, j=Y, l=L, k=Kr],
                r_eh_d_[i, j, l, k] == 
                (1-p.r_c_Helec[l+1, k+1]) * p.r_c_Hfac[l+1, k+1] *
                (p.r_c_H[l+1, k+1] * r_cp_d_[i, j, l, k] + 
                 p.r_rhs_H[l+1, k+1] * y_r[i, j, l, k])
               )
    @constraint(m, r_eh_d_m_i_[i=P, j=Y, l=L, k=Kr],
                r_eh_d_[i, j, l, k] <= p.r_eh_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_eh_s_e_[i=P, j=Y, l=L],
                r_eh[i, j, l] == sum(r_eh_d_[i, j, l, k] for k in Kr)
               )
    # d082923
    # -> fuel required for heat.
    # p.r_c_F in [0,1], (i.e. a fraction, heat by fuel /tot heat)
    @constraint(m, r_ehf_d_e_[i=P, j=Y, l=L, k=Kr, f=Fu],
                r_ehf_d_[i, j, l, k, f] ==
                p.r_c_F[l+1, f+1, k+1] * r_eh_d_[i, j, l, k] 
                + p.r_rhs_F[l+1, f+1, k+1] * y_r[i, j, l, k]
               )
    # r_ehf_bM
    @constraint(m, r_ehf_m_i_[i=P, j=Y, l=L, k=Kr, f=Fu],
                r_ehf_d_[i, j, l, k, f] <= p.r_ehf_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_ehf_s_e_[i=P, j=Y, l=L, f=Fu],
                r_ehf[i, j, l, f] == sum(r_ehf_d_[i, j, l, k, f] for k in Kr)
               )

    # -> electricity requirement
    # r_u and m_ud_, p.Um & p.UmRhs
    @constraint(m, r_u_d_e_[i=P, j=Y, l=L, k=Kr],
                r_u_d_[i, j, l, k] == 
                (1 - p.r_c_UonSite[i+1, j+1, l+1, k+1])* p.r_c_Ufac[l+1, k+1] *
                (p.r_c_U[l+1, k+1] * r_cp_d_[i, j, l, k] 
                 + p.r_rhs_U[l+1,k+1]*y_r[i,j,l,k]) +
                # electrification
                p.r_c_Helec[l+1, k+1] * p.r_c_Hfac[l+1, k+1] * (
                                        p.r_c_H[l+1, k+1] * r_cp_d_[i, j, l, k] 
                                        + p.r_rhs_H[l+1, k+1] * y_r[i, j, l, k])
               )

    # r_u_bM (big-M)
    @constraint(m, r_u_i_[i=P, j=Y, l=L, k=Kr],
                r_u_d_[i, j, l, k] <= p.r_u_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_u_s_e_[i=P, j=Y, l=L],
                r_u[i, j, l] == sum(r_u_d_[i, j, l, k] for k in Kr)
               )

    # -> electricity on-site generation
    @constraint(m, r_u_onsite_d_e_[i=P, j=Y, l=L, k=Kr],
                r_u_onsite_d_[i, j, l, k] == 
                p.r_c_UonSite[i+1, j+1, l+1, k+1] * p.r_c_Ufac[l+1, k+1]*
                (p.r_c_U[l+1, k+1] * r_cp_d_[i, j, l, k] 
                #+ p.r_c_UonSite[i+1, j+1, l+1, k+1]
                * p.r_rhs_U[l+1,k+1]*y_r[i,j,l,k])
               )
    # r_u_bM (big-M)
    @constraint(m, r_u_onsite_i_[i=P, j=Y, l=L, k=Kr],
                r_u_onsite_d_[i, j, l, k] <= p.r_u_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_u_onsite_s_e_[i=P, j=Y, l=L],
                r_u_onsite[i, j, l] == 
                sum(r_u_onsite_d_[i, j, l, k] for k in Kr)
               )

    # -> process (intrinsic) emissions
    # r_cp_e, r_cp_e_d_. p.Cp & p.CpRhs
    @constraint(m, r_cp_e_d_e_[i=P, j=Y, l=L, k=Kr],
                r_cp_e_d_[i, j, l, k] == 
                p.r_c_cp_e[l+1, k+1] * r_cp_d_[i, j, l, k]
                + p.r_rhs_cp_e[l+1, k+1] * y_r[i, j, l, k]
               )
    # r_cp_e_bM
    @constraint(m, r_cp_e_m_i_[i=P, j=Y, l=L, k=Kr],
                r_cp_e_d_[i, j, l, k] <= p.r_cp_e_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_cp_e_s_e_[i=P, j=Y, l=L],
                r_cp_e[i, j, l] == sum(r_cp_e_d_[i, j, l, k] for k in Kr))
    

    # -> electricity emissions
    # GcI : grid carbon intensity
    #@constraint(m, r_u_em_e_[i=P, j=Y, l=L], 
    #            r_u_em[i, j, l] == p.GcI[i, j, l] * r_u[i, j, l]
    #           )
    # -> -> process (disaggregated) emissions

    # -> scope 0 emission
    # c_Fe (fuel emission factor) r_Hr (heat rate)
    @constraint(m, r_ep0_d_e_[i=P, j=Y, l=L, k=Kr],
                r_ep0_d_[i, j, l, k] == 
                # fuel
                sum(p.r_c_Fe[f+1, k+1]*r_ehf_d_[i, j, l, k, f] for f in Fu) 
                # process
                + r_cp_e_d_[i, j, l, k]
                # in-site electricity
                + sum(p.r_c_Fe[f+1, k+1]
                    *p.r_c_Hr[l+1, f+1, k+1]
                    *p.r_c_Fgenf[l+1, f+1, k+1]
                    for f in Fu) * r_u_onsite_d_[i, j, l, k] 
               )
    # r_ep0_bM
    @constraint(m, r_ep0_m_i_[i=P, j=Y, l=L, k=Kr],
                r_ep0_d_[i, j, l, k] <= p.r_ep0_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_ep0_s_e_[i=P, j=Y, l=L],
                r_ep0[i, j, l] == sum(r_ep0_d_[i, j, l, k] for k in Kr)
               )

    # -> scope 1 emitted
    # p.r_chi
    @constraint(m, r_ep1ge_d_e_[i=P, j=Y, l=L, k=Kr],
                r_ep1ge_d_[i, j, l, k] == (1.0 - p.r_chi[l+1, k+1]) 
                * r_ep0_d_[i, j, l, k]
               )
    # r_ep1ge_bM
    @constraint(m, r_ep1ge_m_i_[i=P, j=Y, l=L, k=Kr],
                r_ep1ge_d_[i, j, l, k] <= p.r_ep1ge_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_ep1ge_s_e_[i=P, j=Y, l=L],
                r_ep1ge[i, j, l] == sum(r_ep1ge_d_[i, j, l, k] for k in Kr)
               )

    # -> scope 1 captured
    # p.r_sigma
    @constraint(m, r_ep1gce_d_e_[i=P, j=Y, l=L, k=Kr],
                r_ep1gce_d_[i, j, l, k] == 
                p.r_chi[l+1, k+1] * (1 - p.r_sigma[l+1, k+1]) 
                * r_ep0_d_[i, j, l, k]
               )
    @constraint(m, r_ep1gce_m_i_[i=P, j=Y, l=L, k=Kr],
                r_ep1gce_d_[i, j, l, k] <= p.r_ep1gce_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_ep1gce_s_e_[i=P, j=Y, l=L],
                r_ep1gce[i, j, l] == sum(r_ep1gce_d_[i, j, l, k] for k in Kr)
               )

    # -> scope 1 stored
    # p.sigma ?
    @constraint(m, r_ep1gcs_d_e_[i=P, j=Y, l=L, k=Kr],
                r_ep1gcs_d_[i, j, l, k] ==
                p.r_chi[l+1, k+1] * p.r_sigma[l+1, k+1] * r_ep0_d_[i, j, l, k]
               )
    # ep1gcsm_bM
    @constraint(m, r_ep1gcs_m_i_[i=P, j=Y, l=L, k=Kr],
                r_ep1gcs_d_[i, j, l, k] <= p.r_ep1gcs_bM * y_r[i, j, l, k]
               )

    @constraint(m, r_ep1gcs_s_e_[i=P, j=Y, l=L],
                r_ep1gcs[i, j, l] == sum(r_ep1gcs_d_[i, j, l, k] for k in Kr)
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
    @constraint(m, r_conm_d_e_[i=P, j=Y, l=L, k=Kr],
                r_conm_d_[i, j, l, k] == p.r_c_Onm[l+1, k+1] * r_cp_d_[i, j, l, k]
                + p.r_rhs_Onm[l+1, k+1] * y_r[i, j, l, k]
               )
    @constraint(m, r_conm_m_i_[i=P, j=Y, l=L, k=Kr],
                r_conm_d_[i, j, l, k] <= p.r_conm_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_conm_s_e_[i=P, j=Y, l=L],
                r_conm[i, j, l] == sum(r_conm_d_[i, j, l, k] for k in Kr)
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> retrofit disagg LOAN added amount (ladd), p.r_loanFact
    @constraint(m, r_l_md_e_[i=P, j=Y, l=L, k=Kr], # retrofit disagg
                r_l_md_[i, j, l, k] == 
                p.r_loanFact[i+1, j+1, l+1, k+1]*r_cp_d_[i, j, l, k]
               ) # uses c0+expc
    @constraint(m, r_l_mm_i_[i=P, j=Y, l=L, k=Kr], # big-M
                r_l_md_[i, j, l, k] <= p.r_l_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_l_ms_e_[i=P, j=Y, l=L],
                r_l[i, j, l] == sum(r_l_md_[i, j, l, k] for k in Kr)
               )
    # 76 
    # 76 #######################################################################
    ##
    # -> retrofit disagg payment (associated with the payment)
    @constraint(m, r_ann_md_e_[i=P, j=Y, l=L, k=Kr],
                r_ann_md_[i, j, l, k] == 
                p.r_annf[i+1, j+1, l+1, k+1]*r_cp_d_[i, j, l, k]
               ) # uses c0+expc
    @constraint(m, r_ann_mm_i_[i=P, j=Y, l=L, k=Kr],
                r_ann_md_[i, j, l, k] <= p.r_ann_bM * y_r[i, j, l, k]
               )
    @constraint(m, r_ann_s_e_[i=P, j=Y, l=L],
                r_ann[i, j, l] == sum(r_ann_md_[i, j, l, k] for k in Kr)
               )

    # 76 
    # 76 #######################################################################
    ##
    # --> padd switch (loan-switch), associated with the loan
    @constraint(m, r_ladd_d_e_[i=P, j=Y, l=L], 
                r_ladd_d_[i, j, l, 0] == r_l_pd_[i, j, l, 0]
               )
    # p_add_bM
    @constraint(m, r_ladd_m_0_y_i_[i=P, j=Y, l=L; j<n_years],
                r_ladd_d_[i, j, l, 0] <= # goes from 1 to 0 only
                p.r_ladd_bM * (y_r[i, j, l,0] - y_r[i,j+1,l,0])
               )

    # (zero otw.)
    #@constraint(m, m_ladd_m_i1_[i=P, j=Y, l=L],
    #            r_ladd_d_[i, j, l, 1] <=
    #            maximum(r_ladd_bM[i, j, l, k])*(1-y_r[t-1,l,0]+y_r[t,l,0])
    #           )
    #
    @constraint(m, r_ladd_s_e_[i=P, j=Y, l=L],
                r_ladd[i, j, l] == 
                #sum(r_ladd_d_[i, j, l, k] for k in (0, 1))
                r_ladd_d_[i, j, l, 0]
               )
    # padd switch (loan-retrofit)
    #@constraint(m, r_l_m_0_y_i_[i=P, j=Y, l=L; j<n_years],
    #            r_l_pd_[i, j, l, 0] <=
    #            p.r_l_bM * (y_r[i,j,l,0] - y_r[i,j+1,l,0])
    #           )
    #
    @constraint(m, r_l_m_1_y_i_[i=P, j=Y, l=L; j<n_years],
                r_l_pd_[i, j, l, 1] <=
                #p.r_l_bM * (1.0 - y_r[t-1,l,0] + y_r[t,l,0])
                p.r_l_bM * (1.0 - y_r[i,j,l,0] + y_r[i, j+1,l,0])
               )
    # connection to the retrofit 
    @constraint(m, m_l_s_e_[i=P, j=Y, l=L; j<n_years], #  here
                r_l[i, j+1, l] == sum(r_l_pd_[i, j, l, k] for k in (0, 1))
               )
    
    # 76 
    # 76 #######################################################################
    ##
    # -> loan disaggregation, r_yps = 1 if paid, 0 otw
    @constraint(m, r_loan_ps_s_e_[i=P, j=Y, l=L],
                r_loan[i, j, l] == r_loan_p[i, j, l] - r_loan_n[i, j, l]
               )
    @constraint(m, r_loan_p_m0_i_[i=P, j=Y, l=L],
                r_loan_p[i, j, l] <= p.r_loan_bM * (1 - r_yps[i, j, l])
               )
    @constraint(m, r_loan_n_m0_i_[i=P, j=Y, l=L],
                r_loan_n[i, j, l] <= p.r_loan_bM * r_yps[i, j, l]
               )
    # payment
    @constraint(m, r_pay_s_e_[i=P, j=Y, l=L],
                r_pay[i, j, l] == r_pay_1[i, j, l]
               )
    @constraint(m, r_pay_1_m0_i_[i=P, j=Y, l=L],
                r_pay_1[i, j, l] <= p.r_pay_bM * (1 - r_yps[i, j, l])
               )
    # annuity
    @constraint(m, r_ann_01_s_e_[i=P, j=Y, l=L],
                r_ann[i, j, l] == r_ann_0[i, j, l] + r_ann_1[i, j, l]
               )
    @constraint(m, r_ann_0_m_i_[i=P, j=Y, l=L],
                r_ann_0[i, j, l] <= p.r_ann_bM * r_yps[i, j, l]
               )
    @constraint(m, r_ann_1_m_i_[i=P, j=Y, l=L],
                r_ann_1[i, j, l] <= p.r_ann_bM * (1 - r_yps[i, j, l])
               )
    @constraint(m, r_ann_d0_e_[i=P, j=Y, l=L],
                r_ann_1[i, j, l] == r_pay_1[i, j, l]
               )

    # if the plant becomes retired r_pay_p might still be positive
    # 76 
    # 76 #######################################################################
    ##
    # -> loan balance
    @constraint(m, r_loan_bal_y_e_[i=P, j=Y, l=L; j<n_years],
                r_loan[i, j+1, l] == r_loan[i, j, l] 
                - r_pay[i, j, l] 
                + r_ladd[i, j, l]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> retirement
    # retirement is just r_loan, we just need a way to activate it. 
    # a) have a switch using y_o going from 0 to 1
    # b) take the snapshot of the current value of r_loan and use it as the cost
    #
    @constraint(m, t_ret_c_d_e_[i=P, j=Y, l=L],  # only enforceable at the switch
                t_ret_cost_d_[i, j, l, 0] == t_loan_d_[i, j, l, 0]
               )

    # t_ret_c_bM
    @constraint(m, t_ret_c_bm_0_y_i_[i=P, j=Y, l=L; j<n_years],
                t_ret_cost_d_[i, j, l, 0] <= 
                #p.t_ret_c_bM * (y_o[t-1, l] - y_o[i, j, l])
                p.t_ret_c_bM * (y_o[i, j, l] - y_o[i, j+1, l])
               )

    # m_loan_d_bM
    @constraint(m, r_loan_d_bm_0_y_i_[i=P, j=Y, l=L; j<n_years],
                t_loan_d_[i, j, l, 0] <=  # retired
                #p.t_loan_bM * (y_o[t-1, l] - y_o[i, j, l])
                p.t_loan_bM * (y_o[i, j, l] - y_o[i,j+1, l])
               )
    # m_loan_d_bM
    @constraint(m, r_loan_d_bm_1_y_i_[i=P, j=Y, l=L; j<n_years], # not retired
                t_loan_d_[i, j, l, 1] <= p.t_loan_bM * 
                #(1 + y_o[i, j, l] - y_o[t-1, l])
                (1 + y_o[i,j+1, l] - y_o[i, j, l])
               )
    @constraint(m, t_ret_c_s_e_[i=P, j=Y, l=L],
                t_ret_cost[i, j, l] == t_ret_cost_d_[i, j, l, 0]  
                # only one needed
               )
               
    @constraint(m, r_loan_s_e_[i=P, j=Y, l=L],  # total loan
                r_loan_p[i, j, l] + e_loan_p[i, j, l]
                == t_loan_d_[i, j, l, 0] + t_loan_d_[i, j, l, 1]
               )
    # -> total payment (retrofit/existing + )
    @constraint(m, o_pay_s_e_[i=P, j=Y, l=L],
                o_pay[i, j, l] == o_pay_d_[i, j, l, 1]
               )
    @constraint(m, o_pay_d1_e_[i=P, j=Y, l=L], # online
                # there is pay for retrof but not for expansion
                o_pay_d_[i, j, l, 1] == o_tpay_d_[i, j, l, 1]
               )
    @constraint(m, o_pay_m_1_i_[i=P, j=Y, l=L], # on
                o_pay_d_[i, j, l, 1] <= p.o_pay_bM * y_o[i, j, l]
               )
    # this one includes both expansion and retrof
    @constraint(m, o_tpay_s_e_[i=P, j=Y, l=L],
                r_pay[i, j, l] + e_pay[i, j, l] == 
                o_tpay_d_[i, j, l, 0] + o_tpay_d_[i, j, l, 1]
               )
    #
    @constraint(m, o_tpay_m1_i_[i=P, j=Y, l=L], # on
                o_tpay_d_[i, j, l, 1] <= p.o_pay_bM * y_o[i, j, l]
               )
    #
    @constraint(m, o_tpay_m0_i_[i=P, j=Y, l=L], # off
                o_tpay_d_[i, j, l, 0] <= p.o_pay_bM * (1 - y_o[i, j, l])
               )
    # -> o&m
    @constraint(m, o_conm_s_e_[i=P, j=Y, l=L],
                o_conm[i, j, l] == o_conm_d_[i, j, l, 1]
               )
    @constraint(m, o_conm_d1_e_[i=P, j=Y, l=L],
                o_conm_d_[i, j, l, 1] == o_tconm_d_[i, j, l, 1]
               )

    @constraint(m, o_conm_m_1_i_[i=P, j=Y, l=L], # on
                o_conm_d_[i, j, l, 1] <= p.o_conm_bM * y_o[i, j, l]
               )
    #
    @constraint(m, o_tconm_m1_i_[i=P, j=Y, l=L],
                o_tconm_d_[i, j, l, 1] <= p.o_conm_bM * y_o[i, j, l]
               )
    @constraint(m, o_tconm_m0_i_[i=P, j=Y, l=L],
                o_tconm_d_[i, j, l, 0] <= p.o_conm_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tconm_s_e_[i=P, j=Y, l=L],
                r_conm[i, j, l] == o_tconm_d_[i, j, l, 0] 
                + o_tconm_d_[i, j, l, 1]
               )


    # 76 
    # 76 #######################################################################
    ##
    # -> new plant capacity
    #@constraint(m, n_c0_e_[l=L],
    #            n_c0[l] >= p.c0[l+1]
    #           )
    # -> new plant capacity disaggregation
    @constraint(m, n_cp0_d_e_[i=P, j=Y, l=L],
                n_cp_d_[i, j, l, 0] == 0
               ) # how do we make this 0 at k=0?
    @constraint(m, n_cp_d_e_[i=P, j=Y, l=L, k=Kn; k>0],
                n_cp_d_[i, j, l, k] == n_c0_d_[i, j, l, k]
               )
    #@constraint(m, n_c0_d0_i_[i=P, j=Y, l=L, k=Kn; k>0],
    #            n_c0_d_[i, j, l, k] >= p.c0[l+1] * y_n[i, j, l, k]
    #           )

    @constraint(m, n_cp_m_i_[i=P, j=Y, l=L, k=Kn; k>0], # skip the 0-th
                n_cp_d_[i, j, l, k] <= p.n_cp_bM * y_n[i, j, l, k]
               )

    @constraint(m, n_cap_add_d_lo_i_[i=P, j=Y, l=L, k=Kn],
                n_c0_d_[i, j, l, k] >= p.n_c0_lo[l+1] * y_n[i, j, l, k]
               )
    @constraint(m, n_cap_add_d_m_i_[i=P, j=Y, l=L, k=Kn],
                n_c0_d_[i, j, l, k] <= p.n_c0_bM[l+1] * y_n[i, j, l, k]
               )
    #
    @constraint(m, n_cp_s_e_[i=P, j=Y, l=L],
                n_cp[i, j, l] == sum(n_cp_d_[i, j, l, k] for k in Kn)
               )

    @constraint(m, n_c0_s_e_[i=P, j=Y, l=L],
                n_c0[l] == sum(n_c0_d_[i, j, l, k] for k in Kn)
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> base loan (proportional to the capacity)
    @constraint(m, n_l_d0_e_[i=P, j=Y, l=L, k=Kn],
                n_l_d_[i, j, l, k] == 
                p.n_loanFact[l+1, k+1] * n_cp_d_[i, j, l, k]
               ) # this should be 0 at 0th  
    # perhaps this should be n_c0_d_
    @constraint(m, n_l_m_i0_[i=P, j=Y, l=L, k=Kn],
                n_l_d_[i, j, l, k] <= p.n_l_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_l_s_[i=P, j=Y, l=L],
                n_l[i, j, l] == sum(n_l_d_[i, j, l, k] for k in Kn)
               )
    # -> annuity (how much we pay)
    @constraint(m, n_ann_d0_e_[i=P, j=Y, l=L, k=Kn],
                n_ann_d_[i, j, l, k] == p.n_Ann[l+1, k+1] * n_cp_d_[i, j, l, k]
               )
    @constraint(m, n_ann_m_i0_[i=P, j=Y, l=L, k=Kn],
                n_ann_d_[i, j, l, k] <= p.n_ann_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_ann_s_[i=P, j=Y, l=L],
                n_ann[i, j, l] == sum(n_ann_d_[i, j, l, k] for k in Kn)
               )
    
    # 76 
    # 76 #######################################################################
    ##
    # -> 
    # link the n_cost_ to the ladd in a single time period.
    @constraint(m, n_ladd_m_0_y_i_[i=P, j=Y, l=L; j<n_years], 
                n_ladd_d_[i, j, l, 0] <= 
                p.n_ladd_bM * (y_n[i, j, l, 0] - y_n[i, j+1, l, 0])
               ) # y_n goes from 1 to 0

    #@constraint(m, n_l_pd_m_0_y_i_[i=P, j=Y, l=L; j<n_years],
    #            n_l_pd_[i, j, l, 0] <= 
    #            p.n_l_bM * (y_n[i,j, l, 0] - y_n[i, j+1, l, 0])
    #           )

    @constraint(m, n_l_pd_m_1_y_i_[i=P, j=Y, l=L; j<n_years],
                n_l_pd_[i, j, l, 1] <= 
                #p.n_l_bM * (1 + y_n[i, j, l, 0] - y_n[t-1, l, 0])
                p.n_l_bM * (1 + y_n[i,j+1, l, 0] - y_n[i, j, l, 0])
               )
    #
    @constraint(m, n_ladd_d0_e_[i=P, j=Y, l=L],
                n_ladd_d_[i, j, l, 0] == n_l_pd_[i, j, l, 0]
               )
    @constraint(m, n_ladd_s_e_[i=P, j=Y, l=L],
                n_ladd[i, j, l] == n_ladd_d_[i, j, l, 0]
               )
    @constraint(m, n_l_s_e_[i=P, j=Y, l=L; j<n_years],
                n_l[i, j+1, l] ==
                n_l_pd_[i, j, l, 0] + n_l_pd_[i, j, l, 1]
               )

    # loan
    @constraint(m, n_loan_s_e_[i=P, j=Y, l=L],
                n_loan[i, j, l] == n_loan_p[i, j, l] - n_loan_n[i, j, l]
               )
    @constraint(m, n_loan_p_m0_i_[i=P, j=Y, l=L],
                n_loan_p[i, j, l] <= p.n_loan_bM * (1 - n_yps[i, j, l])
               )
    @constraint(m, n_loan_n_m0_i_[i=P, j=Y, l=L],
                n_loan_n[i, j, l] <= p.n_loan_bM * n_yps[i, j, l]
               )
    # pay
    @constraint(m, n_pay_s_e_[i=P, j=Y, l=L],
                n_pay[i, j, l] == n_pay_1[i, j, l]
               )
    # note: this mechanism seems to have a lag of 1 period.
    @constraint(m, n_pay_n_m0_i_[i=P, j=Y, l=L],
                n_pay_1[i, j, l] <= p.n_pay_bM * (1 - n_yps[i, j, l])
               )
    # annuity
    @constraint(m, n_ann_s_e_[i=P, j=Y, l=L],
                n_ann[i, j, l] == n_ann_0[i, j, l] + n_ann_1[i, j, l]
               )
    @constraint(m, n_ann_0_m_i_[i=P, j=Y, l=L],
                n_ann_0[i, j, l] <= p.n_ann_bM * n_yps[i, j, l]
               )
    @constraint(m, n_ann_1_m_i_[i=P, j=Y, l=L],
                n_ann_1[i, j, l] <= p.n_ann_bM * (1 - n_yps[i, j, l])
               )
    @constraint(m, n_pay_d0_e_[i=P, j=Y, l=L],
                n_pay_1[i, j, l] == n_ann_1[i, j, l]
               )

    # 76 
    # 76 #######################################################################
    ##
    # -> expansion loan balance
    @constraint(m, n_loan_bal_y_e_[i=P, j=Y, l=L; j<n_years],
                n_loan[i, j+1, l] == n_loan[i, j, l]
                - n_pay[i, j, l]
                + n_ladd[i, j, l]
               )
    # 76 
    # 76 #######################################################################
    ##
    #n_c_Helec
    # -> heating requirement.
    # p.Hm (heating factor, i.e. heat / product)
    @constraint(m, n_eh_d_e_[i=P, j=Y, l=L, k=Kn],
                n_eh_d_[i, j, l, k] == 
                (1-p.n_c_Helec[l+1, k+1]) * p.n_c_Hfac[l+1, k+1] * 
                (p.n_c_H[l+1, k+1] * n_cp_d_[i, j, l, k] + 
                 p.n_rhs_h[l+1, k+1] * y_n[i, j, l, k])
               )
    @constraint(m, n_eh_d_m_i_[i=P, j=Y, l=L, k=Kn],
                n_eh_d_[i, j, l, k] <= p.n_eh_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_eh_s_e_[i=P, j=Y, l=L],
                n_eh[i, j, l] == sum(n_eh_d_[i, j, l, k] for k in Kn)
               )
    # -> fuel required for heat.
    # p.n_c_F in [0,1], (i.e. a fraction, heat by fuel /tot heat)
    @constraint(m, n_ehf_d_e_[i=P, j=Y, l=L, k=Kn, f=Fu],
                n_ehf_d_[i, j, l, k, f] ==
                p.n_c_F[l+1, f+1, k+1] * n_eh_d_[i, j, l, k] 
                + p.n_rhs_F[l+1, f+1, k+1] * y_n[i, j, l, k]
               )
    # n_ehf_bM
    @constraint(m, n_ehf_m_i_[i=P, j=Y, l=L, k=Kn, f=Fu],
                n_ehf_d_[i, j, l, k, f] <= p.n_ehf_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_ehf_s_e_[i=P, j=Y, l=L, f=Fu],
                n_ehf[i, j, l, f] == sum(n_ehf_d_[i, j, l, k, f] for k in Kn)
               )

    # -> electricity requirement
    # n_u and m_ud_, p.Um & p.UmRhs
    @constraint(m, n_u_d_e_[i=P, j=Y, l=L, k=Kn],
                n_u_d_[i, j, l, k] == 
                (1 - p.n_c_UonSite[i+1, j+1, l+1, k+1])*p.n_c_Ufac[l+1,k+1]*
                (p.n_c_U[l+1, k+1] * n_cp_d_[i, j, l, k] 
                 + p.n_rhs_U[l+1, k+1] * y_n[i, j, l, k]) +
                # electrification
                p.n_c_Helec[l+1, k+1] * p.n_c_Hfac[l+1, k+1] * (
                                        p.n_c_H[l+1, k+1]*n_cp_d_[i, j, l, k] 
                                        + p.n_rhs_h[l+1, k+1]*y_n[i, j, l, k]
                                       )
               )
    # n_u_bM (big-M)
    @constraint(m, n_u_i_[i=P, j=Y, l=L, k=Kn],
                n_u_d_[i, j, l, k] <= p.n_u_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_u_s_e_[i=P, j=Y, l=L],
                n_u[i, j, l] == sum(n_u_d_[i, j, l, k] for k in Kn)
               )
    # -> electricity on-site 
    # n_u_onsite, n_c_Uonsite viz. the onsite fraction
    @constraint(m, n_u_onsite_d_e_[i=P, j=Y, l=L, k=Kn],
                n_u_onsite_d_[i, j, l, k] == 
                p.n_c_UonSite[i+1,j+1,l+1,k+1]*p.n_c_Ufac[l+1,k+1]*
                (p.n_c_U[l+1, k+1] * n_cp_d_[i, j, l, k] 
                #+ p.n_c_UonSite[i+1,j+1,l+1,k+1]*
                +p.n_rhs_U[l+1, k+1] * y_n[i, j, l, k])
               )
    # n_u_bM (big-M)
    @constraint(m, n_u_onsite_i_[i=P, j=Y, l=L, k=Kn],
                n_u_onsite_d_[i, j, l, k] <= p.n_u_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_u_onsite_s_e_[i=P, j=Y, l=L],
                n_u_onsite[i, j, l] == sum(n_u_onsite_d_[i, j, l, k] for k in Kn)
               )

    
    # -> process (intrinsic) emissions
    # n_cp_e, n_cp_e_d_. p.Cp & p.CpRhs
    @constraint(m, n_cp_e_d_e_[i=P, j=Y, l=L, k=Kn],
                n_cp_e_d_[i, j, l, k] == 
                p.n_c_cp_e[l+1, k+1] * n_cp_d_[i, j, l, k]
                + p.n_rhs_cp_e[l+1, k+1] * y_n[i, j, l, k]
               )
    # n_cp_e_bM
    @constraint(m, n_cp_e_m_i_[i=P, j=Y, l=L, k=Kn],
                n_cp_e_d_[i, j, l, k] <= p.n_cp_e_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_cp_e_s_e_[i=P, j=Y, l=L],
                n_cp_e[i, j, l] == sum(n_cp_e_d_[i, j, l, k] for k in Kn)
               )

    # -> -> process (disaggregated) emissions

    # -> scope 0 emission
    # c_Fe (fuel emission factor), n_Hr (heat rate)
    @constraint(m, n_ep0_d_e_[i=P, j=Y, l=L, k=Kn],
                n_ep0_d_[i, j, l, k] == 
                # fuel
                sum(p.n_c_Fe[f+1, k+1]*n_ehf_d_[i, j, l, k, f] for f in Fu) 
                # process
                + n_cp_e_d_[i, j, l, k]
                # in-site electricity
                + sum(p.n_c_Fe[f+1, k+1]
                      *p.n_c_Hr[l+1, f+1, k+1] 
                      *p.n_c_Fgenf[l+1, f+1, k+1]
                      for f in Fu) * n_u_onsite_d_[i, j, l, k]
               )
    # n_ep0_bM
    @constraint(m, n_ep0_m_i_[i=P, j=Y, l=L, k=Kn],
                n_ep0_d_[i, j, l, k] <= p.n_ep0_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_ep0_s_e_[i=P, j=Y, l=L],
                n_ep0[i, j, l] == sum(n_ep0_d_[i, j, l, k] for k in Kn)
               )
    # -> scope 1 emitted
    # p.n_chi
    @constraint(m, n_ep1ge_d_e_[i=P, j=Y, l=L, k=Kn],
                n_ep1ge_d_[i, j, l, k] == (1.0 - p.n_chi[l+1, k+1]) 
                * n_ep0_d_[i, j, l, k]
               )
    # n_ep1ge_bM
    @constraint(m, n_ep1ge_m_i_[i=P, j=Y, l=L, k=Kn],
                n_ep1ge_d_[i, j, l, k] <= p.n_ep1ge_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_ep1ge_s_e_[i=P, j=Y, l=L],
                n_ep1ge[i, j, l] == sum(n_ep1ge_d_[i, j, l, k] for k in Kn)
               )

    # -> scope 1 captured
    # p.n_sigma
    @constraint(m, n_ep1gce_d_e_[i=P, j=Y, l=L, k=Kn],
                n_ep1gce_d_[i, j, l, k] == 
                p.n_chi[l+1, k+1] * (1 - p.n_sigma[l+1, k+1]) 
                * n_ep0_d_[i, j, l, k]
               )
    @constraint(m, n_ep1gce_m_i_[i=P, j=Y, l=L, k=Kn],
                n_ep1gce_d_[i, j, l, k] <= p.n_ep1gce_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_ep1gce_s_e_[i=P, j=Y, l=L],
                n_ep1gce[i, j, l] == sum(n_ep1gce_d_[i, j, l, k] for k in Kn)
               )

    # -> scope 1 stored
    # p.sigma ?
    @constraint(m, n_ep1gcs_d_e_[i=P, j=Y, l=L, k=Kn],
                n_ep1gcs_d_[i, j, l, k] ==
                p.n_chi[l+1, k+1] * p.n_sigma[l+1, k+1] * n_ep0_d_[i, j, l, k]
               )
    # ep1gcsm_bM
    @constraint(m, n_ep1gcs_m_i_[i=P, j=Y, l=L, k=Kn],
                n_ep1gcs_d_[i, j, l, k] <= p.n_ep1gcs_bM * y_n[i, j, l, k]
               )

    @constraint(m, n_ep1gcs_s_e_[i=P, j=Y, l=L],
                n_ep1gcs[i, j, l] == sum(n_ep1gcs_d_[i, j, l, k] for k in Kn)
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
    @constraint(m, n_conm_d_e_[i=P, j=Y, l=L, k=Kn],
                n_conm_d_[i, j, l, k] == 
                p.n_c_Onm[l+1, k+1] * n_cp_d_[i, j, l, k]
                + p.n_rhs_Onm[l+1, k+1] * y_n[i, j, l, k]
               )
    @constraint(m, n_conm_m_i_[i=P, j=Y, l=L, k=Kn],
                n_conm_d_[i, j, l, k] <= p.n_conm_bM * y_n[i, j, l, k]
               )
    @constraint(m, n_conm_s_e_[i=P, j=Y, l=L],
                n_conm[i, j, l] == sum(n_conm_d_[i, j, l, k] for k in Kn)
               )

    # 76 
    # 76 #######################################################################
    ##
    #if index_p isa Int && index_p == 0
        # -> new alloc logic
        @constraint(m, n_logic_0_init[l=L], 
                    y_n[0, 0, l, 0] == 1  # start with no new
                   )
        @constraint(m, n_logic_1_init[l=L, k=Kn; k>0], 
                    y_n[0, 0, l, k] == 0  # start with no new facility
                   )
    #end
    @constraint(m, n_logic_0_y_i_[i=P, j=Y, l=L; j<n_years],
                y_n[i, j+1, l, 0] <= y_n[i, j, l, 0]  # only build in the future
               )
    @constraint(m, n_logic_1_y_i_[i=P, j=Y, l=L, k=Kn; j<n_years && k>0],
                y_n[i, j+1, l, k] >= y_n[i, j, l, k]  # only build in the future
               )
    @constraint(m, n_logic_s_e[i=P, j=Y, l=L], 
                sum(y_n[i, j, l, k] for k in Kn) == 1
               )
    # for all points besides the initial one
    @constraint(m, n_logic_2[i=P, j=Y, l=L, k=Kn; k>0 && i>0 && j>0],
                y_n[i, j, l, k] + y_o[i, j, l] <= 1
               )
    
    # 76
    # 76 #######################################################################
    ##
    # -> initial loans
    #if index_p isa Int && index_p == 0
        @constraint(m, r_loan_initial_cond[l=L],
                    r_loan[0, 0, l] == p.r_loan0[l+1]
                   )
        @constraint(m, e_loan_init_cond[l=L],
                    e_loan[0, 0, l] == 0.0
                   )
        @constraint(m, n_loan_init_cond[l=L],
                    n_loan[0, 0, l] == 0.0
                   )
    #end

    # 76 
    # 76 #######################################################################
    ##
    # -> output capacity
    @constraint(m, o_cp_d_e[i=P, j=Y, l=L],
                o_cp_d_[i, j, l, 1] == o_tcp_d_[i, j, l, 1]
               )
    @constraint(m, o_cp_m1_i_[i=P, j=Y, l=L], # on
                o_cp_d_[i, j, l, 1] <= p.o_cp_bM * y_o[i, j, l]
               )
    #
    @constraint(m, o_cp_s_e_[i=P, j=Y, l=L],
                o_cp[i, j, l] == o_cp_d_[i, j, l, 1]
               )
    #
    @constraint(m, o_tcp_d_m1_i_[i=P, j=Y, l=L], # on
                o_tcp_d_[i, j, l, 1] <= p.o_cp_bM * y_o[i, j, l]
               )
    @constraint(m, o_tcp_d_m0_i_[i=P, j=Y, l=L], # off
                o_tcp_d_[i, j, l, 0] <= p.o_cp_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tcp_d_s_e_[i=P, j=Y, l=L],
                r_cp[i, j, l] == o_tcp_d_[i, j, l, 0] + o_tcp_d_[i, j, l, 1]
               ) # total cap
    # 76 
    # 76 #######################################################################
    ##
    # -> existing plant electricity consumption
    @constraint(m, o_u_d_e_[i=P, j=Y, l=L], 
                o_u_d_[i, j, l, 1] == o_tu_d_[i, j, l, 1]
               )
    @constraint(m, o_u_m_i_[i=P, j=Y, l=L],
                o_u_d_[i, j, l, 1] <= p.o_u_bM * y_o[i, j, l]
               )
    @constraint(m, o_u_s_[i=P, j=Y, l=L],
                o_u[i, j, l] == o_u_d_[i, j, l, 1]
               )
    @constraint(m, o_tu_d_m1_i_[i=P, j=Y, l=L],
                o_tu_d_[i, j, l, 1] <= p.o_u_bM * y_o[i, j, l]
               )
    @constraint(m, o_tu_d_m0_i_[i=P, j=Y, l=L],
                o_tu_d_[i, j, l, 0] <= p.o_u_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tu_s_e_[i=P, j=Y, l=L],
                r_u[i, j, l] == o_tu_d_[i, j, l, 0] + o_tu_d_[i, j, l, 1]
               )


    #
    # -> output emissions (why?)
    # -> scope 0
    @constraint(m, o_ep0_d_e_[i=P, j=Y, l=L],
                o_ep0_d_[i, j, l, 1] == o_tep0_d_[i, j, l, 1]
               )
    @constraint(m, o_ep0_m_i_[i=P, j=Y, l=L], # on
                o_ep0_d_[i, j, l, 1] <= p.o_ep0_bM * y_o[i, j, l]
               )
    @constraint(m, o_ep0_s_[i=P, j=Y, l=L],
                o_ep0[i, j, l] == o_ep0_d_[i, j, l, 1]
               )
    @constraint(m, o_tep0_d_m1_i_[i=P, j=Y, l=L], # on
                o_tep0_d_[i, j, l, 1] <= p.o_ep0_bM * y_o[i, j, l]
               )
    @constraint(m, o_tep0_d_m0_i_[i=P, j=Y, l=L], # off
                o_tep0_d_[i, j, l, 0] <= p.o_ep0_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tep0_s_e_[i=P, j=Y, l=L],
                r_ep0[i, j, l] == o_tep0_d_[i, j, l, 0] + o_tep0_d_[i, j, l, 1]
               )
    # ->  
    @constraint(m, o_ep1ge_d_e_[i=P, j=Y, l=L],
                o_ep1ge_d_[i, j, l, 1] == o_tep1ge_d_[i, j, l, 1]
               )
    @constraint(m, o_ep1ge_m_i_[i=P, j=Y, l=L], # on
                o_ep1ge_d_[i, j, l, 1] <= p.o_ep1ge_bM * y_o[i, j, l]
               )
    @constraint(m, o_ep1ge_s_[i=P, j=Y, l=L],
                o_ep1ge[i, j, l] == o_ep1ge_d_[i, j, l, 1]
               )
    @constraint(m, o_tep1ge_d_m1_i_[i=P, j=Y, l=L], # on
                o_tep1ge_d_[i, j, l, 1] <= p.o_ep1ge_bM * y_o[i, j, l]
               )
    @constraint(m, o_tep1ge_d_m0_i_[i=P, j=Y, l=L], # off
                o_tep1ge_d_[i, j, l, 0] <= p.o_ep1ge_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tep1ge_s_e_[i=P, j=Y, l=L],
                r_ep1ge[i, j, l] == 
                o_tep1ge_d_[i, j, l, 0] + o_tep1ge_d_[i, j, l, 1]
               )
    # ->
    @constraint(m, o_ep1gcs_d_e_[i=P, j=Y, l=L],
                o_ep1gcs_d_[i, j, l, 1] == o_tep1gcs_d_[i, j, l, 1]
               )
    @constraint(m, o_ep1gcs_m_i_[i=P, j=Y, l=L], # on
                o_ep1gcs_d_[i, j, l, 1] <= p.o_ep1gcs_bM * y_o[i, j, l]
               )
    @constraint(m, o_ep1gcs_s_[i=P, j=Y, l=L],
                o_ep1gcs[i, j, l] == o_ep1gcs_d_[i, j, l, 1]
               )
    @constraint(m, o_tep1gcs_d_m1_i_[i=P, j=Y, l=L], # on
                o_tep1gcs_d_[i, j, l, 1] <= p.o_ep1gcs_bM * y_o[i, j, l]
               )
    @constraint(m, o_tepgcs_d_m0_i_[i=P, j=Y, l=L], # off
                o_tep1gcs_d_[i, j, l, 0] <= p.o_ep1gcs_bM * (1 - y_o[i, j, l])
               )
    @constraint(m, o_tep1gcs_s_e_[i=P, j=Y, l=L],
                r_ep1gcs[i, j, l] == o_tep1gcs_d_[i, j, l, 0] 
                + o_tep1gcs_d_[i, j, l, 1]
               )
    # ->>
    @constraint(m, o_last_loan_s_e_[i=P, l=L], 
                r_loan_p[i, n_years, l] 
                + e_loan_p[i, n_years, l] == 
                o_loan_last[i, l, 0] + o_loan_last[i, l, 1]
               )
    @constraint(m, o_last_loan_d_m1_i_[i=P, l=L],
                o_loan_last[i, l, 1] <= 1e5 * y_o[i, n_years, l]
               )
    @constraint(m, o_last_loan_d_m0_i_[i=P, l=L],
                o_loan_last[i, l, 0] <= 1e5 * (1 - y_o[i, n_years, l])
               )


    # 76 
    # 76 #######################################################################
    ##
    # objective
    # pricing-problem
    return m
end

function attachFullObjectiveBlock(m::JuMP.Model, p::parms, s::sets)

    Y = s.Y
    P = s.P
    L = s.L

    o_pay = m[:o_pay]
    o_conm = m[:o_conm]
    n_pay = m[:n_pay]
    n_conm = m[:n_conm]
    o_loan_last = m[:o_loan_last]
    
    t_ret_cost = m[:t_ret_cost]  # 15
    n_loan_p = m[:n_loan_p]
    @objective(m, Min,
               # loan
               sum(sum(sum(p.discount[i+1, j+1] * o_pay[i, j, l] 
                           for l in L) for j in Y) for i in P)
               # o&m
               + sum(sum(sum(p.discount[i+1, j+1] * o_conm[i, j, l] 
                             for l in L) for j in Y) for i in P)
               # loan new
               + sum(sum(sum(p.discount[i+1, j+1] * n_pay[i, j, l] 
                             for l in L) for j in Y) for i in P)
               # o&m new
               + sum(sum(sum(p.discount[i+1, j+1] * n_conm[i, j, l] 
                             for l in L) for j in Y) for i in P)
               # retirement
               + sum(sum(sum(p.discount[i+1, j+1]*t_ret_cost[i, j, l] 
                             for l in L) for j in Y) for i in P)
               + sum(p.discount[n_periods+1,n_years+1]*
                     o_loan_last[n_periods,l,0] 
                     for l in L)
               + sum(p.discount[n_periods+1,n_years+1]*
                     n_loan_p[n_periods,l,0] 
                     for l in L)

               # if you retire but still have unpayed loan it is gonna be
               # reflected here :()
              )
    @printf "Full objective was be attached\n" 
end



function complicatingMatrix(m::Model, p::parms, s::sets)
    P = s.P
    Y = s.Y
    L = s.L

    Kr = s.Kr
    Kn = s.Kn

    @printf "Size of the periods\t%5i\n" length(P)
    @printf "Size of the years\t%5i\n" length(Y)
    @printf "Size of the locations\t%5i\n" length(L)
    @printf "Size of the retrofits\t%5i\n" length(Kr)
    @printf "Size of the newplants\t%5i\n" length(Kn)

    n_periods = p.n_periods
    n_years = p.n_years

    @printf "n_periods\t%5i\n" n_periods
    @printf "n_years\t%5i\n" n_years

    size_var_vec = 0
    
    # order online, retrofit, expansion, new.
    #
    y_o = m[:y_o]  # 1
    l_yo = length(P)*length(L)*2
    size_var_vec += l_yo

    y_r = m[:y_r]  # 2
    l_yr = length(P)*length(L)*length(Kr)*2
    size_var_vec += l_yr

    y_e = m[:y_e]  # 3
    l_ye = length(P)*length(L)*2
    size_var_vec += l_ye

    y_n = m[:y_n]  # 4
    l_yn = length(P)*length(L)*length(Kn)*2
    size_var_vec += l_yn
    
    r_ladd_d_ = m[:r_ladd_d_]  # 5
    l_rladdd = length(P)*length(L)
    size_var_vec += l_rladdd

    r_l_pd_ = m[:r_l_pd_]  # 6
    l_rlpd = length(P)*length(L)*2
    size_var_vec += l_rlpd

    r_loan = m[:r_loan]  # 7
    l_rloan = length(P)*length(L)*2
    size_var_vec += l_rloan

    r_pay = m[:r_pay]  # 8
    l_rpay = length(P)*length(L)
    size_var_vec += l_rpay

    r_ladd = m[:r_ladd]  # 9
    l_rladd = length(P)*length(L)
    size_var_vec += l_rladd

    e_ladd_d_ = m[:e_ladd_d_]  # 10
    l_eladdd = length(P)*length(L)
    size_var_vec += l_eladdd

    e_l_pd_ = m[:e_l_pd_]  # 11
    l_elpd = length(P)*length(L)
    size_var_vec += l_elpd

    e_loan = m[:e_loan]  # 12
    l_eloan = length(P)*length(L)*2
    size_var_vec += l_eloan

    e_pay = m[:e_pay]  # 13
    l_epay = length(P)*length(L)
    size_var_vec += l_epay

    e_ladd = m[:e_ladd]  # 14
    l_eladd = length(P)*length(L)
    size_var_vec += l_eladd
    
    # (retirement)
    t_ret_cost_d_ = m[:t_ret_cost_d_]  # 15
    l_tretcostd = length(P)*length(L)
    size_var_vec += l_tretcostd

    t_loan_d_ = m[:t_loan_d_]  # 16
    l_tloand = length(P)*length(L)*2
    size_var_vec += l_tloand
    
    n_ladd_d_ = m[:n_ladd_d_]  # 17
    l_nladdd = length(P)*length(L)
    size_var_vec += l_nladdd

    n_l_pd_ = m[:n_l_pd_]  # 18
    l_nlpd = length(P)*length(L)
    size_var_vec += l_nlpd

    n_loan = m[:n_loan]  # 19
    l_nloan = length(P)*length(L)*2
    size_var_vec += l_nloan

    n_pay = m[:n_pay]  # 20
    l_npay = length(P)*length(L)
    size_var_vec += l_npay

    n_ladd = m[:n_ladd]  # 21
    l_nladd = length(P)*length(L)
    size_var_vec += l_nladd
    
    o_cp = m[:o_cp]  # 22
    l_ocp = length(P)*length(Y)*length(L)
    size_var_vec += l_ocp

    n_cp = m[:n_cp]  # 23
    l_ncp = length(P)*length(Y)*length(L)
    size_var_vec += l_ncp

    o_ep0 = m[:o_ep0]  # 24
    l_oep0 = length(P)*length(Y)*length(L)
    size_var_vec += l_oep0

    n_ep0 = m[:n_ep0]  # 25
    l_nep0 = length(P)*length(Y)*length(L)
    size_var_vec += l_nep0

    o_u = m[:o_u]  # 26
    l_ou = length(P)*length(Y)*length(L)
    size_var_vec += l_ou

    n_u = m[:n_u]  # 27
    l_nu = length(P)*length(Y)*length(L)
    size_var_vec += l_nu

    #
    @printf "Size of the variable vector\t %5i\n" size_var_vec
    #
    var_vect = Vector{VariableRef}(undef, size_var_vec)
    var_ptr = Dict{VariableRef, Int64}()

    nvars_block = 1

    j = 1
    for i in P
        for l in L
            # block variables ####################
            var_vect[j] = y_o[i, 0, l]  # 1
            var_ptr[y_o[i,0,l]] = j
            j += 1
            var_vect[j] = y_o[i, n_years, l]  # 2
            var_ptr[y_o[i,n_years,l]] = j
            j += 1
            for k in Kr
                var_vect[j] = y_r[i, 0, l, k]  # 3
                var_ptr[y_r[i,0,l,k]] = j
                j += 1
                var_vect[j] = y_r[i, n_years, l, k]  # 4
                var_ptr[y_r[i,n_years,l,k]] = j
                j += 1
            end
            var_vect[j] = y_e[i, 0, l]  # 5
            var_ptr[y_e[i,0,l]] = j
            j += 1
            var_vect[j] = y_e[i, n_years, l]  # 6
            var_ptr[y_e[i,n_years,l]] = j
            j += 1
            for k in Kn
                var_vect[j] = y_n[i, 0, l, k]  # 7
                var_ptr[y_n[i,0,l,k]] = j
                j += 1
                var_vect[j] = y_n[i, n_years, l, k]  # 8
                var_ptr[y_n[i,n_years,l,k]] = j
                j += 1
            end
            var_vect[j] = r_ladd_d_[i, n_years, l, 0]  # 9
            var_ptr[r_ladd_d_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = r_l_pd_[i, n_years, l, 0]  # 10
            var_ptr[r_l_pd_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = r_l_pd_[i, n_years, l, 1]  # 11
            var_ptr[r_l_pd_[i,n_years,l,1]] = j
            j += 1
            var_vect[j] = r_loan[i, 0, l]  # 12
            var_ptr[r_loan[i,0,l]] = j
            j += 1
            var_vect[j] = r_loan[i, n_years, l]  # 13
            var_ptr[r_loan[i,n_years,l]] = j
            j += 1
            var_vect[j] = r_pay[i, n_years, l]  # 14
            var_ptr[r_pay[i,n_years,l]] = j
            j += 1
            var_vect[j] = r_ladd[i, n_years, l]  # 15
            var_ptr[r_ladd[i,n_years,l]] = j
            j += 1
            #y_e[i, n_years, l]
            #y_e[i+1, 0, l]
            var_vect[j] = e_ladd_d_[i, n_years, l, 0]  # 16
            var_ptr[e_ladd_d_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = e_l_pd_[i, n_years, l, 1]  # 17
            var_ptr[e_l_pd_[i,n_years,l,1]] = j
            j += 1
            var_vect[j] = e_loan[i, 0, l]  # 18
            var_ptr[e_loan[i,0,l]] = j
            j += 1
            var_vect[j] = e_loan[i, n_years, l]  # 19
            var_ptr[e_loan[i,n_years,l]] = j
            j += 1
            var_vect[j] = e_pay[i, n_years, l]  # 20
            var_ptr[e_pay[i,n_years,l]] = j
            j += 1
            var_vect[j] = e_ladd[i, n_years, l]  # 21
            var_ptr[e_ladd[i,n_years,l]] = j
            j += 1
            #y_o[i, n_years, l]
            #y_o[i+1,0, l]
            var_vect[j] = t_ret_cost_d_[i, n_years, l, 0]  # 22
            var_ptr[t_ret_cost_d_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = t_loan_d_[i, n_years, l, 0]  # 23
            var_ptr[t_loan_d_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = t_loan_d_[i, n_years, l, 1]  # 24
            var_ptr[t_loan_d_[i,n_years,l,1]] = j
            j += 1
            #y_n[i, n_years, l, 0]
            #y_n[i+1, 0, l, 0]
            var_vect[j] = n_ladd_d_[i, n_years, l, 0]  # 25
            var_ptr[n_ladd_d_[i,n_years,l,0]] = j
            j += 1
            var_vect[j] = n_l_pd_[i, n_years, l, 1]  # 26
            var_ptr[n_l_pd_[i,n_years,l,1]] = j
            j += 1
            var_vect[j] = n_loan[i, 0, l]  # 27
            var_ptr[n_loan[i,0,l]] = j
            j += 1
            var_vect[j] = n_loan[i, n_years, l]  # 28
            var_ptr[n_loan[i,n_years,l]] = j
            j += 1
            var_vect[j] = n_pay[i, n_years, l]  # 29
            var_ptr[n_pay[i,n_years,l]] = j
            j += 1
            var_vect[j] = n_ladd[i, n_years, l]  # 30
            var_ptr[n_ladd[i,n_years,l]] = j
            j += 1
            for y in Y
                var_vect[j] = o_cp[i, y, l]  # 31
                var_ptr[o_cp[i,y,l]] = j
                j += 1
                var_vect[j] = n_cp[i, y, l]  # 32
                var_ptr[n_cp[i,y,l]] = j
                j += 1
                var_vect[j] = o_ep0[i, y, l]  # 33
                var_ptr[o_ep0[i,y,l]] = j
                j += 1
                var_vect[j] = n_ep0[i, y, l]  # 34
                var_ptr[n_ep0[i,y,l]] = j
                j += 1
                var_vect[j] = o_u[i, y, l]  # 35
                var_ptr[o_u[i,y,l]] = j
                j += 1
                var_vect[j] = n_u[i, y, l]  # 36
                var_ptr[n_u[i,y,l]] = j
                j += 1
            end
            # block variables ####################
            if i == first(P) && l == first(L)
                nvars_block = j-1
                @printf "the block of variables is %i long\t" nvars_block
            end
        end
    end

    #
    size_con_vec = ((length(P) - 1) * length(L)  # 1
    + (length(P) - 1) * length(L) * (length(Kr) - 1)  # 2
    + (length(P) - 1) * length(L) * length(Kr)  # 3
    + (length(P) - 1) * length(L) * length(Kr)  # 4
    + (length(P) - 1) * length(L)  # 5
    + (length(P) - 1) * length(L)  # 6
    + (length(P) - 1) * length(L)  # 7
    + (length(P) - 1) * length(L)  # 8
    + (length(P) - 1) * length(L)  # 9
    + (length(P) - 1) * length(L)  # 10
    + (length(P) - 1) * length(L)  # 11
    + (length(P) - 1) * length(L)  # 12
    + (length(P) - 1) * length(L)  # 13
    + (length(P) - 1) * length(L)  # 14
    + (length(P) - 1) * length(L)  # 15
    + (length(P) - 1) * length(L)  # 16
    + (length(P) - 1) * length(L)  # 17
    + (length(P) - 1) * length(L)  # 18
    + (length(P) - 1) * length(L)  # 19
    + (length(P) - 1) * length(L) * (length(Kn) - 1)  # 20
    + length(P) * length(Y)  # 21
    + length(P) * length(Y)  # 22
   )
    
    @printf "size_con_vec\t%5i\n" size_con_vec
    #
    D = zeros(size_con_vec, size_var_vec)
    rhs = zeros(size_con_vec)
    consense = ones(Int32, size_con_vec)

    #
    ncon = 1
    for i in P
        if i < n_periods
            for l in L
                # 1 o_logic_1_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1, 0, l]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                for k in Kr
                    if k > 0
                        # 2 r_logic_budget_p_link_i_
                        v1 = var_ptr[y_r[i, n_years, l, k]]
                        v2 = var_ptr[y_r[i+1, 0, l, k]]
                        D[ncon, v1] = -1.
                        D[ncon, v2] = 1.
                        rhs[ncon] = 0.0
                        ncon += 1
                    end
                    # 3 r_logic_onoff_1_p_link_i_
                    v1 = var_ptr[y_o[i, n_years, l]]
                    v2 = var_ptr[y_r[i, n_years, l, k]]
                    v3 = var_ptr[y_r[i+1, 0, l, k]]
                    D[ncon, v1] = 1.
                    D[ncon, v2] = -1.
                    D[ncon, v3] = 1.
                    rhs[ncon] = 0.0
                    ncon += 1
                    # 4 r_logic_onoff_2_p_link_i_
                    v1 = var_ptr[y_o[i, n_years, l]]
                    v2 = var_ptr[y_r[i, n_years, l, k]]
                    v3 = var_ptr[y_r[i+1, 0, l, k]]
                    D[ncon, v1] = 1.
                    D[ncon, v2] = 1.
                    D[ncon, v3] = -1.
                    rhs[ncon] = 0.0
                    ncon += 1
                end
                # 5 r_ladd_m_0_p_link_i_
                v1 = var_ptr[y_r[i,n_years,l,0]]
                v2 = var_ptr[y_r[i+1,0,l,0]]
                v3 = var_ptr[r_ladd_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.r_ladd_bM
                D[ncon, v2] = -p.r_ladd_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 6 r_l_m_0_p_link_i_
                v1 = var_ptr[y_r[i,n_years,l,0]]
                v2 = var_ptr[y_r[i+1,0,l,0]]
                v3 = var_ptr[r_l_pd_[i, n_years, l, 0]]
                D[ncon, v1] = p.r_l_bM
                D[ncon, v2] = -p.r_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 7 r_l_m_1_p_link_i_
                v1 = var_ptr[y_r[i,n_years,l,0]]
                v2 = var_ptr[y_r[i+1,0,l,0]]
                v3 = var_ptr[r_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = -p.r_l_bM
                D[ncon, v2] = p.r_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.r_l_bM
                ncon += 1
                # 8 r_loan_bal_p_link_e_
                v1 = var_ptr[r_loan[i, n_years, l]]
                v2 = var_ptr[r_pay[i, n_years, l]]
                v3 = var_ptr[r_ladd[i, n_years, l]]
                v4 = var_ptr[r_loan[i+1, 0, l]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                D[ncon, v3] = 1.
                D[ncon, v4] = -1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 9 e_logic_1_p_link_i_
                v1 = var_ptr[y_e[i, n_years, l]]
                v2 = var_ptr[y_e[i+1, 0, l]]
                D[ncon, v1] = -1
                D[ncon, v2] = 1
                rhs[ncon] = 0.0
                ncon += 1
                # 10 e_ladd_m_0_p_link_i_
                v1 = var_ptr[y_e[i, n_years, l]]
                v2 = var_ptr[y_e[i+1, 0, l]]
                v3 = var_ptr[e_ladd_d_[i, n_years, l, 0]]
                D[ncon, v1] = -p.e_ladd_bM
                D[ncon, v2] = p.e_ladd_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 11 e_l_m_1_p_link_i_
                v1 = var_ptr[y_e[i, n_years, l]]
                v2 = var_ptr[y_e[i+1, 0, l]]
                v3 = var_ptr[e_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = p.e_l_bM
                D[ncon, v2] = -p.e_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.e_l_bM
                ncon += 1
                # 12 e_loan_bal_p_link_e_
                v1 = var_ptr[e_loan[i, n_years, l]]
                v2 = var_ptr[e_pay[i, n_years, l]]
                v3 = var_ptr[e_ladd[i, n_years, l]]
                v4 = var_ptr[e_loan[i+1, 0, l]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                D[ncon, v3] = 1.
                D[ncon, v4] = -1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 13 t_ret_c_bm_0_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1,0, l]]
                v3 = var_ptr[t_ret_cost_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.t_ret_c_bM
                D[ncon, v2] = -p.t_ret_c_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 14 r_loan_d_bm_0_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1,0, l]]
                v3 = var_ptr[t_loan_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.t_loan_bM
                D[ncon, v2] = -p.t_loan_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 15 r_loan_d_bm_1_p_link_i_
                #
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1, 0, l]]
                v3 = var_ptr[t_loan_d_[i, n_years, l, 1]]
                D[ncon, v1] = -p.t_loan_bM
                D[ncon, v2] = p.t_loan_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.t_loan_bM
                ncon += 1
                # 16 n_ladd_m_0_p_link_i_
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1, 0, l, 0]]
                v3 = var_ptr[n_ladd_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.n_ladd_bM
                D[ncon, v2] = -p.n_ladd_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 17 n_ladd_m_1_p_link_i_ 
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1,0, l, 0]]
                v3 = var_ptr[n_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = -p.n_l_bM
                D[ncon, v2] = p.n_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.n_l_bM
                ncon += 1
                # 18 n_loan_bal_p_link_e_
                v1 = var_ptr[n_loan[i, n_years, l]]
                v2 = var_ptr[n_pay[i, n_years, l]]
                v3 = var_ptr[n_ladd[i, n_years, l]]
                v4 = var_ptr[n_loan[i+1, 0, l]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                D[ncon, v3] = 1.
                D[ncon, v4] = -1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 19 n_logic_0_p_link_i_
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1, 0, l, 0]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                for k in Kn
                    if k > 0
                        # 20 n_logic_1_p_link_i_
                        v1 = var_ptr[y_n[i, n_years, l, k]]
                        v2 = var_ptr[y_n[i+1, 0, l, k]]
                        D[ncon, v1] = -1.
                        D[ncon, v2] = 1.
                        rhs[ncon] = 0.0
                        ncon += 1
                    end
                end
            end
        end
        # 21 ag_dem_l_link_i_
        for y in Y
            for l in L
                v1 = var_ptr[o_cp[i, y, l]]
                v2 = var_ptr[n_cp[i, y, l]]
                D[ncon, v1] = 1
                D[ncon, v2] = 1
                rhs[ncon] = p.demand[i+1, y+1]
            end
        ncon += 1
        end
        # 22 ag_co2_l_link_i_
        for y in Y
            for l in L
                v1 = var_ptr[o_ep0[i, y, l]]
                v2 = var_ptr[n_ep0[i, y, l]]
                v3 = var_ptr[o_u[i, y, l]]
                v4 = var_ptr[n_u[i, y, l]]
                D[ncon, v1] = -1.
                D[ncon, v2] = -1.
                D[ncon, v3] = -p.GcI[i+1, y+1, l+1]
                D[ncon, v4] = -p.GcI[i+1, y+1, l+1]
                rhs[ncon] = -p.co2_budget[i+1, l+1]
            end
        ncon += 1
        end
    end
    @printf "ncon=%i\n" ncon
    return D, var_vect, var_ptr, nvars_block, rhs
end


function attachPeriodBlock(m::Model, p::parms, s::sets)
    P = s.P
    Y = s.Y
    L = s.L

    Kr = s.Kr
    Kn = s.Kn

    n_periods = p.n_periods
    n_years = p.n_years
    
    # order online, retrofit, expansion, new.
    #
    y_o = m[:y_o]  # 1
    y_r = m[:y_r]  # 2
    y_e = m[:y_e]  # 3
    y_n = m[:y_n]  # 4
    r_ladd_d_ = m[:r_ladd_d_]  # 5
    r_l_pd_ = m[:r_l_pd_]  # 6
    r_loan = m[:r_loan]  # 7
    r_pay = m[:r_pay]  # 8
    r_ladd = m[:r_ladd]  # 9
    e_ladd_d_ = m[:e_ladd_d_]  # 10
    e_l_pd_ = m[:e_l_pd_]  # 11
    e_loan = m[:e_loan]  # 12
    e_pay = m[:e_pay]  # 13
    e_ladd = m[:e_ladd]  # 14
    # (retirement)
    t_ret_cost_d_ = m[:t_ret_cost_d_]  # 15
    t_loan_d_ = m[:t_loan_d_]  # 16
    n_ladd_d_ = m[:n_ladd_d_]  # 17
    n_l_pd_ = m[:n_l_pd_]  # 18
    n_loan = m[:n_loan]  # 19
    n_pay = m[:n_pay]  # 20
    n_ladd = m[:n_ladd]  # 21
    o_ep0 = m[:o_ep0]  # 22
    n_ep0 = m[:n_ep0]  # 23
    o_u = m[:o_u]  # 24
    n_u = m[:n_u]  # 25

    # 1: o_logic_1
    @constraint(m, o_logic_1_p_link_i_[i=P, l=L; i<n_periods],
                #y_o[i+1, 0, l] <= y_o[i, n_years, l]
                y_o[i, n_years, l] - y_o[i+1, 0, l]  >= 0.0
               )
    # 2: 
    @constraint(m, r_logic_budget_p_link_i_[i=P, l=L, k=Kr; 
                                       k>0 && i<n_periods],
                #y_r[i+1, 0, l, k] >= y_r[i, n_years, l, k]
                -y_r[i, n_years, l, k] + y_r[i+1, 0, l, k] >= 0.0
               )
    # 3: 
    @constraint(m, r_logic_onoff_1_p_link_i_[i=P, l=L, k=Kr;
                                        i<n_periods],
                y_o[i, n_years, l] - y_r[i, n_years, l, k] 
                + y_r[i+1, 0, l, k] >= 0.
               )
    # 4: 
    @constraint(m, r_logic_onoff_2_p_link_i_[i=P, l=L, k=Kr;
                                        i<n_periods],
                #y_o[i, n_years, l] + 1 - y_r[i+1, 0, l, k] 
                #+ y_r[i, n_years, l, k] >= 1
                y_o[i, n_years, l] + y_r[i, n_years, l, k]
                - y_r[i+1, 0, l, k] >= 0
               )
    # 5: 
    # p_add_bM
    @constraint(m, r_ladd_m_0_p_link_i_[i=P, l=L; i<n_periods],
                # r_ladd_d_[i, n_years, l, 0] <= # goes from 1 to 0 only
                # p.r_ladd_bM * (y_r[i,n_years,l,0] - y_r[i+1,0,l,0])
                p.r_ladd_bM * (y_r[i,n_years,l,0] - y_r[i+1,0,l,0])
                - r_ladd_d_[i, n_years, l, 0] 
                >= 0.
               )
    # 6: 
    # padd switch (loan-retrofit)
    @constraint(m, r_l_m_0_p_link_i_[i=P, l=L; i<n_periods],
                #r_l_pd_[i, n_years, l, 0] <=
                #p.r_l_bM * (y_r[i,n_years,l,0] - y_r[i+1,0,l,0])
                p.r_l_bM * (y_r[i,n_years,l,0] - y_r[i+1,0,l,0])
                - r_l_pd_[i, n_years, l, 0]
                >= 0.0
               )
    # 7: 
    @constraint(m, r_l_m_1_p_link_i_[i=P, l=L; i<n_periods],
                # r_l_pd_[i, n_years, l, 1] <=
                # p.r_l_bM * (1.0 - y_r[i,n_years,l,0] + y_r[i+1,0,l,0])
                p.r_l_bM * (-y_r[i,n_years,l,0] + y_r[i+1,0,l,0])
                - r_l_pd_[i, n_years, l, 1]
                >= -p.r_l_bM
               )
    # !!!!!!!!!!!!
    r_l = m[:r_l]
    @constraint(m, m_l_s_link_e_[i=P, l=L; i<n_periods], #  here
                r_l[i+1, 0, l] == 
                sum(r_l_pd_[i, n_years, l, k] for k in (0, 1))
               )
    # 8: 
    @constraint(m, r_loan_bal_p_link_e_[i=P, l=L; i<n_periods],
                #r_loan[i+1, 0, l] == r_loan[i, n_years, l] 
                #- r_pay[i, n_years, l] 
                #+ r_ladd[i, n_years, l]
                r_loan[i, n_years, l] 
                - r_pay[i, n_years, l] 
                + r_ladd[i, n_years, l] - r_loan[i+1, 0, l] == 0.
               )

    # 9: 
    @constraint(m, e_logic_1_p_link_i_[i=P, l=L; i<n_periods],
                #y_e[i+1, 0, l] >= y_e[i, n_years, l]  #
                -y_e[i, n_years, l] + y_e[i+1, 0, l] >= 0 #
               )
    # 10: 
    @constraint(m, e_ladd_m_0_p_link_i_[i=P, l=L; i<n_periods], 
                #e_ladd_d_[i, n_years, l, 0] <= 
                #p.e_ladd_bM * (y_e[i+1, 0, l] - y_e[i, n_years, l])
                p.e_ladd_bM * (-y_e[i, n_years, l] + y_e[i+1, 0, l])
                - e_ladd_d_[i, n_years, l, 0] >= 0.
               )
    # 11: 
    @constraint(m, e_l_m_1_p_link_i_[i=P, l=L; i<n_periods], 
                # need to set this to 0
                #e_l_pd_[i, n_years, l, 1] <= 
                #p.e_l_bM * (1 - y_e[i+1, 0, l] + y_e[i, n_years, l])
                p.e_l_bM * (y_e[i, n_years, l] - y_e[i+1, 0, l])
                - e_l_pd_[i, n_years, l, 1] >= -p.e_l_bM
               )  #  the 0th component is implied by the e_ladd_m constr
    # !!!!!!!!!
    e_l = m[:e_l]
    @constraint(m, e_l_s_link_e_[i=P, l=L; i<n_periods],
                e_l[i+1, 0, l] == 
                e_l_pd_[i, n_years, l, 0] + e_l_pd_[i, n_years, l, 1]
               )
    # 12: 
    @constraint(m, e_loan_bal_p_link_e_[i=P, l=L; i<n_periods],
                #e_loan[i+1, 0, l] == e_loan[i, n_years, l]
                #- e_pay[i, n_years, l]
                #+ e_ladd[i, n_years, l]
                e_loan[i, n_years, l]
                - e_pay[i, n_years, l]
                + e_ladd[i, n_years, l]
                - e_loan[i+1, 0, l]
                == 0.
               )
    ####
    # 13: , 
    # t_ret_c_bM
    @constraint(m, t_ret_c_bm_0_p_link_i_[i=P, l=L; i<n_periods],
                #t_ret_cost_d_[i, n_years, l, 0] <= 
                #p.t_ret_c_bM * (y_o[i, n_years, l] - y_o[i+1,0, l])
                p.t_ret_c_bM * (y_o[i, n_years, l] - y_o[i+1,0, l])
                - t_ret_cost_d_[i, n_years, l, 0] >= 0.
               )
    # 14: 
    # m_loan_d_bM
    @constraint(m, r_loan_d_bm_0_p_link_i_[i=P, l=L; i<n_periods],
                #t_loan_d_[i, n_years, l, 0] <=  # retired
                #p.t_loan_bM * (y_o[i, n_years, l] - y_o[i+1,0, l])
                p.t_loan_bM * (y_o[i, n_years, l] - y_o[i+1,0, l])
                - t_loan_d_[i, n_years, l, 0] >= 0.
               )
    # 15: 
    @constraint(m, r_loan_d_bm_1_p_link_i_[i=P, l=L; i<n_periods], # 
                #t_loan_d_[i, n_years, l, 1] <= p.t_loan_bM * 
                #(1 + y_o[i+1, 0, l] - y_o[i, n_years, l])
                p.t_loan_bM * (-y_o[i, n_years, l] + y_o[i+1, 0, l])
                - t_loan_d_[i, n_years, l, 1] >= -p.t_loan_bM
               )
    ###

    # 16: 
    @constraint(m, n_ladd_m_0_p_link_i_[i=P, l=L; i<n_periods], 
                #n_ladd_d_[i, n_years, l, 0] <= 
                #p.n_ladd_bM * (y_n[i, n_years, l, 0] - y_n[i+1, 0, l, 0])
                p.n_ladd_bM*(y_n[i, n_years, l, 0] - y_n[i+1, 0, l, 0])
                -n_ladd_d_[i, n_years, l, 0] >= 0.
               ) # y_n goes from 1 to 0
    # 17: 
    @constraint(m, n_ladd_m_1_p_link_i_[i=P, l=L;i<n_periods],
                #n_l_pd_[i, n_years, l, 1] <= 
                #p.n_l_bM * (1 + y_n[i+1,0, l, 0] - y_n[i, n_years, l, 0])
                p.n_l_bM * (-y_n[i, n_years, l, 0] + y_n[i+1,0, l, 0])
                -n_l_pd_[i, n_years, l, 1] >= -p.n_l_bM
               )
    n_l = m[:n_l]
    @constraint(m, n_l_s_e_link_i_[i=P, l=L; i<n_periods],
                n_l[i+1, 0, l] ==
                n_l_pd_[i, n_years, l, 0] + n_l_pd_[i, n_years, l, 1]
               )
    # 18: 
    @constraint(m, n_loan_bal_p_link_e_[i=P, l=L; i<n_periods],
                #n_loan[i+1, 0, l] == n_loan[i, n_years, l]
                #- n_pay[i, n_years, l]
                #+ n_ladd[i, n_years, l]
                n_loan[i, n_years, l]
                - n_pay[i, n_years, l]
                + n_ladd[i, n_years, l] - n_loan[i+1, 0, l] == 0.0
               )
    # 19: 
    @constraint(m, n_logic_0_p_link_i_[i=P, l=L; i<n_periods],
                #y_n[i+1, 0, l, 0] <= y_n[i, n_years, l, 0]
                y_n[i, n_years, l, 0] - y_n[i+1, 0, l, 0] >= 0.
               )
    # 20: 
    @constraint(m, n_logic_1_p_link_i_[i=P, l=L, k=Kn; i<n_periods && k>0],
                #y_n[i+1, 0, l, k] >= y_n[i, n_years, l, k]
                -y_n[i, n_years, l, k] + y_n[i+1, 0, l, k]  >= 0.
               )

end

function attachLocationBlock(m::Model, p::parms, s::sets)
    L = s.L

    P = s.P
    Y = s.Y
    
    o_cp = m[:o_cp]
    n_cp = m[:n_cp]
    # 21
    # aggregate demand constraint
    @constraint(m, ag_dem_l_link_i_[i=P, j=Y],
                sum(o_cp[i, j, l] for l in L) + sum(n_cp[i, j, l] for l in L)
                >= p.demand[i+1, j+1]
               )

    o_ep1ge = m[:o_ep1ge]
    n_ep1ge = m[:n_ep1ge]

    o_u = m[:o_u]
    n_u = m[:n_u]
    # 22
    #@constraint(m, ag_co2_l_link_i_[i=P, j=Y],
    #            # existing
    #            -sum(o_ep1ge[i, j, l] for l in L) 
    #            # new
    #            - sum(n_ep1ge[i, j, l] for l in L)
    #            # grid associated emissions
    #            - sum(p.GcI[i+1,j+1,l+1]*0.29329722222222*
    #                  (o_u[i, j, l] + n_u[i, j, l])
    #                  for l in L)
    #            >= -p.co2_budget[i+1, j+1]
    #           )


end

function turnover_con!(m::JuMP.Model, p::parms, s::sets)
    L = s.L

    P = s.P
    Y = s.Y
    
    o_cp = m[:o_cp]
    n_cp = m[:n_cp]
    
    # @constraint(m, turnover[i=P, j=Y],
    #             sum(n_cp[i, j, l] for l in L) <= 
    #             0.5 * sum(o_cp[i, j, l] for l in L)
    #            )



end

function min_ep1ge(m::JuMP.Model, p::parms, s::sets)

    Y = s.Y
    P = s.P
    L = s.L

    o_ep1ge = m[:o_ep1ge]
    n_ep1ge = m[:n_ep1ge]

    o_u = m[:o_u]
    n_u = m[:n_u]

    @objective(m, Min,
               sum(sum(sum(o_ep1ge[i, j, l] for l in L) 
                       # new
                       +sum(n_ep1ge[i, j, l] for l in L)
                       # grid associated emissions
                       +sum(p.GcI[i+1,j+1,l+1]*0.29329722222222*
                            (o_u[i, j, l] + n_u[i, j, l])
                            for l in L) for j in Y) for i in P)
              )

end

function reattachBlockMod!(m::Model, index_l, p::parms,
        s::sets)
    L = index_l 
    
    P = s.P
    Y = s.Y
    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu


end

function attachBlockObjective(m::Model, p::parms, s::sets,
        D_ij::Matrix, vv::Vector, pi_::Vector, i_::Int64, l_::Int64)
    o_pay = m[:o_pay]
    o_conm = m[:o_conm]
    n_pay = m[:n_pay]
    n_conm = m[:n_conm]
    o_loan_last = m[:o_loan_last]
    
    t_ret_cost = m[:t_ret_cost]  # 15
    Y = s.Y


    if l_ == last(s.L)
        clast = p.discount[n_periods+1, n_years+1]
    else
        clast = 0
    end

    @objective(m, Min,
               # loan
               sum(p.discount[i_+1, j+1] * o_pay[i_, j, l_] for j in Y)
               # o&m
               + sum(p.discount[i_+1, j+1] * o_conm[i_, j, l_] for j in Y)
               # loan new
               + sum(p.discount[i_+1, j+1] * n_pay[i_, j, l_] for j in Y)
               # o&m new
               + sum(p.discount[i_+1, j+1] * n_conm[i_, j, l_] for j in Y)
               # retirement
               + sum(p.discount[i_+1, j+1]*t_ret_cost[i_, j, l_]  for j in Y)
               # complicating block
               - pi_'*D_ij*vv

              )
end

function genBlockVarVec(m::Model, p::parms, s::sets, i_, l_)

    Y = s.Y

    Kr = s.Kr
    Kn = s.Kn

    n_periods = p.n_periods
    n_years = p.n_years

    size_var_vec = 0

    y_o = m[:y_o]  # 1
    size_var_vec += 2

    y_r = m[:y_r]  # 2
    size_var_vec += length(Kr)*2

    y_e = m[:y_e]  # 3
    size_var_vec += 2

    y_n = m[:y_n]  # 4
    size_var_vec += length(Kn)*2
    
    r_ladd_d_ = m[:r_ladd_d_]  # 5
    size_var_vec += 1

    r_l_pd_ = m[:r_l_pd_]  # 6
    size_var_vec += 2

    r_loan = m[:r_loan]  # 7
    size_var_vec += 2

    r_pay = m[:r_pay]  # 8
    size_var_vec += 1

    r_ladd = m[:r_ladd]  # 9
    size_var_vec += 1

    e_ladd_d_ = m[:e_ladd_d_]  # 10
    size_var_vec += 1

    e_l_pd_ = m[:e_l_pd_]  # 11
    size_var_vec += 1

    e_loan = m[:e_loan]  # 12
    size_var_vec += 2

    e_pay = m[:e_pay]  # 13
    size_var_vec += 1

    e_ladd = m[:e_ladd]  # 14
    size_var_vec += 1
    
    # (retirement)
    t_ret_cost_d_ = m[:t_ret_cost_d_]  # 15
    size_var_vec += 1

    t_loan_d_ = m[:t_loan_d_]  # 16
    size_var_vec += 2
    
    n_ladd_d_ = m[:n_ladd_d_]  # 17
    size_var_vec += 1

    n_l_pd_ = m[:n_l_pd_]  # 18
    size_var_vec += 1

    n_loan = m[:n_loan]  # 19
    size_var_vec += 2

    n_pay = m[:n_pay]  # 20
    size_var_vec += 1

    n_ladd = m[:n_ladd]  # 21
    size_var_vec += 1
    
    o_cp = m[:o_cp]  # 22
    size_var_vec += length(Y)

    n_cp = m[:n_cp]  # 23
    size_var_vec += length(Y)

    o_ep0 = m[:o_ep0]  # 24
    size_var_vec += length(Y)

    n_ep0 = m[:n_ep0]  # 25
    size_var_vec += length(Y)

    o_u = m[:o_u]  # 26
    size_var_vec += length(Y)

    n_u = m[:n_u]  # 27
    size_var_vec += length(Y)
   

    var_vect = Vector{VariableRef}(undef, size_var_vec)
    var_ptr = Dict{VariableRef, Int64}()

    j = 1


    var_vect[j] = y_o[i_, 0, l_]  # 1
    var_ptr[y_o[i_,0,l_]] = j
    j += 1
    var_vect[j] = y_o[i_, n_years, l_]  # 2
    var_ptr[y_o[i_,n_years,l_]] = j
    j += 1
    for k in Kr
        var_vect[j] = y_r[i_, 0, l_, k]  # 3
        var_ptr[y_r[i_,0,l_,k]] = j
        j += 1
        var_vect[j] = y_r[i_, n_years, l_, k]  # 4
        var_ptr[y_r[i_,n_years,l_,k]] = j
        j += 1
    end
    var_vect[j] = y_e[i_, 0, l_]  # 5
    var_ptr[y_e[i_,0,l_]] = j
    j += 1
    var_vect[j] = y_e[i_, n_years, l_]  # 6
    var_ptr[y_e[i_,n_years,l_]] = j
    j += 1
    for k in Kn
        var_vect[j] = y_n[i_, 0, l_, k]  # 7
        var_ptr[y_n[i_,0,l_,k]] = j
        j += 1
        var_vect[j] = y_n[i_, n_years, l_, k]  # 8
        var_ptr[y_n[i_,n_years,l_,k]] = j
        j += 1
    end
    var_vect[j] = r_ladd_d_[i_, n_years, l_, 0]  # 9
    var_ptr[r_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = r_l_pd_[i_, n_years, l_, 0]  # 10
    var_ptr[r_l_pd_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = r_l_pd_[i_, n_years, l_, 1]  # 11
    var_ptr[r_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = r_loan[i_, 0, l_]  # 12
    var_ptr[r_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = r_loan[i_, n_years, l_]  # 13
    var_ptr[r_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = r_pay[i_, n_years, l_]  # 14
    var_ptr[r_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = r_ladd[i_, n_years, l_]  # 15
    var_ptr[r_ladd[i_,n_years,l_]] = j
    j += 1
    #y_e[i_, n_years, l_]
    #y_e[i_+1, 0, l_]
    var_vect[j] = e_ladd_d_[i_, n_years, l_, 0]  # 16
    var_ptr[e_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = e_l_pd_[i_, n_years, l_, 1]  # 17
    var_ptr[e_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = e_loan[i_, 0, l_]  # 18
    var_ptr[e_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = e_loan[i_, n_years, l_]  # 19
    var_ptr[e_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = e_pay[i_, n_years, l_]  # 20
    var_ptr[e_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = e_ladd[i_, n_years, l_]  # 21
    var_ptr[e_ladd[i_,n_years,l_]] = j
    j += 1
    #y_o[i_, n_years, l_]
    #y_o[i_+1,0, l_]
    var_vect[j] = t_ret_cost_d_[i_, n_years, l_, 0]  # 22
    var_ptr[t_ret_cost_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = t_loan_d_[i_, n_years, l_, 0]  # 23
    var_ptr[t_loan_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = t_loan_d_[i_, n_years, l_, 1]  # 24
    var_ptr[t_loan_d_[i_,n_years,l_,1]] = j
    j += 1
    #y_n[i_, n_years, l_, 0]
    #y_n[i_+1, 0, l_, 0]
    var_vect[j] = n_ladd_d_[i_, n_years, l_, 0]  # 25
    var_ptr[n_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = n_l_pd_[i_, n_years, l_, 1]  # 26
    var_ptr[n_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = n_loan[i_, 0, l_]  # 27
    var_ptr[n_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = n_loan[i_, n_years, l_]  # 28
    var_ptr[n_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = n_pay[i_, n_years, l_]  # 29
    var_ptr[n_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = n_ladd[i_, n_years, l_]  # 30
    var_ptr[n_ladd[i_,n_years,l_]] = j
    j += 1
    for y in Y
        var_vect[j] = o_cp[i_, y, l_]  # 31
        var_ptr[o_cp[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_cp[i_, y, l_]  # 32
        var_ptr[n_cp[i_,y,l_]] = j
        j += 1
        var_vect[j] = o_ep0[i_, y, l_]  # 33
        var_ptr[o_ep0[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_ep0[i_, y, l_]  # 34
        var_ptr[n_ep0[i_,y,l_]] = j
        j += 1
        var_vect[j] = o_u[i_, y, l_]  # 35
        var_ptr[o_u[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_u[i_, y, l_]  # 36
        var_ptr[n_u[i_,y,l_]] = j
        j += 1
    end
    return var_vect, var_ptr
end
