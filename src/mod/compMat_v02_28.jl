using Printf

function complicatingMatrix(m::JuMP.Model, p::parms, s::sets)
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


    # these variables have to follow the same order. 
    # order online, retrofit, expansion, new.
   
    # calculate the size of the vector of variables
    lv, _ = initComplVars(m, s)
    y_o = lv[1]  # 1*
    y_r = lv[2]  # 2*
    y_e = lv[3]  # 3*
    y_n = lv[4]  # 4*
    r_ladd_d_ = lv[5]  # 5
    r_l_pd_ = lv[6]  # 6*
    r_l = lv[7]  # 7
    r_loan = lv[8]  # 8*
    r_pay = lv[9]  # 9
    r_ladd = lv[10]  # 10
    e_ladd_d_ = lv[11]  # 11
    e_l_pd_ = lv[12]  # 12*
    e_l = lv[13]  # 13
    e_loan = lv[14]  # 14*
    e_pay = lv[15]  # 15
    e_ladd = lv[16]  # 16
    t_ret_cost_d_ = lv[17]  # 17
    t_loan_d_ = lv[18]  # 18*
    n_ladd_d_ = lv[19]  # 19
    n_l_pd_ = lv[20]  # 20*
    n_l = lv[21]  # 21
    n_loan = lv[22]  # 22*
    n_pay = lv[23]  # 23
    n_ladd = lv[24]  # 24
    o_cp = lv[25]  # 25
    n_cp = lv[26]  # 26
    o_ep1ge = lv[27]  # 27
    n_ep1ge = lv[28]  # 28
    o_u = lv[29]  # 29
    n_u = lv[30]  # 30
    x = lv[31]  # 31
    n_c0 = lv[32]  # 32

    #
    var_vect = Vector{VariableRef}(undef, 0)
    var_ptr = Dict{VariableRef, Int64}()
    nvb = 0
    j = 0
    for i_ in P
        for l_ in L
            ############## block variables #################
            vv, vp, nvb = genBlockVarVec(m, p, s, i_, l_)
            var_vect = vcat(var_vect, vv) 
            for k in keys(vp)
                var_ptr[k] = vp[k] + nvb * j
            end
            j += 1
        end
    end
    @printf "Variable block size\t%i\n" nvb
    @printf "Size of the variable vector (result)\t%i \n" length(var_vect)
    @printf "Number of blocks\t%i \n" j

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
    + (length(P) - 1) * length(L)  # 20
    + (length(P) - 1) * length(L)  # 21
    + (length(P) - 1) * length(L) * (length(Kn) - 1)  # 22
    + (length(P) - 1) * length(L)  # 23
    + (length(P) - 1) * length(L)  # 24
    + length(P) * length(Y)  # 25
    + length(P) * length(Y)  # 26
   )

    size_var_vec = length(var_vect) 
    @printf "size_con_vec\t%5i\n" size_con_vec
    #
    D = zeros(size_con_vec, size_var_vec)
    rhs = zeros(size_con_vec)
    consense = ones(Int32, size_con_vec)
    # 1 => >=
    # 0 => ==
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
                    end  # if k > 0
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
                end  # for k in Kr
                # 5 r_ladd_m_0_p_link_i_
                v1 = var_ptr[y_r[i,n_years,l,0]]
                v2 = var_ptr[y_r[i+1,0,l,0]]
                v3 = var_ptr[r_ladd_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.r_ladd_bM
                D[ncon, v2] = -p.r_ladd_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 6 r_l_m_1_y_p_link_i_
                v1 = var_ptr[y_r[i,n_years,l,0]]
                v2 = var_ptr[y_r[i+1,0,l,0]]
                v3 = var_ptr[r_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = -p.r_l_bM
                D[ncon, v2] = p.r_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.r_l_bM
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
                # 14 t_ret_c_bm_0_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1,0, l]]
                v3 = var_ptr[t_ret_cost_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.t_ret_c_bM
                D[ncon, v2] = -p.t_ret_c_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 15 r_loan_d_bm_0_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1,0, l]]
                v3 = var_ptr[t_loan_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.t_loan_bM
                D[ncon, v2] = -p.t_loan_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 16 r_loan_d_bm_1_p_link_i_
                v1 = var_ptr[y_o[i, n_years, l]]
                v2 = var_ptr[y_o[i+1, 0, l]]
                v3 = var_ptr[t_loan_d_[i, n_years, l, 1]]
                D[ncon, v1] = -p.t_loan_bM
                D[ncon, v2] = p.t_loan_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.t_loan_bM
                ncon += 1
                # 17 n_ladd_m_0_p_link_i_
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1, 0, l, 0]]
                v3 = var_ptr[n_ladd_d_[i, n_years, l, 0]]
                D[ncon, v1] = p.n_ladd_bM
                D[ncon, v2] = -p.n_ladd_bM
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                # 18 n_ladd_m_1_p_link_i_ 
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1,0, l, 0]]
                v3 = var_ptr[n_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = -p.n_l_bM
                D[ncon, v2] = p.n_l_bM
                D[ncon, v3] = -1.
                rhs[ncon] = -p.n_l_bM
                ncon += 1
                # 21 n_logic_0_p_link_i_
                v1 = var_ptr[y_n[i, n_years, l, 0]]
                v2 = var_ptr[y_n[i+1, 0, l, 0]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                rhs[ncon] = 0.0
                ncon += 1
                for k in Kn
                    if k > 0
                        # 22 n_logic_1_p_link_i_
                        v1 = var_ptr[y_n[i, n_years, l, k]]
                        v2 = var_ptr[y_n[i+1, 0, l, k]]
                        D[ncon, v1] = -1.
                        D[ncon, v2] = 1.
                        rhs[ncon] = 0.0
                        ncon += 1
                    end # if k > 0
                end # for k in Kn
            end  # for l in L
        end  # if i < n_periods
        #
        # 26 ag_dem_l_link_i_
        for y in Y
            for l in L
                v1 = var_ptr[o_cp[i, y, l]]
                v2 = var_ptr[n_cp[i, y, l]]
                D[ncon, v1] = 1
                D[ncon, v2] = 1
                rhs[ncon] = p.demand[i+1, y+1]
            end  # for l in L
        ncon += 1  # the previous one is a single constraint
        end  # for y in Y
        # 27 ag_co2_l_link_i_
        for y in Y
            for l in L
                v1 = var_ptr[o_ep1ge[i, y, l]]
                v2 = var_ptr[n_ep1ge[i, y, l]]
                v3 = var_ptr[o_u[i, y, l]]
                v4 = var_ptr[n_u[i, y, l]]
                D[ncon, v1] = -1.
                D[ncon, v2] = -1.
                D[ncon, v3] = -p.GcI[i+1, y+1, l+1]
                D[ncon, v4] = -p.GcI[i+1, y+1, l+1]
                rhs[ncon] = -p.co2_budget[i+1, l+1]
            end  # for l in L
        ncon += 1  # the previous one is a single constraint
        end  # for y in Y
    end  # for i in P
    ##
    for i in P
        if i < n_periods
            for l in L
                # 7 r_l_s_p_link_e_ ! eq
                v1 = var_ptr[r_l[i+1, 0, l]]
                v2 = var_ptr[r_l_pd_[i, n_years, l, 0]]
                v3 = var_ptr[r_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 8 r_loan_bal_p_link_e_ eq
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
                # 12 e_l_s_p_link_e_ ! eq
                v1 = var_ptr[e_l_pd_[i, n_years, l, 0]]
                v2 = var_ptr[e_l_pd_[i, n_years, l, 1]] # should this go ealier
                v3 = var_ptr[e_l[i+1, 0, l]]
                D[ncon, v1] = -1.
                D[ncon, v2] = -1.
                D[ncon, v3] = 1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 13 e_loan_bal_p_link_e_ ! eq
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
                # 19 n_l_s_e_link_i_ ! eq
                v1 = var_ptr[n_l[i+1, 0, l]]
                v2 = var_ptr[n_l_pd_[i, n_years, l, 0]]
                v3 = var_ptr[n_l_pd_[i, n_years, l, 1]]
                D[ncon, v1] = 1.
                D[ncon, v2] = -1.
                D[ncon, v3] = -1.
                rhs[ncon] = 0.0
                consense[ncon] = 0 # equality
                ncon += 1
                # 20 n_loan_bal_p_link_e_ eq
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
                # 23 x_p_link_e_
                v1 = var_ptr[x[i, l]]
                v2 = var_ptr[x[i+1, l]]
                D[ncon, v1] = 1.0
                D[ncon, v2] = -1.0
                rhs[ncon] = 0
                consense[ncon] = 0 # equality
                ncon += 1
                # 24 n_c0_p_link_e_
                v1 = var_ptr[n_c0[i, l]]
                v2 = var_ptr[n_c0[i+1, l]]
                D[ncon, v1] = 1.0
                D[ncon, v2] = -1.0
                rhs[ncon] = 0
                consense[ncon] = 0 # equality
                ncon += 1
            end  # for l in L

        end  # if i < n_periods
    end  # for i in P

    @printf "ncon=%i\n" ncon
    return D, var_vect, var_ptr, nvb, rhs, consense
end



function genBlockVarVec(m::JuMP.Model, p::parms, s::sets, i_, l_)

    Y = s.Y
    Kr = s.Kr
    Kn = s.Kn
    n_periods = p.n_periods
    n_years = p.n_years

    lv, size_var_vec = initComplVars(m, s)
    # assign variables as required
    y_o = lv[1]  # 1*
    y_r = lv[2]  # 2*
    y_e = lv[3]  # 3*
    y_n = lv[4]  # 4*
    r_ladd_d_ = lv[5]  # 5
    r_l_pd_ = lv[6]  # 6*
    r_l = lv[7]  # 7
    r_loan = lv[8]  # 8*
    r_pay = lv[9]  # 9
    r_ladd = lv[10]  # 10
    e_ladd_d_ = lv[11]  # 11
    e_l_pd_ = lv[12]  # 12*
    e_l = lv[13]  # 13
    e_loan = lv[14]  # 14*
    e_pay = lv[15]  # 15
    e_ladd = lv[16]  # 16
    t_ret_cost_d_ = lv[17]  # 17
    t_loan_d_ = lv[18]  # 18*
    n_ladd_d_ = lv[19]  # 19
    n_l_pd_ = lv[20]  # 20*
    n_l = lv[21]  # 21
    n_loan = lv[22]  # 22*
    n_pay = lv[23]  # 23
    n_ladd = lv[24]  # 24
    o_cp = lv[25]  # 25
    n_cp = lv[26]  # 26
    o_ep1ge = lv[27]  # 27
    n_ep1ge = lv[28]  # 28
    o_u = lv[29]  # 29
    n_u = lv[30]  # 30
    x = lv[31]  # 31
    n_c0 = lv[32]  # 32

   

    var_vect = Vector{VariableRef}(undef, size_var_vec)
    var_ptr = Dict{VariableRef, Int64}()
    
    # populate the vector elements with the variable references
    j = 1
    var_vect[j] = y_o[i_, 0, l_]  # 1-0
    var_ptr[y_o[i_,0,l_]] = j
    j += 1
    var_vect[j] = y_o[i_, n_years, l_]  # 1-1
    var_ptr[y_o[i_,n_years,l_]] = j
    j += 1
    for k in Kr
        var_vect[j] = y_r[i_, 0, l_, k]  # 2-0
        var_ptr[y_r[i_,0,l_,k]] = j
        j += 1
        var_vect[j] = y_r[i_, n_years, l_, k]  # 2-1
        var_ptr[y_r[i_,n_years,l_,k]] = j
        j += 1
    end
    var_vect[j] = y_e[i_, 0, l_]  # 3-0
    var_ptr[y_e[i_,0,l_]] = j
    j += 1
    var_vect[j] = y_e[i_, n_years, l_]  # 3-1
    var_ptr[y_e[i_,n_years,l_]] = j
    j += 1
    for k in Kn
        var_vect[j] = y_n[i_, 0, l_, k]  # 4-0
        var_ptr[y_n[i_,0,l_,k]] = j
        j += 1
        var_vect[j] = y_n[i_, n_years, l_, k]  # 4-1
        var_ptr[y_n[i_,n_years,l_,k]] = j
        j += 1
    end
    var_vect[j] = r_ladd_d_[i_, n_years, l_, 0]  # 5
    var_ptr[r_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = r_l_pd_[i_, n_years, l_, 0]  # 6-0
    var_ptr[r_l_pd_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = r_l_pd_[i_, n_years, l_, 1]  # 6-1
    var_ptr[r_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = r_l[i_, 0, l_]  # 7
    var_ptr[r_l[i_, 0, l_]] = j
    j += 1
    var_vect[j] = r_loan[i_, 0, l_]  # 8-0
    var_ptr[r_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = r_loan[i_, n_years, l_]  # 8-1
    var_ptr[r_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = r_pay[i_, n_years, l_]  # 9
    var_ptr[r_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = r_ladd[i_, n_years, l_]  # 10
    var_ptr[r_ladd[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = e_ladd_d_[i_, n_years, l_, 0]  # 11
    var_ptr[e_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = e_l_pd_[i_, n_years, l_, 0]  # 12-0
    var_ptr[e_l_pd_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = e_l_pd_[i_, n_years, l_, 1]  # 12-1
    var_ptr[e_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = e_l[i_, 0, l_]  # 13
    var_ptr[e_l[i_, 0, l_]] = j
    j += 1
    var_vect[j] = e_loan[i_, 0, l_]  # 14-0
    var_ptr[e_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = e_loan[i_, n_years, l_]  # 14-1
    var_ptr[e_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = e_pay[i_, n_years, l_]  # 15
    var_ptr[e_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = e_ladd[i_, n_years, l_]  # 16
    var_ptr[e_ladd[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = t_ret_cost_d_[i_, n_years, l_, 0]  # 17
    var_ptr[t_ret_cost_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = t_loan_d_[i_, n_years, l_, 0]  # 18-0
    var_ptr[t_loan_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = t_loan_d_[i_, n_years, l_, 1]  # 18-1
    var_ptr[t_loan_d_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = n_ladd_d_[i_, n_years, l_, 0]  # 19
    var_ptr[n_ladd_d_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = n_l_pd_[i_, n_years, l_, 0]  # 20-0
    var_ptr[n_l_pd_[i_,n_years,l_,0]] = j
    j += 1
    var_vect[j] = n_l_pd_[i_, n_years, l_, 1]  # 20-1
    var_ptr[n_l_pd_[i_,n_years,l_,1]] = j
    j += 1
    var_vect[j] = n_l[i_, 0, l_]  # 21
    var_ptr[n_l[i_, 0, l_]] = j
    j += 1

    var_vect[j] = n_loan[i_, 0, l_]  # 22-0
    var_ptr[n_loan[i_,0,l_]] = j
    j += 1
    var_vect[j] = n_loan[i_, n_years, l_]  # 22-1
    var_ptr[n_loan[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = n_pay[i_, n_years, l_]  # 23
    var_ptr[n_pay[i_,n_years,l_]] = j
    j += 1
    var_vect[j] = n_ladd[i_, n_years, l_]  # 24
    var_ptr[n_ladd[i_,n_years,l_]] = j
    j += 1
    for y in Y
        var_vect[j] = o_cp[i_, y, l_]  # 25
        var_ptr[o_cp[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_cp[i_, y, l_]  # 26
        var_ptr[n_cp[i_,y,l_]] = j
        j += 1
        var_vect[j] = o_ep1ge[i_, y, l_]  # 27
        var_ptr[o_ep1ge[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_ep1ge[i_, y, l_]  # 28
        var_ptr[n_ep1ge[i_,y,l_]] = j
        j += 1
        var_vect[j] = o_u[i_, y, l_]  # 29
        var_ptr[o_u[i_,y,l_]] = j
        j += 1
        var_vect[j] = n_u[i_, y, l_]  # 30
        var_ptr[n_u[i_,y,l_]] = j
        j += 1
    end
    var_vect[j] = x[i_, l_]  # 31
    var_ptr[x[i_, l_]] = j
    j += 1
    var_vect[j] = n_c0[i_, l_]  # 32
    var_ptr[n_c0[i_, l_]] = j
    j += 1

    @printf "Size of the blockvar(computed)\t%i\n" size_var_vec
    @printf "Size of the blockvar(resulting)\t%i\n" j-1
    return var_vect, var_ptr, (j-1)
end



function initComplVars(m::JuMP.Model, s::sets)
    Kr = s.Kr
    Kn = s.Kn
    Y = s.Y
    
    # calculate the size of the vector of variables
    size_bvar_v = 0

    list_cont = []
    push!(list_cont, m[:y_o])  # 1*
    size_bvar_v += 2

    push!(list_cont, m[:y_r])  # 2*
    size_bvar_v += length(Kr)*2

    push!(list_cont, m[:y_e])  # 3*
    size_bvar_v += 2

    push!(list_cont, m[:y_n])  # 4*
    size_bvar_v += length(Kn)*2
    
    push!(list_cont, m[:r_ladd_d_])  # 5
    size_bvar_v += 1

    push!(list_cont, m[:r_l_pd_])  # 6*
    size_bvar_v += 2

    push!(list_cont, m[:r_l])  # 7
    size_bvar_v += 1

    push!(list_cont, m[:r_loan])  # 8*
    size_bvar_v += 2

    push!(list_cont, m[:r_pay])  # 9
    size_bvar_v += 1

    push!(list_cont, m[:r_ladd])  # 10
    size_bvar_v += 1

    push!(list_cont, m[:e_ladd_d_])  # 11
    size_bvar_v += 1

    push!(list_cont, m[:e_l_pd_])  # 12*
    size_bvar_v += 2

    push!(list_cont, m[:e_l])  # 13
    size_bvar_v += 1

    push!(list_cont, m[:e_loan])  # 14*
    size_bvar_v += 2

    push!(list_cont, m[:e_pay])  # 15
    size_bvar_v += 1

    push!(list_cont, m[:e_ladd])  # 16
    size_bvar_v += 1
    
    # (retirement)
    push!(list_cont, m[:t_ret_cost_d_])  # 17
    size_bvar_v += 1

    push!(list_cont, m[:t_loan_d_])  # 18*
    size_bvar_v += 2
    
    push!(list_cont, m[:n_ladd_d_])  # 19
    size_bvar_v += 1

    push!(list_cont, m[:n_l_pd_])  # 20*
    size_bvar_v += 2

    push!(list_cont, m[:n_l])  # 21
    size_bvar_v += 1

    push!(list_cont, m[:n_loan])  # 22*
    size_bvar_v += 2

    push!(list_cont, m[:n_pay])  # 23
    size_bvar_v += 1

    push!(list_cont, m[:n_ladd])  # 24
    size_bvar_v += 1
    
    #
    push!(list_cont, m[:o_cp])  # 25
    size_bvar_v += length(Y)

    push!(list_cont, m[:n_cp])  # 26
    size_bvar_v += length(Y)

    push!(list_cont, m[:o_ep1ge])  # 27
    size_bvar_v += length(Y)

    push!(list_cont, m[:n_ep1ge])  # 28
    size_bvar_v += length(Y)

    push!(list_cont, m[:o_u])  # 29
    size_bvar_v += length(Y)

    push!(list_cont, m[:n_u])  # 30
    size_bvar_v += length(Y)

    push!(list_cont, m[:x])
    size_bvar_v += 1

    push!(list_cont, m[:n_c0])
    size_bvar_v += 1


    return list_cont, size_bvar_v   
end


function primalCompValueVec!(m::JuMP.Model, p::parms, s::sets)
    lv, size_var_vec = initComplVars(m, s)
    # assign variables as required
    y_o = lv[1]  # 1*
    y_r = lv[2]  # 2*
    y_e = lv[3]  # 3*
    y_n = lv[4]  # 4*
    r_ladd_d_ = lv[5]  # 5
    r_l_pd_ = lv[6]  # 6*
    r_l = lv[7]  # 7
    r_loan = lv[8]  # 8*
    r_pay = lv[9]  # 9
    r_ladd = lv[10]  # 10
    e_ladd_d_ = lv[11]  # 11
    e_l_pd_ = lv[12]  # 12*
    e_l = lv[13]  # 13
    e_loan = lv[14]  # 14*
    e_pay = lv[15]  # 15
    e_ladd = lv[16]  # 16
    t_ret_cost_d_ = lv[17]  # 17
    t_loan_d_ = lv[18]  # 18*
    n_ladd_d_ = lv[19]  # 19
    n_l_pd_ = lv[20]  # 20*
    n_l = lv[21]  # 21
    n_loan = lv[22]  # 22*
    n_pay = lv[23]  # 23
    n_ladd = lv[24]  # 24
    o_cp = lv[25]  # 25
    n_cp = lv[26]  # 26
    o_ep1ge = lv[27]  # 27
    n_ep1ge = lv[28]  # 28
    o_u = lv[29]  # 29
    n_u = lv[30]  # 30
    x = lv[31]
    n_c0 = lv[32]
   

    var_vect = Vector{Float64}(undef, size_var_vec)

    j = 1
    var_vect[j] = value(y_o[i_, 0, l_])  # 1-0
    j += 1
    var_vect[j] = value(y_o[i_, n_years, l_])  # 1-1
    j += 1
    for k in Kr
        var_vect[j] = value(y_r[i_, 0, l_, k])  # 2-0
        j += 1
        var_vect[j] = value(y_r[i_, n_years, l_, k])  # 2-1
        j += 1
    end
    var_vect[j] = value(y_e[i_, 0, l_])  # 3-0
    j += 1
    var_vect[j] = value(y_e[i_, n_years, l_])  # 3-1
    j += 1
    for k in Kn
        var_vect[j] = value(y_n[i_, 0, l_, k])  # 4-0
        j += 1
        var_vect[j] = value(y_n[i_, n_years, l_, k])  # 4-1
        j += 1
    end
    var_vect[j] = value(r_ladd_d_[i_, n_years, l_, 0])  # 5
    j += 1
    var_vect[j] = value(r_l_pd_[i_, n_years, l_, 0])  # 6-0
    j += 1
    var_vect[j] = value(r_l_pd_[i_, n_years, l_, 1])  # 6-1
    j += 1
    var_vect[j] = value(r_l[i_, 0, l_])  # 7
    j += 1
    var_vect[j] = value(r_loan[i_, 0, l_])  # 8-0
    j += 1
    var_vect[j] = value(r_loan[i_, n_years, l_])  # 8-1
    j += 1
    var_vect[j] = value(r_pay[i_, n_years, l_])  # 9
    j += 1
    var_vect[j] = value(r_ladd[i_, n_years, l_])  # 10
    j += 1
    var_vect[j] = value(e_ladd_d_[i_, n_years, l_, 0])  # 11
    j += 1
    var_vect[j] = value(e_l_pd_[i_, n_years, l_, 0])  # 12-0
    j += 1
    var_vect[j] = value(e_l_pd_[i_, n_years, l_, 1])  # 12-1
    j += 1
    var_vect[j] = value(e_l[i_, 0, l_])  # 13
    j += 1
    var_vect[j] = value(e_loan[i_, 0, l_])  # 14-0
    j += 1
    var_vect[j] = value(e_loan[i_, n_years, l_])  # 14-1
    j += 1
    var_vect[j] = value(e_pay[i_, n_years, l_])  # 15
    j += 1
    var_vect[j] = value(e_ladd[i_, n_years, l_])  # 16
    j += 1
    var_vect[j] = value(t_ret_cost_d_[i_, n_years, l_, 0])  # 17
    j += 1
    var_vect[j] = value(t_loan_d_[i_, n_years, l_, 0])  # 18-0
    j += 1
    var_vect[j] = value(t_loan_d_[i_, n_years, l_, 1])  # 18-1
    j += 1
    var_vect[j] = value(n_ladd_d_[i_, n_years, l_, 0])  # 19
    j += 1
    var_vect[j] = value(n_l_pd_[i_, n_years, l_, 0])  # 20-0
    j += 1
    var_vect[j] = value(n_l_pd_[i_, n_years, l_, 1])  # 20-1
    j += 1
    var_vect[j] = value(n_l[i_, 0, l_])  # 21
    j += 1

    var_vect[j] = value(n_loan[i_, 0, l_])  # 22-0
    j += 1
    var_vect[j] = value(n_loan[i_, n_years, l_])  # 22-1
    j += 1
    var_vect[j] = value(n_pay[i_, n_years, l_])  # 23
    j += 1
    var_vect[j] = value(n_ladd[i_, n_years, l_])  # 24
    j += 1
    for y in Y
        var_vect[j] = value(o_cp[i_, y, l_])  # 25
        j += 1
        var_vect[j] = value(n_cp[i_, y, l_])  # 26
        j += 1
        var_vect[j] = value(o_ep1ge[i_, y, l_])  # 27
        j += 1
        var_vect[j] = value(n_ep1ge[i_, y, l_])  # 28
        j += 1
        var_vect[j] = value(o_u[i_, y, l_])  # 29
        j += 1
        var_vect[j] = value(n_u[i_, y, l_])  # 30
        j += 1
    end

    return var_vect
end


function primalValVec(mv::Vector{JuMP.Model}, p::parms, s::sets)
    # A single column vector for the primal values is considered.
    # Iteration by block is required.
    P = s.P
    Y = s.Y
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu
    #
    size_v = (
              length(Y)	  #	1
              + length(Y)*length(Kr)	  #	2
              + length(Y)	  #	3
              + length(Y)*2	  #	4
              + length(Y)	  #	5
              + length(Y)	  #	6
              + length(Y)*2	  #	7
              + 1	  #	8
              + length(Y)*2	  #	9
              + length(Y)	  #	10
              + length(Y)*length(Kr)	  #	11
              + length(Y)	  #	12
              + length(Y)*length(Kr)	  #	13
              + length(Y)	  #	14
              + length(Y)*length(Kr)	  #	15
              + length(Y)*length(Fu)	  #	16
              + length(Y)*length(Kr)*length(Fu)	  #	17
              + length(Y)	  #	18
              + length(Y)*length(Kr)	  #	19
              + length(Y)	  #	20
              + length(Y)*length(Kr)	  #	21
              + length(Y)	  #	22
              + length(Y)*length(Kr)	  #	23
              + length(Y)	  #	24
              + length(Y)*length(Kr)	  #	25
              + length(Y)	  #	26
              + length(Y)*length(Kr)	  #	27
              + length(Y)	  #	28
              + length(Y)*length(Kr)	  #	29
              + length(Y)	  #	30
              + length(Y)*length(Kr)	  #	31
              + length(Y)*length(Kr)	  #	32
              + length(Y)	  #	33
              + length(Y)	  #	34
              + length(Y)	  #	35
              + length(Y)	  #	36
              + length(Y)	  #	37
              + length(Y)	  #	38
              + length(Y)	  #	39
              + length(Y)	  #	40
              + length(Y)	  #	41
              + length(Y)	  #	42
              + length(Y)*length(Kr)	  #	43
              + length(Y)*2	  #	44
              + length(Y)	  #	45
              + length(Y)	  #	46
              + length(Y)	  #	47
              + length(Y)*length(Kr)	  #	48
              + length(Y)	  #	49
              + length(Y)	  #	50
              + length(Y)*2	  #	51
              + length(Y)	  #	52
              + length(Y)	  #	53
              + length(Y)	  #	54
              + length(Y)	  #	55
              + length(Y)	  #	56
              + length(Y)	  #	57
              + length(Y)	  #	58
              + length(Y)	  #	59
              + length(Y)	  #	60
              + length(Y)	  #	61
              + length(Y)*2	  #	62
              + length(Y)	  #	63
              + length(Y)	  #	64
              + length(Y)*length(Kn)	  #	65
              + 1	  #	66
              + length(Y)*length(Kn)	  #	67
              + length(Y)	  #	68
              + length(Y)*length(Kn)	  #	69
              + length(Y)*length(Kn)	  #	70
              + length(Y)	  #	71
              + length(Y)*length(Kn)	  #	72
              + length(Y)	  #	73
              + length(Y)	  #	74
              + length(Y)	  #	75
              + length(Y)	  #	76
              + length(Y)	  #	77
              + length(Y)*2	  #	78
              + length(Y)	  #	79
              + length(Y)	  #	80
              + length(Y)	  #	81
              + length(Y)	  #	82
              + length(Y)	  #	83
              + length(Y)	  #	84
              + length(Y)	  #	85
              + length(Y)*length(Kn)	  #	86
              + length(Y)*length(Fu)	  #	87
              + length(Y)*length(Kn)*length(Fu)	  #	88
              + length(Y)	  #	89
              + length(Y)*length(Kn)	  #	90
              + length(Y)	  #	91
              + length(Y)*length(Kn)	  #	92
              + length(Y)	  #	93
              + length(Y)*length(Kn)	  #	94
              + length(Y)	  #	95
              + length(Y)*length(Kn)	  #	96
              + length(Y)	  #	97
              + length(Y)*length(Kn)	  #	98
              + length(Y)	  #	99
              + length(Y)*length(Kn)	  #	100
              + length(Y)	  #	101
              + length(Y)*length(Kn)	  #	102
              + length(Y)*length(Kn)	  #	103
              + length(Y)	  #	104
              + length(Y)	  #	105
              + length(Y)	  #	106
              + length(Y)	  #	107
              + length(Y)	  #	108
              + length(Y)*2	  #	109
              + length(Y)	  #	110
              + length(Y)	  #	111
              + length(Y)*2	  #	112
              + length(Y)*2	  #	113
              + length(Y)	  #	114
              + length(Y)	  #	115
              + length(Y)*2	  #	116
              + length(Y)	  #	117
              + length(Y)	  #	118
              + length(Y)*2	  #	119
              + length(Y)	  #	120
              + length(Y)	  #	121
              + length(Y)*2	  #	122
              + length(Y)	  #	123
              + length(Y)	  #	124
              + length(Y)*2	  #	125
              + 2	  #	126
             )
    @printf "Size v = \t %i\n" size_v
    xk = Array{Float64, 3}(undef, size_v, length(mv), 1)


    K = 1:length(mv)
    for k in K
        nz = 1
        m = mv[k]
        (i_, l_) = m[:_blk_ij]
        @printf "i = %i, j = %i \n" i_ l_
        # complicating variables
        vc, _, _ = genBlockVarVec(m, p, s, i_, l_)
        # other variables
        v = all_variables(m)
        # filter out the complicating variables
        v = filter(y->!in(y, vc), v)
        @printf "size of vc = %i\tsize of v = %i\n" length(vc) length(v)
        xk[:, k, 1] = value.(vcat(vc, v))
    end

    return xk
end

function blockObjAttach!(mv::Vector{JuMP.Model}, p::parms, s::sets,
        D::Matrix{Float64}, pi_::Vector{Float64})
    # generate the variable reference vector
    k = 1
    for i_ in s.P
        for l_ in s.L
            m = mv[k]
            vv, vp, nvb = genBlockVarVec(m, p, s, i_, l_)
            c0 = nvb * (k-1) + 1
            c1 = nvb * k
            println(c0, "\t", c1)
            D_ij = D[:, c0:c1]
            attachBlockObjective(m, p, s, D_ij, vv, pi_, i_, l_)
            k += 1
        end
    end
end


function primValVecLabels(mv::Vector{JuMP.Model}, p::parms, s::sets)
    # A single column vector for the primal values is considered.
    # Iteration by block is required.
    P = s.P
    Y = s.Y
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu = s.Fu
    #
    size_v = (
              length(Y)	  #	1
              + length(Y)*length(Kr)	  #	2
              + length(Y)	  #	3
              + length(Y)*2	  #	4
              + length(Y)	  #	5
              + length(Y)	  #	6
              + length(Y)*2	  #	7
              + 1	  #	8
              + length(Y)*2	  #	9
              + length(Y)	  #	10
              + length(Y)*length(Kr)	  #	11
              + length(Y)	  #	12
              + length(Y)*length(Kr)	  #	13
              + length(Y)	  #	14
              + length(Y)*length(Kr)	  #	15
              + length(Y)*length(Fu)	  #	16
              + length(Y)*length(Kr)*length(Fu)	  #	17
              + length(Y)	  #	18
              + length(Y)*length(Kr)	  #	19
              + length(Y)	  #	20
              + length(Y)*length(Kr)	  #	21
              + length(Y)	  #	22
              + length(Y)*length(Kr)	  #	23
              + length(Y)	  #	24
              + length(Y)*length(Kr)	  #	25
              + length(Y)	  #	26
              + length(Y)*length(Kr)	  #	27
              + length(Y)	  #	28
              + length(Y)*length(Kr)	  #	29
              + length(Y)	  #	30
              + length(Y)*length(Kr)	  #	31
              + length(Y)*length(Kr)	  #	32
              + length(Y)	  #	33
              + length(Y)	  #	34
              + length(Y)	  #	35
              + length(Y)	  #	36
              + length(Y)	  #	37
              + length(Y)	  #	38
              + length(Y)	  #	39
              + length(Y)	  #	40
              + length(Y)	  #	41
              + length(Y)	  #	42
              + length(Y)*length(Kr)	  #	43
              + length(Y)*2	  #	44
              + length(Y)	  #	45
              + length(Y)	  #	46
              + length(Y)	  #	47
              + length(Y)*length(Kr)	  #	48
              + length(Y)	  #	49
              + length(Y)	  #	50
              + length(Y)*2	  #	51
              + length(Y)	  #	52
              + length(Y)	  #	53
              + length(Y)	  #	54
              + length(Y)	  #	55
              + length(Y)	  #	56
              + length(Y)	  #	57
              + length(Y)	  #	58
              + length(Y)	  #	59
              + length(Y)	  #	60
              + length(Y)	  #	61
              + length(Y)*2	  #	62
              + length(Y)	  #	63
              + length(Y)	  #	64
              + length(Y)*length(Kn)	  #	65
              + 1	  #	66
              + length(Y)*length(Kn)	  #	67
              + length(Y)	  #	68
              + length(Y)*length(Kn)	  #	69
              + length(Y)*length(Kn)	  #	70
              + length(Y)	  #	71
              + length(Y)*length(Kn)	  #	72
              + length(Y)	  #	73
              + length(Y)	  #	74
              + length(Y)	  #	75
              + length(Y)	  #	76
              + length(Y)	  #	77
              + length(Y)*2	  #	78
              + length(Y)	  #	79
              + length(Y)	  #	80
              + length(Y)	  #	81
              + length(Y)	  #	82
              + length(Y)	  #	83
              + length(Y)	  #	84
              + length(Y)	  #	85
              + length(Y)*length(Kn)	  #	86
              + length(Y)*length(Fu)	  #	87
              + length(Y)*length(Kn)*length(Fu)	  #	88
              + length(Y)	  #	89
              + length(Y)*length(Kn)	  #	90
              + length(Y)	  #	91
              + length(Y)*length(Kn)	  #	92
              + length(Y)	  #	93
              + length(Y)*length(Kn)	  #	94
              + length(Y)	  #	95
              + length(Y)*length(Kn)	  #	96
              + length(Y)	  #	97
              + length(Y)*length(Kn)	  #	98
              + length(Y)	  #	99
              + length(Y)*length(Kn)	  #	100
              + length(Y)	  #	101
              + length(Y)*length(Kn)	  #	102
              + length(Y)*length(Kn)	  #	103
              + length(Y)	  #	104
              + length(Y)	  #	105
              + length(Y)	  #	106
              + length(Y)	  #	107
              + length(Y)	  #	108
              + length(Y)*2	  #	109
              + length(Y)	  #	110
              + length(Y)	  #	111
              + length(Y)*2	  #	112
              + length(Y)*2	  #	113
              + length(Y)	  #	114
              + length(Y)	  #	115
              + length(Y)*2	  #	116
              + length(Y)	  #	117
              + length(Y)	  #	118
              + length(Y)*2	  #	119
              + length(Y)	  #	120
              + length(Y)	  #	121
              + length(Y)*2	  #	122
              + length(Y)	  #	123
              + length(Y)	  #	124
              + length(Y)*2	  #	125
              + 2	  #	126
             )
    @printf "Size v = \t %i\n" size_v
    xk = Matrix{String}(undef, size_v, length(mv))


    K = 1:length(mv)
    for k in K
        nz = 1
        m = mv[k]
        (i_, l_) = m[:_blk_ij]
        @printf "i = %i, j = %i \n" i_ l_
        # complicating variables
        vc, _, _ = genBlockVarVec(m, p, s, i_, l_)
        # other variables
        v = all_variables(m)
        # filter out the complicating variables
        v = filter(y->!in(y, vc), v)
        @printf "size of vc = %i\tsize of v = %i\n" length(vc) length(v)
        xk[:, k] = name.(vcat(vc, v))
    end

    return xk
end

# rmp model.
# There are several situations in which this model is interacted with.
# 1. creation of the model from scratch
# 2. adding columns
# 3. creation of the model from a set of existing columns 


