mutable struct parms
    n_periods::Int64 # 0
    n_years::Int64 # 
    n_location::Int64 # 1
    n_rtft::Int64 # 2
    n_new::Int64 # 3
    n_fu::Int64 # 4
    
    # 76 
    # 76 #######################################################################
    ##
    e_C::Vector # 5 capacity factor for expansion [l]

    e_c_bM::Float64 # 6 *capacity* expansion big-M (capacity units)
    x_bM::Float64  # 7 *allocation* big-M 

    e_loanFact::Vector # 8 capacity expansion cost fact (cost units)
    e_l_bM::Float64  # 9 capacity cost big-M (cost u)

    e_Ann::Vector # 10 capacity expansion annuity [l] (cost units)
    e_ann_bM::Float64 # 11 capacity expansion annuity big-M (cost units)

    e_ladd_bM::Float64 # 12 ladd big-M (cost u)

    e_loan_bM::Float64 # 13 expansion loan big-M

    e_pay_bM::Float64 # 14

    # 76 
    # 76 #######################################################################
    ##

    c0::Vector # 15 initial capacity [l]

    r_c_C::Matrix # 16 *mode* capacity factor [l, k]
    r_rhs_C::Matrix # 17 mode capacity rhs

    r_cp_bM::Float64 # 18 mode capacity big-M (capacity units)
    r_cpb_bM::Float64 # 19 mode-base capacity big-M (capacity units)

    r_c_H::Matrix  # 20 mode heat factor [l, k] (heat u/cap u)
    r_rhs_H::Matrix  # 21 mode heat factor rhs (heat u)

    r_eh_bM::Float64  # 22 mode heat big-M (heat u)

    r_c_F::Array  # 23 mode fuel factor [l, k, f] (fuel u/heat u)
    r_rhs_F::Array  # 24 mode fuel rhs

    r_ehf_bM::Float64 # 25 mode fuel big-M

    r_c_U::Matrix  # 26 mode electricity requirement [l, k] (elec u/cap u)
    r_rhs_U::Matrix # 27
    
    r_c_UonSite::Array  # on-site capacity factor (0-1)

    r_u_bM # 28 mode elec big-M

    r_c_cp_e::Matrix # 29 process intrinsic emission factor [l, k] (em u/cap u)
    r_rhs_cp_e::Matrix # 30 process intrinsic emission rhs

    r_cp_e_bM::Float64 # 31 process intrinsic emission big-M

    r_c_Fe::Array  # 32 fuel emission factor [l, k, f] (em u/fu u)
    r_c_Fgenf::Array  # generation by fuel factor (0,1)
    r_c_Hr::Array  # fuel heat rate

    r_ep0_bM::Float64 # 33

    r_chi::Matrix # 34 captured [l, k]
    r_ep1ge_bM::Float64 # 35 fml

    r_sigma::Matrix # 36 stored emissions [l, k]

    r_ep1gce_bM::Float64 # 37
    r_ep1gcs_bM::Float64 # 38

    r_c_Onm::Matrix # 39
    r_rhs_Onm::Matrix # 40
    r_conm_bM::Float64 # 41

    r_loanFact::Array # 42
    r_l_bM::Float64  # 43 mode loan big-M

    r_annf::Array # 44 annuity factor for mode  [t, l, k]
    r_ann_bM::Float64 # 45 mode annuity big-M

    r_ladd_bM::Float64 # 46 mode loan-add big-M

    r_loan_bM::Float64 # 47
    r_pay_bM::Float64 # 48 mode payment big-M

    t_ret_c_bM::Float64 # 49 total retirement cost big-M

    t_loan_bM::Float64 # 50 total loan (for retirement) big-M

    o_pay_bM::Float64 # 51
    o_conm_bM::Float64 # 52

    # 76 
    # 76 #######################################################################
    ##
    
    n_cp_bM # 53
    n_c0_bM # 54
    n_loanFact::Matrix # 55
    n_l_bM # 56

    n_Ann::Matrix # 57
    n_ann_bM # 58
    n_ladd_bM # 59
    n_loan_bM # 60
    n_pay_bM # 61

    n_c_H::Matrix # 62 [l, k]
    n_rhs_h::Matrix # 
    n_eh_bM # 63

    n_c_F::Array # 64 [l, k, f]
    n_rhs_F::Array # 65 [l, k, f]
    n_ehf_bM::Float64 # 66

    n_c_U::Matrix # 67
    n_rhs_U::Matrix # 68
    
    n_c_UonSite::Array

    n_u_bM::Float64 # 69

    n_c_cp_e::Matrix # 70
    n_rhs_cp_e::Matrix # 71
    n_cp_e_bM # 72

    n_c_Fe::Array # 73 [l, k, f]
    n_c_Fgenf::Array  # [l, k, f]
    n_c_Hr::Array  # [l, k, f]

    n_ep0_bM # 74
    n_chi::Matrix # 75

    n_ep1ge_bM # 76
    n_sigma::Matrix # 77
    n_ep1gce_bM # 78

    n_ep1gcs_bM # 79

    n_c_Onm::Matrix # 80
    n_rhs_Onm::Matrix # 81

    n_conm_bM # 82
    
    # 76 
    # 76 #######################################################################
    ##

    o_cp_bM::Float64 # 83
    o_cp_e_bM::Float64 # 83
    o_u_bM::Float64
    o_ep0_bM::Float64 # 83

    o_ep1ge_bM::Float64
    o_ep1gcs_bM::Float64


    # 76 
    # 76 #######################################################################
    ##

    r_loan0::Vector # 83
    
    discount::Matrix # 84

    demand::Matrix # 85
    co2_budget::Matrix # 86
    GcI::Array  # 87

end
