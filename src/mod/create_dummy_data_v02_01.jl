
include("./params_v02_01.jl")

n_periods = 4 # 0
n_years = 4
n_loc = 2 # 1
n_rtft = 5 # 2
n_new = 3 # 3
n_fu = 10 # 4


# 80 
# 80 ###########################################################################
##

#e_C = LinRange(1, 10., n_loc+1) # 5 capacity factor
e_C = LinRange(1, 10., n_loc+2) # 5 capacity factor

e_c_bM = 1e5 # 6
x_bM = 1e2 # 7

e_loanFact = ones(n_loc+1) # 8

for i in 1:n_loc+1
    e_loanFact[i] = 1.1^i
end

e_l_bM = 1e5 # 9

e_Ann = ones(n_loc+1) # 10

for i in 1:n_loc+1
    e_Ann[i] = 1.1^i
end

e_ann_bM = 1e5 # 11

e_ladd_bM = 1e5 # 12

e_loan_bM = 1e5 # 13
e_pay_bM = 1e5 # 14

# 80 
# 80 ###########################################################################
##

#c0 = LinRange(1, 10, n_loc+1) # 15
c0 = LinRange(1, 10, n_loc+2) # 15

#r_c_C = ones(n_loc+1, n_rtft+1) # 16
r_c_C = ones(n_loc+2, n_rtft+1) # 16

for i in 1:n_rtft+1
    #r_c_C[:, i] = LinRange(1.0, 1.0 + i^1.1, n_loc+1)
    r_c_C[:, i] = LinRange(1.0, 1.0 + i^1.1, n_loc+2)
end

r_rhs_C = zeros(n_loc+1, n_rtft+1) # 17

r_cp_bM = 1e5 # 18
r_cpb_bM = 1e5 # 19

r_c_H = ones(n_loc+1, n_rtft+1)*10. # 20
r_rhs_H = zeros(n_loc+1, n_rtft+1) # 21

r_eh_bM = 1e5 # 22

r_c_F = ones(n_loc+1, n_rtft+1, n_fu+1) # 23
r_rhs_F = zeros(n_loc+1, n_rtft+1, n_fu+1) # 24

r_ehf_bM = 1e5 # 25

r_c_U = ones(n_loc+1, n_rtft+1) # 26
r_rhs_U = zeros(n_loc+1, n_rtft+1) # 27

r_c_UonSite = ones(n_periods+1,
                   n_years+1,
                   n_loc+1,
                   n_rtft+1).*(5/100)

r_u_bM = 1e5 # 28

r_c_cp_e = ones(n_loc+1, n_rtft+1) # 29
r_rhs_cp_e = zeros(n_loc+1, n_rtft+1) # 30

r_cp_e_bM = 1e5 # 31

r_c_Fe = ones(n_loc+1, n_rtft+1, n_fu+1) # 32

fracs = ones(n_fu+1).*(1/(n_fu+1))
r_c_Fgenf = ones(n_loc+1, n_rtft+1, n_fu+1)


for i in 1:n_loc+1
    for j in 1:n_rtft+1
        r_c_Fgenf[i, j, :] .= fracs
    end
end

r_c_Hr = ones(n_loc+1, n_rtft+1, n_fu+1).*0.5

r_ep0_bM = 1e5 # 33

r_chi = ones(n_loc+1, n_rtft+1)*0.5 # 34

r_ep1ge_bM = 1e5 # 35

r_sigma = ones(n_loc+1, n_rtft+1)*0.5 # 36

r_ep1gce_bM = 1e5 # 37
r_ep1gcs_bM = 1e5 # 38

#r_c_Onm = ones(n_loc+1, n_rtft+1)
r_c_Onm = ones(n_loc+2, n_rtft+1)

for i in 1:n_rtft+1
    # r_c_Onm[:, i] = LinRange(1.0, 1+i^4, n_loc+1)
    r_c_Onm[:, i] = LinRange(1.0, 1+i^4, n_loc+2)
end


r_rhs_Onm = ones(n_loc+1, n_rtft+1) # 40
r_conm_bM = 1e6 # 41

r_loanFact = ones(n_periods+1, n_years+1, n_loc+1, n_rtft+1)*100 # 42

r_l_bM = 1e5 # 43

r_annf = ones(n_periods+1, n_years+1, n_loc+1, n_rtft+1) # 44
r_ann_bM = 1e6 # 45

r_ladd_bM = 1e5 # 46

r_loan_bM = 1e5 # 47
r_pay_bM = 1e5 # 48

t_ret_c_bM = 1e5 # 49
t_loan_bM = 1e5 # 50

o_pay_bM = 1e5 # 51
o_conm_bM = 1e5 # 52

# 80 
# 80 ###########################################################################
##
#
n_cp_bM = 1e5 # 53
n_c0_bM = 1e5 # 54
#n_loanFact = ones(n_loc+1, n_new+1) # 55
n_loanFact = ones(n_loc+2, n_new+1) # 55
for i in 1:n_new+1
    #n_loanFact[:, i] = LinRange(1.0, 1+i^1.1, n_loc+1)
    n_loanFact[:, i] = LinRange(1.0, 1+i^1.1, n_loc+2)
end

n_l_bM = 1e5 # 56



#n_Ann = ones(n_loc+1, n_new+1) # 57
n_Ann = ones(n_loc+2, n_new+1) # 57

for i in 1:n_new+1
    #n_Ann[:, i] = collect(LinRange(1, 10, length=n_loc+1)).^i
    n_Ann[:, i] = collect(LinRange(1, 10, n_loc+2)).^i
end


n_ann_bM = 1e7 # 58

n_ladd_bM = 1e5 # 59
n_loan_bM = 1e5 # 60
n_pay_bM = 1e5 # 61

n_c_H = ones(n_loc+1, n_new+1) # 62
n_rhs_h = ones(n_loc+1, n_new+1)
n_rhs_h[:, 1] .= 0.0
n_eh_bM = 1e5 # 63

n_c_F = ones(n_loc+1, n_new+1, n_fu+1) # 64
n_rhs_F = ones(n_loc+1, n_new+1, n_fu+1) # 65
n_rhs_F[:, 1, :] .= 0.0

n_ehf_bM = 1e5 # 66

n_c_U = ones(n_loc+1, n_new+1) # 67
n_rhs_U = ones(n_loc+1, n_new+1) # 68
n_rhs_U[:, 1] .= 0.0

n_c_UonSite = ones(n_periods+1, 
                   n_years+1, 
                   n_loc+1, 
                   n_new+1).*(5/100)

n_u_bM = 1e5 # 69

n_c_cp_e = ones(n_loc+1, n_new+1) # 70
n_rhs_cp_e = ones(n_loc+1, n_new+1) # 71
n_rhs_cp_e[:, 1] .= 0.0

n_cp_e_bM = 1e5 # 72

n_c_Fe = ones(n_loc+1, n_new+1, n_fu+1) # 73

fracs = ones(n_fu+1).*(1/(n_fu+1))
n_c_Fgenf = ones(n_loc+1, n_new+1, n_fu+1) # 73

for i in 1:n_loc+1
    for j in 1:n_new+1
        n_c_Fgenf[i, j, :] .= fracs
    end
end

n_c_Hr = ones(n_loc+1, n_new+1, n_fu+1).*0.5 # 73

n_ep0_bM = 1e5 # 74

n_chi = ones(n_loc+1, n_new+1)*0.5 # 75

n_ep1ge_bM = 1e5 # 76
n_sigma = ones(n_loc+1, n_new+1)*0.5 # 77
n_ep1gce_bM = 1e5 # 78

n_ep1gcs_bM = 1e5 # 79

#n_c_Onm = ones(n_loc+1, n_new+1) # 80
n_c_Onm = ones(n_loc+2, n_new+1) # 80

for i in 1:n_new+1
    #n_c_Onm[:, i] = LinRange(1., 1.0+i^1.1, n_loc+1)
    n_c_Onm[:, i] = LinRange(1., 1.0+i^1.1, n_loc+2)
end

n_rhs_Onm = ones(n_loc+1, n_new+1) # 81
n_rhs_Onm[:, 1] .= 0.0

n_conm_bM = 1e5 # 82

# 80
# 80 ###########################################################################
##

o_cp_bM = 1e5 # 82
o_cp_e_bM = 1e5 # 83
o_u_bM = 1e5 # 83
o_ep0_bM = 1e5

o_cp_bM = 1e5 # 82
o_cp_e_bM = 1e5 # 83
o_ep0_bM = 1e5

o_cp_bM = 1e5 # 82
o_cp_e_bM = 1e5 # 83
o_ep0_bM = 1e5
o_ep1ge_bM = 1e5
o_ep1gcs_bM = 1e5
# 80
# 80 ###########################################################################
##

r_loan0 = ones(n_loc+1) * 20 # 83
#r_loan0[2] = 1.
r_loan0[1] = 1.
#r_loan0[3] = 3e2

discount = ones(n_periods+1, n_years+1)

for i in 0:n_periods
    for j in 0:n_years
        year = j + i * n_years
        discount[i+1, j+1] = (1/(1+0.05)^(year))
    end
end

#discount = [1/(1+0.05)^t for t in 0:n_periods]
#discount = ones(n_periods+1) # 84


demand = ones(n_periods+1, n_years+1)
d0 = range(n_loc, n_loc+10, length=(n_periods+1)*(n_years+1))
d0 = collect(d0)

for i in 0:n_periods
    for j in 0:n_years
        year = j + i * n_years
        demand[i+1, j+1] = d0[year+1]
    end
end

co2_budget = ones(n_periods+1, n_loc+1).*1e5


GcI = ones(n_periods+1, n_years+1, n_loc+1)

p = parms(n_periods, # 0 
          n_years,
          n_loc, # 1
          n_rtft, # 2
          n_new, # 3
          n_fu,  # 4 
          e_C,  # 5
          e_c_bM,  # 6
          x_bM,  # 7
          e_loanFact,  # 8
          e_l_bM,  # 9
          e_Ann,  # 10
          e_ann_bM,   # 11
          e_ladd_bM,   # 12
          e_loan_bM,   # 13
          e_pay_bM,   # 14
          c0,   # 15
          r_c_C,   # 16
          r_rhs_C,   # 17
          r_cp_bM,   # 18
          r_cpb_bM,  # 19
          r_c_H,  # 20
          r_rhs_H,  # 21
          r_eh_bM,  # 22
          r_c_F,  # 23
          r_rhs_F,  # 24
          r_ehf_bM,  # 25
          r_c_U,  # 25
          r_rhs_U,  # 26 
          r_c_UonSite,  # 27
          r_u_bM,  # 28
          r_c_cp_e,  # 29
          r_rhs_cp_e,  # 30
          r_cp_e_bM,  # 31
          r_c_Fe,  # 32
          r_c_Fgenf,  # 33
          r_c_Hr,  # 34
          r_ep0_bM,  # 35
          r_chi,  # 36
          r_ep1ge_bM,  # 37
          r_sigma,  # 38
          r_ep1gce_bM, # 39
          r_ep1gcs_bM, # 40
          r_c_Onm, # 41
          r_rhs_Onm, # 42
          r_conm_bM, # 43
          r_loanFact, # 44
          r_l_bM, # 45
          r_annf, # 46
          r_ann_bM, # 47
          r_ladd_bM, # 48
          r_loan_bM, # 49
          r_pay_bM, # 50
          t_ret_c_bM, # 50
          t_loan_bM, # 51
          o_pay_bM, # 52
          o_conm_bM, # 53
          n_cp_bM, # 54
          n_c0_bM, # 55
          n_loanFact, # 56
          n_l_bM, # 57
          n_Ann, # 58
          n_ann_bM, # 59
          n_ladd_bM, # 60
          n_loan_bM, # 61
          n_pay_bM, # 62
          n_c_H, # 63
          n_rhs_h, # 
          n_eh_bM, # 64
          n_c_F, # 65
          n_rhs_F, # 66
          n_ehf_bM, # 67
          n_c_U, # 68
          n_rhs_U, # 69
          n_c_UonSite,
          n_u_bM, # 70
          n_c_cp_e, # 71
          n_rhs_cp_e, # 72
          n_cp_e_bM, # 73
          n_c_Fe, # 74
          n_c_Fgenf,
          n_c_Hr,
          n_ep0_bM, # 75
          n_chi, # 76
          n_ep1ge_bM, # 77
          n_sigma, # 78
          n_ep1gce_bM, # 79
          n_ep1gcs_bM, # 80
          n_c_Onm, # 81
          n_rhs_Onm, # 82
          n_conm_bM, # 83
          o_cp_bM,
          o_cp_e_bM,
          o_u_bM,
          o_ep0_bM,
          o_ep1ge_bM,
          o_ep1gcs_bM,
          r_loan0, # 84
          discount, # 85
          demand,
          co2_budget,
          GcI
         )
