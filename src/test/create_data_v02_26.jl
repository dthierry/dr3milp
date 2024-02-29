
using DataFrames
using CSV
using Printf

include("../mod/params_v02_10.jl")

# 1
fac_cap_f = "/Users/dthierry/Projects/dr3milp/src/data/state_cap.csv"
# 2
heat_f = "/Users/dthierry/Projects/dr3milp/src/data/heat_intensity.csv"
# 3
fuel_f = "/Users/dthierry/Projects/dr3milp/src/data/fuel_frac.csv"
# 4
elec_f = "/Users/dthierry/Projects/dr3milp/src/data/elec_int.csv"
# 5
em_f = "/Users/dthierry/Projects/dr3milp/src/data/fuel_emr.csv"
# 6
grid_f = "/Users/dthierry/Projects/dr3milp/src/data/processed_grid_co2.csv"

# 80 ###########################################################################
cap_sf = 1e3
d_cap = DataFrame(CSV.File(fac_cap_f))
d_cap[:, "ClinkerCap"] = d_cap[:, "ClinkerCap"]./cap_sf  # scaled


d_heat = DataFrame(CSV.File(heat_f))
d_heat[:, 2] = d_heat[:, 2].*cap_sf  # scaled


d_fuel = DataFrame(CSV.File(fuel_f))  # fractions

d_elec = DataFrame(CSV.File(elec_f))
d_elec[:, 2] = d_elec[:, 2].*cap_sf  # scaled

d_em = DataFrame(CSV.File(em_f))

d_grd = DataFrame(CSV.File(grid_f))

g_grd = groupby(d_grd, "State")

# number of locations
n_loc = 2 # nrow(d_cap) - 1
n_loc -= 1
loc_range = 1:n_loc+1
# number of fuels
n_fu = ncol(d_fuel) - 2
# number of year per period
n_years = 5 - 1
# number of periods
#n_periods = 6 - 1 # 0
n_periods = 2 - 1 # 0
# we need a couple retrofits
# 1: normal plant
# 2: carbon capture 30%
# 3: heat rate efficiency 4%
# 4: fuel switch coal to gas
# 5: cc and heat rate
n_rtft = 5-1 # 2

# fuel-switch merges bituminous, coalcoke, coke, and subbituminous into
# naturalgas

# we need a couple new plants
# 1: no plant 
# 2: carbon capture 90%
# 3: partial electrification
n_new = 3-1 # 3


r_loan0 = ones(n_loc+1) * 20 # 83
#r_loan0[2] = 1.
r_loan0[1] = 1.
#r_loan0[3] = 3e2
# 80 
# 80 ###########################################################################
##

#e_C = LinRange(1, 10., n_loc+1) # 5 capacity factor
e_C = floor.(d_cap[loc_range, "ClinkerCap"]/10) # 5 capacity factor

x_bM = 11. # 7
e_c_bM = maximum(e_C) * x_bM # 6

e_loanFact = ones(n_loc+1)*1e3 # 8


e_l_bM = e_c_bM * maximum(e_loanFact) # 9

e_Ann = e_loanFact./30 # 10

for i in 1:n_loc+1
    e_Ann[i] = 1.1^i
end

e_ann_bM = e_c_bM * maximum(e_Ann) # 11

e_ladd_bM = e_l_bM # 12

e_loan_bM = e_l_bM + 0.0 # 13
e_pay_bM = e_c_bM * maximum(e_Ann) # 14

# 80 
# 80 ###########################################################################
##

#c0 = LinRange(1, 10, n_loc+1) # 15
c0 = d_cap[loc_range, "ClinkerCap"] # 15

#r_c_C = ones(n_loc+1, n_rtft+1) # 16
r_c_C = ones(n_loc+1, n_rtft+1) # 16 # these do not change capacity

r_rhs_C = zeros(n_loc+1, n_rtft+1) # 17

r_cp_bM = e_c_bM + maximum(d_cap[loc_range, "ClinkerCap"]) # 18
r_cpb_bM = r_cp_bM # 19

# heat intensity
r_c_H = ones(n_loc+1, n_rtft+1) # 20
for k in 1:n_rtft+1
    r_c_H[:, k] .= d_heat[loc_range, 2]
end
r_rhs_H = zeros(n_loc+1, n_rtft+1) # 21


# hi factor 
r_c_Hfac = ones(n_loc+1, n_rtft+1)

r_c_Hfac[:, 2] .= 1.129  # increases
r_c_Hfac[:, 3] .= 0.96  # decreases
r_c_Hfac[:, 4] .= 1.05  # increases
r_c_Hfac[:, 5] .= 1.129 - 0.04

r_eh_bM = r_cp_bM * maximum(d_heat[loc_range, 2]) * maximum(r_c_Hfac) # 22

# electrification factor
r_c_Helec = zeros(n_loc+1, n_rtft+1)

# fuel fraction
r_c_F = ones(n_loc+1, n_fu+1, n_rtft+1) # 23

for i in 1:n_rtft+1
    r_c_F[:, :, i] .= d_fuel[loc_range, 2:end]
end
drf_f = copy(d_fuel)
# merge all fractions into rf natural gas
drf_f[:, "naturalgas"] = (drf_f[:, "naturalgas"] + 
                             drf_f[:, "bituminous"] +
                             drf_f[:, "coalcoke"] +
                             drf_f[:, "coke"] +
                             drf_f[:, "subbituminous"])

# zero out the rest
drf_f[:, "bituminous"] .= 0
drf_f[:, "coalcoke"] .= 0
drf_f[:, "coke"] .= 0
drf_f[:, "subbituminous"] .= 0

r_c_F[:, :, 4] .= drf_f[loc_range, 2:end]

r_rhs_F = zeros(n_loc+1, n_fu+1, n_rtft+1) # 24

r_ehf_bM = r_eh_bM # 25

# electricity intensity mmbtu/ton
r_c_U = ones(n_loc+1, n_rtft+1) # 26
for i in 1:n_rtft+1
    r_c_U[:, i] = d_elec[loc_range, 2]
end
r_rhs_U = zeros(n_loc+1, n_rtft+1) # 27

# five percent
r_c_UonSite = ones(n_periods+1,
                   n_years+1,
                   n_loc+1,
                   n_rtft+1).*(5/100)

r_c_Ufac = ones(n_loc+1, n_rtft+1)

r_c_Ufac[:, 2] .= 1.05  # increases
r_c_Ufac[:, 5] .= 1.05

r_u_bM = r_cp_bM * maximum(d_elec[loc_range, 2])*maximum(r_c_Ufac) # 28

# process emissions
r_c_cp_e = ones(n_loc+1, n_rtft+1).*cap_sf*501674/(1000)# 29
r_rhs_cp_e = zeros(n_loc+1, n_rtft+1) # 30
r_cp_e_bM = r_cp_bM*cap_sf*501674/(1000)  # 31

# fuel emission factor
r_c_Fe = zeros(n_fu+1, n_rtft+1) # 32
for i in 1:n_rtft+1
    r_c_Fe[:, i] .= d_em[:, 3]
end

# in-site elec fuel fraction
r_c_Fgenf = zeros(n_loc+1, n_fu+1, n_rtft+1)
for i in 1:n_rtft+1
    r_c_Fgenf[:, :, i] .= d_fuel[loc_range, 2:end]
end

# in-site generation heat rate, one mmbtu per mwh
# 1 MWh = 3.412142 MMBTU
r_c_Hr = ones(n_loc+1, n_fu+1, n_rtft+1)


r_ep0_bM = r_ehf_bM + r_cp_e_bM + r_u_bM*(5/100)*maximum(d_em[:, 3]) # 33

# capture fraction
r_chi = zeros(n_loc+1, n_rtft+1) # 34

r_chi[:, 2] .= 0.3
r_chi[:, 5] .= 0.3

r_ep1ge_bM = r_ep0_bM # 35

r_sigma = ones(n_loc+1, n_rtft+1)*0.5 # 36

r_ep1gce_bM = r_ep0_bM # 37
r_ep1gcs_bM = r_ep0_bM # 38

#r_c_Onm = ones(n_loc+1, n_rtft+1)
r_c_Onm = ones(n_loc+1, n_rtft+1)

if n_loc > 0
    for i in 1:n_rtft+1
        r_c_Onm[:, i] = LinRange(1.0, 1+i^4, n_loc+1)
    end
end

r_rhs_Onm = zeros(n_loc+1, n_rtft+1) # 40
r_conm_bM = r_cp_bM * maximum(r_c_Onm) # 41

r_loanFact = ones(n_periods+1, n_years+1, n_loc+1, n_rtft+1)*1e3 # 42
r_loanFact[:,:,:,2] .= 1e1
r_loanFact[:,:,:,3] .= 1e1
r_loanFact[:,:,:,4] .= 5e1
r_loanFact[:,:,:,5] .= 5e1


r_l_bM = r_cp_bM * maximum(r_loanFact) # 43

r_annf = ones(n_periods+1, n_years+1, n_loc+1, n_rtft+1) # 44
r_annf = r_loanFact./30

r_ann_bM = r_cp_bM * maximum(r_annf) # 45

r_ladd_bM = r_l_bM # 46

r_loan_bM = r_l_bM + maximum(r_loan0) # 47
r_pay_bM = r_ann_bM # 48

t_loan_bM = r_loan_bM + e_loan_bM # 50
t_ret_c_bM = t_loan_bM # 49

o_pay_bM = r_ann_bM + e_ann_bM  # 51
o_conm_bM = r_conm_bM  # 52

# 80 
# 80 ###########################################################################
##
#
n_c0_bM = c0.*2 # 54
n_c0_lo = c0 # 54
n_cp_bM = maximum(n_c0_bM) # 53

#n_loanFact = ones(n_loc+1, n_new+1) # 55
n_loanFact = ones(n_loc+1, n_new+1) # 55
n_loanFact[:, 2] .= 1e3
n_loanFact[:, 3] .= 1e4

n_l_bM = n_cp_bM * maximum(n_loanFact) # 56


#n_Ann = ones(n_loc+1, n_new+1) # 57
n_Ann = ones(n_loc+1, n_new+1) # 57
n_Ann = n_loanFact./30

if n_loc > 0
    for i in 1:n_new+1
        n_Ann[:, i] = collect(LinRange(1, 10, n_loc+1)).^i
    end
end

n_ann_bM = n_cp_bM * maximum(n_Ann) # 58

n_ladd_bM = n_l_bM # 59
n_loan_bM = n_l_bM # 60
n_pay_bM = n_ann_bM # 61

n_c_H = ones(n_loc+1, n_new+1) # 62
for k in 1:n_new+1
    n_c_H[:, k] .= d_heat[loc_range, 2]
end
n_rhs_h = zeros(n_loc+1, n_new+1)

# hi factor 
n_c_Hfac = ones(n_loc+1, n_new+1)
n_c_Hfac[:, 2] .= 1.448

n_eh_bM = n_cp_bM * maximum(d_heat[loc_range, 2])*maximum(n_c_Hfac)  # 63
# electrification
n_c_Helec = zeros(n_loc+1, n_new+1)
n_c_Helec[:, 3] .= 0.4


n_c_F = ones(n_loc+1, n_fu+1, n_new+1) # 64
for i in 1:n_new+1
    n_c_F[:, :, i] .= d_fuel[loc_range, 2:end]
end

n_rhs_F = zeros(n_loc+1, n_fu+1, n_new+1) # 65

n_ehf_bM = n_eh_bM + maximum(n_rhs_F)  # 66

n_c_U = ones(n_loc+1, n_new+1) # 67
for i in 1:n_new+1
    n_c_U[:, i] = d_elec[loc_range, 2]
end
n_rhs_U = zeros(n_loc+1, n_new+1) # 68

n_c_UonSite = ones(n_periods+1, 
                   n_years+1, 
                   n_loc+1, 
                   n_new+1).*(5/100)

n_c_Ufac = ones(n_loc+1, n_new+1)

n_c_Ufac[:, 2] .= 1.05  # increases
n_c_Ufac[:, 3] .= 1.05
n_u_bM = n_cp_bM * maximum(n_c_U) * maximum(n_c_Ufac) # 69

# 501674 gCO2/ton
n_c_cp_e = ones(n_loc+1, n_new+1).*cap_sf*501674/(1000) # 70
n_rhs_cp_e = zeros(n_loc+1, n_new+1) # 71

n_cp_e_bM = n_cp_bM * maximum(n_c_cp_e) # 72

n_c_Fe = ones(n_fu+1, n_new+1) # 73
for i in 1:n_new+1
    n_c_Fe[:, i] .= d_em[:, 3]
end

n_c_Fgenf = zeros(n_loc+1, n_fu+1, n_new+1) # 73
for i in 1:n_new+1
    n_c_Fgenf[:, :, i] .= d_fuel[loc_range, 2:end]
end

n_c_Hr = ones(n_loc+1, n_fu+1, n_new+1) # 73

n_ep0_bM = n_ehf_bM + n_cp_e_bM + n_u_bM*(5/100)*maximum(d_em[:, 3]) # 33

n_chi = zeros(n_loc+1, n_new+1)  # 75
n_chi[:, 2] .= 0.9 # 75

n_ep1ge_bM = n_ep0_bM # 76
n_sigma = ones(n_loc+1, n_new+1)*0.5 # 77

n_ep1gce_bM = n_ep0_bM  # 78
n_ep1gcs_bM = n_ep0_bM  # 79

#n_c_Onm = ones(n_loc+1, n_new+1) # 80
n_c_Onm = ones(n_loc+1, n_new+1) # 80
if n_loc > 0
    for i in 1:n_new+1
        n_c_Onm[:, i] = LinRange(1., 1.0+i^1.1, n_loc+1)
    end
end
n_rhs_Onm = zeros(n_loc+1, n_new+1) # 81

n_conm_bM = n_cp_bM * maximum(n_c_Onm) # 82

# 80
# 80 ###########################################################################
##

o_cp_bM = r_cp_bM # 82
o_cp_e_bM = 1e5 # 83 this one does not exists
o_u_bM = r_u_bM # 83
o_ep0_bM = r_ep0_bM


o_ep1ge_bM = r_ep1ge_bM
o_ep1gcs_bM = r_ep1gcs_bM
# 80
# 80 ###########################################################################
##


discount = ones(n_periods+1, n_years+1)

for i in 0:n_periods
    for j in 0:n_years
        year = j + i * (n_years+1)
        discount[i+1, j+1] = (1/(1+0.05)^(year))
    end
end

#discount = [1/(1+0.05)^t for t in 0:n_periods]
#discount = ones(n_periods+1) # 84


demand = ones(n_periods+1, n_years+1)

agg_dem0 = sum(d_cap[loc_range, :ClinkerCap])
d_growth = 1.2

demand_proj = LinRange(agg_dem0, agg_dem0 * d_growth,
                       (n_periods+1)*(n_years+1))

for i in 0:n_periods
    for j in 0:n_years
        year = j + i * (n_years+1)
        demand[i+1, j+1] = demand_proj[year+1]
    end
end

efuel = sum(c0.*r_c_H[:, 1].*r_c_F[:, :, 1].*r_c_Fe[: ,1]', dims=2)
eproc = c0.*r_c_cp_e[:, 1]
# onsite electricity mmbtu
el_ons = r_c_UonSite[1,1,:,1].*r_c_U[:, 1].*c0
eelins =sum(el_ons.*r_c_Hr[:,:,1].*r_c_Fgenf[:,:, 1].*r_c_Fe[:, 1]', dims=2)

el_grd = (1 .- r_c_UonSite[1,1,:,1]).*r_c_U[:, 1].*c0


GcI = zeros(n_periods+1, n_years+1, n_loc+1)
for i in 1:n_loc+1
    state = d_cap[i, "State"]
    for j in 0:n_periods
        for k in 0:n_years
            y = (n_years-1)*j + k
            GcI[j+1, k+1, i] = g_grd[(State=state,)][y+1, :Co2KgMwh]
        end
    end
end

# 1 MMBtu = 0.29329722222222 MWh
e_grd = el_grd.*GcI[1, 1, :].*0.29329722222222
co2_total0 = efuel .+ eproc .+ eelins .+ e_grd
#co2_total0 = 1.149872065065227e9
co2_reduction = 0.71
co2_endpoint = sum(co2_total0) * (1-co2_reduction)
co2_budget_v = collect(LinRange(sum(co2_total0), co2_endpoint, 30))

co2_budget = ones(n_periods+1, n_years+1)

for i in 0:n_periods
    for j in 0:n_years
        year = j + i * (n_years+1)
        co2_budget[i+1, j+1] = co2_budget_v[year+1]
    end
end


# 80 ###########################################################################
##

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
          r_c_Hfac,
          r_c_Helec,
          r_c_F,  # 23
          r_rhs_F,  # 24
          r_ehf_bM,  # 25
          r_c_U,  # 25
          r_rhs_U,  # 26 
          r_c_UonSite,  # 27
          r_c_Ufac,
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
          n_c0_lo,
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
          n_c_Hfac,
          n_c_Helec,
          n_c_F, # 65
          n_rhs_F, # 66
          n_ehf_bM, # 67
          n_c_U, # 68
          n_rhs_U, # 69
          n_c_UonSite,
          n_c_Ufac,
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

