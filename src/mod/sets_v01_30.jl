
struct sets
    P::UnitRange  # periods
    Y::UnitRange  # years per period
    L::UnitRange
    Kr::UnitRange
    Kn::UnitRange
    Fu::UnitRange
    function sets(
            n_periods::Int64,
            n_years::Int64,  # years per period
            #t_horizon::Int64, 
            n_locations::Int64, 
            n_rtrf::Int64, n_new::Int64, n_fu::Int64)
        P = 0:n_periods
        Y = 0:n_years
        L = 0:n_locations
        Kr = 0:n_rtrf
        Kn = 0:n_new
        Fu = 0:n_fu
        new(P, Y, L, Kr, Kn, Fu)
    end
end

