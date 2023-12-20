struct sets
    T::UnitRange
    L::UnitRange
    Kr::UnitRange
    Kn::UnitRange
    Fu::UnitRange
    function sets(t_horizon::Int64, n_locations::Int64, 
            n_rtrf::Int64, n_new::Int64, n_fu::Int64)
        T = 0:t_horizon
        L = 0:n_locations
        Kr = 0:n_rtrf
        Kn = 0:n_new
        Fu = 0:n_fu
        new(T, L, Kr, Kn, Fu)
    end
end

