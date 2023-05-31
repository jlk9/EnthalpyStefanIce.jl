module EnthalpyStefanIce

using Plots, Roots, SparseArrays, SpecialFunctions

# Write your package code here.

#= Runs a mushy sea ice thermodynamic model using the Stefan number

Inputs:
stefan  stefan number defined as Latent heat / sensible heat (L/cp*DeltaT), usually 10.0

=#
function runStefanModel(stefan)

    # quick finite volume
    L       = 2                         # box size 
    N       = 2^8
    dy      = L/N;                      # number of cells and cell size
    y       = range(start = dy/2, stop=L-(dy/2), step=dy) # location of cell centers
    y_edges = LinRange(0,L,N+1)         # location of cell edges
    k       = ones(length(y))           # thermal conductivity

    # TODO: can use collect on y, possible y_edges, i.e. y = collect(range(start = dy/2, stop=L-(dy/2), step=dy))

    # time information
    totaltime = 20              # total time
    dt        = (1/2.5)*(dy.^2) # timestep
    Nt        = round(Int, totaltime/dt)    # total timesteps

    kup    = k[2:N]                              # thermal conductivity on i+1
    kdown  = k[1:(N-1)]                          # thermal conductivity on i-1
    Trans  = (2/dy) ./ ((1 ./ kdown)+(1 ./ kup)) # transmissability as harmonic average 
    Trans1 = (2/dy)*k[1]
    TransN = (2/dy)*k[N]                         # transmissability at 1 and N
    Tup    = push!(deepcopy(Trans), TransN)                # i+1 transmissability to build matrix
    Tdown  = pushfirst!(deepcopy(Trans), Trans1)           # i-1 transmissability to build matrix

    # derivative matrix
    M_i = append!(collect(2:N), collect(1:N), collect(1:(N-1)))
    M_j = append!(collect(1:(N-1)), collect(1:N), collect(2:N))
    M_v = append!(deepcopy(Trans), deepcopy(-(Tup+Tdown)), deepcopy(Trans))
    M   = sparse(M_i, M_j, M_v);

    # boundary conditions vector
    bcv    = zeros(size(y))
    bcv[1] = -1*Trans1
    bcv[N] = 0*TransN

    # solve
    initial_enthalpy = stefan .* ones(size(y)); # start with enthalpy = Stefan (all liquid)
    # forward euler
    enthalpy = initial_enthalpy
    a = zeros(Nt)
    t = zeros(Nt) # initialize
    
    for i in 1:Nt
        t[i]     = i*dt; # time 
        theta    = min.(enthalpy,0); # temperature
        enthalpy = enthalpy+(dt./dy)*(M*theta+bcv); # step equation
        if minimum(enthalpy) <= 0
            a[i] = maximum(y[enthalpy .<= 0]); # determine location of ice-water interface
        end
    end

    return a, t, L
end

function plotStefanModel!(a, t, L, stefan)

    plot(t, a ./ L, label="enthalpy fvm", lc="black")
    plot!(t, sqrt.((2/stefan).*t)./L, label="large St scaling", lc=:blue)

    fλ(x) = (1 / stefan) - sqrt(π).*x.*exp.(x.^2).*erf.(x)
    λ     = find_zero(fλ, 1.0)

    plot!(t, (2λ/L).*sqrt.(t), label="stefan problem", lc=:red)

    title!("Enthalpy Stefan Comparison")
    xlabel!("Time, t")
    ylabel!("Ice thickness, a/L")

end

end