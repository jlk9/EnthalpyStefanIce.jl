module EnthalpyStefanIce

using LinearAlgebra, Plots, Roots, SpecialFunctions

# Write your package code here.

#= Builds a mushy sea ice thermodynamic model using the Stefan number

Inputs:
L           box size, usually 2
N           number of cells, usually 2^8
totaltime   total time usually 20

=#
function build_stefan_model(L, N, totaltime)

    # quick finite volume
    dy      = L/N;                      # cell size
    y       = range(start = dy/2, stop=L-(dy/2), step=dy) # location of cell centers
    y_edges = LinRange(0,L,N+1)         # location of cell edges
    k       = ones(length(y))           # thermal conductivity

    # TODO: can use collect on y, possible y_edges, i.e. y = collect(range(start = dy/2, stop=L-(dy/2), step=dy))

    # time information
    dt = (1/2.5)*(dy.^2) # timestep
    Nt = round(Int, totaltime/dt)    # total timesteps

    kup    = k[2:N]                              # thermal conductivity on i+1
    kdown  = k[1:(N-1)]                          # thermal conductivity on i-1
    Trans  = (2/dy) ./ ((1 ./ kdown)+(1 ./ kup)) # transmissability as harmonic average 
    Trans1 = (2/dy)*k[1]
    TransN = (2/dy)*k[N]                         # transmissability at 1 and N
    Tup    = push!(deepcopy(Trans), TransN)                # i+1 transmissability to build matrix
    Tdown  = pushfirst!(deepcopy(Trans), Trans1)           # i-1 transmissability to build matrix
    
    M   = SymTridiagonal(-(Tup+Tdown), Trans)
    M .*= dt/dy
    

    # boundary conditions pair
    bcv_start = -1*Trans1 * (dt/dy)
    bcv_end   = 0*TransN * (dt/dy)

    # solve using forward euler
    enthalpy = ones(size(y));

    θ = zeros(size(y)) # temperature
    a = zeros(Nt)
    t = collect(dt:dt:(Nt*dt)) #time

    return Nt, y, M, enthalpy, bcv_start, bcv_end, θ, a, t
end

#= Runs a mushy sea ice thermodynamic model using the Stefan number

Inputs:
Nt
M
enthalpy
bcv_start
bcv_end
θ
a
stefan

=#
function run_stefan_model!(N, Nt, y, M, enthalpy, bcv_start, bcv_end, θ, a, stefan)

    enthalpy .*= stefan # start with enthalpy = Stefan (all liquid)
    
    for i in 1:Nt
        θ .= min.(enthalpy,0); # temperature

        # step equation
        enthalpy    += M * θ;
        enthalpy[1] += bcv_start # add boundary conditions
        enthalpy[N] += bcv_end
        
        # determine location of ice-water interface
        #=
        y_temp = y[enthalpy .< 0]
        if length(y_temp) > 0
            a[i] = maximum(y_temp);
        end
        =#

        # This approach is less concise but faster, can swap to array filtering / maximum above if desired
        largest_negative = typemin(eltype(enthalpy))
        largest_negative_index = 0
        for j in 1:length(y)
            if enthalpy[j] < 0.0 && enthalpy[j] > largest_negative
                largest_negative = enthalpy[j]
                largest_negative_index = j
            end
        end

        if largest_negative_index != 0
            a[i] = y[largest_negative_index]
        end
        
    end
end

function plot_stefan_model(a, t, L, stefan)

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