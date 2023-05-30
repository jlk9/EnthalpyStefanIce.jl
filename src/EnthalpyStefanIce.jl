module EnthalpyStefanIce

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
    Nt        = totaltime/dt    # total timesteps

    kup    = k[2:N] # thermal conductivity on i+1
    kdown  = k[1:(N-1)] # thermal conductivity on i-1
    Trans  = (2/dy) ./ ((1./kdown)+(1./kup)) # transmissability as harmonic average 
    Trans1 = (2/dy)*k[1]
    TransN = (2/dy)*k[N] # transmissability at 1 and N
    Tup    = [Trans, TransN] # i+1 transmissability to build matrix
    Tdown  = [Trans1, Trans] # i-1 transmissability to build matrix


    #M = sparse([2:N, 1:N, 1:(N-1)],[1:(N-1), 1:N, 2:N],[Trans, -(Tup+Tdown), Trans]); % derivative matrix
    #bcv = zeros(size(y))'; bcv(1) = -1*Trans1; bcv(N) = 0*TransN; % boundary conditions vectory

end

end
