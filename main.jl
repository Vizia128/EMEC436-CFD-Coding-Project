using LinearAlgebra, CairoMakie, Observables, Interpolations, JLD, Dates

struct Parameters
    L # length in the x and y directions
    t_end # end time
    t_step_max # maximum time steps
    t_record # recording time scale

    U_dim # non-dimensional x-velocity
    μ # dynamic viscosity
    ρ # density
    Re # Reynolds number

    τ # stability factor for timestep change
    ϵ
    ω # system of overrelaxation step
    itermax
    BC # boundary condition switch

    n # x-divisions
    m # y-divisions
    dx
    dy
    dxi
    dyi
    dxi2
    dyi2
end

function Parameters(;L = 1.0, t_end = 35, t_step_max = 128, U_dim = 1.0, 
        μ = 1.0, ρ = 1000.0, τ = 0.5, ϵ = 0.001, ω = 1.7,
        itermax = 100, BC = 1, n = 25, m = 25, t_record = 1)
    Re = ρ*U_dim*L/μ

    dx, dy = L/n, L/m
    dxi, dyi = 1/dx, 1/dy
    dxi2, dyi2 = dxi^2, dyi^2
    return Parameters(L, t_end,t_step_max, t_record, U_dim, μ, ρ, Re, τ, ϵ, ω, itermax, BC, n, m, dx, dy, dxi, dyi, dxi2, dyi2)
end

function comp_delt(U, V, divg, params)
    (; τ, Re, ϵ, dx, dy, dxi2, dyi2) = params
    return τ*max(1.0e-9, min(0.5*Re/(dxi2 + dyi2), dx/maximum(abs.(U)), dy/maximum(abs.(V))))/max(2log10(divg/ϵ)+1, 1)
end

function placeBCs(U, V, params)
    (; n, m, U_dim) = params
    # i == 1, left wall
    V[1,:] = -V[2,:] # stationary
    U[1,:] = zeros(size(U, 2)) # no permeation
    # i == n+1, right wall
    U[n+1,:] = zeros(size(U, 2)) # no permeation
    # i == n+2, right wall
    V[n+2,:] = -V[n+1,:] # stationary

    # j == 1, bottom wall
    U[:,1] = -U[:,2] # stationary
    V[:,1] = zeros(size(V, 1)) # no permeation
    # j == m+1, top wall
    V[:,m+1] = zeros(size(V, 1)) # no permeation
    # j == m+2, top wall
    U[:,m+2] = 2*U_dim .- U[:,m+1]; # moving in the x-direction

    return U, V
end

function solveFG(U, V, Us, Vs, dt, params)
    (; n, m, dx, dy, μ, ρ, dxi, dyi, dxi2, dyi2) = params

    for i in 2:n, j in 2:m+1
        d2udx2 = dxi2*(U[i-1,j] - 2*U[i,j] + U[i+1,j]) # second derivative of u in x
        d2udy2 = dyi2*(U[i,j-1] - 2*U[i,j] + U[i,j+1]) # second derivative of u in y
        du2dx = dxi*(((U[i,j] + U[i+1,j])/2)^2 - ((U[i-1,j] + U[i,j])/2)^2) # first derivative of u^2 in x
        duvdy = dyi*((U[i,j] + U[i,j+1])*(V[i,j] + V[i+1,j])/4 - (U[i,j-1] + U[i,j])*(V[i,j-1] + V[i+1,j-1])/4) # first derivative of uv in y
        
        Fij = μ/ρ*(d2udx2 + d2udy2) - du2dx - duvdy
        Us[i,j] = U[i,j] + Fij*dt # predict the updated u velocity
    end

    for i in 2:n+1, j in 2:m
        d2vdx2 = dxi2*(V[i-1,j] - 2*V[i,j] + V[i+1,j]) # second derivative of v in x
        d2vdy2 = dyi2*(V[i,j-1] - 2*V[i,j] + V[i,j+1]) # second derivative of v in y
        dv2dx = dyi*(((V[i,j] + V[i,j+1])/2)^2 - ((V[i,j-1] + V[i,j])/2)^2) # first derivative of u^2 in x
        dvudy = dxi*((U[i,j] + U[i,j+1])*(V[i,j] + V[i+1,j])/4 - (U[i-1,j] + U[i-1,j+1])*(V[i-1,j] + V[i,j])/4) # first derivative of uv in y
        
        Gij = μ/ρ*(d2vdx2 + d2vdy2) - dv2dx - dvudy
        Vs[i,j] = V[i,j] + Gij*dt # predict the updated u velocity  
    end

    # update the boundary cells of Us and Vs, having arbitrarily chosen the velocities

    Us[1,2:m+1] = U[1,2:m+1]
    Us[n+1,2:m+1] = U[n+1,2:m+1]
    Vs[2:n+1,1] = V[2:n+1,1]
    Vs[2:n+1,m+1] = V[2:n+1,m+1]

    # Us[1,:] = U[1,:]
    # Us[n+1,:] = U[n+1,:]
    # Vs[:,1] = V[:,1]
    # Vs[:,m+1] = V[:,m+1]

    return Us, Vs
end

function solvePoisson(Us, Vs, P, dt, params)
    (; ρ, n, m, dxi, dyi, dxi2, dyi2, ϵ, ω, itermax) = params
    dti = 1/dt
    rhs::Matrix{Float64} = zeros(size(P))
    for i in 2:n+1, j in 2:m+1
        rhs[i-1, j-1] = ρ*dti*(dxi*(Us[i,j] - Us[i-1,j]) + dyi*(Vs[i,j] - Vs[i,j-1]))
    end

    iter = 0
    divg = 1
    while iter ≤ itermax && divg ≥ ϵ
        for i in 2:n+1, j in 2:m+1
            P[i,j] = (1-ω)*P[i,j] + ω/(2*(dxi2+dyi2))*(dxi2*(P[i-1,j] + P[i+1,j]) + dyi2*(P[i,j-1] + P[i,j+1]) - rhs[i-1,j-1])
        end

        # update the pressure boundary conditions in the ghost cells
        P[1,:] = P[2,:]
        P[n+2,:] = P[n+1,:]
        P[:,1] = P[:,2]
        P[:,m+2] = P[:,m+1]

        # check for convergence
        divg_ij_sum = 0
        for i in 2:n+1, j in 2:m+1
            divg_ij_sum += (dxi2*(P[i-1,j] - 2P[i,j] + P[i+1,j]) + dxi2*(P[i,j-1] - 2P[i,j] + P[i,j+1]) - rhs[i-1,j-1])^2
        end
        divg = divg_ij_sum / (n*m)
        iter += 1
    end
    return P, divg
end

function solveVel(U, V, Us, Vs, P, dt, params)
    (; n, m, ρ, dxi, dyi) = params
    for i in 2:n, j in 2:m+1
        U[i,j] = Us[i,j] - dt/ρ*dxi*(P[i+1,j] - P[i,j])
    end
    for i in 2:n+1, j in 2:m
        V[i,j] = Vs[i,j] - dt/ρ*dyi*(P[i,j+1] - P[i,j])
    end
    return U, V
end

function create_velocity_func(Ui, Vi, t)
    # Uinterp = interpolate(Ut[:,:,t], BSpline(Linear()))
    # Vinterp = interpolate(Vt[:,:,t], BSpline(Linear()))

    function velocity(x,y)
        return Point2(Ui[x,y,t], Vi[x,y,t])
    end
    return velocity
end

function postProc(Ut, Vt, Pt, params)
    (;t_end, n, m) = params
    t_index = Observable(1)
    sliceU = @lift(Ut[:,:, $t_index])
    sliceV = @lift(Vt[:,:, $t_index])
    sliceP = @lift(Pt[:,:, $t_index])
    sliceUV = @lift(sqrt.(Ut[:,1:end-1, $t_index].^2 + Vt[1:end-1,:, $t_index].^2))

    Ui = interpolate(Ut, BSpline(Linear()))
    Vi = interpolate(Vt, BSpline(Linear()))

    f = Figure()

    _, au = heatmap(f[1,1], sliceU; axis = (; title = "U"))
    Colorbar(f[1,2], au; vertical=true)
    _, av = heatmap(f[1,3], sliceV; axis = (; title = "V"))
    Colorbar(f[1,4], av; vertical=true)
    _, ap = heatmap(f[2,1], sliceP; axis = (; title = "P"))
    Colorbar(f[2,2], ap; vertical=true)
    _, auv = streamplot(f[2,3], @lift(create_velocity_func(Ui, Vi, $t_index)), 1:n+1, 1:m+1; axis = (; title = "(U, V)"))

    sl = Slider(f[4, 1:4], horizontal = true, range = 1:t_end)
    connect!(t_index, sl.value)

    return f
end

function postProc(dir::String)
    dict = load(dir)
    fig = postProc(dict["Ut"], dict["Vt"], dict["Pt"], dict["params"])
    return fig
end

function fluidsolve(params::Parameters)
    (; n, m, t_end, t_record, Re) = params

    datetime = replace("$(now())", ":" => ".")
    isdir("output") ? true : mkdir("output")
    dir = mkdir("output/n=$n , Re=$Re , t_end=$t_end , $datetime")

    U::Matrix{Float64} = zeros(n + 1, m + 2)
    V::Matrix{Float64} = zeros(n + 2, m + 1)
    P::Matrix{Float64} = zeros(n + 2, m + 2)

    Us::Matrix{Float64} = zeros(n + 1, m + 2)
    Vs::Matrix{Float64} = zeros(n + 2, m + 1)

    Ut::Array{Float64} = zeros(n + 1, m + 2, Int(ceil(t_end/t_record)))
    Vt::Array{Float64} = zeros(n + 2, m + 1, Int(ceil(t_end/t_record)))
    Pt::Array{Float64} = zeros(n + 2, m + 2, Int(ceil(t_end/t_record)))
    
    t = 0.
    t_step::Int = 0
    t_slice::Int = 1
    divg = 1
    while t < t_end && divg < 1e100
        dt = comp_delt(U, V, divg, params)
        U, V = placeBCs(U, V, params)
        Us, Vs = solveFG(U, V, Us, Vs, dt, params)
        P, divg = solvePoisson(Us, Vs, P, dt, params)
        U, V = solveVel(U, V, Us, Vs, P, dt, params)

        t += dt
        t_step += 1
        if abs(t_slice*t_record - t) < t_record/10
            Ut[:,:,t_slice] = U
            Vt[:,:,t_slice] = V
            Pt[:,:,t_slice] = P
            t_slice += 1   
            @save "$dir/Ut.jld" Ut Vt Pt params
        end
        println("t_step: $t_step , t: $(round(t; digits=6)) , divg: $(round(divg; digits=6)) , maxP: $(round(maximum(P); digits=6)) , maxU: $(round(maximum(U); digits=6)) , maxV: $(round(maximum(V); digits=6))")
    end
    fig = postProc(Ut, Vt, Pt, params)

    return fig
end

function run(;L = 1.0, t_end = 35, t_step_max = 128, U_dim = 1.0, 
        μ = 1.0, ρ = 1000.0, τ = 0.5, ϵ = 0.001, ω = 1.7, 
        mu = 1.0, rho = 1000.0, tau = 0.5, epsilon = 0.001, omega = 1.7, 
        itermax = 100, BC = 1, n = 25, m = 25, t_record = 1)
    μ, ρ, τ, ϵ, ω = mu, rho, tau, epsilon, omega
    params = Parameters(;L, t_end, t_step_max, U_dim, μ, ρ, τ, ϵ, ω, itermax, BC, n, m, t_record)
    return fluidsolve(params)
end