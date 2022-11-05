using CairoMakie, JLD, CSV, DataFrames, Interpolations

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

function loadsolution(dir::String)
    dict = load(dir)
    return dict["Ut"], dict["Vt"], dict["Pt"], dict["params"]
end

function Xt_to_center_velocity(dir::String)
    Ut, Vt, Pt, params = loadsolution(dir)
    U_center = Ut[Int(floor(end/2)),:,end]
    V_center = Vt[:,Int(floor(end/2)),end]
    return U_center, V_center
end

function CSV_to_center_velocity(dir::String)
    df = CSV.read(dir, DataFrame)
    return df[:,2]
end

function center_velocity_compare(zia_dir::String, star_dir_U::String, star_dir_V::String, 
                                ghia_dir_U::String, ghia_dir_V::String, Re::Real)
    U_zia, V_zia = Xt_to_center_velocity(zia_dir)

    U_star = CSV_to_center_velocity(star_dir_U)
    V_star = CSV_to_center_velocity(star_dir_V)

    U_ghia = CSV_to_center_velocity(ghia_dir_U)
    V_ghia = CSV_to_center_velocity(ghia_dir_V)

    ghia_index_U = [1, 8, 9, 10, 14, 23, 37, 59, 65, 80, 95, 110, 123, 124, 125, 126, 129]
    ghia_index_V = [1, 9, 10, 11, 13, 21, 30, 31, 65, 104, 111, 117, 122, 123, 124, 125, 129]

    U_ghia = [ghia_index_U U_ghia]
    V_ghia = [ghia_index_V V_ghia]

    fig = Figure(resolution = (800, 450))

    lines(fig[1,1], U_zia; label = "Zia", color = :blue, with = 2, axis = (; title = "U", ylabel = "Velocity", xlabel = "Y Cell Position"))
    lines!(fig[1,1], U_star; label = "Star", color = :red, linestyle = :dash, width = 2)
    scatter!(fig[1,1], U_ghia; label = "Ghia", color = :transparent, marker = :circle, markersize = 12, strokewidth = 1, strokecolor = :green)
    # axislegend()

    lines(fig[1,2], V_zia; label = "Zia", color = :blue, with = 2, axis = (; title = "V", xlabel = "X Cell Position"))
    lines!(fig[1,2], V_star; label = "Star", color = :red, linestyle = :dash, width = 2)
    scatter!(fig[1,2], V_ghia; label = "Ghia", color = :transparent, marker = :circle, markersize = 12, strokewidth = 1, strokecolor = :green)
    axislegend()

    return fig
end

function time_slice_plot(zia_dir::String, t_slice::Int)
    Ut, Vt, Pt, params = loadsolution(zia_dir)
    (;t_record, n, m, Re) = params
    time = t_slice * t_record

    # U = Ut[4:end-3,4:end-3,t_slice]
    # V = Vt[4:end-3,4:end-3,t_slice]
    # P = Pt[4:end-3,4:end-3,t_slice]    
    U = Ut[:,:,t_slice]
    V = Vt[:,:,t_slice]
    P = Pt[:,:,t_slice]

    f = Figure(resolution = (1200, 1000))

    # Label(f[0, 1:3], "t = $(round(time, digits=1))", textsize=40)

    _, au = heatmap(f[1,1], U; axis = (; title = "U"), colorrange = (minimum(U),1))
    Colorbar(f[1,2], au; vertical=true)
    _, av = heatmap(f[1,3], V; axis = (; title = "V"))
    Colorbar(f[1,4], av; vertical=true)
    _, ap = heatmap(f[2,1], P; axis = (; title = "P"))
    Colorbar(f[2,2], ap; vertical=true)

    Ui = interpolate(U, BSpline(Linear()))
    Vi = interpolate(V, BSpline(Linear()))
    UV(x,y) = Point2(Ui[x,y], Vi[x,y])

    _, auv = streamplot(f[2,3], UV, 1:size(U,1), 1:size(V,2); axis = (; title = "(U, V)"))
    # mkdir("output/Re=$Re")
    save("output/Re=$Re/UVPstreamline_ts=$t_slice.png", f)
    return f
end

# for t in 1:160
#     time_slice_plot("output/n=128 , Re=1000.0 , t_end=16 , 2022-11-03T12.31.40.002/Ut.jld", t)
# end
# for t in [1,2,4,8,16,32]
#     time_slice_plot("output/n=128 , Re=400.0 , t_end=35 , 2022-11-02T19.11.14.253/Ut.jld", t)
# end
# ffig = center_velocity_compare("output/n=128 , Re=400.0 , t_end=35 , 2022-11-02T19.11.14.253/Ut.jld",
#     "StarCCM/Re400 U.csv", "StarCCM/Re400 V.csv", "Ghia/Re400 U.csv", "Ghia/Re400 V.csv", 100)
