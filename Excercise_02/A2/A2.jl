using Plots

include("fdtd.jl")

function main()
    # Spatial discretization
    n = 2000 + 1
    Δx = 1.0e-3

    # Time discretization
    Δt = Δx

    # spatial coordinates (not needed for calculation)
    x0 = (n - 1.0) / 2.0 * Δx # artificial midpoint
    x_e = collect(range(-x0 / 2.0         , length = n,     step = Δx))
    x_h = collect(range(-x0 / 2.0 + Δx / 2, length = n - 1, step = Δx))

    #initial fields: gaussian packet
    e = exp.(-0.5e2 * (n - 1) * Δx .* x_e .^ 2)
    h = exp.(-0.5e2 * (n - 2) * Δx .* x_h .^ 2)
    invϵ = map(x -> if (x > 1.0) 1.0 / 4.0 else 1.0 end, x_e)

    #simulation
    plotmod = 5
    for i = 1:2000
        fdtd_step!(e, h, invϵ, Δx, Δt)
        if (i % plotmod) == 0
            p1 = plot(x_e, e, ylims = [-1.0, 1.0], size = [1200, 600]);
            plot!(p1, x_h, h, show = true)
        end
        #sleep(0.1)
    end

end

main()
