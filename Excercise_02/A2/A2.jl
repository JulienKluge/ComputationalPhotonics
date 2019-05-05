using Plots
using FFTW
using LinearAlgebra

include("fdtd.jl")

function A2b()
    # Spatial discretization
    n = 2000 + 1
    Δx = 1.0e-3

    # Time discretization
    Δt = Δx

    #steps to take
    steps = 2000

    #result at first and last step
    snapshots = Int.(round.(collect(range(1, stop = steps, length = 20))))

    #simulate
    result = FDTD1D(steps, snapshots, n, Δx, Δt)

    #save
    SaveFDTDResult("A2b.csv", result)

    #plot result
    ps = []
    for s = 1:length(snapshots)
        p = plot(result.xe, result.efield[s,:], ylims = [-1.5, 1.5])
        plot!(p, result.xh, result.hfield[s,:])
        push!(ps, p)
    end
    plot(ps..., layout = (length(snapshots), 1), show = true, size = [1200, 600])
end

function A2c()
    # Spatial discretization
    n = 2000 + 1
    Δx = 1.0e-3

    # Time discretization
    Δt = 1.01 * Δx
    steps = 130
    snapshots = [130, 140]
    result = FDTD1D(steps, snapshots, n, Δx, Δt)
    SaveFDTDResult("A2c.csv", result)
    p = plot(result.xe, result.efield[1,:], ylims = [-1.5, 1.5])
    plot!(p, result.xh, result.hfield[1,:])
end

function A2d()
    n = 2000 + 1
    Δx = 1.0e-3
    Δt = Δx
    steps = 2000
    snapshots = collect(1:steps)
    result = FDTD1D(steps, snapshots, n, Δx, Δt, 1.0)
    e = result.efield[:, 1000]
    h = result.hfield[:, 1000]
    ef = fft(e .* ((-1.0).^collect(1:length(e))))
    hf = fft(h .* ((-1.0).^collect(1:length(h))))
    s = 0.5 .* real.(ef .* hf)
    plot(abs.(s[950:1050]))
end

function A2f()
    n = 2000 + 1
    Δx = 1.0e-3
    Δt = Δx
    steps = 1500
    snapshots = Int.(round.(collect(range(1, stop = steps, length = 20))))
    result = FDTD1D(steps, snapshots, n, Δx, Δt, 1.0)
    SaveFDTDResult("A2f.csv", result)
end

function main()
    #A2b()
    #A2c()
    A2d()
    #A2f()
end

main()
