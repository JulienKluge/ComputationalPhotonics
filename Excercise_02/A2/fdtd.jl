#fdtd.jl
using DelimitedFiles

#
# single 1d fdtd step in-place
#
function FDTDStep1D!(e::Array{T}, h::Array{T}, invϵ::Array{T}, Δx::T, Δt::T) where T
    ne = length(e)
    nh = length(h)
    dtdxFactor = Δt / Δx

    #e-field pass
    e[1] = e[end] = zero(T) #ensure boundary condition
    @inbounds for i = 2:(ne - 1) #spare first and last element
        e[i] += dtdxFactor * invϵ[i] * (h[i - 1] - h[i])
    end
    #h-field pass
    @inbounds for i = 1:nh
        h[i] += dtdxFactor * (e[i] - e[i + 1])
    end
end

#
# fdtd 1d simulation
#
function FDTD1D(steps::Int, snapshots::Array{Int}, n::Int, Δx::Float64, Δt::Float64, ϵBorder::Float64 = Inf)

    # spatial coordinates (not needed for calculation)
    x0 = (n - 1.0) / 2.0 * Δx # artificial midpoint
    x_e = collect(range(-x0 / 2.0         , length = n,     step = Δx))
    x_h = collect(range(-x0 / 2.0 + Δx / 2, length = n - 1, step = Δx))

    #initial fields: gaussian packet
    e = exp.(-0.5e2 * (n - 1) * Δx .* x_e .^ 2)
    h = exp.(-0.5e2 * (n - 2) * Δx .* x_h .^ 2)
    invϵ = map(x -> if (x > ϵBorder) 1.0 / 4.0 else 1.0 end, x_e)

    #initialize the result variables and snapshot-track-variables
    snapshotIndex = 1
    nextSnapshotIndex = snapshots[1]
    eResult = Array{Float64}(undef, length(snapshots), n)
    hResult = Array{Float64}(undef, length(snapshots), n - 1)

    #lets go
    for i = 1:steps
        FDTDStep1D!(e, h, invϵ, Δx, Δt) #magic

        if i == nextSnapshotIndex #if we reached a snapshot point:
            #copy the current fields
            for m = 1:n
                eResult[snapshotIndex, m] = e[m]
            end
            for m = 1:(n - 1)
                hResult[snapshotIndex, m] = h[m]
            end
            if snapshotIndex < length(snapshots) #and advance the snapshot-track-variables
                snapshotIndex += 1
                nextSnapshotIndex = snapshots[snapshotIndex]
            end
        end
    end

    #return out the results
    (efield = eResult, xe = x_e, hfield = hResult, xh = x_h, snaps = snapshots)
end


#
# save result in csv file
#
function SaveFDTDResult(resultFile, result)
    snapshots = result.snaps
    open(resultFile, "w") do io
        writedlm(io, [result.snaps, result.xe, result.xh], ',')
        writedlm(io, result.efield, ',')
        writedlm(io, result.hfield, ',')
    end
end