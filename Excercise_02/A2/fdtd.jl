#fdtd.jl

function fdtd_step!(e::Array{T}, h::Array{T}, invϵ::Array{T}, Δx::T, Δt::T) where T
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
