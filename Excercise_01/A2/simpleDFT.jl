function simpleDFT(hk::Union{Array{T}, SubArray{T}}) where T
    N = length(hk)
    if N == 2
        return [Complex(hk[1] + hk[2]), Complex(hk[1] - hk[2])]
    end
    ev = simpleDFT(@view hk[1:2:end]) #even
    od = simpleDFT(@view hk[2:2:end]) #odd

    twiddle = exp(-2 * π * im / N)
    w = one(Complex{T}) #twiddle^0

    Hk = Array{Complex{T}}(undef, N)
    NHalf = Int(N / 2)
    @inbounds @simd for i = 1:NHalf
        Hk[i] = ev[i] + w * od[i]
        Hk[i + NHalf] = ev[i] - w * od[i] #took me longer than i care to admit - had to look it up
        w = w * twiddle
    end

    return Hk
end
