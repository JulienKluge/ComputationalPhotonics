function simpleDFT(hk::Array{T}) where T
    N = length(hk)
    if N == 1
        return hk # * exp(-2iπn*0/N) = exp(0) = 1
    end
    ev = simpleDFT(hk[1:2:end]) #even
    od = simpleDFT(hk[2:2:end]) #odd

    twiddle = exp(-2 * π * im / N)
    w = 1 #twiddle^0

    Hk = Array{Complex{T}}(undef, N)
    NHalf = Int(N / 2)
    for i = 1:NHalf
        Hk[i] = ev[i] + w * od[i]
        Hk[i + NHalf] = ev[i] - w * od[i] #took me longer than i care to admit - had to look it up
        w = w * twiddle
    end

    return Hk
end
