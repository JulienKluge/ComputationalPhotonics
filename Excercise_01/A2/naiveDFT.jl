function naiveDFT(hk::Array{T}) where T
    N = length(hk);
    Wnk = Array{Complex{T}}(undef, N, N)

    prefactor = -2 * π * im / N

    for n = 0:(N - 1), k = 0:(N - 1)
        Wnk[k + 1, n + 1] = exp(prefactor * n * k)
    end

    return Wnk * hk
end

function naiveDFTIterative(hk::Array{T}) where T
    N = length(hk);
    Hk = zeros(Complex{T}, N)

    prefactor = -2 * π * im / N

    for k = 1:N
        s = zero(Complex{T})
        km1 = k - 1
        for j = 1:N
            s += exp(prefactor * (j - 1) * km1) * hk[j]
        end
        Hk[k] = s
    end

    return Hk
end
