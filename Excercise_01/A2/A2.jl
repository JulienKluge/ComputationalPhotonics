using FFTW
using Plots

include("naiveDFT.jl")
include("simpleDFT.jl")
include("benchmarkFunction.jl")

function main()
    #
    # Create test set
    #
    N = 2^8 # 2^8=256
    tmax = 30.0
    Δ = tmax / N

    t = collect(range(0.0, tmax, length = N))
    ω = collect(range(0, stop = N - 1))
    hk = sin.(t) .+ 0.5 .* sin.(5.0 .* t .+ 1.0) .+ sin.(4.0 .* t .- 1.5)


    #
    # Calculations
    #
    resNaive = naiveDFT(hk)                 #naive approach
    resNaiveIter = naiveDFTIterative(hk)    #naive approach iterative
    resSimple = simpleDFT(hk)               #simple recursive dft
    resFFTW = fft(hk)                       #fftw
    resRealFFTW = rfft(hk)                  #real fftw

    
    println(sum(abs.(resFFTW .- resNaive))) #differences to naive<->fftw
    println(sum(abs.(resFFTW .- resNaiveIter))) #differences to naive<->fftw
    println(sum(abs.(resFFTW .- resSimple))) #differences to simple<->fftw

    
    #
    # Benchmarking
    #
    benchmarkDFTFunction(fft, 1.0, 5, "FFTW.csv")
    benchmarkDFTFunction(simpleDFT, 1.0, 5, "simpleDFT.csv")
    benchmarkDFTFunction(naiveDFT, 1.0, 5, "naiveDFT.csv")
    benchmarkDFTFunction(naiveDFTIterative, 1.0, 5, "naiveDFTIter.csv")

    #
    # plotting
    #
    gr()
    p1re = plot(ω, real.(resFFTW), title = "FFTW Real", label = "");
    p1im = plot(ω, imag.(resFFTW), title = "FFTW Imag", label = "");
    p1abs = plot(ω, abs.(resFFTW) .^ 2 , title = "FFTW Abs^2", label = "");

    p2re = plot(ω, real.(resNaive), title = "Naive Real", label = "");
    p2im = plot(ω, imag.(resNaive), title = "Naive Imag", label = "");
    p2abs = plot(ω, abs.(resNaive) .^ 2 , title = "Naive Abs^2", label = "");

    p3re = plot(ω, real.(resSimple), title = "Simple Real", label = "");
    p3im = plot(ω, imag.(resSimple), title = "Simple Imag", label = "");
    p3abs = plot(ω, abs.(resSimple) .^ 2 , title = "Simple Abs^2", label = "");

    p = plot(p1re, p1im, p1abs, p2re, p2im, p2abs, p3re, p3im, p3abs, layout = (3, 3), size = (1024, 800))
    savefig(p, "fftPlot.png")

end

main()
