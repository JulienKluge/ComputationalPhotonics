using Random
using DelimitedFiles

function benchmarkDFTFunction(f::Function, maxTime, repititions::Int, resultFile)
    n = 2
    times = Array{Float64, 1}(undef, 0)
    sizes = Array{Int, 1}(undef, 0)

    Random.seed!(30011995) #my birthday...
    duration = 0.0
    while duration <= maxTime
        n *= 2
        inputArr = rand(n)

        startTime = time_ns()
        for i = 1:repititions
            f(inputArr)
        end
        endTime = time_ns()

        duration = (endTime - startTime) / repititions * 1.0e-9
        push!(sizes, n)
        push!(times, duration)
    end

    open(resultFile, "w") do io
        writedlm(io, [sizes, times], ',')
    end
end
