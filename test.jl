using JLD

A = 1.0*[1,2,3,4,5,6,7,8]
B = 1.0*[1,2,3,4,5,6,7,8]
function adder(A)
    for i in 1:8
        for j in 1:100000000
            A[i] += sin(j)
        end
    end
end

function thread_adder(A)
    Threads.@threads for i in 1:8
        for j in 1:100000000
            A[i] += sin(j)
        end
    end
end

AA = load("test_A.jld")