using Distributions

function test1()
    b = Beta(0.9, 0.5)
    for i in 1:9000
        cdf(b, 0.2222222222222222)
    end
end

function test2()
    b = Beta(0.9, 0.5)
    c = Beta(2, 7)
    for i in 1:9000
        cdf(b, mean(c))
    end
end

function test3()
    b = Beta(0.9, 0.5)
    c = Beta(2, 7)
    for i in 1:9000
        cdf(b, median(c))
    end
end


function test4()
    b = Beta(0.9, 0.5)
    c = Beta(2, 7)
    for i in 1:9000
        sum(rand(b, 9000) .> rand(c, 9000)) / 9000
    end
end