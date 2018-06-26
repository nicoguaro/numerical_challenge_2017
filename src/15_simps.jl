function simps(y; x=nothing, dx=1.0)
    n = length(y)
    if x == nothing
        inte = sum(y[1:2:n-2] + 4*y[2:2:n-1] + y[3:2:n])*dx
    else
        dx = x[2:2:end] - x[1:2:end-1]
        inte = (y[1:2:n-2] + 4*y[2:2:n-1] + y[3:2:n])' * dx
    end
    return inte/3
end


n = 21
x = linspace(0, 10, n)
y = sin.(x)
println(simps(y, x=x))
println(1 - cos(10))

