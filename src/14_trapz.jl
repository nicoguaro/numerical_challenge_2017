function trapz(y; x=nothing, dx=1.0)
    if x == nothing
        inte = 0.5*dx*(2*sum(y[2:end-1]) + y[1] + y[end])
    else
        dx = x[2:end] - x[1:end-1]
        inte = 0.5* (y[1:end-1] + y[2:end])' * dx
    end
    return inte
end


dx = 0.001
x = 0:dx:10
y = sin.(x)
println(trapz(y, dx=dx))
println(1 - cos(10))

