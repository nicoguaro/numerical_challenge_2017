using Distributions


function monte_carlo_int(fun, N, low, high; args=())
    ndims = length(low)
    pts = rand(Uniform(0, 1), N, ndims)
    for cont = 1:ndims
        pts[:, cont] = pts[:, cont]*(high[cont] - low[cont]) + low[cont]
    end
    V = prod(high - low)
    return V*sum(fun(pts, args...))/N
end


function circ(x, rad)
    return 0.5*(1 - sign.(x[:, 1].^2 + x[:, 2].^2 - rad^2))
end


N = 1000000
low = [-1, -1]
high = [1, 1]
rad = 1
inte = monte_carlo_int(circ, N, low, high, args=(rad,))
