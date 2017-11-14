using PyPlot

function lagrange(x_int, y_int, x_new)
    y_new = zeros(x_new)
    for (xi, yi) in zip(x_int, y_int)
        prod = ones(x_new)
        for xj in x_int
            if xi != xj
                prod = prod.* (x_new - xj)/(xi - xj)
            end
        end
        y_new += yi*prod
    end
    return y_new
end


x_int = linspace(-5, 5, 11)
y_int = 1./(1 + exp.(-x_int))
x_new = linspace(-5, 5, 1000)
y_new = lagrange(x_int, y_int, x_new)
plot(x_int, y_int, "ok")
plot(x_new, y_new)

