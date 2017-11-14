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


function gauss_lobatto(N; tol=1e-15)
    x = -cos.(linspace(0, pi, N))
    P = zeros(N, N)  # Vandermonde Matrix
    x_old = 2
    while maximum(abs.(x - x_old)) > tol
        x_old = x
        P[:, 1] = 1
        P[:, 2] = x
        for k = 3:N
            P[:, k] = ((2 * k - 1) * x .* P[:, k - 1] -
                       (k - 1) * P[:, k - 2]) / k
        end
        x = x_old - (x .* P[:, N] - P[:, N - 1]) ./ (N* P[:, N])
    end
    return x
end


runge(x) =  1./(1 + 25*x.^2)
x = linspace(-1, 1, 100)
x_int = linspace(-1, 1, 11)
x_int2 = gauss_lobatto(11)
x_new = linspace(-1, 1, 1000)
y_new = lagrange(x_int, runge(x_int), x_new)
y_new2 = lagrange(x_int2, runge(x_int2), x_new)
plot(x, runge(x), "k")
plot(x_new, y_new)
plot(x_new, y_new2)
legend(["Runge function",
            "Uniform interpolation",
            "Lobatto-sampling interpolation"])
xlabel("x")
ylabel("y")
