using PyPlot


function spline_coeff(x, y)
    h = x[2:end] - x[1:end-1]
#    println(h)
    d_up = copy(h)
    d_up[1] = 0.0
    d_down = copy(h)
    d_down[end-1] = 0.0
    d_cent = ones(x)
    d_cent[2:end-1] = 2*(h[1:end-1] + h[2:end])
    mat = diagm(d_cent) + diagm(d_up, 1) + diagm(d_down, -1)
    alpha = zeros(x)
    alpha[2:end-1] = 3./h[2:end].*(y[3:end] - y[2:end-1]) - 3./h[1:end-1].*(y[2:end-1] - y[1:end-2])
    c = mat \ alpha
    b = zeros(x)
    d = zeros(x)
    b[1:end-1] = (y[2:end] - y[1:end-1])./h - h./3.*(c[2:end] + 2*c[1:end-1])
    d[1:end-1] = (c[2:end] - c[1:end-1])./(3*h)
    return y, b, c, d
end


n = 20
x = collect(linspace(0, 2*pi, n))
y = sin.(x)
a, b, c, d = spline_coeff(x, y)
for cont = 1:n - 1
    x_loc = linspace(x[cont], x[cont + 1], 100)
    y_loc = a[cont] + b[cont]*(x_loc - x[cont]) +
            c[cont]*(x_loc - x[cont]).^2 +
            d[cont]*(x_loc - x[cont]).^3
    plot(x_loc, y_loc, "red")
end
plot(x, y, "o")
show()
