using PyPlot


function euler(dydt, y0, t; args=())
    ndof = length(y0)
    ntimes = length(t)
    y = zeros(ndof, ntimes)
    y[:, 1] = y0
    for cont = 2:ntimes
        h = t[cont] - t[cont - 1]
        y[:, cont] = y[:, cont - 1] + h*dydt(y[:, cont - 1], t[cont], args...)
    end
    return y
end


function pend(y, t, b, c)
    theta, omega = y
    dydt = [omega, -b*omega - c*sin(theta)]
    return dydt
end


b = 0.25
c = 5.0
y0 = [pi - 0.1, 0.0]
t = linspace(0, 10, 1001)
y = euler(pend, y0, t, args=(b, c))
plot(t, y[1, :])
plot(t, y[2, :])
xlabel(L"$t$")
legend([L"$\theta(t)$", L"$\omega(t)$"])
show()
