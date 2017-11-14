using PyPlot


function RK4(dydt, y0, t; args=())
    ndof = length(y0)
    ntimes = length(t)
    y = zeros(ndof, ntimes)
    y[:, 1] = y0
    for cont = 2:ntimes
        h = t[cont] - t[cont - 1]
        k1 = dydt(y[:, cont - 1], t[cont], args...)
        k2 = dydt(y[:, cont - 1] + 0.5*h*k1, t[cont] + 0.5*h, args...)
        k3 = dydt(y[:, cont - 1] + 0.5*h*k2, t[cont] + 0.5*h, args...)
        k4 = dydt(y[:, cont - 1] + h*k3, t[cont] + h, args...)
        y[:, cont] = y[:, cont - 1] + h/6*(k1 + 2*k2 + 2*k3 + k4)
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
y = RK4(pend, y0, t, args=(b, c))
plot(t, y[1, :])
plot(t, y[2, :])
xlabel(L"$t$")
legend([L"$\theta(t)$", L"$\omega(t)$"])
show()
