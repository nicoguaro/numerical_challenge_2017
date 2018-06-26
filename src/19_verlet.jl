using PyPlot


function verlet(force, x0, v0, m, t; args=())
    ndof = length(x0)
    ntimes = length(t)
    x = zeros(ndof, ntimes)
    dt = t[2] - t[1]
    x[:, 1] = x0
    F = force(x0, v0, m, t, args...)
    x[1:2:end, 2] = x0[1:2:end] + v0[1:2:end]*dt + 0.5*F[1:2:end]*dt^2./m
    x[2:2:end, 2] = x0[2:2:end] + v0[2:2:end]*dt + 0.5*F[2:2:end]*dt^2./m
    for cont = 3:ntimes
        dt = t[cont] - t[cont - 1]
        v = (x[:, cont - 1] - x[:, cont - 2])/dt
        acel = force(x[:, cont - 1], v, m, t, args...)
        acel[1:2:end] = acel[1:2:end]./m
        acel[2:2:end] = acel[2:2:end]./m
        x[:, cont] = 2*x[:, cont - 1] - x[:, cont - 2] + acel*dt^2
    end
    return x
end


spring_force(x, v, m, t, k) = -k*x
x0 = [1.0, 0.0]
v0 = [0.0, 1.0]
m = [1.0]
t = linspace(0, 10.0, 1000)
k = 1.0
x = verlet(spring_force, x0, v0, m, t, args=(k,))

#%% Plot
figure(figsize=(6, 3))
subplot(121)
plot(t, x[1, :])
plot(t, x[2, :])
subplot(122)
plot(x[1, :], x[2, :])
show()
