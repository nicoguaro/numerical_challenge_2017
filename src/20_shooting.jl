using PyPlot, DifferentialEquations, Roots


function shooting(dydx, x, y0, yf; shoot=nothing)
    if shoot == nothing
        shoot = rand(-20:0.1:20)
    end
    function F(s)
        prob = ODEProblem(dydx, [y0, s], (x[1], x[end]))
        sol = solve(prob)
        return sol(x[end])[1] - yf
    end
    shoot = fzero(F, shoot)
    prob = ODEProblem(dydx, [y0, shoot], (x[1], x[end]))
    sol = solve(prob, solveat=x)
    return sol(x)[1, :], shoot
end


func(x, y) = [y[2], 1.5*y[1]^2]
x = linspace(0, 1, 1000)
y, shoot = shooting(func, x , 4.0, 1.0, shoot=-5.0)
plot(x, y)
xlabel(L"$x$")
ylabel(L"$y$")
show()
