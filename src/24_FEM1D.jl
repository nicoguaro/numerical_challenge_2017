using PyPlot


function FEM1D(coords, source)
    N = length(coords)
    stiff_loc = [2.0 -2.0; -2.0 2.0]
    eles = [[cont, cont + 1] for cont in 1:N-1]
    stiff = zeros(N, N)
    rhs = zeros(N)
    for ele in eles
        jaco = coords[ele[2]] - coords[ele[1]]
        rhs[ele] = rhs[ele] + jaco*source(coords[ele])
        stiff[ele, ele] = stiff[ele, ele] +  stiff_loc/jaco
    end
    return stiff, rhs
end


N = 100
fun(x) = x.^3
x = linspace(0, 1, N)
stiff, rhs = FEM1D(x, fun)
sol = zeros(N)
sol[2:end-1] = -stiff[2:end-1, 2:end-1]\rhs[2:end-1]


#%% Plotting
figure(figsize=(4, 3))
plot(x, sol)
plot(x, x.*(x.^4 - 1)/20, linestyle="dashed")
xlabel(L"$x$")
ylabel(L"$y$")
legend(["FEM solution", "Exact solution"])
tight_layout()
show()



