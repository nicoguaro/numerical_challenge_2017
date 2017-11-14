using PyPlot


function ritz(N, source)
    stiff_mat = zeros(N, N)
    rhs = zeros(N)
    for row in 0:N-1
        for col in 0:N-1
            numer = (2 + 2*row + 2*col + 2*row*col)
            denom = (row + col + 1) * (row + col + 2) * (row + col + 3)
            stiff_mat[row + 1, col + 1] = numer/denom
        end
        fun(x) = x^(row + 1)*(1 - x)*source(x)
        rhs[row + 1], _  = quadgk(fun, 0, 1)
    end
    return stiff_mat, rhs
end


N = 2
source(x) = x^3
mat, rhs = ritz(N, source)
c = -mat\rhs
x = linspace(0, 1, 100)
y = zeros(x)
for cont in 0:N - 1
    y += c[cont + 1]*x.^(cont + 1).*(1 - x)
end

#%% Plotting
figure(figsize=(4, 3))
plot(x, y)
plot(x, x.*(x.^4 - 1)/20, linestyle="dashed")
xlabel(L"$x$")
ylabel(L"$y$")
legend(["Ritz solution", "Exact solution"])
tight_layout()
show()
