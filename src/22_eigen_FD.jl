using PyPlot


function beam_FDM_eigs(L, N)
    x = linspace(0, L, N)
    dx = x[2] - x[1]
    stiff = diagm(6*ones(N - 2)) -
            diagm(4*ones(N - 3), -1) - diagm(4*ones(N - 3), 1) +
            diagm(1*ones(N - 4), 2) + diagm(1*ones(N - 4), -2)
    vals, vecs = eig(stiff/dx^4)     
    return vals, vecs, x
end


N = 1001
nvals = 20
nvecs = 4
vals, vecs, x = beam_FDM_eigs(1.0, N)

#%% Plotting
num = 1:nvals
# Missing line for setting the math font
figure(figsize=(8, 3))
subplot(1, 2, 1)
plot(num, sqrt.(vals[1:nvals]), "o")
xlabel(L"$N$")
ylabel(L"$\omega\sqrt{\frac{\lambda}{EI}}$")
subplot(1, 2 ,2)
for k in 1:nvecs
    vec = zeros(N)
    vec[2:end-1] = vecs[:, k]
    plot(x, vec, label="n=$(k)")
end

xlabel(L"$x$")
ylabel(L"$w$")
legend(ncol=2, framealpha=0.8, loc=1)
tight_layout()
show()
