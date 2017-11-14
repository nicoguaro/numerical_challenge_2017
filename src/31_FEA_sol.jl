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


function conj_grad(A, b, x; tol=1e-8)
    r = b - A * x
    p = r
    rsq_old = dot(r, r)
    niter = 1
    for cont = 1:length(b)
        Ap = A * p
        alpha = rsq_old / dot(p, Ap)
        x = x + alpha*p
        r = r - alpha*Ap
        rsq_new = dot(r, r)
        if sqrt(rsq_new) < tol
            break
        end
        p = r + (rsq_new / rsq_old) * p
        rsq_old = rsq_new
        niter += 1
    end
    return x, niter, norm(r)
end



fun(x) = x.^3
N_vec = round.(logspace(0.5, 3, 50))
err_vec = zeros(N_vec)
niter_vec = zeros(N_vec)
figure(figsize=(8, 3))
for (cont, N) in enumerate(N_vec)
    x = linspace(0.0, 1.0,N)
    stiff, rhs = FEM1D(x, fun)
    sol = zeros(N)
    sol[2:end-1], niter, _ = conj_grad(stiff[2:end-1, 2:end-1],
                                -rhs[2:end-1], rhs[2:end-1])
    err = norm(sol - x.*(x.^4 - 1)/20)/norm(x.*(x.^4 - 1)/20)
    err_vec[cont] = err
    niter_vec[cont] = niter
end
subplot(121)
loglog(N_vec, err_vec)
xlabel("Number of nodes")
ylabel("Relative error")
subplot(122)
loglog(N_vec, niter_vec)
xlabel("Number of nodes")
ylabel("Number of iterations")
tight_layout()
show()



