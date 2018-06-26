using PyPlot


function jacobi(u, update; tol=1e-6, niter=500)
    num = niter
    for n = 1:niter
        u_new = update(u)
        if norm(u_new - u) < tol
            num = n
            break
        else
            u = copy(u_new)
        end
    end
    return u, num
end


function heat_FDM(u)
    u_new = copy(u)
    u_new[2:end-1, 2:end-1] = 0.25*(u_new[1:end-2, 2:end-1] +
        u_new[3:end, 2:end-1] + u_new[2:end-1, 1:end-2] + u_new[2:end-1, 3:end])
    return u_new
end

    
nx = 50
ny = 50
x_vec = linspace(-0.5, 0.5, nx)
y_vec = linspace(-0.5, 0.5, ny)
u0 = zeros(nx, ny)
u0[:, 1] = 1 - x_vec
u0[1, :] = 1 - y_vec
nvec = [100, 1000, 10000, 100000]
for (num, niter) = enumerate(nvec)
    u, n = jacobi(u0, heat_FDM, tol=1e-12, niter=niter)
    subplot(2, 2, num)
    contourf(u, cmap="hot")
    xlabel(L"$x$")
    ylabel(L"$y$")
    title("$(n) iterations")
    axis("image")
end
tight_layout()
show()
