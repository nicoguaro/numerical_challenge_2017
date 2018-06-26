function grad_descent(fun, grad, x; niter=50, gtol=1e-8, verbose=false)
    msg = "Maximum number of iterations reached."
    g_old = grad(x)
    gamma = 0.1
    for cont = 1:niter
        if verbose
            println("n: $(cont), x: $(x), g: $(g_old)")
        end
        dx = - gamma*g_old
        x = x + dx
        g = grad(x)
        dg = g - g_old
        g_old = g
        gamma = dx' * dg / (dg' * dg)
        if norm(g) < gtol
            msg = "Extremum found with desired accuracy."
            break
        end
    end
    return x, fun(x), msg
end


function rosen(x)
    return (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
end


function rosen_grad(x)
    return [-2*(1 - x[1]) - 400*x[1]*(x[2] - x[1]^2);
            200*(x[2] - x[1]^2)]
end


println(grad_descent(rosen, rosen_grad, [2.0, 1.0], verbose=true))

