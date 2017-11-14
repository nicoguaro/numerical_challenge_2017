function newton_opt(fun, grad, hess, x; niter=50, gtol=1e-8, verbose=false)
    msg = "Maximum number of iterations reached."
    g = grad(x)
    for cont = 1:niter
        if verbose
            println("n: $(cont), x: $(x), g: $(g)")
        end
        x = x - hess(x)\g
        g = grad(x)
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


function rosen_hess(x)
    return [-400*(x[2] - x[1]^2) + 800*x[1]^2 + 2 -400*x[1];
            -400*x[1] 200]
end



println(newton_opt(rosen, rosen_grad, rosen_hess, [2.0, 1.0], verbose=true))

