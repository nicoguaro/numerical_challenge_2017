function newton(fun, grad, x, niter=50, ftol=1e-12, verbose=false)
    msg = "Maximum number of iterations reached."
    for cont = 1:niter
        if abs(grad(x)) < ftol
            x = nothing
            msg = "Derivative near to zero."
            break
        end
        if verbose
            println("n: $(cont), x: $(x)")
        end
        x = x - fun(x)/grad(x)
        if abs(fun(x)) < ftol
            msg = "Root found with desired accuracy."
            break
        end
    end
    return x, msg
end


function fun(x)
    return cos(x) - x
end


function grad(x)
    return -sin(x) - 1.0
end


println(newton(fun, grad, 1.0))

