function newton(fun, jaco, x, niter=50, ftol=1e-12, verbose=false)
    msg = "Maximum number of iterations reached."
    for cont = 1:niter
        J = jaco(x)
        f = fun(x)
        if det(J) < ftol
            x = nothing
            msg = "Derivative near to zero."
            break
        end
        if verbose
            println("n: $(cont), x: $(x)")
        end
        x = x - J\f
        if norm(f) < ftol
            msg = "Root found with desired accuracy."
            break
        end
    end
    return x, msg
end


function fun(x)
    return [x[1] + 2*x[2] - 2, x[1]^2 + 4*x[2]^2 - 4]
end


function jaco(x)
    return [1.0 2.0;
            2*x[1] 8*x[2]]
end


println(newton(fun, jaco, [1.0, 10.0]))

