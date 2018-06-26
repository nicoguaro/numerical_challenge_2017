function broyden(fun, jaco, x, niter=50, ftol=1e-12, verbose=false)
    msg = "Maximum number of iterations reached."
    J = jaco(x)
    for cont = 1:niter
        if det(J) < ftol
            x = nothing
            msg = "Derivative near to zero."
            break
        end
        if verbose
            println("n: $(cont), x: $(x)")
        end
        f_old = fun(x)
        dx = -J\f_old
        x = x + dx
        f = fun(x)
        df = f - f_old
        J = J + (df - J*dx) * dx'/ (dx' * dx)
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
    return [1 2;
           2*x[1] 8*x[2]]
end


println(broyden(fun, jaco, [1.0, 2.0]))

