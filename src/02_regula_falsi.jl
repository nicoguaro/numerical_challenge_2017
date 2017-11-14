function regula_falsi(fun, a, b, niter=50, ftol=1e-12, verbose=false)
    if fun(a) * fun(b) > 0
        c = nothing
        msg = "The function should have a sign change in the interval."
    else
        for cont = 1:niter
            qa = fun(a)
            qb = fun(b)
            c = (a*qb - b*qa)/(qb - qa)
            qc = fun(c)
            if verbose
                println("n: $(cont), c: $(c)")
            end
            if abs(fun(c)) < ftol
                msg = "Root found with desired accuracy."
                break
            elseif qa * qc < 0
                b = c
            elseif qb * qc < 0
                a = c
            end
            msg = "Maximum number of iterations reached."
        end
    end
    return c, msg
end

function fun(x)
    return cos(x) - x
end

println(regula_falsi(fun, 0.5, 0.25*pi))

