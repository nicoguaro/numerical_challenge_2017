function bisection(fun, a, b, xtol=1e-6, ftol=1e-12)
    if fun(a) * fun(b) > 0
        c = nothing
        msg = "The function should have a sign change in the interval."
    else
        nmax = ceil(log2((b - a)/xtol))
        for cont = 1:nmax
            c = 0.5*(a + b)
            if abs(fun(c)) < ftol
                msg = "Root found with desired accuracy."
                break
            elseif fun(a) * fun(c) < 0
                b = c
            elseif fun(b) * fun(c) < 0
                a = c
            end
            msg = "Maximum number of iterations reached."
        end
    end
    return c, msg
end

function fun(x)
    return cos(x) - x^2
end

println(bisection(fun, 0.0, 1.0))

