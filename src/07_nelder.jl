function nelder_mead_step(fun, verts; alpha=1, gamma=2, rho=0.5,
                     sigma=0.5)
    nverts, _ = size(verts)
    f = [fun(verts[row, :]) for row in 1:nverts]
    # 1. Order
    order = sortperm(f)
    verts = verts[order, :]
    f = f[order]
    # 2. Calculate xo, the centroid
    xo = mean(verts[1:end - 1, :], 1)
    # 3. Reflection
    xr = xo + alpha*(xo - verts[end, :]')
    fr = fun(xr)
    if f[1]<=fr && fr<f[end - 1]
        new_verts = [verts[1:end-1, :]; xr]
    # 4. Expansion
    elseif fr<f[1]
        xe = xo + gamma*(xr - xo)
        fe = fun(xe)
        if fe < fr
            new_verts = [verts[1:end-1, :]; xe]
        else
            new_verts = [verts[1:end-1, :]; xr]
        end
    # 5. Contraction
    else
        xc = xo + rho*(verts[end, :]' - xo)
        fc = fun(xc)
        if fc < f[end]
            new_verts = [verts[1:end-1, :]; xc]
    # 6. Shrink
        else
            new_verts = zeros(verts)
            new_verts[1, :] = verts[1, :]
            for k =  1:nverts
                new_verts[k, :] = sigma*(verts[k,:] - verts[1,:])
            end
        end
    end

    return new_verts
end


function nelder_mead(fun, x; niter=50, atol=1e-8, verbose=false)
    msg = "Maximum number of iterations reached."
    f_old = fun(mean(x, 1))
    for cont = 1:niter
        if verbose
            println("n: $(cont), x: $(mean(x, 1))")
        end
        x = nelder_mead_step(fun, x)
        df = fun(mean(x, 1)) - f_old
        f_old = fun(mean(x, 1))
        if norm(df) < atol
            msg = "Extremum found with desired accuracy."
            break
        end
    end
    return mean(x, 1), f_old, msg
end


function rosen(x)
    return (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2
end


x = [1 0;
    1 1;
    2 0]
println(nelder_mead(rosen, x, verbose=false))
