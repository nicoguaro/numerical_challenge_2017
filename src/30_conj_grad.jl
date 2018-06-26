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


N = 1000
A = -diagm(2*ones(N)) + diagm(ones(N-1), -1) + diagm(ones(N-1), 1)
b = ones(N)
x0 = ones(N)
x, niter, accu = conj_grad(A, b, x0)
