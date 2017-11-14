function clenshaw_curtis(fun; N=50, a=-1, b=1)
    nodes = linspace(-1, 1, N + 1)*pi
    jac = 0.5*(b - a)
    tfun(x) = fun(a + jac*(x + 1))
    inte = 0
    for k = 0:2:N
        coef = 1/N*(tfun(1) + tfun(-1)*(-1)^k +
            2*sum(tfun(cos.(nodes[2:end-1])).*cos.(k*nodes[2:end-1])))
        if k == 0
            inte += coef
        elseif k == N
            inte += coef/(1 - N^2)
        else
            inte += 2*coef/(1 - k^2)
        end
    end
    return inte*jac
end


N = 100
fun(x) = exp.(-x.^2)
print(clenshaw_curtis(fun, N=N, a=0, b=3))

