using PyPlot


function vander_mat(x)
    n = length(x)
    van = zeros(n, n)
    power = 0:n-1
    for row = 1:n
        van[row, :] = x[row].^power
    end
    return van
end


function inter_coef(x)
    vand_mat = vander_mat(x)
    coef = vand_mat \ eye(length(x))
    return coef
end


function compute_interp(x, f, x_eval)
    n = length(x)
    coef = inter_coef(x)
    f_eval = zeros(x_eval)
    for row = 1:n
        for col = 1:n
            f_eval += coef[row, col]*x_eval.^(row - 1)*f[col]
        end
    end
    return f_eval
end


n = 11
x = - cos.(linspace(0, pi, n))
f = 1./(1 + 25*x.^2)
x_eval = linspace(-1, 1, 500)
interp_f = compute_interp(x, f, x_eval)
plot(x_eval, 1./(1 + 25*x_eval.^2))
plot(x_eval, interp_f)
plot(x, f, ".")
ylim(0, 1.2)
show()
