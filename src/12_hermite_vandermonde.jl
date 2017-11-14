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


function conf_vander_mat(x)
    n = length(x)
    conf_van = zeros(2*n, 2*n)
    power = 0:2*n-1
    for row = 1:n
        conf_van[row, :] = x[row].^power
        conf_van[row + n, :] = power.*x[row].^(power - 1)
    end
    return conf_van
end


function inter_coef(x; inter_type="lagrange")
    if inter_type == "lagrange"
        vand_mat = vander_mat(x)
    elseif inter_type == "hermite"
        vand_mat = conf_vander_mat(x)
    end
    coef = vand_mat \ eye(size(vand_mat)[1])
    return coef
end


function compute_interp(x, f, x_eval; df=nothing)
    n = length(x)
    if df == nothing
        coef = inter_coef(x, inter_type="lagrange")
    else
        coef = inter_coef(x, inter_type="hermite")
    end
    f_eval = zeros(x_eval)
    nmat = size(coef)[1]
    for row = 1:nmat
        for col = 1:nmat
            if col <= n || nmat == n
                f_eval += coef[row, col]*x_eval.^(row - 1)*f[col]
            else
                f_eval += coef[row, col]*x_eval.^(row - 1)*df[col - n]
            end
        end
    end
    return f_eval
end


n = 7
x = -cos.(linspace(0, pi, n))
f(x) = 1./(1 + 25*x.^2)
df(x) = -50*x./(1 + 25*x.^2).^2
x_eval = linspace(-1, 1, 500)
interp_f = compute_interp(x, f(x), x_eval, df=df(x))
plot(x_eval, f(x_eval))
plot(x_eval, interp_f)
plot(x, f(x), ".")
ylim(0, 1.2)
show()
