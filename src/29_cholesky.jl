function cholesky(mat)
    m, _ = size(mat)
    mat = copy(mat)
    for col = 1:m
        mat[col, col] = sqrt(mat[col, col] -
            dot(mat[col, 1:col-1], mat[col, 1:col-1]))
        for row = col + 1:m
            mat[row, col] = (mat[row, col] -
               dot(mat[row, 1:col-1], mat[col, 1:col-1]))/mat[col, col]
        end
    end
    for row = 2:m
        mat[1:row-1, row] = 0
    end
    return mat
end


A = [4 -1 1;
    -1 4.25 2.75;
    1 2.75 3.5]
B = cholesky(A)
