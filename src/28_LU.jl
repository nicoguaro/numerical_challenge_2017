function LU(mat)
    m, _ = size(mat)
    mat = copy(mat)
    for col = 1:m - 2
        for row = col + 1:m
            if mat[row, col] != 0.0
                lam = mat[row, col]/mat[col, col]
                mat[row, col + 1:m] = mat[row, col + 1:m] -
                                      lam * mat[col, col + 1:m]
                mat[row, col] = lam
            end
        end
    end
    return mat
end


A = [1.0 1.0 0.0 3.0;
    2.0 1.0 -1.0 1.0;
    3.0 -1.0 -1.0 2.0;
    -1.0 2.0 3.0 -1.0]
B = LU(A)
