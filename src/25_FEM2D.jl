using PyPlot


function assem(coords, elems, source)
    ncoords = size(coords)[1]
    nelems = size(elems)[1]
    stiff = zeros(ncoords, ncoords)
    rhs = zeros(ncoords)
    for el_cont = 1:nelems
        elem = elems[el_cont, :]
        stiff_loc, jaco = local_stiff(coords[elem, :])
        rhs[elem] += jaco*mean(source[elem])
        for row = 1:3
            for col = 1:3
                row_glob = elem[row]
                col_glob = elem[col]
                stiff[row_glob, col_glob] += stiff_loc[row, col]
            end
        end
    end
    return stiff, rhs
end


function local_stiff(coords)
    dHdr = [-1 1 0; -1 0 1]
    jaco = dHdr * coords
    stiff = 0.5*[2 -1 -1; -1 1 0; -1 0 1]
    return stiff/det(jaco), det(jaco)
end


sq3 = sqrt(3)
coords =[sq3 -1;
        0 0;
        2*sq3 0;
        0 2;
        2*sq3 2;
        sq3 3;
        sq3 1]
elems =[2 1 7;
        1 3 7;
        3 5 7;
        5 6 7;
        6 4 7;
        4 2 7]
source = -ones(7)
stiff, rhs = assem(coords, elems, source)
free = 1:6
sol = stiff[free, free] \ rhs[free]
sol_c = zeros(size(coords)[1])
sol_c[free] = sol
tricontourf(coords[:, 1], coords[:, 2], sol_c, cmap="hot")
colorbar()
axis("image")
show()
