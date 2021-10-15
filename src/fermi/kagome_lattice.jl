"""
kagome lattice的相关功能
"""



"""
上面那个能带
"""
macro upperband(kx, ky)
    return esc(quote
        xval = $kx / 4
        yval = √3 * ($ky) / 4
        sqr = 2 * cos(2*xval - 2*yval)
        sqr += 2 * cos(2*xval + 2*yval)
        sqr += 2 * cos(4*xval) + 3
        if sqr < 0
            if isapprox(sqr, 0., atol=1e-12)
                sqr = 0
            else
                throw(error("能带有错误"))
            end
        end
        eng = - 1 + sqrt(sqr)
    end)
end


"""
上半个能带的kagome lattice
"""
function upperband_kagome_lattice(μ)
    disp(x, y) = (@upperband x y) + μ
    return TriangularSystem{:KAGOME}(
        [disp]
    )
end

