function calc_int_ps(nm, ps_pot :: Vector{Float64})
    ## N = 2s+1, get V[i,j,k,l]
    int_el = zeros(nm, nm, nm)
    s = .5 * (nm - 1)
    for m1 in 1 : nm 
        m1r = m1 - s - 1
        for m2 in 1 : nm 
            m2r = m2 - s - 1
            for m3 = 1 : nm
                m3r = m3 - s - 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm)
                    continue
                end
                m4r = m4 - s - 1
                for l in 1 : length(ps_pot)
                    if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l)
                        break
                    end 
                    int_el[m1, m2, m3] += ps_pot[l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r)
                end 
            end 
        end 
    end
    return int_el
end