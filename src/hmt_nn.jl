function calc_int_ps_fl(nm, nf, mat_v, ps_pot)
    int_el = zeros(no, no, no, nf)
    s = .5 * (nm - 1)
    for o1 in 1 : no 
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        m1r = m1 - s - 1
        for o2 in 1 : no 
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            m2r = m2 - s - 1
            for o3 in 1 : no 
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                m3r = m3 - s - 1
                for f4 in 1 : nf
                    m4 = m1 + m2 - m3 
                    if (m4 <= 0 || m4 > nm)
                        continue
                    end
                    m4r = m4 - s - 1
                    for i = 1 : length(mat_v)
                        for l in 1 : length(ps_pot[i])
                            if (abs(m1r + m2r) > nm - l || abs(m3r + m4r) > nm - l)
                                break
                            end 
                            int_el[o1, o2, o3, f4] += ps_pot[i][l] * (2 * nm - 2 * l + 1) * wigner3j(s, s, nm - l, m1r, m2r, -m1r - m2r) * wigner3j(s, s, nm - l, m4r, m3r, -m3r - m4r) * mat_v[i][f1, f4] * mat_v[i][f2, f3]
                        end
                    end
                end
            end
        end
    end
    return int_el
end

function generate_hmt_ops_fuzzy_nn(nm :: Int, nf :: Int, mat_v, ps_pot, mat_h ; m0 = 0)
    os = OpSum()
    int_el = calc_int_ps_fl(nm, nf, mat_v, ps_pot)
    no = nm * nf
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for f2 = 1 : nf 
            m2 = m1 
            o2 = (m2 - 1) * nf + f2
            if (abs(mat_h[f1, f2]) < 1E-8) continue end 
            val = mat_h[f1, f2]
            os += val, "Cdag", o1 + m0, "C", o2 + m0
        end
    end
    for o1 = 1 : no
        m1 = div(o1 - 1, nf) + 1
        f1 = mod(o1 - 1, nf) + 1
        for o2 = o1 + 1 : no # consider only o1 < o2, o3 > o4
            m2 = div(o2 - 1, nf) + 1
            f2 = mod(o2 - 1, nf) + 1
            for o3 = 1 : no 
                m3 = div(o3 - 1, nf) + 1
                f3 = mod(o3 - 1, nf) + 1
                m4 = m1 + m2 - m3 
                if (m4 <= 0 || m4 > nm) continue end
                for f4 = 1 : nf
                    o4 = (m4 - 1) * nf + f4
                    if (o3 <= o4) continue end
                    if (abs(int_el[o1, o2, o3, f4]) < 1E-8 && abs(int_el[o2, o1, o3, f4]) < 1E-8) continue end
                    val = 2. * (int_el[o1, o2, o3, f4] - int_el[o2, o1, o3, f4])
                    os += val, "Cdag", o1 + m0, "Cdag", o2 + m0, "C", o3 + m0, "C", o4 + m0
                end
            end
        end
    end
    return os
end
