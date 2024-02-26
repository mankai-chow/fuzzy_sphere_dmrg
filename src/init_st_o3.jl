function generate_init_st_o3(no, sites) # ; lz = 0, z2 = 1)
    conf = [(mod(o, 4) == 1 || mod(o, 4) == 2 ? "Occ" : "Emp") for o = 1 : no]
    # if (lz >= 0)
    #     conf[1] = "Emp"
    #     conf[2 + 2 * lz] = "Occ"
    # else
    #     conf[2] = "Occ"
    #     conf[1 - 2 * lz] = "Emp"
    # end
    # if (z2 == 1)
    #     conf[no    ] = "Occ"
    #     conf[no - 1] = "Emp"
    # end
    st0 = MPS(sites, conf)
    return st0 
end