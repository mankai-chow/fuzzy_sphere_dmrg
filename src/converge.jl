mutable struct EntanglementObserver <: AbstractObserver
    nb :: Int
    e_tol :: Float64
    e_last :: Float64
    EntanglementObserver(nb, e_tol = 0.0) = new(nb, e_tol, 1000.0)
end

function ITensors.measure!(o::EntanglementObserver; bond, psi, half_sweep, kwargs...)
    if (bond != o.nb)
        return
    end
    wf_center, other = half_sweep==1 ? (psi[bond+1],psi[bond]) : (psi[bond],psi[bond+1])
    U,S,V = svd(wf_center, uniqueinds(wf_center,other))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    println("  Entanglement across bond $bond = $SvN")
end

function ITensors.checkdone!(o :: EntanglementObserver; kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    if abs(energy - o.e_last) < o.e_tol
        println("Stopping DMRG after sweep $sw")
        return true
    end
    # Otherwise, update last_energy and keep going
    o.e_last = energy
    return false
end