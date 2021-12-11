function find_largest_α(Q_m::SortedDict{T,T}, pval::T) where T <: Real
    q_sum = Q_sum(Q_m);
    largest_k = nothing;
    for (k,v) in Q_m
        if q_sum ≥ pval
            largest_k = k;
        else
            return k
        end
        q_sum -= v;
    end
    return largest_k
end

function pval_w_Qm(Qm::SortedDict{T,T}, α::Real) where T <: Real
    pval = 0;
    for (k,v) in Qm
        if k ≥ α
            pval += v;
        end
    end 
    return pval
end

function find_δ(Q_m::SortedDict{T,T}, pval_ϵ::Real, pval::Real) where T <: Real
    q_sum_plus_pval_ϵ = Q_sum(Q_m)+pval_ϵ;
    largest_δ = nothing;
    for (k,v) in Q_m
        if q_sum_plus_pval_ϵ ≥ pval
            largest_δ = k;
        else
            return k
        end
        q_sum_plus_pval_ϵ -= v;
    end
    return largest_δ
end


"""
    pval2score(pwm, pval, ϵ=1e-1, k=10, bg=[.25,.25,.25,.25])
Returns the highest score(M,pval) of a `pwm` such that p-value is greater or equal to `pval`.

Input:
* `pwm`: a 4 x m matrix
* `pval`: a p-value; e.g. pval = 1e-3
* `ϵ`: initial granularity  (optional) 
* `k`: Refinement parameter (optional)
* `bg`: multinomial background (optional)

Output    
* `α`: the highest score-threshold
"""
function pvalue2score(pwm::Matrix{T}, pval::Real, ϵ=1e-2, k=10, bg=[.25,.25,.25,.25]) where T <: Real
    @assert 0 ≤ pval ≤ 1 "pvalue must be in [0,1]"
    @assert size(pwm,1) == 4 "The input matrix must have 4 and only 4 rows"
    pwm_double64 = Double64.(pwm);
    pval_double64 = Double64(pval);
    bg_double64 = Double64.(bg);
    mpwm = min_score_range(pwm_double64);
    m = size(pwm, 2);
    pwm_ϵ = round_pwm(mpwm, ϵ);
    E = calc_E(mpwm, pwm_ϵ);
    Q = create_Q(m);
    Q = score_distribution(pwm_ϵ,worst_score(pwm_ϵ),Inf,bg_double64);
    α = find_largest_α(Q[m], pval_double64);

    @inbounds while !(pval_w_Qm(Q[m], α-E) == pval_w_Qm(Q[m], α))
        # println("err: ", pval_w_Qm(Q[m], α-E) - pval_w_Qm(Q[m], α));
        ϵ = ϵ/k; 
        pwm_ϵ = round_pwm(mpwm, ϵ);
        E = calc_E(mpwm, pwm_ϵ);
        Q = create_Q(m);
        Q = score_distribution(pwm_ϵ,α-E,α+E,bg_double64);
        #=
        note:
        Sometimes the score range is simply 
        too small and hence the score_distribution subroutine 
        cannot find any scores within the range [α-E, α+E].
        This happens when the error of the round matrix, E, is very 
        small (see definition 3 in Touzet and Varre's paper: 
        https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15)
        When this happens we return α.
        =#
        if isempty(Q[m])
            # α = α + E;
            # break;
            return α
        end
        pval_ϵ = fast_pvalue(pwm_ϵ,α+E);
        δ = find_δ(Q[m],pval_ϵ,pval);
        α = δ;
    end

    return α;
end