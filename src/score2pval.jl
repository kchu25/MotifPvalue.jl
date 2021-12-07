#=
Helper function for score2pvalue. Find the lowest 
score s such that P-value(M_ϵ,s) = P-value(M_ϵ,s-E).

Input:
    Q_m: a sorted dict
    E: the max column-wise error between matrix M and M_ϵ
Output:
=#
function find_s(Q_m::SortedDict{Float64,Float64}, E::Real)
    keys_ = collect(keys(Q_m)); ℓ = length(keys_);
    if ℓ > 1
        for i = 1:(ℓ-1)
            @inbounds if keys_[i]+E < keys_[i+1]
                return keys_[i];
            end
        end
    end
    return nothing
end

"""
    score2pvalue(pwm, α, ϵ=1e-1, k=100, bg=[.25,.25,.25,.25])
Returns P-value(M,α) of a `pwm` with a given threshold `α`.

Input:
* `pwm`: a 4 x m matrix
* `α`: the score 
* `ϵ`: initial granularity  (optional) 
* `k`: Refinement parameter (optional)
* `bg`: multinomial background (optional)

Output:
* `pval`: p-value 
"""
function score2pvalue(pwm::Matrix{T}, α::Real, ϵ=1e-1, k=100, bg=[.25,.25,.25,.25]) where T <: Real
    @assert size(pwm,1) == 4 "The input matrix must have 4 and only 4 rows"
    mpwm = min_score_range(pwm);
    m = size(mpwm,2);    
    β = best_score(mpwm)+1;
    pval = 0; 
    s = 0; 
    i = 1;
    @inbounds while !(α ≈ s)
        # println(α-s)
        ϵ = i == 1 ? ϵ : ϵ/k; i+=1;
        # isinf(ϵ) || isnan(ϵ) && (println(i); break);
        pwm_ϵ = round_pwm(mpwm, ϵ);
        # any(isnan.(pwm_ϵ)) && (println(i); break);
        E = calc_E(mpwm, pwm_ϵ);
        Q = create_Q(m);
        Q = score_distribution(pwm_ϵ, α-E, β, bg);
        #=
        note: 
        Sometimes the score range is small 
        (i.e. β-(α-E) is very  small) and hence we won't be 
        able to find such s. In such case, we return the 
        calculated p-value. 

        Note that this p-value is an underestimate. To get a 
        more accurate result, set k to be a larger number, 
        e.g. k=100.
        =#
        s = find_s(Q[m],E);
        if !isnothing(s)
            for (k,v) in Q[m]
                pval += k ≥ s ? v : 0;
            end
            β = s;
        else
            break;
        end
    end    
    return pval
end