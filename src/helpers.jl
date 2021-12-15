#=
Best possible score of a PWM

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    the best possible score this matrix can get
=#
function best_score(pwm::Matrix{T}) where T <: Real
    sum(maximum(pwm[:,i]) for i = 1:size(pwm,2));
end

function best_score(pwm_col::Vector{T}) where T <: Real
    maximum(pwm_col);
end

#=
Worst possible score of a PWM

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    the worst possible score this matrix can get
=#
function worst_score(pwm::Matrix{T}) where T <: Real
    sum(minimum(pwm[:,i]) for i = 1:size(pwm,2));
end

function worst_score(pwm_col::Vector{T}) where T <: Real
    minimum(pwm_col);
end
#=
Return a column-permuted PWM that minimize the score range so that δ₁ ≥ δ₂ ≥ … ≥ δₘ
where δᵢ = best_score(pwm[:,i])-worst_score(pwm[:,i]). 

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
# Output
    a column-permuted pwm 
=#
function min_score_range(pwm::Matrix{T}) where T <: Real
    pwm[:,sortperm([best_score(pwm[:,i])-worst_score(pwm[:,i]) for i = 1:size(pwm,2)],rev=true)];
end

#=
"Round the PWM"

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
    `granularity`: a small positive real number e.g. 0.01 or 0.001, etc.
# Output
    a rounded pwm of the input pwm
=#
function round_pwm(pwm::Matrix{T}, granularity::Real) where T <: Real
    floor.(pwm ./ granularity) * granularity;
end

#=
The maximum error induced by the rounded pwm M_ϵ 
(see definition 3 in https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15; this is the quantity E) 

# Input
    `pwm::Matrix{Real}`: a 4 x m matrix
    `granularity`:
# Output
    A positve real number that's the maximum error induced by the rounded pwm M_ϵ 
=#
function calc_E(pwm, pwm_rounded) 
    sum(maximum(pwm[:,i]-pwm_rounded[:,i]) for i = 1:size(pwm,2));
end

#=
Note: Use a nested dictionary to represent the distribution Q
    Since Q is used to reference the probability of (M[1…i],score),
    the keys in the first layer is i, and the value in the first layer 
    are dictionaries with scores as keys and probability as values

    call create_Q(m) to initialize such a distribution Q
        where m is the "width" of the PWM
=#
create_Q(m::Integer) = Dict{Int16,SortedDict{Double64,Double64}}(i==0 ? i=>SortedDict(0=>1) : i=>SortedDict() for i=0:m);

#=
Input: 
    pwm: a 4 x m matrix
    α, β: score interval [α, β]
    bg: 4 x 1 vector that specifies the multinomial genomic background; default to flat background.    
Output:
    Q: a probability mass table
        e.g. Q[m] shows all the weights of P[pwm_score = η] for α ≤ η ≤ β
=#
function score_distribution(pwm_::Matrix{T}, α::Real, β::Real, bg=[.25,.25,.25,.25]) where T <: Real    
    m = size(pwm_,2);
    Q = create_Q(m);
    @inbounds for i = 1:m        
        bs = i+1 > m ? 0 : best_score(pwm_[:,i+1:m]);
        ws = i+1 > m ? 0 : worst_score(pwm_[:,i+1:m]);
        for score in keys(Q[i-1])
            for j = 1:4
                t = score + pwm_[j,i];
                if α - bs ≤ t ≤ β - ws
                    if haskey(Q[i], t)
                        Q[i][t] += Q[i-1][score]*bg[j];
                    else
                        Q[i][t] = Q[i-1][score]*bg[j];
                    end
                end
            end
        end
    end
    return Q
end

#=
Return the probability of the background model that has score ≥ α with respect to the input pwm

Input: 
    pwm: a 4 x m matrix
    α: a score threshold
    bg: 4 x 1 vector that specifies the multinomial genomic background; default to flat background.    
Output:
    pval: the probability that the background can have score ≥ α with respect to pwm
=#
function fast_pvalue(pwm::Matrix{T}, α::Real, bg=[.25,.25,.25,.25]) where T <: Real
    m = size(pwm,2);
    Q = create_Q(m);
    pval = 0f0;
    @inbounds for i = 1:m
        bs = i+1 > m ? 0 : best_score(pwm[:,i+1:m]);
        ws = i+1 > m ? 0 : worst_score(pwm[:,i+1:m]);
        for (score,_) in Q[i-1]
            for j = 1:4
                t = score + pwm[j,i];
                if α - ws ≤ t
                    pval = pval + Q[i-1][score]*bg[j];
                elseif α - bs ≤ t
                    if haskey(Q[i], t)
                        Q[i][t] += Q[i-1][score]*bg[j];
                    else
                        Q[i][t] = Q[i-1][score]*bg[j];
                    end
                end
            end
        end
    end
    return pval
end

# return the sum of all the weights 
Q_sum(Q_m::SortedDict{Double64,Double64}) = sum(v for v in values(Q_m));
