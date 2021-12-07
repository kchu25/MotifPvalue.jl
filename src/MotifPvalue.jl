module MotifPvalue

#### Dependencies ##################
using DataStructures
####################################

#### Exported methods and types ####
export score2pvalue, pval2score
####################################

#### Load files ####################
include("helpers.jl")
include("score2pval.jl")
include("pval2score.jl")
####################################

end
