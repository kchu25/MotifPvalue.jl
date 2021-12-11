module MotifPvalue

#### Dependencies ##################
using DataStructures, DoubleFloats
####################################

#### Exported methods and types ####
export score2pvalue, pvalue2score
####################################

#### Load files ####################
include("helpers.jl")
include("score2pval.jl")
include("pval2score.jl")
####################################

end
