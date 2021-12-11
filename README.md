# MotifPvalue

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/MotifPvalue.jl/dev)
[![Build Status](https://github.com/kchu25/MotifPvalue.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/MotifPvalue.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/MotifPvalue.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/MotifPvalue.jl)


# Introduction
This package provides two useful subroutines:
*  Given a score threshold `α`, `score2pvalue` calculates the p-value, i.e. the probability of a set of strings generated from a specfied background model such that each string in this set can attain a score higher than or equal to `α` with respect to the input `pwm`.
* Given a p-value `pval`, `pval2score` calculates the highest score `α` such that `score2pvalue(pwm,α)` is larger than or equal to `pval`.

By default, the background model is specified as i.i.d discrete uniform.

This is an implementation of:

Efficient and accurate P-value computation for Position Weight Matrices by Touzet et al.
https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15

# Basic examples

    using MotifPvalue

    # An example PWM
    pwm = [-2.86995   1.3814    1.36906  -4.88479   -4.19164  -4.88479   1.36607  -5.57793   -0.336375    1.38238;
            1.36707  -5.57793  -2.86995  -0.294869  -4.88479   1.34492  -3.49849   1.26497   -2.1768     -4.88479;
           -5.57793  -4.47899  -5.57793   1.11923   -3.78584  -5.57793  -4.88479  -5.57793    0.798703   -5.57793;
           -4.19164  -5.57793  -4.88479  -1.68576    1.375    -1.88885  -3.17976  -0.798616  -0.0526305  -5.57793];
    
    # Compute pvalue with pwm and score = 4
    score2pvalue(pwm, 4)
    > 0.00021457672119140625

    # Compute the score-threshold for p-value 1e-4
    pvalue2score(pwm, 1e-4)
    > 5.57