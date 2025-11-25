# Robust score (Rao) tests with finite-sample correction

The default behavior is to do no finite sample correction (for the
covariance of the score) for correlated data, but to do it for
uncorrelated data. This choice performed best under a small simulation
study. See Guo et al for details; the proposed modification is in
equation 20.

## Usage

``` r
robust_score_test(glm_object, call_to_model, param = 1, id = NA)
```

## Arguments

- glm_object:

  The fitted glm under the alternative.

- call_to_model:

  The call used to fit the model. Used internally.

- param:

  which parameter do you want to test? Used internally.

- id:

  observations with the same id are in the same cluster

## References

Guo, X., Pan, W., Connett, J. E., Hannan, P. J., & French, S. A. (2005).
Small-sample performance of the robust score test and its modifications
in generalized estimating equations. *Statistics in Medicine, 24*(22),
3479–3495. Wiley Online Library. <doi:10.1002/sim.2161>
