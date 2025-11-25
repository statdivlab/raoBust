# Combine independent p-values that test the same null hypothesis into an overall summary

For a vector of p-values \\(p_1...p_n)\\ that independently test a
specific null hypothesis \$H_0\$, find the relevant quantile of a chi-sq
distribution for each p-value, then sum them up. That's a
\\\chi^2_n\\-distributed test statistic under \$H_0\$, and the relevant
p-value can be found. Useful for combining analyses after stratification
(to deal with potentially complex experimental designs).

## Usage

``` r
combine_independent_p_values(ps)
```

## Arguments

- ps:

  A vector of independent p-values (that test the same null hypothesis)

## Value

A single p-value that combines the evidence from the vector

## Author

Amy Willis
