# Jackknife standard errors

Jackknife standard errors

## Usage

``` r
jackknife_se(object, dat, id = NULL)
```

## Arguments

- object:

  The fitted object under the alternative.

- dat:

  The data used to fit the model

- id:

  Observations with the same id are in the same cluster. If not
  included, independence between observations is assumed.
