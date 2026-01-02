# Spectral entropy of a time series

Computes spectral entropy from a univariate normalized spectral density,
estimated using an AR model.

## Usage

``` r
entropy(x)
```

## Arguments

- x:

  a univariate time series

## Value

A non-negative real value for the spectral entropy \\H_s(x_t)\\.

## Details

The *spectral entropy* equals the Shannon entropy of the spectral
density \\f_x(\lambda)\\ of a stationary process \\x_t\\: \$\$ H_s(x_t)
= - \int\_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda, \$\$
where the density is normalized such that \\\int\_{-\pi}^{\pi}
f_x(\lambda) d \lambda = 1\\. An estimate of \\f(\lambda)\\ can be
obtained using [`spec.ar`](https://rdrr.io/r/stats/spec.ar.html) with
the `burg` method.

## References

Jerry D. Gibson and Jaewoo Jung (2006). “The Interpretation of Spectral
Entropy Based Upon Rate Distortion Functions”. IEEE International
Symposium on Information Theory, pp. 277-281.

Goerg, G. M. (2013). “Forecastable Component Analysis”. Proceedings of
the 30th International Conference on Machine Learning (PMLR) 28 (2):
64-72, 2013. Available at
<https://proceedings.mlr.press/v28/goerg13.html>.

## See also

[`spec.ar`](https://rdrr.io/r/stats/spec.ar.html)

## Author

Rob J Hyndman

## Examples

``` r
entropy(rnorm(1000))
#> entropy 
#>       1 
entropy(lynx)
#>   entropy 
#> 0.7331515 
entropy(sin(1:20))
#>     entropy 
#> 0.003481715 
```
