# Overview

`distributions` is a Javascript library to compute some functions
associated with statistical distributions.

## available functions

name | density | distributon | log-density | random | quantile
-----|---------|-------------|-------------|--------|---------
beta | ✅ | ✅ | ✅ | ✅ | ✅ |
binomial | ✅ | ✅ | ✅ | ✅ | ✅ |
exponential | ✅ | ✅ | ✅ | ✅ | ✅ |
F | ❌ | ❌ | ❌ | ❌| ❌ |
gamma | ✅ | ✅ | ✅ | ✅ | ✅ |
geometric | ❌ | ❌ | ❌ | ❌| ❌ |
negative binomial | ❌ | ❌ | ❌ | ❌| ❌ |
normal | ✅ | ✅ | ✅ |  ✅ | ✅ |
poisson | ✅ | ✅ | ✅ | ✅ | ✅ |
rayleigh | ❌ | ❌ | ❌ | ❌| ❌ |
T | ❌ | ❌ | ❌ | ❌| ❌ |
uniform | ✅ | ✅ | ✅ | ✅ | ✅ |

## internals

This package makes use of the following

* [regularized incomplete beta function](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function)
* [regularized lower incomplete gamma function](https://en.wikipedia.org/wiki/Incomplete_gamma_function#Regularized_gamma_functions_and_Poisson_random_variables)
* [Cornish-Fisher expansion](https://en.wikipedia.org/wiki/Cornish–Fisher_expansion#cite_note-Abramowitz935-3)
* the root finding algorithm named after [Richard Brent](https://en.wikipedia.org/wiki/Richard_P._Brent), [Brent's method](https://en.wikipedia.org/wiki/Brent%27s_method)
* ...

See also the file `references.bib`.
