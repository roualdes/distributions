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
T | ✅ | ✅ | ✅ | ✅ | ✅ |
uniform | ✅ | ✅ | ✅ | ✅ | ✅ |

The central Student's t distribution is available as `t_density(x, df)`,
`t_distribution(x, df)`, `t_log_density(x, df)`, `t_quantile(p, df)`, and
`t_random(N, df)`. Here `df` is a positive number of degrees of freedom,
including fractional values; `Infinity` gives the standard normal limit.
These functions use zero location and unit scale. The quantile accepts
probabilities in `[0, 1]`, returning infinities at the endpoints. Invalid
arguments return `NaN` (an array of `NaN`s for invalid `df` in `t_random`);
an invalid sample count throws `RangeError`.

## build

To build the package, first clone this repository.  Then run

```
npm i
npm run build
```

This will make a `dist` folder with `distributions.es.js` and
`distributions.cjs.js`.

## test

Tests are run with

```
npm run test
```

## internals

This package makes use of the following

* [regularized incomplete beta function](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function)
* [regularized lower incomplete gamma function](https://en.wikipedia.org/wiki/Incomplete_gamma_function#Regularized_gamma_functions_and_Poisson_random_variables)
* [Cornish-Fisher expansion](https://en.wikipedia.org/wiki/Cornish–Fisher_expansion#cite_note-Abramowitz935-3)
* the root finding algorithm named after [Richard Brent](https://en.wikipedia.org/wiki/Richard_P._Brent), [Brent's method](https://en.wikipedia.org/wiki/Brent%27s_method)
* ...

See also the file `references.bib`.
