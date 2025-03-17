import { gamma } from "./gamma_fn";
import { LOG_SQRT_2PI } from "./constants";
import { brent } from "./brent";
import { normal_quantile } from "./normal";

/**
 * Poisson(lambda) distribution's density function
 *
 * @param x argument to density function
 * @param lambda rate parameter
 * @returns poisson density at x
 *
 * @category poisson
 */
export function poisson_density(x: number, lambda: number): number {
  return Math.exp(poisson_log_density(x, lambda));
}

/**
 * Poisson(lambda) distribution function
 *
 * @param x argument to distribution function
 * @param lambda rate parameter
 * @returns poisson distribution at x
 *
 * @category poisson
 */
export function poisson_distribution(x: number, lambda: number): number {
  const k = Math.floor(x + 1);
  if (k <= 0) {
    return 0;
  }
  return 1 - gamma(lambda, k);
}

/**
 * Poisson(lambda) distribution's log-density function
 *
 * @param x argument to log-density function
 * @param lambda rate parameter
 * @returns poisson log-density at x
 *
 * @category poisson
 */
export function poisson_log_density(x: number, lambda: number): number {
  const k = Math.floor(x);
  if (x < 0) {
    return 0.0;
  }

  if (lambda == 0.0) {
    return k == 0 ? 1.0 : 0.0;
  }

  let num = k * Math.log(lambda) - lambda;

  let den = 0;
  for (let xi = 2; xi <= k; xi++) {
    den += Math.log(xi);
  }

  return num - den;
}

function _quantile_cornish_fisher_poisson(u: number, lambda: number): number {
  // https://en.wikipedia.org/wiki/Cornish–Fisher_expansion
  let x = normal_quantile(u, 0.0, 1.0);
  let m = lambda;
  let v = m;
  let kappa3 = m;
  let gamma1 = kappa3 / v ** 1.5;
  let He2 = x * x - 1.0;
  let h1 = He2 / 6.0;
  let wp = x + gamma1 * h1;
  return m + Math.sqrt(v) * wp;
}

function _quantile_poisson_root(x: number, u: number, lambda: number): number {
  return poisson_distribution(x, lambda) - u;
}

/**
 * Poisson(lambda) distribution's quantile function
 *
 * @param p probability
 * @param lambda rate parameter
 * @returns poisson distribution at p
 *
 * @category poisson
 */
export function poisson_quantile(p: number, lambda: number): number {
  if (p == 0.0) {
    return 0.0;
  }
  if (lambda < 0.0) {
    return Number.NaN;
  }
  let qhat = _quantile_cornish_fisher_poisson(p, lambda);
  let sd = Math.sqrt(lambda);
  let lower = Math.floor(qhat - sd);
  if (lower < 0.0) {
    lower = 0.0;
  }
  let upper = Math.ceil(qhat + sd);
  let q = brent(lower, upper, p, 1e-8, _quantile_poisson_root, lambda);
  return Math.round(q);
}

/**
 * Poisson(lambda) random numbers
 *
 * @param N quantity of random numbers
 * @param lambda rate parameter
 * @returns N poisson random numbers
 *
 * @category poisson
 */
export function poisson_random(N: number, lambda: number): number[] {
  let out: number[] = new Array(N).fill(0);

  if (lambda < 10) {
    // Devroye X.3 Poisson generator based upon the inversion by
    // sequential search
    let u = 0.0;
    let expnl = Math.exp(-lambda);
    let p = expnl;
    for (let n = 0; n < N; n++) {
      out[n] = 0;
      p = expnl;
      u = Math.random();
      while (u > p) {
        u -= p;
        out[n] += 1;
        p *= lambda / out[n];
      }
    }
    return out;
  } else {
    // don't yet know.  Same place as Binomial rng?
    // Option: if can't find source of below, use
    // IGNPOI from https://people.sc.fsu.edu/~jburkardt/f_src/ranlib/ranlib.f90
    // which is MIT licensed and the same algorithm from R
    const logfactorial = [
      0.0, 0.0, 0.693147180559945, 1.791759469228055, 3.178053830347946,
      4.787491742782046, 6.579251212010101, 8.525161361065415,
      10.60460290274525, 12.801827480081469,
    ];

    let smu = Math.sqrt(lambda);
    let b = 0.931 + 2.53 * smu;
    let a = -0.059 + 0.02483 * b;
    let invalpha = 1.1239 + 1.1328 / (b - 3.4);
    let vr = 0.9227 - 3.6224 / (b - 2.0);
    let v = 0.0;
    let u = 0.0;
    let us = 0.0;
    let k = 0.0;

    for (let n = 0; n < N; n++) {
      while (1) {
        v = Math.random();
        if (v <= 0.86 * vr) {
          u = v / vr - 0.43;
          out[n] = Math.floor(
            ((2 * a) / (0.5 - Math.abs(u)) + b) * u + lambda + 0.445,
          );
          break;
        }

        if (v >= vr) {
          u = Math.random() - 0.5;
        } else {
          u = v / vr - 0.93;
          u = (u < 0 ? -0.5 : 0.5) - u;
          v = Math.random() * vr;
        }

        us = 0.5 - Math.abs(u);
        if (us < 0.013 && v > us) {
          continue;
        }

        k = Math.floor(((2 * a) / us + b) * u + lambda + 0.445);
        v = (v * invalpha) / (a / (us * us) + b);
        if (k >= 10.0) {
          if (
            Math.log(v * smu) <=
            (k + 0.5) * Math.log(lambda / k) -
              lambda -
              LOG_SQRT_2PI +
              k -
              (1.0 / 12.0 - (1.0 / 360.0 - 1.0 / (1260.0 * k * k)) / (k * k)) /
                k
          ) {
            out[n] = k;
            break;
          }
        } else if (k >= 0.0) {
          if (Math.log(v) <= k * Math.log(lambda) - lambda - logfactorial[k]) {
            out[n] = k;
            break;
          }
        }
      }
    }
  }
  return out;
}
