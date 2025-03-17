import { bratio } from "./beta_inc_fn";
import { brent } from "./brent";
import { log_gamma } from "./log_gamma_fn";
import { normal_quantile } from "./normal";
import { PI2 } from "./constants";

/**
 * Binomial(K, p) distribution's density function
 *
 * @param x argument to density function
 * @param K number of trials
 * @param p probability of success for each trial
 * @returns binomial density at x
 *
 * @category binomial
 */
export function binomial_density(x: number, K: number, p: number): number {
  return Math.exp(binomial_log_density(x, K, p));
}

/**
 * Binomial(K, p) distribution function
 *
 * @param x argument to distribution function
 * @param K number of trials
 * @param p probability of success for each trial
 * @returns binomial distribution at x
 *
 * @category binomial
 */
export function binomial_distribution(x: number, K: number, p: number): number {
  const k = Math.floor(x);
  if (k <= 0) {
    return 0;
  }
  if (k >= K) {
    return 1;
  }

  return bratio(K - k, k + 1, 1 - p);
}

// Derived from
// Catherine Loader (2000). _Fast and Accurate Computation of Binomial
// Probabilities
// obtained at http://cm.bell-labs.com/stat/catherine/research.html
// via the Internet Archive's (https://archive.org) Wayback Machine in 2024.

function stirling_error(n: number): number {
  if (n < 16) {
    return log_gamma(n + 1.0) + n - n * Math.log(n) - 0.5 * Math.log(PI2 * n);
  }

  const nn = n * n;
  const S0 = 0.083333333333333333333; // 1/12
  const S1 = 0.00277777777777777777778; // 1/360
  const S2 = 0.00079365079365079365079365; // 1/1260
  const S3 = 0.000595238095238095238095238; // 1/1680
  const S4 = 0.0008417508417508417508417508; // 1/1188

  if (n > 500) {
    return (S0 - S1 / nn) / n;
  }

  if (n > 80) {
    return (S0 - (S1 - S2 / nn) / nn) / n;
  }

  if (n > 35) {
    return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  }

  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

function bd0(x: number, np: number): number {
  let ej, s, s1, v, j;

  if (Math.abs(x - np) < 0.1 * (x + np)) {
    s = ((x - np) * (x - np)) / (x + np);
    v = (x - np) / (x + np);
    ej = 2 * x * v;
    v = v * v;
    for (j = 1; ; j++) {
      ej *= v;
      s1 = s + ej / ((j << 1) + 1);
      if (s1 == s) {
        return s1;
      }
      s = s1;
    }
  }
  return x * Math.log(x / np) + np - x;
}

/**
 * Binomial(K, p) distribution's log-density function
 *
 * @param x argument to log-density function
 * @param K number of trials
 * @param p probability of success for each trial
 * @returns binomial log-density at x
 *
 * @category binomial
 */
export function binomial_log_density(x: number, K: number, p: number): number {
  if (x < 0 || x > K) {
    return -Infinity;
  }

  if (p == 0.0) {
    return x == 0 ? 0 : -Infinity;
  }

  if (p == 1.0) {
    return x == K ? 0 : -Infinity;
  }

  if (x == 0) {
    return K * Math.log1p(-p);
  }

  if (x == K) {
    return K * Math.log(p);
  }

  let lc = stirling_error(K);
  lc -= stirling_error(x);
  lc -= stirling_error(K - x);
  lc -= bd0(x, K * p) + bd0(K - x, K * (1.0 - p));

  let out = lc + 0.5 * Math.log(K / (PI2 * x * (K - x)));
  return out;
}

function _quantile_cornish_fisher_binomial(
  u: number,
  K: number,
  p: number,
): number {
  // https://en.wikipedia.org/wiki/Cornish–Fisher_expansion
  let x = normal_quantile(u, 0, 1);
  let m = K * p;
  let v = m * (1 - p);
  let kappa3 = v - 2 * v * p; // from https://arxiv.org/abs/2012.06270
  let gamma1 = kappa3 / v ** 1.5;
  let He2 = x * x - 1;
  let h1 = He2 / 6;
  let wp = x + gamma1 * h1;
  return m + Math.sqrt(v) * wp;
}

function _quantile_binomial_root(
  x: number,
  u: number,
  K: number,
  p: number,
): number {
  return binomial_distribution(x, K, p) - u;
}

/**
 * Binomial(K, p) distribution's quantile function
 *
 * @param u probability
 * @param K number of trials
 * @param p probability of success for each trial
 * @returns binomial distribution at u
 *
 * @category binomial
 */
export function binomial_quantile(u: number, K: number, p: number): number {
  if (u == 0.0) {
    return 0.0;
  }

  if (u == 1.0) {
    return K;
  }

  let qhat = _quantile_cornish_fisher_binomial(u, K, p);
  let sd = Math.sqrt(K * p * (1 - p));
  let lower = Math.floor(qhat - sd);

  if (lower < 0) {
    lower = 0.0;
  }
  let upper = Math.ceil(qhat + sd);
  if (upper > K) {
    upper = K;
  }

  let q = brent(lower, upper, u, 1e-8, _quantile_binomial_root, K, p);
  return Math.round(q);
}

function fc(k: number): number {
  const invkp1 = 1.0 / (k + 1.0);
  const invkp1_sq = invkp1 * invkp1;
  return (
    (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0) * invkp1_sq) * invkp1_sq) *
    invkp1
  );
}

/**
 * Binomial(K, p) random numbers
 *
 * @param N quantity of random numbers
 * @param K number of trials
 * @param p probability of success for each trial
 * @returns N binomial random numbers
 *
 * @category binomial
 */
export function binomial_random(N: number, K: number, p: number): number[] {
  let out: number[] = new Array(N).fill(0);

  let E = K * p;
  let p_big = p > 0.5;
  let p_ = p_big ? 1.0 - p : p;

  if (E < 10) {
    // Devroy, X.4 Binomial Distribution: first waiting time algorithm
    let ot = 0;
    let sum = 0;
    let g = 0;
    for (let n = 0; n < N; n++) {
      ot = -1.0;
      sum = 0.0;
      while (sum <= K) {
        g = Math.ceil(Math.log(Math.random()) / Math.log1p(-p_));
        sum += g;
        ot += 1;
      }
      out[n] = p_big ? K - ot : ot;
    }
  } else {
    // don't yet know where I got this. That's annoying.
    // Check ~/Documents/research/dists_literature
    // Option: An MIT licensed version of what R uses in FORTRAN90
    // https://people.sc.fsu.edu/~jburkardt/f_src/ranlib/ranlib.f90
    // to replace what's below. which has citation
    // Voratas Kachitvichyanukul and Bruce
    // W. Schmeiser. 1988. Binomial random variate
    // generation. Commun. ACM 31, 2 (Feb. 1988),
    // 216–222. https://doi.org/10.1145/42372.42381

    const binomial_table = [
      0.08106146679532726, 0.04134069595540929, 0.02767792568499834,
      0.02079067210376509, 0.01664469118982119, 0.01387612882307075,
      0.01189670994589177, 0.01041126526197209, 0.009255462182712733,
      0.008330563433362871,
    ];

    let m = Math.floor((K + 1) * p_);
    let r = p_ / (1.0 - p_);
    let kr = (K + 1.0) * r;
    let kpq = K * p_ * (1.0 - p_);
    let sqrt_kpq = Math.sqrt(kpq);
    let b = 1.15 + 2.53 * sqrt_kpq;
    let a = -0.0873 + 0.0248 * b + 0.01 * p_;
    let c = K * p_ + 0.5;
    let alpha = (2.83 + 5.1 / b) * sqrt_kpq;
    let vr = 0.92 - 4.2 / b;
    let urvr = 0.86 * vr;

    let v = 0.0;
    let u = 0.0;
    let us = 0.0;
    let k = 0.0;
    let km = 0.0;
    let rho = 0.0;
    let t = 0.0;
    let f = 0.0;
    let h = 0.0;
    let Kk = 0.0;
    let fc1 = 0.0;
    let fc2 = 0.0;
    let ot = 0.0;

    for (let n = 0; n < N; n++) {
      while (true) {
        v = Math.random();
        if (v <= urvr) {
          u = v / vr - 0.43;
          ot = Math.floor(((2 * a) / (0.5 - Math.abs(u)) + b) * u + c);
          out[n] = p_big ? K - ot : ot;
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
        k = Math.floor(((2 * a) / us + b) * u + c);
        if (k < 0 || k > K) {
          continue;
        }

        v = (v * alpha) / (a / (us * us) + b);
        km = Math.abs(k - m);
        if (km <= 15) {
          // 3.1
          f = 1.0;
          let i = 0;
          if (m < k) {
            i = m;
            while (i != k) {
              i += 1;
              f = f * (kr / i - r);
            }
          } else {
            i = k;
            while (i != m) {
              i += 1;
              v = v * (kr / i - r);
            }
          }

          if (v <= f) {
            out[n] = p_big ? K - k : k;
            break;
          } else {
            continue;
          }
        } else {
          // 3.2
          v = Math.log(v);
          rho =
            (km / kpq) * (((km / 3.0 + 0.625) * km + 1.0 / 6.0) / kpq + 0.5);
          t = (-km * km) / (2.0 * kpq);
          if (v < t - rho) {
            out[n] = p_big ? K - k : k;
          }
          if (v > t + rho) {
            continue;
          }

          km = K - m + 1.0;
          fc1 = m > 9.9 ? fc(m) : binomial_table[Math.floor(m)];
          fc2 = K - m > 9.0 ? fc(K - m) : binomial_table[Math.floor(K - m)];
          h = (m + 0.5) * Math.log((m + 1.0) / (r * km)) + fc1 + fc2;

          Kk = K - k + 1.0;
          fc1 = k > 9.0 ? fc(k) : binomial_table[Math.floor(k)];
          fc2 = K - k > 9.0 ? fc(K - k) : binomial_table[Math.floor(K - k)];
          if (
            v <=
            h +
              (K + 1.0) * Math.log(km / Kk) +
              (k + 0.5) * Math.log((Kk * r) / (k + 1.0)) -
              fc1 -
              fc2
          ) {
            out[n] = p_big ? K - k : k;
          } else {
            continue;
          }
        }
      }
    }
  }
  return out;
}
