import { bratio } from "./beta_inc_fn";
import { xinbeta } from "./beta_inv_fn";
import { log_gamma } from "./log_gamma_fn";
import { binomial_log_density } from "./binomial";

/**
 * Beta(a, b) distribution's density function
 *
 * @param x argument to density function
 * @param a measure of shape
 * @param b measure of shape
 * @returns beta density at x
 *
 * @category beta
 */
export function beta_density(x: number, a: number, b: number): number {
  return Math.exp(beta_log_density(x, a, b));
}

/**
 * Beta(a, b) distribution function
 *
 * @param x argument to distribution function
 * @param a measure of shape
 * @param b measure of shape
 * @returns beta distribution at x
 *
 * @category beta
 */
export function beta_distribution(x: number, a: number, b: number): number {
  return bratio(a, b, x);
}

/**
 * Beta(a, b) distribution's log-density function
 *
 * @param x argument to log-density function
 * @param a measure of shape
 * @param b measure of scale
 * @returns beta log-density at x
 *
 * @category beta
 */
export function beta_log_density(x: number, a: number, b: number): number {
  let ld = 0.0;
  if (a > 1.0 || b > 1.0) {
    ld = binomial_log_density(a - 1.0, a + b - 1.0, x);
    ld += Math.log(b) - Math.log1p(-x);
  } else {
    ld = log_gamma(a + b) - log_gamma(a) - log_gamma(b);
    ld += (a - 1.0) * Math.log(x) + (b - 1.0) * Math.log1p(-x);
  }
  return ld;
}

/**
 * Beta(a, b) distribution's quantile function
 *
 * @param p probability
 * @param a measure of shape
 * @param b measure of shape
 * @returns beta quantile at p
 *
 * @category beta
 */
export function beta_quantile(p: number, a: number, b: number): number {
  return xinbeta(a, b, p);
}

// R. C. H. Cheng. 1978. Generating beta variates with nonintegral
// shape parameters. Commun. ACM 21, 4 (April 1978),
// 317–322. https://doi.org/10.1145/359460.359482
// an alternative implementation at
// https://people.math.sc.edu/Burkardt/c_src/ranlib/ranlib.html
// see GENBET

/**
 * Beta(a, b) random numbers
 *
 * @param N quantity of random numbers
 * @param a measure of shape
 * @param b measure of shape
 * @returns N beta random numbers
 *
 * @category beta
 */
export function beta_random(N: number, a: number, b: number): number[] {
  const log4 = 1.3862944;
  const log5p1 = 2.609438;

  let a0 = a;
  let b0 = b;
  let alpha;
  let beta;
  let gamma;
  let delta;
  let k1;
  let k2;
  let r;
  let s;
  let t;
  let u1;
  let u2;
  let v;
  let w;
  let y;
  let z;

  let out: number[] = new Array(N).fill(0);

  for (let n = 0; n < N; ++n) {
    a0 = a;
    b0 = b;

    if (a > 1.0 && b > 1.0) {
      if (a > b) {
        a0 = b;
        b0 = a;
      }

      alpha = a0 + b0;
      beta = Math.sqrt((alpha - 2.0) / (2 * a0 * b0 - alpha));
      gamma = a0 + 1.0 / beta;

      while (true) {
        u1 = Math.random();
        u2 = Math.random();

        v = beta * Math.log(u1 / (1.0 - u1));
        w = a0 * Math.exp(v);
        z = u1 * u1 * u2;
        r = gamma * v - log4;
        s = a0 + r - w;

        if (s + log5p1 >= 5.0 * z) {
          break;
        }

        t = Math.log(z);
        if (s >= t) {
          break;
        }

        if (r + alpha * Math.log(alpha / (b0 + w)) >= t) {
          break;
        }
      }
    } else {
      if (b > a) {
        a0 = b;
        b0 = a;
      }

      alpha = a0 + b0;
      beta = 1.0 / b0;
      delta = 1.0 + a0 - b0;
      k1 = (delta * (1.0 / 72.0 + b0 / 24.0)) / (a0 * beta - 7.0 / 9.0);
      k2 = 0.25 + (0.5 + 0.25 / delta) * b0;

      while (true) {
        u1 = Math.random();
        u2 = Math.random();

        if (u1 < 0.5) {
          y = u1 * u2;
          z = u1 * y;

          if (k1 <= 0.25 * u2 + z - y) {
            continue;
          }
        } else {
          z = u1 * u1 * u2;

          if (z <= 0.25) {
            v = beta * Math.log(u1 / (1.0 - u1));
            w = a0 * Math.exp(v);

            if (a0 == a) {
              out[n] = w / (b0 + w);
            } else {
              out[n] = b0 / (b0 + w);
            }
            break;
          }

          if (k2 < z) {
            continue;
          }
        }

        v = beta * Math.log(u1 / (1.0 - u1));
        w = a0 * Math.exp(v);

        if (Math.log(z) <= alpha * (Math.log(alpha / (b0 + w)) + v) - log4) {
          break;
        }
      }
    }

    if (a0 == a) {
      out[n] = w / (b0 + w);
    } else {
      out[n] = b0 / (b0 + w);
    }
  }

  return out;
}
