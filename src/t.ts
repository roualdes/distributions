import { betaln, bratio, bratio_complement } from "./beta_inc_fn";
import { xinbeta } from "./beta_inv_fn";
import { gamma_random } from "./gamma";
import {
  normal_distribution,
  normal_log_density,
  normal_quantile,
  normal_random,
} from "./normal";

/**
 * Central Student's t density, with zero location and unit scale.
 * @param x argument to the density function
 * @param df positive degrees of freedom (Infinity gives the standard normal)
 * @returns density at x, or NaN for invalid arguments
 * @category t
 */
export function t_density(x: number, df: number): number {
  return Math.exp(t_log_density(x, df));
}

/**
 * Central Student's t log-density, with zero location and unit scale.
 * @param x argument to the log-density function
 * @param df positive degrees of freedom (Infinity gives the standard normal)
 * @returns log-density at x, or NaN for invalid arguments
 * @category t
 */
export function t_log_density(x: number, df: number): number {
  if (!(df > 0) || Number.isNaN(x)) return Number.NaN;
  if (!Number.isFinite(x)) return -Infinity;
  if (df === Infinity) return normal_log_density(x, 0, 1);

  const a = df / 2;
  const normalizer =
    a === 0
      ? 0.5 * Math.log(df) - Math.LN2
      : -0.5 * Math.log(df) - betaln(a, 0.5);
  if (x === 0) return normalizer;

  // Evaluate log(1 + x^2/df) without overflowing x^2 or x^2/df.
  const logRatio = 2 * Math.log(Math.abs(x)) - Math.log(df);
  if (logRatio < -36) {
    // Avoid underflow of x^2/df before multiplication by a large df.
    return normalizer - 0.5 * x * x - 0.5 * x * (x / df);
  }
  const kernel =
    logRatio > 0
      ? logRatio + Math.log1p(Math.exp(-logRatio))
      : Math.log1p(Math.exp(logRatio));
  return normalizer - (0.5 * df + 0.5) * kernel;
}

/**
 * Central Student's t cumulative distribution P(T <= x).
 * @param x argument to the distribution function
 * @param df positive degrees of freedom (Infinity gives the standard normal)
 * @returns cumulative probability at x, or NaN for invalid arguments
 * @category t
 */
export function t_distribution(x: number, df: number): number {
  if (!(df > 0) || Number.isNaN(x)) return Number.NaN;
  if (x === -Infinity) return 0;
  if (x === Infinity) return 1;
  if (x === 0) return 0.5;
  if (df === Infinity) return normal_distribution(x, 0, 1);

  const a = df / 2;
  if (a === 0) return 0.5;
  const logRatio = 2 * Math.log(Math.abs(x)) - Math.log(df);
  if (logRatio < 0) {
    // Near the median, use x^2/(df+x^2); df/(df+x^2) can round to one.
    const ratio = Math.exp(logRatio);
    if (ratio === 0) return normal_distribution(x, 0, 1);
    const center = 0.5 * bratio(0.5, a, ratio / (1 + ratio));
    if (center > 0.25) {
      const tail = 0.5 * bratio_complement(0.5, a, ratio / (1 + ratio));
      return x < 0 ? tail : 1 - tail;
    }
    return x < 0 ? 0.5 - center : 0.5 + center;
  }

  const inverseRatio = Math.exp(-logRatio);
  const z = inverseRatio / (1 + inverseRatio);
  // If z underflows, the leading beta power can still give a finite t tail.
  const tail =
    z === 0
      ? Math.exp(-Math.LN2 - a * logRatio - Math.log(a) - betaln(a, 0.5))
      : 0.5 * bratio(a, 0.5, z);
  return x < 0 ? tail : 1 - tail;
}

/**
 * Central Student's t quantile (inverse cumulative distribution).
 * @param p probability in [0, 1]
 * @param df positive degrees of freedom (Infinity gives the standard normal)
 * @returns quantile at p; endpoints give infinities, invalid arguments give NaN
 * @category t
 */
export function t_quantile(p: number, df: number): number {
  if (!(df > 0) || !(p >= 0 && p <= 1)) return Number.NaN;
  if (p === 0) return -Infinity;
  if (p === 1) return Infinity;
  if (p === 0.5) return 0;
  if (df === Infinity) return normal_quantile(p, 0, 1);
  if (df >= 1e8) {
    // Cornish-Fisher expansion through df^-2. Here even the largest normal
    // quantile representable by p has a negligible df^-3 correction. This
    // avoids losing precision when the inverse-beta result rounds near one.
    const z = normal_quantile(p, 0, 1);
    const z2 = z * z;
    return (
      z +
      (z * (z2 + 1)) / (4 * df) +
      (z * (5 * z2 * z2 + 16 * z2 + 3)) / (96 * df) / df
    );
  }

  const sign = p < 0.5 ? -1 : 1;
  const a = df / 2;
  if (a === 0) return sign * Infinity;
  const tail = Math.min(p, 1 - p);
  let magnitude: number;
  const center = 2 * Math.abs(p - 0.5);
  const y = tail < 0.25 ? 1 : xinbeta(0.5, a, center);
  if (y > 0.5) {
    const z = xinbeta(a, 0.5, 2 * tail);
    if (z === 0) {
      const logZ = (Math.log(2 * tail) + Math.log(a) + betaln(a, 0.5)) / a;
      magnitude = Math.exp(0.5 * (Math.log(df) - logZ));
    } else {
      magnitude = (Math.sqrt(df) * Math.sqrt(1 - z)) / Math.sqrt(z);
    }
  } else {
    // Invert the central probability directly, avoiding 1 - z near zero.
    if (y === 0) {
      const logY = 2 * (Math.log(center) - Math.LN2 + betaln(0.5, a));
      magnitude = Math.exp(0.5 * (Math.log(df) + logY));
    } else {
      magnitude = (Math.sqrt(df) * Math.sqrt(y)) / Math.sqrt(1 - y);
    }
  }
  // Transforming a beta quantile near one loses a few bits in 1-z. Refine in
  // t coordinates, using the small CDF tail and the log-density directly.
  if (Number.isFinite(magnitude)) {
    for (let iteration = 0; iteration < 2; iteration++) {
      const error = t_distribution(-magnitude, df) - tail;
      if (error === 0) break;
      const step =
        Math.sign(error) *
        Math.exp(Math.log(Math.abs(error)) - t_log_density(magnitude, df));
      const next = magnitude + step;
      if (!(next > 0) || !Number.isFinite(next)) break;
      magnitude = next;
    }
  }
  return sign * magnitude;
}

/**
 * Central Student's t random numbers, with zero location and unit scale.
 * @param N nonnegative integer quantity of random numbers
 * @param df positive degrees of freedom (Infinity gives the standard normal)
 * @returns N random numbers, or N NaNs for invalid degrees of freedom
 * @throws RangeError if N is not a nonnegative integer
 * @category t
 */
export function t_random(N: number, df: number): number[] {
  if (!Number.isInteger(N) || N < 0) {
    throw new RangeError("N must be a nonnegative integer.");
  }
  if (!(df > 0)) return new Array(N).fill(Number.NaN);
  if (df === Infinity) return normal_random(N, 0, 1);
  if (df / 2 === 0) {
    return Array.from({ length: N }, () => t_quantile(Math.random(), df));
  }

  const normal = normal_random(N, 0, 1);
  // Gamma(df/2, rate=df/2) is chi-squared(df)/df. Normalize at generation
  // time to avoid overflowing a chi-squared draw when df is very large.
  const scaledChiSquared = gamma_random(N, df / 2, df / 2);
  return normal.map((z, i) => z / Math.sqrt(scaledChiSquared[i]));
}
