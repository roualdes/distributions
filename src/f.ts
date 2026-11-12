import { betaln, bratio, bratio_complement, rlog1 } from "./beta_inc_fn";
import { xinbeta, xinbeta_complement } from "./beta_inv_fn";
import { gamma_random } from "./gamma";
import { stirling_error } from "./stirling_error";
import { LOG_2PI } from "./constants";

function valid_degrees(df1: number, df2: number): boolean {
  return Number.isFinite(df1) && df1 > 0 && Number.isFinite(df2) && df2 > 0;
}

function log_beta_degrees(df1: number, df2: number): number {
  const a = df1 / 2;
  const b = df2 / 2;
  if (a === 0 && b === 0) return Math.log(df1) - 2 * (Math.log(df1) - Math.LN2);
  if (a === 0) return Math.LN2 - Math.log(df1);
  if (b === 0) return Math.LN2 - Math.log(df2);
  return betaln(a, b);
}

// Both beta coordinates are computed from the smaller odds, avoiding overflow
// and preserving the small coordinate when the other rounds to one.
function beta_coordinates(logOdds: number): [number, number, number, number] {
  if (logOdds <= 0) {
    const odds = Math.exp(logOdds);
    const logY = -Math.log1p(odds);
    return [odds / (1 + odds), 1 / (1 + odds), logOdds + logY, logY];
  }
  const odds = Math.exp(-logOdds);
  const logZ = -Math.log1p(odds);
  return [1 / (1 + odds), odds / (1 + odds), logZ, -logOdds + logZ];
}

/**
 * Central F density on [0, Infinity).
 * @param x argument to the density function
 * @param df1 positive, finite numerator degrees of freedom
 * @param df2 positive, finite denominator degrees of freedom
 * @returns density at x, or NaN for invalid arguments
 * @category F
 */
export function f_density(x: number, df1: number, df2: number): number {
  return Math.exp(f_log_density(x, df1, df2));
}

/**
 * Central F log-density on [0, Infinity).
 * @param x argument to the log-density function
 * @param df1 positive, finite numerator degrees of freedom
 * @param df2 positive, finite denominator degrees of freedom
 * @returns log-density at x, or NaN for invalid arguments
 * @category F
 */
export function f_log_density(x: number, df1: number, df2: number): number {
  if (!valid_degrees(df1, df2) || Number.isNaN(x)) return Number.NaN;
  if (x < 0 || x === Infinity) return -Infinity;
  if (x === 0) return df1 < 2 ? Infinity : df1 === 2 ? 0 : -Infinity;

  const a = df1 / 2;
  const b = df2 / 2;
  const logX = Math.log(x);
  const [z, y, logZ, logY] = beta_coordinates(
    Math.log(df1) - Math.log(df2) + logX,
  );
  if (a === 1) {
    return Math.log(df1) - Math.log(df2) + logX < -36
      ? -x - x / b
      : (b + 1) * logY;
  }
  if (a >= 8 && b >= 8) {
    // Stirling corrections and beta deviances avoid cancellation between
    // terms of order df when the density is concentrated near x=1.
    const sum = a + b;
    const p = a / sum;
    const q = b / sum;
    const e1 =
      Math.abs(logX) < 0.5 ? (x - 1) * y : Math.expm1(logZ - Math.log(p));
    const e2 =
      Math.abs(logX) < 0.5
        ? -((x - 1) / x) * z
        : Math.expm1(logY - Math.log(q));
    const d1 = Math.abs(e1) <= 0.6 ? rlog1(e1) : e1 - logZ + Math.log(p);
    const d2 = Math.abs(e2) <= 0.6 ? rlog1(e2) : e2 - logY + Math.log(q);
    return (
      0.5 * (Math.log(a) + Math.log(b) - Math.log(sum) - LOG_2PI) +
      stirling_error(sum) -
      stirling_error(a) -
      stirling_error(b) -
      a * d1 -
      b * d2 -
      logX
    );
  }
  return a * logZ + b * logY - logX - log_beta_degrees(df1, df2);
}

/**
 * Central F cumulative distribution P(F <= x).
 * @param x argument to the distribution function
 * @param df1 positive, finite numerator degrees of freedom
 * @param df2 positive, finite denominator degrees of freedom
 * @returns cumulative probability at x, or NaN for invalid arguments
 * @category F
 */
export function f_distribution(x: number, df1: number, df2: number): number {
  if (!valid_degrees(df1, df2) || Number.isNaN(x)) return Number.NaN;
  if (x <= 0) return 0;
  if (x === Infinity) return 1;
  const a = df1 / 2;
  const b = df2 / 2;
  if (a === 0 || b === 0) return a === b ? 0.5 : a === 0 ? 1 : 0;
  const logOdds = Math.log(df1) - Math.log(df2) + Math.log(x);
  const [z, y, logZ, logY] = beta_coordinates(logOdds);
  if (logOdds <= 0) {
    if (z === 0) return Math.exp(a * logZ - Math.log(a) - betaln(a, b));
    return bratio(a, b, z);
  }
  if (y === 0) return -Math.expm1(b * logY - Math.log(b) - betaln(a, b));
  return bratio_complement(b, a, y);
}

/**
 * Central F quantile (inverse cumulative distribution).
 * @param p probability in [0, 1]
 * @param df1 positive, finite numerator degrees of freedom
 * @param df2 positive, finite denominator degrees of freedom
 * @returns quantile at p; endpoints give 0 and Infinity, invalid arguments give NaN
 * @category F
 */
export function f_quantile(p: number, df1: number, df2: number): number {
  if (!valid_degrees(df1, df2) || !(p >= 0 && p <= 1)) return Number.NaN;
  if (p === 0) return 0;
  if (p === 1) return Infinity;

  // Reciprocal symmetry keeps the supplied probability in the smaller tail.
  const reflect = p > 0.5;
  const a = (reflect ? df2 : df1) / 2;
  const b = (reflect ? df1 : df2) / 2;
  if (a === 0 || b === 0) {
    if (a === b) return p === 0.5 ? 1 : reflect ? Infinity : 0;
    return (a === 0) !== reflect ? 0 : Infinity;
  }
  const probability = reflect ? 1 - p : p;
  const logScale = reflect
    ? Math.log(df1) - Math.log(df2)
    : Math.log(df2) - Math.log(df1);
  const z = xinbeta(a, b, probability);
  let logOdds: number;
  if (z <= 0.5) {
    const logZ =
      z === 0
        ? (Math.log(probability) + Math.log(a) + betaln(a, b)) / a
        : Math.log(z);
    logOdds = logZ - Math.log1p(-z);
  } else {
    const y = xinbeta_complement(b, a, probability);
    const logY =
      y === 0
        ? (Math.log1p(-probability) + Math.log(b) + betaln(b, a)) / b
        : Math.log(y);
    logOdds = Math.log1p(-y) - logY;
  }
  return Math.exp((reflect ? -1 : 1) * (logScale + logOdds));
}

/**
 * Central F random numbers from independent normalized chi-squared draws.
 * @param N nonnegative integer quantity of random numbers
 * @param df1 positive, finite numerator degrees of freedom
 * @param df2 positive, finite denominator degrees of freedom
 * @returns N random numbers, or N NaNs for invalid degrees of freedom
 * @throws RangeError if N is not a nonnegative integer
 * @category F
 */
export function f_random(N: number, df1: number, df2: number): number[] {
  if (!Number.isInteger(N) || N < 0)
    throw new RangeError("N must be a nonnegative integer.");
  if (!valid_degrees(df1, df2)) return new Array(N).fill(Number.NaN);
  // Tiny gamma shapes can both underflow to zero before their ratio is taken.
  if (df1 < 1 || df2 < 1) {
    return Array.from({ length: N }, () => f_quantile(Math.random(), df1, df2));
  }
  const numerator = gamma_random(N, df1 / 2, df1 / 2);
  const denominator = gamma_random(N, df2 / 2, df2 / 2);
  return numerator.map((x, i) => x / denominator[i]);
}
