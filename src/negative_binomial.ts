import { betaln, bratio, bratio_complement, rlog1 } from "./beta_inc_fn";
import { binomial_log_density } from "./binomial";
import {
  geometric_distribution,
  geometric_log_density,
  geometric_quantile,
  geometric_random,
} from "./geometric";
import { gamma_random } from "./gamma";
import { poisson_random } from "./poisson";
import { normal_quantile } from "./normal";
import { stirling_error } from "./stirling_error";
import { LOG_2PI } from "./constants";

function valid_parameters(size: number, p: number): boolean {
  return Number.isFinite(size) && size > 0 && p > 0 && p <= 1;
}

/**
 * Negative binomial probability mass: failures before size successes.
 * @param x nonnegative integer failure count
 * @param size positive, finite success target (fractional values are allowed)
 * @param p success probability, 0 < p <= 1
 * @returns probability mass at x, or NaN for invalid parameters
 * @category negative binomial
 */
export function negative_binomial_density(
  x: number,
  size: number,
  p: number,
): number {
  return Math.exp(negative_binomial_log_density(x, size, p));
}

/**
 * Negative binomial log-probability mass: failures before size successes.
 * @param x nonnegative integer failure count
 * @param size positive, finite success target (fractional values are allowed)
 * @param p success probability, 0 < p <= 1
 * @returns log-probability mass at x, or NaN for invalid parameters
 * @category negative binomial
 */
export function negative_binomial_log_density(
  x: number,
  size: number,
  p: number,
): number {
  if (!valid_parameters(size, p) || Number.isNaN(x)) return Number.NaN;
  if (!Number.isInteger(x) || x < 0) return -Infinity;
  if (p === 1) return x === 0 ? 0 : -Infinity;
  if (size === 1) return geometric_log_density(x, p);
  if (x === 0) return size * Math.log(p);

  const maximum = Math.max(size, x);
  const logSum = Math.log(maximum) + Math.log1p(Math.min(size, x) / maximum);
  if (size >= 8 && x >= 8) {
    if (maximum < 1e150) {
      // Supply the smaller probability explicitly: constructing 1-p and
      // subtracting it again can lose many bits when p is small.
      const mass =
        p <= 0.5
          ? binomial_log_density(size, size + x, p)
          : binomial_log_density(x, size + x, 1 - p);
      return mass + Math.log(size) - logSum;
    }
    // Logarithmic saddlepoint form avoids products and sums overflowing at
    // very large counts, including when size+x is beyond Number.MAX_VALUE.
    const scaledSum = size / maximum + x / maximum;
    const fraction = size / maximum / scaledSum;
    const complement = x / maximum / scaledSum;
    const d1 = Math.log(p) - Math.log(fraction);
    const d2 = Math.log1p(-p) - Math.log(complement);
    const e1 = Math.expm1(d1);
    const e2 = Math.expm1(d2);
    const deviance1 = Math.abs(e1) <= 0.6 ? rlog1(e1) : e1 - d1;
    const deviance2 = Math.abs(e2) <= 0.6 ? rlog1(e2) : e2 - d2;
    const correction = Number.isFinite(size + x) ? stirling_error(size + x) : 0;
    return (
      0.5 * (Math.log(size) - Math.log(x) - logSum - LOG_2PI) +
      correction -
      stirling_error(size) -
      stirling_error(x) -
      size * deviance1 -
      x * deviance2
    );
  }
  return size * Math.log(p) + x * Math.log1p(-p) - logSum - betaln(size, x + 1);
}

/**
 * Negative binomial cumulative probability of at most floor(x) failures.
 * @param x failure count (nonintegers are rounded down)
 * @param size positive, finite success target (fractional values are allowed)
 * @param p success probability, 0 < p <= 1
 * @returns cumulative probability at x, or NaN for invalid parameters
 * @category negative binomial
 */
export function negative_binomial_distribution(
  x: number,
  size: number,
  p: number,
): number {
  if (!valid_parameters(size, p) || Number.isNaN(x)) return Number.NaN;
  if (x < 0) return 0;
  if (x === Infinity || p === 1) return 1;
  if (size === 1) return geometric_distribution(x, p);
  if (x < 1) return Math.exp(size * Math.log(p));
  return bratio(size, Math.floor(x) + 1, p);
}

/**
 * Smallest failure count whose negative binomial CDF is at least u.
 * @param u cumulative probability in [0, 1]
 * @param size positive, finite success target (fractional values are allowed)
 * @param p success probability, 0 < p <= 1
 * @returns failure count; u=1 gives Infinity unless p=1
 * @category negative binomial
 */
export function negative_binomial_quantile(
  u: number,
  size: number,
  p: number,
): number {
  if (!valid_parameters(size, p) || !(u >= 0 && u <= 1)) return Number.NaN;
  if (u === 0 || p === 1) return 0;
  if (u === 1) return Infinity;
  if (size === 1) return geometric_quantile(u, p);

  // Compare small survival probabilities directly in the upper tail.
  // Allow a few ulps of beta evaluation error at exact CDF jumps; otherwise
  // a value such as F(2)=0.6875 can spuriously select the next integer.
  const tolerance = 8 * Number.EPSILON * Math.min(u, 1 - u);
  const residual = (k: number): number =>
    u > 0.5
      ? 1 - u - bratio_complement(size, k + 1, p) + tolerance
      : negative_binomial_distribution(k, size, p) - u + tolerance;
  const zero = residual(0);
  if (!Number.isFinite(zero)) return Number.NaN;
  if (zero >= 0) return 0;

  const z = normal_quantile(u, 0, 1);
  const mean = (size * (1 - p)) / p;
  const sd = Math.sqrt(size * (1 - p)) / p;
  const guess = mean + sd * z + ((2 - p) * (z * z - 1)) / (6 * p);
  let lower = 0;
  let upper = Number.isFinite(guess) ? Math.max(1, Math.ceil(guess)) : 1;
  let value = residual(upper);
  for (let iteration = 0; value < 0 && iteration < 2048; iteration++) {
    lower = upper;
    if (upper === Number.MAX_VALUE) return Infinity;
    upper = upper >= Number.MAX_VALUE / 2 ? Number.MAX_VALUE : 2 * upper + 1;
    value = residual(upper);
  }
  if (!Number.isFinite(value)) return Number.NaN;
  if (value < 0) return Infinity;

  // Search integers rather than rounding a root of a discontinuous CDF.
  for (let iteration = 0; iteration < 2048 && upper - lower > 1; iteration++) {
    const middle = Math.floor(lower + (upper - lower) / 2);
    if (middle === lower || middle === upper) break;
    value = residual(middle);
    if (!Number.isFinite(value)) return Number.NaN;
    if (value >= 0) upper = middle;
    else lower = middle;
  }
  return upper;
}

/**
 * Negative binomial random failure counts from a gamma-Poisson mixture.
 * @param N nonnegative integer quantity of random numbers
 * @param size positive, finite success target (fractional values are allowed)
 * @param p success probability, 0 < p <= 1
 * @returns N failure counts, or N NaNs for invalid parameters
 * @throws RangeError if N is not a nonnegative integer
 * @category negative binomial
 */
export function negative_binomial_random(
  N: number,
  size: number,
  p: number,
): number[] {
  if (!Number.isInteger(N) || N < 0)
    throw new RangeError("N must be a nonnegative integer.");
  if (!valid_parameters(size, p)) return new Array(N).fill(Number.NaN);
  if (p === 1) return new Array(N).fill(0);
  if (size === 1) return geometric_random(N, p);
  // Beyond this shape, gamma fluctuations are far below floating-point
  // resolution; avoid overflowing the gamma generator's 2*size expression.
  const means =
    size > Number.MAX_VALUE / 2
      ? new Array(N).fill((size * (1 - p)) / p)
      : gamma_random(N, size, p / (1 - p));
  return means.map(lambda =>
    lambda === 0
      ? 0
      : lambda === Infinity
        ? Infinity
        : poisson_random(1, lambda)[0],
  );
}
