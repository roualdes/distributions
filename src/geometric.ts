/**
 * Geometric probability mass: failures before the first success.
 * @param x nonnegative integer failure count
 * @param p success probability, 0 < p <= 1
 * @returns probability mass at x, or NaN for invalid parameters
 * @category geometric
 */
export function geometric_density(x: number, p: number): number {
  return Math.exp(geometric_log_density(x, p));
}

/**
 * Geometric log-probability mass: failures before the first success.
 * @param x nonnegative integer failure count
 * @param p success probability, 0 < p <= 1
 * @returns log-probability mass at x, or NaN for invalid parameters
 * @category geometric
 */
export function geometric_log_density(x: number, p: number): number {
  if (!(p > 0 && p <= 1) || Number.isNaN(x)) return Number.NaN;
  if (!Number.isInteger(x) || x < 0) return -Infinity;
  if (p === 1) return x === 0 ? 0 : -Infinity;
  return Math.log(p) + x * Math.log1p(-p);
}

/**
 * Geometric cumulative probability of at most floor(x) failures.
 * @param x failure count (nonintegers are rounded down)
 * @param p success probability, 0 < p <= 1
 * @returns cumulative probability at x, or NaN for invalid parameters
 * @category geometric
 */
export function geometric_distribution(x: number, p: number): number {
  if (!(p > 0 && p <= 1) || Number.isNaN(x)) return Number.NaN;
  if (x < 0) return 0;
  if (x === Infinity || p === 1) return 1;
  return -Math.expm1((Math.floor(x) + 1) * Math.log1p(-p));
}

/**
 * Smallest failure count whose geometric CDF is at least u.
 * @param u cumulative probability in [0, 1]
 * @param p success probability, 0 < p <= 1
 * @returns failure count; u=1 gives Infinity unless p=1
 * @category geometric
 */
export function geometric_quantile(u: number, p: number): number {
  if (!(p > 0 && p <= 1) || !(u >= 0 && u <= 1)) return Number.NaN;
  if (u === 0 || p === 1) return 0;
  if (u === 1) return Infinity;
  const logFailure = Math.log1p(-p);
  let k = Math.max(0, Math.ceil(Math.log1p(-u) / logFailure) - 1);
  if (!Number.isFinite(k)) return Infinity;
  // Correct a possible rounding error at a CDF jump. Above the exact integer
  // range of doubles, adjacent counts cannot be distinguished.
  if (k < Number.MAX_SAFE_INTEGER) {
    const reached = (count: number): boolean =>
      u > 0.5
        ? (count + 1) * logFailure <= Math.log1p(-u)
        : geometric_distribution(count, p) >= u;
    if (k > 0 && reached(k - 1)) k--;
    if (!reached(k)) k++;
  }
  return k;
}

/**
 * Geometric random failure counts before the first success.
 * @param N nonnegative integer quantity of random numbers
 * @param p success probability, 0 < p <= 1
 * @returns N failure counts, or N NaNs for invalid p
 * @throws RangeError if N is not a nonnegative integer
 * @category geometric
 */
export function geometric_random(N: number, p: number): number[] {
  if (!Number.isInteger(N) || N < 0)
    throw new RangeError("N must be a nonnegative integer.");
  if (!(p > 0 && p <= 1)) return new Array(N).fill(Number.NaN);
  return Array.from({ length: N }, () => geometric_quantile(Math.random(), p));
}
