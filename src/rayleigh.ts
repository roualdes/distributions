/**
 * Rayleigh density, with zero location and scale s.
 * @param x argument to the density function
 * @param s positive, finite scale parameter
 * @returns density at x, or NaN for invalid arguments
 * @category rayleigh
 */
export function rayleigh_density(x: number, s: number): number {
  return Math.exp(rayleigh_log_density(x, s));
}

/**
 * Rayleigh log-density, with zero location and scale s.
 * @param x argument to the log-density function
 * @param s positive, finite scale parameter
 * @returns log-density at x, or NaN for invalid arguments
 * @category rayleigh
 */
export function rayleigh_log_density(x: number, s: number): number {
  if (!Number.isFinite(s) || s <= 0 || Number.isNaN(x)) return Number.NaN;
  if (x <= 0 || x === Infinity) return -Infinity;
  const z = x / s;
  // Divide before squaring to avoid overflowing x^2 or s^2. Halving first
  // also retains finite log-densities when z^2 alone would overflow.
  return Math.log(x) - 2 * Math.log(s) - 0.5 * z * z;
}

/**
 * Rayleigh cumulative distribution P(X <= x), with zero location and scale s.
 * @param x argument to the distribution function
 * @param s positive, finite scale parameter
 * @returns cumulative probability at x, or NaN for invalid arguments
 * @category rayleigh
 */
export function rayleigh_distribution(x: number, s: number): number {
  if (!Number.isFinite(s) || s <= 0 || Number.isNaN(x)) return Number.NaN;
  if (x <= 0) return 0;
  if (x === Infinity) return 1;
  const z = x / s;
  return -Math.expm1(-(0.5 * z) * z);
}

/**
 * Rayleigh quantile, with zero location and scale s.
 * @param p probability in [0, 1]
 * @param s positive, finite scale parameter
 * @returns quantile at p; endpoints give 0 and Infinity, invalid arguments give NaN
 * @category rayleigh
 */
export function rayleigh_quantile(p: number, s: number): number {
  if (!Number.isFinite(s) || s <= 0 || !(p >= 0 && p <= 1)) return Number.NaN;
  if (p === 0) return 0;
  if (p === 1) return Infinity;
  return s * Math.sqrt(-2 * Math.log1p(-p));
}

/**
 * Rayleigh random numbers, with zero location and scale s.
 * @param N nonnegative integer quantity of random numbers
 * @param s positive, finite scale parameter
 * @returns N random numbers, or N NaNs for invalid scale
 * @throws RangeError if N is not a nonnegative integer
 * @category rayleigh
 */
export function rayleigh_random(N: number, s: number): number[] {
  if (!Number.isInteger(N) || N < 0) {
    throw new RangeError("N must be a nonnegative integer.");
  }
  if (!Number.isFinite(s) || s <= 0) return new Array(N).fill(Number.NaN);
  return Array.from(
    { length: N },
    () => s * Math.sqrt(-2 * Math.log1p(-Math.random())),
  );
}
