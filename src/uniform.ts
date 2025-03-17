/**
 * Continuous Unfirom(a, b) distribution's density function
 *
 * @param x - argument to density function
 * @param a - lower bound of support
 * @param b - upper bound of support
 * @returns uniform density at x
 *
 * @category uniform
 */
export function uniform_density(x: number, a: number, b: number): number {
  return Math.exp(uniform_log_density(x, a, b));
}

/**
 * Continuous Unfirom(a, b) distribution function
 *
 * @param x - argument to distribution function
 * @param a - lower bound of support
 * @param b - upper bound of support
 * @returns uniform distribution at x
 *
 * @category uniform
 */
export function uniform_distribution(x: number, a: number, b: number): number {
  if (b < a) {
    console.log("Error: a can't be larger than b.");
    return -Infinity;
  }
  if (x <= a) return 0.0;
  if (x >= b) return 1.0;
  return (x - a) / (b - a);
}

/**
 * Continuous Unfirom(a, b) distribution's log-density function
 *
 * @param x - argument to log-density function
 * @param a - lower bound of support
 * @param b - upper bound of support
 * @returns uniform log-density at x
 *
 * @category uniform
 */
export function uniform_log_density(x: number, a: number, b: number): number {
  if (b < a) {
    console.log("Error: a can't be bigger than b.");
    return -Infinity;
  }
  if (x < a) return -Infinity;
  if (x > b) return -Infinity;
  return -Math.log(b - a);
}

/**
 * Continuous Unfirom(a, b) distribution's quantilefunction
 *
 * @param p probability
 * @param a lower bound of support
 * @param b upper bound of support
 * @returns uniform quantile at p
 *
 * @category uniform
 */
export function uniform_quantile(p: number, a: number, b: number): number {
  if (p < 0.0 || p > 1.0) {
    console.log("Error: p must be between 0 and 1.");
    return -Infinity;
  }

  if (b < a) {
    console.log("Error: a can't be bigger than b.");
    return -Infinity;
  }

  return a + p * (b - a);
}

/**
 * Continuous Unfirom(a, b) random numbers
 *
 * @param N quantity of numbers to generate
 * @param a lower bound of support
 * @param b upper bound of support
 * @returns N uniform random numbers
 *
 * @category uniform
 */
export function uniform_random(N: number, a: number, b: number): number[] {
  if (b < a) {
    console.log("Error: a can't be bigger than b.");
  }

  const out: number[] = new Array(N).fill(0);
  for (let n = 0; n < N; n++) {
    const u = Math.random();
    out[n] = a + (b - a) * u;
  }

  return out;
}
