/**
 * Exponential(lambda) distribution's density function
 *
 * @param x argument to density funciton
 * @param lambda rate parameter
 * @returns exponential density at x
 *
 * @category exponential
 */
export function exponential_density(x: number, lambda: number): number {
  return Math.exp(exponential_log_density(x, lambda));
}

/**
 * Exponential(lambda) distribution function
 *
 * @param x argument to density funciton
 * @param lambda rate parameter
 * @returns exponential distribution at x
 *
 * @category exponential
 */
export function exponential_distribution(x: number, lambda: number): number {
  if (lambda <= 0) {
    console.log("lambda must be positive.");
    return -Infinity;
  }
  return 1 - Math.exp(-lambda * x);
}

/**
 * Exponential(lambda) distribution's log-density function
 *
 * @param x argument to log-density funciton
 * @param lambda rate parameter
 * @returns exponential log-density at x
 *
 * @category exponential
 */
export function exponential_log_density(x: number, lambda: number): number {
  if (lambda <= 0) {
    console.log("a and b must be greater than 0.");
    return -Infinity;
  }
  return Math.log(lambda) - lambda * x;
}

/**
 * Exponential(lambda) distribution's quantile function
 *
 * @param p probability
 * @param lambda rate parameter
 * @returns exponential quantile at p
 *
 * @category exponential
 */
export function exponential_quantile(p: number, lambda: number): number {
  if (p < 0 || p > 1) {
    console.log("p must be between 0 and 1");
    return -Infinity;
  }
  if (lambda <= 0) {
    console.log("lambda must be positive");
    return -Infinity;
  }
  return -Math.log1p(-p) / lambda;
}

/**
 * Exponential(lambda) random numbers
 *
 * @param N quantity of numbers to generate
 * @param lambda rate parameter
 * @returns N exponential random numbers
 *
 * @category exponential
 */
export function exponential_random(N: number, lambda: number): number[] {
  if (lambda <= 0) {
    console.log("lambda must be positive.");
  }

  let out: number[] = new Array(N).fill(0);
  for (let n = 0; n < N; n++) {
    out[n] = -Math.log(Math.random()) / lambda;
  }
  return out;
}
