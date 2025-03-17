import { erf } from "./erf";
import { LOG_2PI } from "./constants";
import { phi_inv } from "./phi_inv";

/**
 * Normal(m, s) distribution's density function
 *
 * @param x argument to density function
 * @param m measure of center
 * @param s measure of scale
 * @returns normal density at x
 *
 * @category normal
 */
export function normal_density(x: number, m: number, s: number): number {
  return Math.exp(normal_log_density(x, m, s));
}

/**
 * Normal(m, s) distribution function
 *
 * @param x argument to distribution function
 * @param m measure of center
 * @param s measure of scale
 * @returns normal distribution at x
 *
 * @category normal
 */
export function normal_distribution(x: number, m: number, s: number): number {
  const z = (x - m) / s;
  return 0.5 * (1.0 + erf(z / Math.sqrt(2.0)));
}

/**
 * Normal(m, s) distribution's log-density function
 *
 * @param x argument to log-density function
 * @param m measure of center
 * @param s measure of scale
 * @returns normal log-density at x
 *
 * @category normal
 */
export function normal_log_density(x: number, m: number, s: number): number {
  const z = (x - m) / s;
  return -0.5 * LOG_2PI - Math.log(s) - 0.5 * z * z;
}

/**
 * Normal(m, s) distribution's quantile function
 *
 * @param p probability
 * @param m measure of center
 * @param s measure of scale
 * @returns normal quantile at p
 *
 * @category normal
 */
export function normal_quantile(p: number, m: number, s: number): number {
  return phi_inv(p) * s + m;
}

/**
 * Normal(m, s) random numbers
 *
 * @param N quantity of random numbers
 * @param m measure of center
 * @param s measure of scale
 * @returns N normal random numbers
 *
 * @category normal
 */
export function normal_random(N: number, m: number, s: number): number[] {
  if (s <= 0) {
    console.log("scale must be positive.");
  }

  let out: number[] = new Array(N).fill(0);
  for (let n = 0; n < N; n++) {
    out[n] = phi_inv(Math.random()) * s + m;
  }
  return out;
}
