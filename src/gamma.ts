import { PI } from "./constants";

import { log_gamma } from "./log_gamma_fn";
import { gamma_inv } from "./gamma_inv";
import { normal_quantile } from "./normal";

/**
 * Gamma(a, b) distribution's density function
 *
 * @param x argument to density function
 * @param a measure of shape
 * @param b measure of scale
 * @returns normal density at x
 *
 * @category gamma
 */
export function gamma_density(x: number, a: number, b: number): number {
  return Math.exp(gamma_log_density(x, a, b));
}

/**
 * Gamma(a, b) distribution function
 *
 * @param x argument to distribution function
 * @param a measure of shape
 * @param b measure of scale
 * @returns normal distribution at x
 *
 * @category gamma
 */
export function gamma_distribution(x: number, a: number, b: number): number {
  if (a <= 0.0 || b <= 0.0) {
    console.log("Error: a and b must be greater than 0.");
    return -Infinity;
  }
  return gamma_inv(b * x, a);
}

/**
 * Gamma(a, b) distribution's log-density function
 *
 * @param x argument to log-density function
 * @param a measure of shape
 * @param b measure of scale
 * @returns gamma log-density at x
 *
 * @category gamma
 */
export function gamma_log_density(x: number, a: number, b: number): number {
  if (a <= 0.0 || b <= 0.0) {
    console.log("Error: a and b must be greater than 0.");
    return -Infinity;
  }
  return a * Math.log(b) - log_gamma(a) + (a - 1.0) * Math.log(x) - b * x;
}

// adapted by Edward A. Roualdes on 2024-06-18 from
// https://people.sc.fsu.edu/~jburkardt/c_src/asa091/asa091.html
// authored by John Burkardt

/******************************************************************************/
/*
  Purpose:

  ppchi2() evaluates the percentage points of the Chi-squared PDF.

  Discussion

  Incorporates the suggested changes in AS R85 (vol.40(1),
  pages 233-5, 1991) which should eliminate the need for the limited
  range for P, though these limits have not been removed
  from the routine.

  Licensing:

  This code is distributed under the MIT license.

  Modified:

  05 June 2013

  Author:

  Original FORTRAN77 version by Donald Best, DE Roberts.
  C version by John Burkardt.

  Reference:

  Donald Best, DE Roberts,
  Algorithm AS 91:
  The Percentage Points of the Chi-Squared Distribution,
  Applied Statistics,
  Volume 24, Number 3, 1975, pages 385-390.

  Parameters:

  Input, double P,  value of the chi-squared cumulative
  probability density function.
  0.000002 <= P <= 0.999998.

  Input, double V, the parameter of the chi-squared probability
  density function.
  0 < V.

  Input, double G, the value of log ( Gamma ( V / 2 ) ).

  Output, int *IFAULT, is nonzero if an error occurred.
  0, no error.
  1, P is outside the legal range.
  2, V is not positive.
  3, an error occurred in GAMMAD.
  4, the result is probably as accurate as the machine will allow.

  Output, double PPCHI2, the value of the chi-squared random
  deviate with the property that the probability that a chi-squared random
  deviate with parameter V is less than or equal to PPCHI2 is P.
 */

function ppchi2(p: number, v: number): number {
  let a;
  let aa = 0.6931471806;
  let b;
  let c;
  let c1 = 0.01;
  let c2 = 0.222222;
  let c3 = 0.32;
  let c4 = 0.4;
  let c5 = 1.24;
  let c6 = 2.2;
  let c7 = 4.67;
  let c8 = 6.66;
  let c9 = 6.73;
  let c10 = 13.32;
  let c11 = 60.0;
  let c12 = 70.0;
  let c13 = 84.0;
  let c14 = 105.0;
  let c15 = 120.0;
  let c16 = 127.0;
  let c17 = 140.0;
  let c18 = 175.0;
  let c19 = 210.0;
  let c20 = 252.0;
  let c21 = 264.0;
  let c22 = 294.0;
  let c23 = 346.0;
  let c24 = 420.0;
  let c25 = 462.0;
  let c26 = 606.0;
  let c27 = 672.0;
  let c28 = 707.0;
  let c29 = 735.0;
  let c30 = 889.0;
  let c31 = 932.0;
  let c32 = 966.0;
  let c33 = 1141.0;
  let c34 = 1182.0;
  let c35 = 1278.0;
  let c36 = 1740.0;
  let c37 = 2520.0;
  let c38 = 5040.0;
  let ch;
  let e = 0.5e-6;
  let g = log_gamma(0.5 * v);
  let i;
  let maxit = 20;
  let pmax = 0.999998;
  let pmin = 0.000002;
  let p1;
  let p2;
  let q;
  let s1;
  let s2;
  let s3;
  let s4;
  let s5;
  let s6;
  let t;
  let value = -1.0;
  let x;
  let xx;

  if (p < pmin || pmax < p) {
    console.log("Error");
    return value;
  }

  if (v <= 0.0) {
    console.log("Error: v must be positive");
    return value;
  }

  xx = 0.5 * v;
  c = xx - 1.0;

  // Starting approximation for small chi-squared
  if (v < -c5 * Math.log(p)) {
    ch = Math.pow(p * xx * Math.exp(g + xx * aa), 1.0 / xx);

    if (ch < e) {
      value = ch;
      return value;
    }
  } else if (v <= c3) {
    // Starting approximation for V less than or equal to 0.32
    ch = c4;
    a = Math.log1p(-p);

    while (1) {
      q = ch;
      p1 = 1.0 + ch * (c7 + ch);
      p2 = ch * (c9 + ch * (c8 + ch));

      t = -0.5 + (c7 + 2.0 * ch) / p1 - (c9 + ch * (c10 + 3.0 * ch)) / p2;

      ch = ch - (1.0 - (Math.exp(a + g + 0.5 * ch + c * aa) * p2) / p1) / t;

      if (Math.abs(q / ch - 1.0) <= c1) {
        break;
      }
    }
  } else {
    x = normal_quantile(p, 0.0, 1.0);

    // Starting approximation using Wilson and Hilferty estimate
    p1 = c2 / v;
    ch = v * Math.pow(x * Math.sqrt(p1) + 1.0 - p1, 3.0);

    // Starting approximation for P tending to 1.
    if (c6 * v + 6.0 < ch) {
      ch = -2.0 * (Math.log1p(-p) - c * Math.log(0.5 * ch) + g);
    }
  }

  // Call to algorithm AS 239 and calculation of seven term Taylor series
  for (i = 1; i <= maxit; i++) {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p - gamma_inv(p1, xx);

    t = p2 * Math.exp(xx * aa + g + p1 - c * Math.log(ch));
    b = t / ch;
    a = 0.5 * t - b * c;
    s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 + c11 * a))))) / c24;
    s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 * a)))) / c37;
    s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37;
    s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 + c36 * a))) / c38;
    s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37;
    s6 = (c15 + c * (c23 + c16 * c)) / c38;
    ch +=
      t *
      (1.0 +
        0.5 * t * s1 -
        b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));

    if (e < Math.abs(q / ch - 1.0)) {
      value = ch;
      return value;
    }
  }

  value = ch;

  return value;
}

/**
 * Gamma(a, b) distribution's quantile function
 *
 * @param p probability
 * @param a measure of shape
 * @param b measure of rate
 * @returns gamma quantile at p
 *
 * @category gamma
 */
export function gamma_quantile(p: number, a: number, b: number): number {
  if (p < 0.0 || p > 1.0) {
    console.log("Error: p must be between 0 and 1.");
    return -Infinity;
  }

  if (a <= 0.0 || b <= 0.0) {
    console.log("Error: a and b must be positive.");
    return -Infinity;
  }

  return (ppchi2(p, 2 * a) * 0.5) / b;
}

// derived from Boost implementation
// https://www.boost.org/doc/libs/1_87_0/boost/random/gamma_distribution.hpp

// Boost Software License - Version 1.0 - August 17th, 2003
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

/**
 * Gamma(a, b) random numbers
 *
 * @param N quantity of random numbers
 * @param a measure of shape
 * @param b measure of rate
 * @returns N gamma random numbers
 *
 * @category gamma
 */
export function gamma_random(N: number, a: number, b: number): number[] {
  let out: number[] = new Array(N).fill(0);
  if (a < 0.0 || b <= 0.0) {
    console.log("Error: a and b must be positive.");
    return out;
  }

  if (a == 1.0) {
    let e = 0.0;
    for (let n = 0; n < N; ++n) {
      e = -Math.log(Math.random());
      out[n] = e / b;
    }
  } else if (a > 1.0) {
    for (let n = 0; n < N; ++n) {
      while (1) {
        let y = Math.tan(PI * Math.random());
        let x = Math.sqrt(2 * a - 1.0) * y + a - 1.0;

        if (x <= 0.0) {
          continue;
        }

        if (
          Math.random() >
          (1.0 + y * y) *
            Math.exp(
              (a - 1.0) * Math.log(x / (a - 1.0)) -
                Math.sqrt(2.0 * a - 1.0) * y,
            )
        ) {
          continue;
        }
        out[n] = x / b;
        break;
      }
    }
  } else {
    for (let n = 0; n < N; ++n) {
      let p = Math.exp(1.0) / (a + Math.exp(1.0));
      while (1) {
        let u = Math.random();
        let y = -Math.log(Math.random());
        let x = 0.0;
        let q = 0.0;
        if (u < p) {
          x = Math.exp(-y / a);
          q = p * Math.exp(-x);
        } else {
          x = 1.0 + y;
          q = p + (1.0 - p) * Math.pow(x, a - 1.0);
        }

        if (u >= q) {
          continue;
        }
        out[n] = x / b;
        break;
      }
    }
  }
  return out;
}
