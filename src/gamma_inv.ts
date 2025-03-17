import { log_gamma } from "./log_gamma_fn";
import { normal_distribution } from "./normal";

// adapted by Edward A. Roualdes on 2024-06-17 from
// https://people.sc.fsu.edu/~jburkardt/c_src/asa239/asa239.c
export function gamma_inv(x: number, p: number): number {
  let a;
  let an;
  let arg;
  let b;
  let c;
  let elimit = -88.0;
  let oflo = 1.0e37;
  let plimit = 1000.0;
  let pn1;
  let pn2;
  let pn3;
  let pn4;
  let pn5;
  let pn6;
  let rn;
  let tol = 1.0e-14;
  let value = 0.0;
  let xbig = 1.0e8;

  if (x < 0.0) {
    console.log("Error: x must be non-negative.");
    return value;
  }

  if (p <= 0.0) {
    console.log("Error: p must be positive.");
    return value;
  }

  if (x == 0.0) {
    value = 0.0;
    return value;
  }

  // If P is large, use a normal approximation.
  if (p > plimit) {
    pn1 =
      3.0 * Math.sqrt(p) * (Math.pow(x / p, 1.0 / 3.0) + 1.0 / (9.0 * p) - 1.0);
    return normal_distribution(pn1, 0.0, 1.0);
  }

  // If X is large set value = 1.
  if (x > xbig) {
    value = 1.0;
    return value;
  }

  // Use Pearson's series expansion.
  if (x <= 1.0 || x < p) {
    arg = p * Math.log(x) - x - log_gamma(p + 1.0);
    c = 1.0;
    value = 1.0;
    a = p;

    while (1) {
      a = a + 1.0;
      c = (c * x) / a;
      value = value + c;

      if (c <= tol) {
        break;
      }
    }

    arg = arg + Math.log(value);

    if (elimit <= arg) {
      value = Math.exp(arg);
    } else {
      value = 0.0;
    }
  } else {
    // Use a continued fraction expansion.
    arg = p * Math.log(x) - x - log_gamma(p);
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    while (1) {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if (pn6 != 0.0) {
        rn = pn5 / pn6;
        const ttol = tol < tol * rn ? tol : tol * rn;
        if (Math.abs(value - rn) <= ttol) {
          break;
        }
        value = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;

      // Re-scale terms in continued fraction if terms are large.
      if (Math.abs(pn5) >= oflo) {
        pn1 /= oflo;
        pn2 /= oflo;
        pn3 /= oflo;
        pn4 /= oflo;
      }
    }

    arg += Math.log(value);

    if (arg >= elimit) {
      value = 1.0 - Math.exp(arg);
    } else {
      value = 1.0;
    }
  }

  return value;
}
