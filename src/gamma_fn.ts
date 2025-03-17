import { log_gamma } from "./log_gamma_fn";
import { normal_distribution } from "./normal";

// lower regularized incomplete gamma function
export function gamma(x: number, p: number) {
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
  let value;
  let xbig = 1.0e8;

  value = 0.0;

  if (x <= 0.0) {
    return value;
  }

  if (p <= 0.0) {
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
        let ttol = tol < tol * rn ? tol : tol * rn;
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

// adapted from https://www.johndcook.com/blog/cpp_gamma/
// TODO regularized or not?
// export function gamma(x: number): number {
//   if (x <= 0.0) {
//     console.log("Error: x must be positive"); // TOOD better error handling
//   }

//   // Split the function domain into three intervals:
//   // (0, 0.001), [0.001, 12), and (12, infinity)

//   // First interval: (0, 0.001)
//   // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
//   // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
//   // The relative error over this interval is less than 6e-7.

//   if (x < 0.001) {
//     return 1.0 / (x * (1.0 + GAMMA * x));
//   }

//   // Second interval: [0.001, 12)
//   if (x < 12.0) {
//     // The algorithm directly approximates gamma over (1,2) and uses
//     // reduction identities to reduce other arguments to this interval.

//     let y = x;
//     let n = 0;
//     let arg_was_less_than_one = y < 1.0;

//     // Add or subtract integers as necessary to bring y into (1,2)
//     // Will correct for this below
//     if (arg_was_less_than_one) {
//       y += 1.0;
//     } else {
//       n = Math.floor(y) - 1; // will use n later
//       y -= n;
//     }

//     // numerator coefficients for approximation over the interval (1,2)
//     const p = [
//       -1.71618513886549492533811, 2.47656508055759199108314e1,
//       -3.79804256470945635097577e2, 6.29331155312818442661052e2,
//       8.66966202790413211295064e2, -3.14512729688483675254357e4,
//       -3.61444134186911729807069e4, 6.64561438202405440627855e4,
//     ];

//     // denominator coefficients for approximation over the interval (1,2)
//     const q = [
//       -3.08402300119738975254353e1, 3.15350626979604161529144e2,
//       -1.01515636749021914166146e3, -3.10777167157231109440444e3,
//       2.25381184209801510330112e4, 4.75584627752788110767815e3,
//       -1.34659959864969306392456e5, -1.15132259675553483497211e5,
//     ];

//     let num = 0.0;
//     let den = 1.0;
//     let z = y - 1;

//     for (let i = 0; i < 8; i++) {
//       num = (num + p[i]) * z;
//       den = den * z + q[i];
//     }
//     let result = num / den + 1.0;

//     // Apply correction if argument was not initially in (1,2)
//     if (arg_was_less_than_one) {
//       // Use identity gamma(z) = gamma(z+1)/z
//       // The variable "result" now holds gamma of the original y + 1
//       // Thus we use y-1 to get back the orginal y.
//       result /= y - 1.0;
//     } else {
//       // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
//       for (let i = 0; i < n; i++) result *= y++;
//     }

//     return result;
//   }

//   ///////////////////////////////////////////////////////////////////////////
//   // Third interval: [12, infinity)

//   if (x > 171.624) {
//     // Correct answer too large to display. Force +infinity.
//     return Infinity;
//   }

//   return Math.exp(log_gamma(x));
// }
