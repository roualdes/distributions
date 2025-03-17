import { GAMMA, LOG_SQRT_2PI } from "./constants";
import { erf } from "./erf";

// derived from a Claude translation of
// https://dl.acm.org/doi/10.1145/131766.131776

function ipmpar(i: number): number {
  const ipmpar = [2, 31, 2147483647, 2, 24, -125, 128, 53, -1021, 1024];
  return ipmpar[i - 1];
}

function exparg(l: number): number {
  const b = ipmpar(4);
  let lnb;
  if (b == 2) {
    lnb = 0.69314718055995;
  } else if (b == 8) {
    lnb = 2.0794415416798;
  } else if (b == 16) {
    lnb = 2.7725887222398;
  } else {
    lnb = Math.log(b);
  }

  let m;
  if (l == 0) {
    m = ipmpar(7);
  } else {
    m = ipmpar(6) - 1;
  }

  return 0.99999 * (m * lnb);
}

function fpser(a: number, b: number, x: number): number {
  // beta_distribution(0.3, 15.0, 1e-16) partially tests this fn

  // b < min(eps, eps * a) && x <= 0.5
  let an, t, s, c, tol;
  let fpser = 1;

  if (a > 1e-3 * Number.EPSILON) {
    fpser = 0;
    t = a * Math.log(x);
    if (t < exparg(1)) {
      return 0;
    }
    fpser = Math.exp(t);
  }

  /* Note that 1/B(a,b) = b */
  fpser = (b / a) * fpser;
  tol = Number.EPSILON / a;
  an = a + 1;
  t = x;
  s = t / an;

  do {
    an += 1;
    t = x * t;
    c = t / an;
    s += c;
  } while (Math.abs(c) > tol);

  fpser *= 1 + a * s;
  return fpser;
}

function apser(a: number, b: number, x: number): number {
  // beta_distribution(1e-16, 1e-16, 1) partially tests this fn

  const eps = Number.EPSILON;
  const g = GAMMA;
  let bx, t, c, tol, j, s, aj;

  bx = b * x;
  t = x - bx;

  if (b * eps <= 2e-2) {
    c = Math.log(x) + psi(b) + g + t;
  } else {
    c = Math.log(bx) + g + t;
  }

  tol = 5.0 * eps * Math.abs(c);
  j = 1.0;
  s = 0.0;

  do {
    j += 1;
    t *= x - bx / j;
    aj = t / j;
    s += aj;
  } while (Math.abs(aj) > tol);

  return -a * (c + s);
}

function psi(xx: number): number {
  /* Constants */
  const piov4 = 0.785398163397448;
  const dx0 = 1.461632144968362341262659542325721325;

  /* Coefficients for rational approximation of
  psi(x) / (x - x0),  0.5 <= x <= 3.0 */
  const p1 = [
    0.89538502298197e-2, 0.477762828042627e1, 0.142441585084029e3,
    0.118645200713425e4, 0.363351846806499e4, 0.413810161269013e4,
    0.130560269827897e4,
  ];

  const q1 = [
    0.448452573429826e2, 0.520752771467162e3, 0.22100079924783e4,
    0.364127349079381e4, 0.1908310765963e4, 0.691091682714533e-5,
  ];

  /* Coefficients for rational approximation of
  psi(x) - ln(x) + 1/(2*x),  x > 3.0 */
  const p2 = [
    -0.212940445131011e1, -0.701677227766759e1, -0.448616543918019e1,
    -0.648157123766197,
  ];

  const q2 = [
    0.322703493791143e2, 0.892920700481861e2, 0.546117738103215e2,
    0.777788548522962e1,
  ];

  let x = xx;
  let aug = 0.0;

  /* Machine dependent constants */
  let xmax1 = Math.min(ipmpar(3), 1 / Number.EPSILON);
  const xsmall = 1e-9;

  /* Handle case where x < 0.5 */
  if (x < 0.5) {
    /* Check for error conditions first */
    if (x == 0.0) {
      return 0.0;
    }

    if (Math.abs(x) <= xsmall) {
      /* For very small |x|, use 1/x as substitute for pi*cotan(pi*x) */
      aug = -1.0 / x;
    } else {
      /* Reduction of argument for cotan */
      let w = -x;
      let sgn = piov4;

      if (w <= 0.0) {
        w = -w;
        sgn = -sgn;
      }

      /* Check for error condition */
      if (w >= xmax1) {
        return 0.0;
      }

      /* Compute fractional part */
      let nq = w;
      w = w - nq;
      nq = w * 4.0;
      w = 4.0 * (w - nq * 0.25);

      /* Adjust argument and determine sign */
      let n = nq / 2;
      if (n + n != nq) {
        w = 1.0 - w;
      }

      let z = piov4 * w;
      let m = n / 2;
      if (m + m != n) {
        sgn = -sgn;
      }

      n = (nq + 1) / 2;
      m = n / 2;
      m = m + m;

      if (m == n) {
        /* Check for singularity */
        if (z == 0.0) {
          return 0.0;
        }
        /* Use cos/sin as substitute for cotan */
        aug = sgn * ((Math.cos(z) / Math.sin(z)) * 4.0);
      } else {
        /* Use sin/cos as substitute for tan */
        aug = sgn * ((Math.sin(z) / Math.cos(z)) * 4.0);
      }
    }
    x = 1.0 - x;
  }

  /* Now handle x >= 0.5 cases */
  if (x <= 3.0) {
    /* 0.5 <= x <= 3.0 */
    let den = x;
    let upper = p1[0] * x;

    for (let i = 0; i < 5; i++) {
      den = (den + q1[i]) * x;
      upper = (upper + p1[i + 1]) * x;
    }

    den = (upper + p1[6]) / (den + q1[5]);
    let xmx0 = x - dx0;
    return den * xmx0 + aug;
  } else {
    /* x > 3.0 */
    if (x >= xmax1) {
      return aug + Math.log(x);
    }

    /* 3.0 < x < xmax1 */
    let w = 1.0 / (x * x);
    let den = w;
    let upper = p2[0] * w;

    for (let i = 0; i < 3; i++) {
      den = (den + q2[i]) * w;
      upper = (upper + p2[i + 1]) * w;
    }

    aug = upper / (den + q2[3]) - 0.5 / x + aug;
    return aug + Math.log(x);
  }
}

/**
 * @brief Power series expansion for evaluating Ix(a,b)
 *
 * Evaluates the incomplete beta function Ix(a,b) when b <= 1 or b*x <= 0.7
 *
 * @param a First parameter of incomplete beta function
 * @param b Second parameter of incomplete beta function
 * @param x Upper limit of integration
 * @param eps Tolerance used for computation
 * @return Computed value of the incomplete beta function
 */
function bpser(a: number, b: number, x: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) partially tests this fn

  const eps = Number.EPSILON;
  let a0, b0, z, apb, u, c, t, w, n, sum, tol;
  let m, i, result;

  /* Initial check */
  if (x == 0.0) {
    return 0.0;
  }

  /* Compute the factor x^a/(a*beta(a,b)) */
  a0 = Math.min(a, b);
  if (a0 >= 1.0) {
    z = a * Math.log(x) - betaln(a, b);
    result = Math.exp(z) / a;
    return result;
  }

  b0 = Math.max(a, b);
  if (b0 >= 8.0) {
    /* Procedure for a0 < 1 and b0 >= 8 */
    u = gamln1(a0) + algdiv(a0, b0);
    z = a * Math.log(x) - u;
    result = (a0 / a) * Math.exp(z);

    if (result == 0.0 || a <= 0.1 * eps) {
      return result;
    }

    /* Compute the series */
    sum = 0.0;
    n = 0.0;
    c = 1.0;
    tol = eps / a;

    do {
      n = n + 1.0;
      c = c * (0.5 + (0.5 - b / n)) * x;
      w = c / (a + n);
      sum = sum + w;
    } while (Math.abs(w) > tol);

    return result * (1.0 + a * sum);
  }

  if (b0 > 1.0) {
    /* Procedure for a0 < 1 and 1 < b0 < 8 */
    u = gamln1(a0);
    m = b0 - 1.0;

    if (m >= 1) {
      c = 1.0;
      for (i = 1; i <= m; ++i) {
        b0 = b0 - 1.0;
        c = c * (b0 / (a0 + b0));
      }
      u = Math.log(c) + u;
    }

    z = a * Math.log(x) - u;
    b0 = b0 - 1.0;
    apb = a0 + b0;

    if (apb > 1.0) {
      u = a0 + b0 - 1.0;
      t = (1.0 + gam1(u)) / apb;
    } else {
      t = 1.0 + gam1(apb);
    }

    result = (Math.exp(z) * (a0 / a) * (1.0 + gam1(b0))) / t;

    if (result == 0.0 || a <= 0.1 * eps) {
      return result;
    }

    /* Compute the series */
    sum = 0.0;
    n = 0.0;
    c = 1.0;
    tol = eps / a;

    do {
      n = n + 1.0;
      c = c * (0.5 + (0.5 - b / n)) * x;
      w = c / (a + n);
      sum = sum + w;
    } while (Math.abs(w) > tol);

    return result * (1.0 + a * sum);
  }

  /* Procedure for a0 < 1 and b0 <= 1 */
  result = Math.pow(x, a);
  if (result == 0.0) {
    return 0.0;
  }

  apb = a + b;
  if (apb > 1.0) {
    u = a + b - 1.0;
    z = (1.0 + gam1(u)) / apb;
  } else {
    z = 1.0 + gam1(apb);
  }

  c = ((1.0 + gam1(a)) * (1.0 + gam1(b))) / z;
  result = result * c * (b / apb);

  if (result == 0.0 || a <= 0.1 * eps) {
    return result;
  }

  /* Compute the series */
  sum = 0.0;
  n = 0.0;
  c = 1.0;
  tol = eps / a;

  do {
    n = n + 1.0;
    c = c * (0.5 + (0.5 - b / n)) * x;
    w = c / (a + n);
    sum = sum + w;
  } while (Math.abs(w) > tol);

  return result * (1.0 + a * sum);
}

/**
 * @brief Evaluates Ix(a,b) - Ix(a+n,b)
 *
 * Computes the difference between two incomplete beta functions
 * where n is a positive integer.
 *
 * @param a First parameter of incomplete beta function
 * @param b Second parameter of incomplete beta function
 * @param x First upper limit of integration
 * @param y Complement of x (1-x)
 * @param n Positive integer increment for second term
 * @param eps Tolerance used for computation
 * @return Computed value of the difference Ix(a,b) - Ix(a+n,b)
 */
function bup(a: number, b: number, x: number, y: number, n: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) partially tests this fn

  const eps = Number.EPSILON;
  let apb, ap1, mu, d, t, bup_val, r, w, l;
  let k, nm1, kp1, i;

  /* Obtain the scaling factor exp(-mu) and
  exp(mu)*(x^a*y^b/beta(a,b))/a */
  apb = a + b;
  ap1 = a + 1.0;
  mu = 0;
  d = 1.0;

  if (n != 1 && a >= 1.0) {
    if (apb >= 1.1 * ap1) {
      mu = Math.abs(exparg(1));
      k = exparg(0);
      if (k < mu) mu = k;
      t = mu;
      d = Math.exp(-t);
    }
  }

  bup_val = brcmp1(mu, a, b, x, y) / a;

  if (n == 1 || bup_val == 0.0) {
    return bup_val;
  }

  nm1 = n - 1;
  w = d;

  /* Let k be the index of the maximum term */
  k = 0;
  if (b > 1.0) {
    if (y > 1.0e-4) {
      r = ((b - 1.0) * x) / y - a;
      if (r >= 1.0) {
        k = nm1;
        t = nm1;
        if (r < t) k = r;
      }
    } else {
      k = nm1;
    }

    /* Add the increasing terms of the series */
    for (i = 1; i <= k; i++) {
      l = i - 1;
      d = ((apb + l) / (ap1 + l)) * x * d;
      w = w + d;
    }

    if (k == nm1) {
      return bup_val * w;
    }
  }

  /* Add the remaining terms of the series */
  kp1 = k + 1;
  for (i = kp1; i <= nm1; i++) {
    l = i - 1;
    d = ((apb + l) / (ap1 + l)) * x * d;
    w = w + d;
    if (d <= eps * w) {
      break;
    }
  }

  /* Terminate the procedure */
  return bup_val * w;
}

/**
 * Continued fraction expansion for IX(A,B) when A,B > 1.
 * It is assumed that LAMBDA = (A + B)*Y - B.
 *
 * @param {number} a - First parameter (assumed > 1)
 * @param {number} b - Second parameter (assumed > 1)
 * @param {number} x - X value
 * @param {number} y - Y value
 * @param {number} lambda - Precomputed value (A + B)*Y - B
 * @param {number} eps - Epsilon for convergence testing
 * @returns {number} - The result of the continued fraction expansion
 */
function bfrac(
  a: number,
  b: number,
  x: number,
  y: number,
  lambda: number,
  eps: number,
): number {
  // binomial_distribution(40, 130.0, 0.3) partially tests this fn

  let result = brcomp(a, b, x, y);
  if (result === 0.0) return 0.0;

  let c = 1.0 + lambda;
  let c0 = b / a;
  let c1 = 1.0 + 1.0 / a;
  let yp1 = y + 1.0;

  let n = 0.0;
  let p = 1.0;
  let s = a + 1.0;
  let an = 0.0;
  let bn = 1.0;
  let anp1 = 1.0;
  let bnp1 = c / c1;
  let r = c1 / c;

  // Continued fraction calculation
  while (true) {
    n = n + 1.0;
    let t = n / a;
    let w = n * (b - n) * x;
    let e = a / s;
    let alpha = p * (p + c0) * e * e * (w * x);
    e = (1.0 + t) / (c1 + t + t);
    let beta = n + w / s + e * (c + n * yp1);
    p = 1.0 + t;
    s = s + 2.0;

    t = alpha * an + beta * anp1;
    an = anp1;
    anp1 = t;
    t = alpha * bn + beta * bnp1;
    bn = bnp1;
    bnp1 = t;

    let r0 = r;
    r = anp1 / bnp1;
    if (Math.abs(r - r0) <= eps * r) break;

    an = an / bnp1;
    bn = bn / bnp1;
    anp1 = r;
    bnp1 = 1.0;
  }

  return result * r;
}

/**
 * Evaluation of X^A*Y^B/Beta(A,B)
 *
 * @param {number} a - First parameter
 * @param {number} b - Second parameter
 * @param {number} x - X value
 * @param {number} y - Y value
 * @returns {number} - The result of X^A*Y^B/Beta(A,B)
 */
function brcomp(a: number, b: number, x: number, y: number): number {
  // binomial_distribution(40, 130.0, 0.3) partially tests this fn

  // Constant = 1/SQRT(2*PI)
  const CONST = 0.398942280401433;

  if (x === 0.0 || y === 0.0) return 0.0;

  let a0 = Math.min(a, b);
  if (a0 >= 8.0) {
    return algorithmForLargeAB();
  }

  let lnx, lny;
  if (x <= 0.375) {
    lnx = Math.log(x);
    lny = Math.log1p(-x); // alnrel(-x);
  } else if (y <= 0.375) {
    lnx = Math.log1p(-y); // alnrel(-y);
    lny = Math.log(y);
  } else {
    lnx = Math.log(x);
    lny = Math.log(y);
  }

  let z = a * lnx + b * lny;

  if (a0 >= 1.0) {
    z = z - betaln(a, b);
    return Math.exp(z);
  }

  // Procedure for a < 1 or b < 1
  let b0 = Math.max(a, b);

  if (b0 >= 8.0) {
    return algorithmForB0GE8();
  }

  if (b0 > 1.0) {
    return algorithmForB0Between1And8();
  }

  // Algorithm for b0 <= 1
  let result = Math.exp(z);
  if (result === 0.0) return 0.0;

  let apb = a + b;
  let zValue;

  if (apb > 1.0) {
    let u = a + b - 1.0;
    zValue = (1.0 + gam1(u)) / apb;
  } else {
    zValue = 1.0 + gam1(apb);
  }

  let c = ((1.0 + gam1(a)) * (1.0 + gam1(b))) / zValue;
  result = (result * (a0 * c)) / (1.0 + a0 / b0);

  return result;

  // Algorithm for 1 < b0 < 8
  function algorithmForB0Between1And8() {
    let u = gamln1(a0);
    let n = Math.floor(b0 - 1.0);
    let c = 1.0;

    let tempB0 = b0;
    for (let i = 1; i <= n; i++) {
      tempB0 = tempB0 - 1.0;
      c = c * (tempB0 / (a0 + tempB0));
    }

    u = Math.log(c) + u;
    z = z - u;
    tempB0 = b0 - 1.0;
    let apb = a0 + tempB0;
    let t;

    if (apb > 1.0) {
      let u = a0 + tempB0 - 1.0;
      t = (1.0 + gam1(u)) / apb;
    } else {
      t = 1.0 + gam1(apb);
    }

    return (a0 * Math.exp(z) * (1.0 + gam1(tempB0))) / t;
  }

  // Algorithm for b0 >= 8
  function algorithmForB0GE8() {
    let u = gamln1(a0) + algdiv(a0, b0);
    return a0 * Math.exp(z - u);
  }

  // Procedure for a >= 8 and b >= 8
  function algorithmForLargeAB() {
    let h, x0, y0, lambda;

    if (a > b) {
      h = b / a;
      x0 = 1.0 / (1.0 + h);
      y0 = h / (1.0 + h);
      lambda = (a + b) * y - b;
    } else {
      h = a / b;
      x0 = h / (1.0 + h);
      y0 = 1.0 / (1.0 + h);
      lambda = a - (a + b) * x;
    }

    let e = -lambda / a;
    let u;

    if (Math.abs(e) <= 0.6) {
      u = rlog1(e);
    } else {
      u = e - Math.log(x / x0);
    }

    e = lambda / b;
    let v;

    if (Math.abs(e) <= 0.6) {
      v = rlog1(e);
    } else {
      v = e - Math.log(y / y0);
    }

    let z = Math.exp(-(a * u + b * v));
    return CONST * Math.sqrt(b * x0) * z * Math.exp(-bcorr(a, b));
  }
}

/**
 * Evaluation of the function x - ln(1 + x)
 *
 * @param {number} x - Input value
 * @returns {number} - Result of x - ln(1 + x)
 */
function rlog1(x: number): number {
  // binomial_distribution(40, 130.0, 0.3) partially tests this fn

  const A = 0.0566749439387324;
  const B = 0.0456512608815524;

  const P0 = 0.333333333333333;
  const P1 = -0.224696413112536;
  const P2 = 0.00620886815375787;

  const Q1 = -1.27408923933623;
  const Q2 = 0.354508718369557;

  if (x < -0.39 || x > 0.57) {
    // Direct calculation for values outside the series expansion range
    const w = x + 0.5 + 0.5;
    return x - Math.log(w);
  }

  let h,
    w1 = 0.0;

  // Argument reduction based on the value of x
  if (x < -0.18) {
    h = x + 0.3;
    h = h / 0.7;
    w1 = A - h * 0.3;
  } else if (x > 0.18) {
    h = 0.75 * x - 0.25;
    w1 = B + h / 3.0;
  } else {
    h = x;
    w1 = 0.0;
  }

  // Series expansion
  const r = h / (h + 2.0);
  const t = r * r;
  const w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);

  return 2.0 * t * (1.0 / (1.0 - r) - r * w) + w1;
}

/**
 * Evaluation of the logarithm of the Beta function
 * @param {number} a0 - First parameter
 * @param {number} b0 - Second parameter
 * @returns {number} Logarithm of the Beta function
 */
function betaln(a0: number, b0: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) partially tests this fn

  const e = LOG_SQRT_2PI;

  let a = Math.min(a0, b0);
  let b = Math.max(a0, b0);

  // Case: a >= 8.0
  if (a >= 8.0) {
    return largeProcedure(a, b, e);
  }

  // Case: 1.0 <= a < 8.0
  if (a >= 1.0) {
    // Subcase: 1 <= a <= 2
    if (a <= 2.0) {
      if (b <= 2.0) {
        return gamln(a) + gamln(b) - gsumln(a, b);
      }

      if (b >= 8.0) {
        return gamln(a) + algdiv(a, b);
      }
    }

    // Subcase: 2 < a < 8
    // Reduction of a when b <= 1000
    if (b <= 1000.0) {
      const n = Math.floor(a - 1.0);
      let w = 1.0;

      for (let i = 1; i <= n; i++) {
        a = a - 1.0;
        const h = a / b;
        w = w * (h / (1.0 + h));
      }

      w = Math.log(w);

      if (b >= 8.0) {
        return w + gamln(a) + algdiv(a, b);
      }

      // Reduction of b when b < 8
      const n2 = Math.floor(b - 1.0);
      let z = 1.0;

      for (let i = 1; i <= n2; i++) {
        b = b - 1.0;
        z = z * (b / (a + b));
      }

      return w + Math.log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
    } else {
      // Reduction of a when b > 1000
      const n = Math.floor(a - 1.0);
      let w = 1.0;

      for (let i = 1; i <= n; i++) {
        a = a - 1.0;
        w = w * (a / (1.0 + a / b));
      }

      return Math.log(w) - n * Math.log(b) + (gamln(a) + algdiv(a, b));
    }
  }

  // Case: a < 1
  if (b >= 8.0) {
    return gamln(a) + algdiv(a, b);
  } else {
    return gamln(a) + (gamln(b) - gamln(a + b));
  }
}

/**
 * Helper function for large a and b case
 * @param {number} a
 * @param {number} b
 * @param {number} e
 * @returns {number}
 */
function largeProcedure(a: number, b: number, e: number): number {
  const w = bcorr(a, b);
  const h = a / b;
  const c = h / (1.0 + h);
  const u = -(a - 0.5) * Math.log(c);
  const v = b * Math.log1p(h); // alnrel(h);

  if (u <= v) {
    return -0.5 * Math.log(b) + e + w - u - v;
  } else {
    return -0.5 * Math.log(b) + e + w - v - u;
  }
}

/**
 * Evaluation of ln(Gamma(a)) for positive a
 * Originally written by Alfred H. Morris, Naval Surface Warfare Center, Dahlgren, Virginia
 * @param {number} a - Input value (must be positive)
 * @returns {number} Natural logarithm of Gamma(a)
 */
function gamln(a: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) partially tests this fn

  // Constant: 0.5*(ln(2*PI) - 1)
  const d = 0.418938533204673;

  // Constants for the series expansion
  const c0 = 0.0833333333333333;
  const c1 = -0.00277777777760991;
  const c2 = 0.00079365066682539;
  const c3 = -0.00059520293135187;
  const c4 = 0.000837308034031215;
  const c5 = -0.00165322962780713;

  // Case 1: a <= 0.8
  if (a <= 0.8) {
    return gamln1(a) - Math.log(a);
  }

  // Case 2: 0.8 < a <= 2.25
  if (a <= 2.25) {
    const t = a - 0.5 - 0.5;
    return gamln1(t);
  }

  // Case 3: 2.25 < a < 10.0
  if (a < 10.0) {
    const n = Math.floor(a - 1.25);
    let t = a;
    let w = 1.0;

    for (let i = 1; i <= n; i++) {
      t = t - 1.0;
      w = t * w;
    }

    return gamln1(t - 1.0) + Math.log(w);
  }

  // Case 4: a >= 10.0
  const t = (1.0 / a) ** 2;
  const w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;

  return d + w + (a - 0.5) * (Math.log(a) - 1.0);
}

/**
 * Evaluation of ln(Gamma(1 + a)) for -0.2 <= a <= 1.25
 * @param {number} a - Input value
 * @returns {number} Natural logarithm of Gamma(1 + a)
 */
function gamln1(a: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) partially tests this fn

  const p0 = 0.577215664901533;
  const p1 = 0.844203922187225;
  const p2 = -0.168860593646662;
  const p3 = -0.780427615533591;
  const p4 = -0.402055799310489;
  const p5 = -0.0673562214325671;
  const p6 = -0.00271935708322958;

  const q1 = 2.88743195473681;
  const q2 = 3.12755088914843;
  const q3 = 1.56875193295039;
  const q4 = 0.361951990101499;
  const q5 = 0.0325038868253937;
  const q6 = 0.000667465618796164;

  // Constants for polynomial approximation (second set)
  const r0 = 0.422784335098467;
  const r1 = 0.848044614534529;
  const r2 = 0.565221050691933;
  const r3 = 0.156513060486551;
  const r4 = 0.017050248402265;
  const r5 = 0.000497958207639485;

  const s1 = 1.24313399877507;
  const s2 = 0.548042109832463;
  const s3 = 0.10155218743983;
  const s4 = 0.00713309612391;
  const s5 = 0.000116165475989616;

  // For a < 0.6, use first set of coefficients
  if (a < 0.6) {
    const w =
      ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a + p0) /
      ((((((q6 * a + q5) * a + q4) * a + q3) * a + q2) * a + q1) * a + 1.0);
    return -a * w;
  }
  // For a >= 0.6, use second set of coefficients
  else {
    const x = a - 0.5 - 0.5;
    const w =
      (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
      (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0);
    return x * w;
  }
}

/**
 * Computation of 1/Gamma(a+1) - 1 for -0.5 <= a <= 1.5
 * @param {number} a - Input value
 * @returns {number} 1/Gamma(a+1) - 1
 */
function gam1(a: number): number {
  // beta_distribution(0.33, 1.5, 7.24) partially tests this fun

  const p = [
    0.577215664901533, -0.409078193005776, -0.230975380857675,
    0.0597275330452234, 0.0076696818164949, -0.00514889771323592,
    0.000589597428611429,
  ];

  const q = [
    1.0, 0.427569613095214, 0.158451672430138, 0.0261132021441447,
    0.00423244297896961,
  ];

  const r = [
    -0.422784335098468, -0.771330383816272, -0.244757765222226,
    0.118378989872749, 0.000930357293360349, -0.0118290993445146,
    0.00223047661158249, 0.000266505979058923, -0.000132674909766242,
  ];

  const s1 = 0.273076135303957;
  const s2 = 0.0559398236957378;

  let t = a;
  const d = a - 0.5;

  if (d > 0.0) {
    t = d - 0.5;
  }

  if (t === 0.0) {
    return 0.0;
  }

  if (t > 0.0) {
    const top =
      (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]) *
        t +
      p[0];
    const bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + q[0];
    const w = top / bot;

    if (d <= 0.0) {
      return a * w;
    } else {
      return (t / a) * (w - 0.5 - 0.5);
    }
  }

  const top =
    (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]) * t + r[3]) *
      t +
      r[2]) *
      t +
      r[1]) *
      t +
    r[0];
  const bot = (s2 * t + s1) * t + 1.0;
  const w = top / bot;

  if (d <= 0.0) {
    return a * (w + 0.5 + 0.5);
  } else {
    return (t * w) / a;
  }
}

/**
 * Computation of ln(Gamma(b)/Gamma(a+b)) when b >= 8
 *
 * In this algorithm, del(x) is the function defined by
 * ln(Gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*PI) + del(x)
 *
 * @param {number} a - First parameter
 * @param {number} b - Second parameter (must be >= 8)
 * @returns {number} ln(Gamma(b)/Gamma(a+b))
 */
function algdiv(a: number, b: number): number {
  // binomial_distribution(10, 11.0, 0.9) partially tests this fn

  const c0 = 0.0833333333333333;
  const c1 = -0.00277777777760991;
  const c2 = 0.00079365066682539;
  const c3 = -0.00059520293135187;
  const c4 = 0.000837308034031215;
  const c5 = -0.00165322962780713;

  let h, c, x, d;

  if (a > b) {
    h = b / a;
    c = 1.0 / (1.0 + h);
    x = h / (1.0 + h);
    d = a + (b - 0.5);
  } else {
    h = a / b;
    c = h / (1.0 + h);
    x = 1.0 / (1.0 + h);
    d = b + (a - 0.5);
  }

  // Set sn = (1 - x^n)/(1 - x)
  const x2 = x * x;
  const s3 = 1.0 + (x + x2);
  const s5 = 1.0 + (x + x2 * s3);
  const s7 = 1.0 + (x + x2 * s5);
  const s9 = 1.0 + (x + x2 * s7);
  const s11 = 1.0 + (x + x2 * s9);

  const t = (1.0 / b) ** 2;
  const w =
    ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) *
      t +
    c0;
  const w_adjusted = w * (c / b);

  const u = d * Math.log1p(a / b); // alnrel(a / b);
  const v = a * (Math.log(b) - 1.0);

  if (u <= v) {
    return w_adjusted - u - v;
  } else {
    return w_adjusted - v - u;
  }
}

/**
 * Evaluation of del(a0) + del(b0) - del(a0 + b0) where
 * ln(Gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*PI) + del(a).
 * It is assumed that a0 >= 8 and b0 >= 8.
 *
 * @param {number} a0 - First parameter (must be >= 8)
 * @param {number} b0 - Second parameter (must be >= 8)
 * @returns {number} del(a0) + del(b0) - del(a0 + b0)
 */
function bcorr(a0: number, b0: number): number {
  // binomial_distribution(40, 130.0, 0.3) partially tests this fn

  const c0 = 0.0833333333333333;
  const c1 = -0.00277777777760991;
  const c2 = 0.00079365066682539;
  const c3 = -0.00059520293135187;
  const c4 = 0.000837308034031215;
  const c5 = -0.00165322962780713;

  // Ensure a <= b
  const a = Math.min(a0, b0);
  const b = Math.max(a0, b0);

  const h = a / b;
  const c = h / (1.0 + h);
  const x = 1.0 / (1.0 + h);
  const x2 = x * x;

  // Set sn = (1 - x^n)/(1 - x)
  const s3 = 1.0 + (x + x2);
  const s5 = 1.0 + (x + x2 * s3);
  const s7 = 1.0 + (x + x2 * s5);
  const s9 = 1.0 + (x + x2 * s7);
  const s11 = 1.0 + (x + x2 * s9);

  // Set w = del(b) - del(a + b)
  const t = (1.0 / b) ** 2;
  const w =
    ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) *
      t +
    c0;
  const w_adjusted = w * (c / b);

  // Compute del(a) + w
  const t_a = (1.0 / a) ** 2;
  const result =
    (((((c5 * t_a + c4) * t_a + c3) * t_a + c2) * t_a + c1) * t_a + c0) / a +
    w_adjusted;

  return result;
}

/**
 * Asymptotic expansion for the incomplete beta function IX(A,B) when A is larger than B.
 * This function calculates an asymptotic expansion and adds the result to the input value w.
 * It assumes that A >= 15 and B <= 1.
 *
 * @param {number} a - First parameter of the incomplete beta function, must be >= 15
 * @param {number} b - Second parameter of the incomplete beta function, must be <= 1
 * @param {number} x - Value between 0 and 1
 * @param {number} y - Equal to 1-x
 * @param {number} w - Value to which the result of the expansion is added
 * @param {number} eps - Tolerance used for convergence
 * @returns {number} The updated w value after adding the expansion result
 */
function bgrat(
  a: number,
  b: number,
  x: number,
  y: number,
  w: number,
  eps: number,
): number {
  // binomial_distribution(10, 11.0, 0.9) partially tests this fn

  const c = new Array(30).fill(0);
  const d = new Array(30).fill(0);

  const bm1 = b - 0.5 - 0.5;
  const nu = a + 0.5 * bm1;

  let lnx;
  if (y > 0.375) {
    lnx = Math.log(x);
  } else {
    lnx = Math.log1p(-y); // alnrel(-y);
  }

  const z = -nu * lnx;

  if (b * z === 0.0) {
    return w;
  }

  // Computation of the expansion
  // Set r = exp(-z)*z^b/gamma(b)
  let r = b * (1.0 + gam1(b)) * Math.exp(b * Math.log(z));
  r = r * Math.exp(a * lnx) * Math.exp(0.5 * bm1 * lnx);

  const u = algdiv(b, a) + b * Math.log(nu);
  const uExp = r * Math.exp(-u);

  if (uExp === 0.0) {
    return w;
  }

  const q = grat1(b, z, r, eps);

  const v = 0.25 * (1.0 / nu) ** 2;
  const t2 = 0.25 * lnx * lnx;
  const l = w / uExp;
  let j = q / r;
  let sum = j;
  let t = 1.0;
  let cn = 1.0;
  let n2 = 0.0;

  for (let n = 1; n <= 30; n++) {
    const bp2n = b + n2;
    j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
    n2 = n2 + 2.0;
    t = t * t2;
    cn = cn / (n2 * (n2 + 1.0));
    c[n - 1] = cn;

    let s = 0.0;
    if (n > 1) {
      const nm1 = n - 1;
      let coef = b - n;

      for (let i = 1; i <= nm1; i++) {
        s = s + coef * c[i - 1] * d[n - i - 1];
        coef = coef + b;
      }
    }

    d[n - 1] = bm1 * cn + s / n;
    const dj = d[n - 1] * j;
    sum = sum + dj;

    if (sum <= 0.0) {
      return w;
    }

    if (Math.abs(dj) <= eps * (sum + l)) {
      return w + uExp * sum;
    }
  }
  return w + uExp * sum;
}

/**
 * Evaluation of the incomplete gamma ratio functions P(a,x) and Q(a,x)
 *
 * It is assumed that a <= 1. eps is the tolerance to be used.
 * The input argument r has the value e^(-x)*x^a/gamma(a).
 *
 * @param {number} a - Parameter a, assumed <= 1
 * @param {number} x - Parameter x
 * @param {number} r - Value e^(-x)*x^a/gamma(a)
 * @param {number} eps - Tolerance to be used
 * @returns {Object} Object with properties p and q containing the computed values
 */
function grat1(a: number, x: number, r: number, eps: number) {
  // binomial_distribution(10, 11.0, 0.9) partially tests this fn

  let j,
    l,
    an,
    c,
    sum,
    tol,
    t,
    z,
    h,
    g,
    w,
    a2nm1,
    a2n,
    b2nm1,
    b2n,
    am0,
    cma,
    an0,
    p,
    q;

  if (a * x === 0.0) {
    return x <= a ? 1.0 : 0.0;
  }

  if (a === 0.5) {
    if (x < 0.25) {
      return 0.5 + (0.5 - erf(Math.sqrt(x)));
    } else {
      return erfc1(0, Math.sqrt(x));
    }
  }

  if (x < 1.1) {
    // Taylor series for P(a,x)/x**a
    an = 3.0;
    c = x;
    sum = x / (a + 3.0);
    tol = (0.1 * eps) / (a + 1.0);

    while (true) {
      an = an + 1.0;
      c = -c * (x / an);
      t = c / (a + an);
      sum = sum + t;
      if (Math.abs(t) <= tol) {
        break;
      }
    }

    j = a * x * ((sum / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));
    z = a * Math.log(x);
    h = gam1(a);
    g = 1.0 + h;

    if (x < 0.25) {
      if (z <= -0.13394) {
        w = Math.exp(z);
        p = w * g * (0.5 + (0.5 - j));
        return 0.5 + (0.5 - p);
      }
    } else if (a >= x / 2.59) {
      w = Math.exp(z);
      p = w * g * (0.5 + (0.5 - j));
      return 0.5 + (0.5 - p);
    }

    l = Math.expm1(z);
    w = 0.5 + (0.5 + l);
    q = (w * j - l) * g - h;

    if (q < 0.0) {
      return 0.0;
    }

    return 0.5 + (0.5 - q);
  } else {
    // Continued fraction expansion
    a2nm1 = 1.0;
    a2n = 1.0;
    b2nm1 = x;
    b2n = x + (1.0 - a);
    c = 1.0;

    while (true) {
      a2nm1 = x * a2n + c * a2nm1;
      b2nm1 = x * b2n + c * b2nm1;
      am0 = a2nm1 / b2nm1;
      c = c + 1.0;
      cma = c - a;
      a2n = a2nm1 + cma * a2n;
      b2n = b2nm1 + cma * b2n;
      an0 = a2n / b2n;

      if (Math.abs(an0 - am0) < eps * an0) {
        break;
      }
    }

    return r * an0;
  }
}

/**
 * Evaluation of the complementary error function
 *
 * erfc1(ind, x) = erfc(x)            if ind = 0
 * erfc1(ind, x) = exp(x*x)*erfc(x)   otherwise
 *
 * @param {number} ind - Mode indicator
 * @param {number} x - Input value
 * @returns {number} Complementary error function value
 */
function erfc1(ind: number, x: number): number {
  // beta_distribution(0.5, 180, 190) tests this fn

  const c = 0.564189583547756;

  const a = [
    0.77105849500132e-4, -0.133733772997339e-2, 0.323076579225834e-1,
    0.479137145607681e-1, 0.128379167095513,
  ];

  const b = [0.301048631703895e-2, 0.538971687740286e-1, 0.375795757275549];

  const p = [
    -1.36864857382717e-7, 5.64195517478974e-1, 7.21175825088309,
    4.31622272220567e1, 1.5298928504694e2, 3.39320816734344e2,
    4.51918953711873e2, 3.00459261020162e2,
  ];

  const q = [
    1.0, 1.27827273196294e1, 7.70001529352295e1, 2.77585444743988e2,
    6.38980264465631e2, 9.3135409485061e2, 7.90950925327898e2,
    3.00459260956983e2,
  ];

  const r = [
    2.10144126479064, 2.62370141675169e1, 2.13688200555087e1, 4.6580782871847,
    2.82094791773523e-1,
  ];

  const s = [
    9.4153775055546e1, 1.8711481179959e2, 9.90191814623914e1,
    1.80124575948747e1,
  ];

  let top, bot, t, e, w;
  let result;

  // ABS(X) <= 0.5
  const ax = Math.abs(x);
  if (ax <= 0.5) {
    t = x * x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
    result = 0.5 + (0.5 - x * (top / bot));

    if (ind !== 0) {
      result = Math.exp(t) * result;
    }

    return result;
  }

  // 0.5 < ABS(X) <= 4
  if (ax <= 4.0) {
    top =
      ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax +
        p[5]) *
        ax +
        p[6]) *
        ax +
      p[7];
    bot =
      ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax +
        q[5]) *
        ax +
        q[6]) *
        ax +
      q[7];
    result = top / bot;
  } else {
    // ABS(X) > 4
    if (x <= -5.6) {
      result = 2.0;
      if (ind !== 0) {
        result = 2.0 * Math.exp(x * x);
      }
      return result;
    }

    if (ind === 0) {
      if (x > 100.0) {
        return 0.0;
      }

      // Check if x*x would cause underflow in exp function
      if (x * x > 745) {
        return 0.0;
      }
    }

    t = (1.0 / x) ** 2;
    top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    result = (c - (t * top) / bot) / ax;
  }

  // Final assembly
  if (ind === 0) {
    w = x * x;
    t = Math.floor(w); // Integer part of w
    e = w - t; // Fractional part of w
    result = (0.5 + (0.5 - e)) * Math.exp(-t) * result;

    if (x < 0.0) {
      result = 2.0 - result;
    }
  } else {
    if (x < 0.0) {
      result = 2.0 * Math.exp(x * x) - result;
    }
  }

  return result;
}

/**
 * Asymptotic expansion for IX(A,B) for large A and B.
 * LAMBDA = (A + B)*Y - B and EPS is the tolerance used.
 * It is assumed that LAMBDA is nonnegative and that
 * A and B are greater than or equal to 15.
 *
 * @param {number} a - Parameter a, assumed >= 15
 * @param {number} b - Parameter b, assumed >= 15
 * @param {number} lambda - Parameter lambda, (a + b)*y - b, assumed >= 0
 * @param {number} eps - Tolerance to be used
 * @returns {number} Asymptotic expansion result
 */
function basym(a: number, b: number, lambda: number, eps: number): number {
  // beta_distribution(0.5, 180, 190) tests this fn

  // NUM is the maximum value that n can take in the loop.
  // It is required that NUM be even.
  const NUM = 20;

  // Constants
  // e0 = 2/sqrt(pi)
  // e1 = 2**(-3/2)
  const e0 = 1.12837916709551;
  const e1 = 0.353553390593274;

  const a0 = new Array(NUM + 1).fill(0);
  const b0 = new Array(NUM + 1).fill(0);
  const c = new Array(NUM + 1).fill(0);
  const d = new Array(NUM + 1).fill(0);

  let h, r0, r1, w0, f, t, z0, z, z2, j0, j1, sum;
  let s, h2, hn, w, znm1, zn, t0, t1, u;

  let result = 0.0;

  if (a >= b) {
    h = b / a;
    r0 = 1.0 / (1.0 + h);
    r1 = (b - a) / a;
    w0 = 1.0 / Math.sqrt(b * (1.0 + h));
  } else {
    h = a / b;
    r0 = 1.0 / (1.0 + h);
    r1 = (b - a) / b;
    w0 = 1.0 / Math.sqrt(a * (1.0 + h));
  }

  f = a * rlog1(-lambda / a) + b * rlog1(lambda / b);
  t = Math.exp(-f);

  if (t === 0.0) {
    return result;
  }

  z0 = Math.sqrt(f);
  z = 0.5 * (z0 / e1);
  z2 = f + f;

  a0[1] = (2.0 / 3.0) * r1;
  c[1] = -0.5 * a0[1];
  d[1] = -c[1];
  j0 = (0.5 / e0) * erfc1(1, z0);
  j1 = e1;
  sum = j0 + d[1] * w0 * j1;

  s = 1.0;
  h2 = h * h;
  hn = 1.0;
  w = w0;
  znm1 = z;
  zn = z2;

  for (let n = 2; n <= NUM; n += 2) {
    hn = h2 * hn;
    a0[n] = (2.0 * r0 * (1.0 + h * hn)) / (n + 2.0);
    let np1 = n + 1;
    s = s + hn;
    a0[np1] = (2.0 * r1 * s) / (n + 3.0);

    for (let i = n; i <= np1; i++) {
      let r = -0.5 * (i + 1.0);
      b0[1] = r * a0[1];

      for (let m = 2; m <= i; m++) {
        let bsum = 0.0;
        let mm1 = m - 1;

        for (let j = 1; j <= mm1; j++) {
          let mmj = m - j;
          bsum = bsum + (j * r - mmj) * a0[j] * b0[mmj];
        }

        b0[m] = r * a0[m] + bsum / m;
      }

      c[i] = b0[i] / (i + 1.0);

      let dsum = 0.0;
      let im1 = i - 1;

      for (let j = 1; j <= im1; j++) {
        let imj = i - j;
        dsum = dsum + d[imj] * c[j];
      }

      d[i] = -(dsum + c[i]);
    }

    j0 = e1 * znm1 + (n - 1.0) * j0;
    j1 = e1 * zn + n * j1;
    znm1 = z2 * znm1;
    zn = z2 * zn;
    w = w0 * w;
    t0 = d[n] * w * j0;
    w = w0 * w;
    t1 = d[np1] * w * j1;
    sum = sum + (t0 + t1);

    if (Math.abs(t0) + Math.abs(t1) <= eps * sum) {
      break;
    }
  }

  u = Math.exp(-bcorr(a, b)); // External function call
  result = e0 * t * u * sum;

  return result;
}

/**
 * Evaluation of the function ln(gamma(a + b))
 * For 1 <= a <= 2 and 1 <= b <= 2
 *
 * @param {number} a - First parameter (1 <= a <= 2)
 * @param {number} b - Second parameter (1 <= b <= 2)
 * @returns {number} Natural logarithm of gamma(a + b)
 */
function gsumln(a: number, b: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) tests this fn

  const x = a + b - 2.0;

  if (x <= 0.25) {
    return gamln1(1.0 + x);
  } else if (x <= 1.25) {
    return gamln1(x) + Math.log1p(x); // alnrel(x);
  } else {
    return gamln1(x - 1.0) + Math.log(x * (1.0 + x));
  }
}

/**
 * Evaluation of exp(mu) * (x**a*y**b/beta(a,b))
 *
 * @param {number} mu
 * @param {number} a
 * @param {number} b
 * @param {number} x
 * @param {number} y
 * @returns {number}
 */
function brcmp1(
  mu: number,
  a: number,
  b: number,
  x: number,
  y: number,
): number {
  // binomial_distribution(3.0, 10.0, 0.7) tests this fn

  // CONST = 1/SQRT(2*PI)
  const CONST = 0.398942280401433;

  let a0 = Math.min(a, b);

  if (a0 >= 8.0) {
    let h, x0, y0, lambda;

    if (a > b) {
      h = b / a;
      x0 = 1.0 / (1.0 + h);
      y0 = h / (1.0 + h);
      lambda = (a + b) * y - b;
    } else {
      h = a / b;
      x0 = h / (1.0 + h);
      y0 = 1.0 / (1.0 + h);
      lambda = a - (a + b) * x;
    }

    let u, v;
    let e = -lambda / a;
    if (Math.abs(e) > 0.6) {
      u = e - Math.log(x / x0);
    } else {
      u = rlog1(e);
    }

    e = lambda / b;
    if (Math.abs(e) > 0.6) {
      v = e - Math.log(y / y0);
    } else {
      v = rlog1(e);
    }

    let z = esum(mu, -(a * u + b * v));
    return CONST * Math.sqrt(b * x0) * z * Math.exp(-bcorr(a, b));
  }

  // Process for a0 < 8.0
  let lnx, lny;

  if (x > 0.375) {
    if (y > 0.375) {
      lnx = Math.log(x);
      lny = Math.log(y);
    } else {
      lnx = Math.log1p(-y); // alnrel(-y);
      lny = Math.log(y);
    }
  } else {
    lnx = Math.log(x);
    lny = Math.log1p(-x); // alnrel(-x);
  }

  let z = a * lnx + b * lny;

  if (a0 >= 1.0) {
    z = z - betaln(a, b);
    return esum(mu, z);
  }

  // Procedure for a < 1 or b < 1
  let b0 = Math.max(a, b);

  if (b0 >= 8.0) {
    let u = gamln1(a0) + algdiv(a0, b0);
    return a0 * esum(mu, z - u);
  }

  if (b0 > 1.0) {
    let u = gamln1(a0);
    let n = Math.floor(b0 - 1.0);
    let c = 1.0;

    if (n >= 1) {
      for (let i = 1; i <= n; i++) {
        b0 = b0 - 1.0;
        c = c * (b0 / (a0 + b0));
      }
      u = Math.log(c) + u;
    }

    z = z - u;
    b0 = b0 - 1.0;
    let apb = a0 + b0;
    let t;

    if (apb > 1.0) {
      let dblU = a0 + b0 - 1.0;
      t = (1.0 + gam1(dblU)) / apb;
    } else {
      t = 1.0 + gam1(apb);
    }

    return (a0 * esum(mu, z) * (1.0 + gam1(b0))) / t;
  }

  // Algorithm for b0 <= 1
  let result = esum(mu, z);
  if (result === 0.0) {
    return result;
  }

  let apb = a + b;
  let z1;

  if (apb > 1.0) {
    let u = a + b - 1.0;
    z1 = (1.0 + gam1(u)) / apb;
  } else {
    z1 = 1.0 + gam1(apb);
  }

  let c = ((1.0 + gam1(a)) * (1.0 + gam1(b))) / z1;
  return (result * (a0 * c)) / (1.0 + a0 / b0);
}

/**
 * Evaluation of exp(mu + x)
 *
 * This function calculates exp(mu + x) in a way that attempts to avoid
 * numerical issues when adding mu and x directly.
 *
 * @param {number} mu
 * @param {number} x
 * @returns {number} The value of exp(mu + x)
 */
function esum(mu: number, x: number): number {
  // binomial_distribution(3.0, 10.0, 0.7) tests this fn

  if (x <= 0.0) {
    if (mu >= 0) {
      let w = mu + x;
      if (w <= 0.0) {
        return Math.exp(w);
      }
    }
  } else {
    if (mu <= 0) {
      let w = mu + x;
      if (w >= 0.0) {
        return Math.exp(w);
      }
    }
  }

  let w = mu;
  return Math.exp(w) * Math.exp(x);
}

/**
 * Evaluation of the incomplete beta function IX(A,B)
 *
 * It is assumed that A and B are nonnegative, and that X <= 1
 * and Y = 1 - X. This function returns the value of IX(A,B).
 *
 * @param {number} a - First parameter of the beta function
 * @param {number} b - Second parameter of the beta function
 * @param {number} x - Value between 0 and 1
 * @param {number} y - Equal to 1 - x
 * @returns {number} The value of the incomplete beta function IX(A,B)
 */
export function bratio(a: number, b: number, x: number): number {
  const y = 1 - x;
  const eps = Number.EPSILON;

  if (a < 0.0 || b < 0.0) {
    return 0.0;
  }

  if (a === 0.0 && b === 0.0) {
    return 0.0;
  }

  if (x < 0.0 || x > 1.0) {
    return 0.0;
  }

  if (y < 0.0 || y > 1.0) {
    return 0.0;
  }

  const z = x + y - 0.5 - 0.5;
  if (Math.abs(z) > 3.0 * eps) {
    return 0.0;
  }

  // Special cases
  if (x === 0.0) {
    if (a === 0.0) {
      return 0.0;
    }
    return 0.0;
  }

  if (y === 0.0) {
    if (b === 0.0) {
      return 0.0;
    }
    return 1.0;
  }

  if (a === 0.0) {
    return 1.0;
  }

  if (b === 0.0) {
    return 0.0;
  }

  // Procedure for A and B < 1.E-3*EPS
  const epsMax = Math.max(eps, 1.0e-15);
  if (Math.max(a, b) < 1.0e-3 * epsMax) {
    return b / (a + b);
  }

  // Main algorithm
  let ind = 0;
  let a0 = a;
  let b0 = b;
  let x0 = x;
  let y0 = y;
  let lambda = 0;
  let w = 0;
  let w1 = 0;
  let n = 0;

  if (Math.min(a0, b0) <= 1.0) {
    if (x <= 0.5) {
      if (b0 < Math.min(epsMax, epsMax * a0)) {
        w = fpser(a0, b0, x0);
      } else if (a0 < Math.min(epsMax, epsMax * b0) && b0 * x0 <= 1.0) {
        w1 = apser(a0, b0, x0);
        w = 0.5 + (0.5 - w1);
      } else if (Math.max(a0, b0) > 1.0) {
        if (b0 <= 1.0) {
          w = bpser(a0, b0, x0);
        } else if (x0 >= 0.3) {
          w1 = bpser(b0, a0, y0);
          w = 0.5 + (0.5 - w1);
        } else if (x0 >= 0.1 || (x0 * b0) ** a0 <= 0.7) {
          w = bpser(a0, b0, x0);
        } else if (b0 > 15.0) {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        } else {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          w1 = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        }
      } else {
        if (a0 >= Math.min(0.2, b0)) {
          w = bpser(a0, b0, x0);
        } else if (x0 ** a0 <= 0.9) {
          w = bpser(a0, b0, x0);
        } else if (x0 >= 0.3) {
          w1 = bpser(b0, a0, y0);
          w = 0.5 + (0.5 - w1);
        } else {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          w1 = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        }
      }
    } else {
      // Case: x > 0.5, swap parameters
      ind = 1;
      a0 = b;
      b0 = a;
      x0 = y;
      y0 = x;

      if (b0 < Math.min(epsMax, epsMax * a0)) {
        w = fpser(a0, b0, x0);
      } else if (a0 < Math.min(epsMax, epsMax * b0) && b0 * x0 <= 1.0) {
        w1 = apser(a0, b0, x0);
        w = 0.5 + (0.5 - w1);
      } else if (Math.max(a0, b0) > 1.0) {
        if (b0 <= 1.0) {
          w = bpser(a0, b0, x0);
        } else if (x0 >= 0.3) {
          w1 = bpser(b0, a0, y0);
          w = 0.5 + (0.5 - w1);
        } else if (x0 >= 0.1 || (x0 * b0) ** a0 <= 0.7) {
          w = bpser(a0, b0, x0);
        } else if (b0 > 15.0) {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          w1 = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        } else {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          w1 = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        }
      } else {
        if (a0 >= Math.min(0.2, b0)) {
          w = bpser(a0, b0, x0);
        } else if (x0 ** a0 <= 0.9) {
          w = bpser(a0, b0, x0);
        } else if (x0 >= 0.3) {
          w1 = bpser(b0, a0, y0);
          w = 0.5 + (0.5 - w1);
        } else {
          n = 20;
          w1 = bup(b0, a0, y0, x0, n);
          b0 = b0 + n;
          w1 = bgrat(b0, a0, y0, x0, w1, 15.0 * eps);
          w = 0.5 + (0.5 - w1);
        }
      }
    }
  } else {
    // Procedure for a0 > 1 and b0 > 1

    if (a > b) {
      lambda = (a + b) * y - b;
    } else {
      lambda = a - (a + b) * x;
    }

    if (lambda < 0.0) {
      ind = 1;
      a0 = b;
      b0 = a;
      x0 = y;
      y0 = x;
      lambda = Math.abs(lambda);
    }

    if (b0 < 40.0 && b0 * x0 <= 0.7) {
      w = bpser(a0, b0, x0);
    } else if (b0 < 40.0) {
      n = Math.floor(b0);
      b0 = b0 - n;
      if (b0 === 0.0) {
        n = n - 1;
        b0 = 1.0;
      }
      w = bup(b0, a0, y0, x0, n);

      if (x0 > 0.7) {
        if (a0 > 15.0) {
          w = bgrat(a0, b0, x0, y0, w, 15.0 * eps);
        } else {
          n = 20;
          w = w + bup(a0, b0, x0, y0, n);
          a0 = a0 + n;
          w = bgrat(a0, b0, x0, y0, w, 15.0 * eps);
        }
      } else {
        w = w + bpser(a0, b0, x0);
      }
    } else {
      let minab = Math.min(a0, b0);

      if (minab <= 100.0 || lambda > 0.03 * minab) {
        w = bfrac(a0, b0, x0, y0, lambda, 15.0 * eps);
      } else {
        w = basym(a0, b0, lambda, 100.0 * eps);
      }
    }
  }

  if (ind === 1) {
    w = 0.5 + (0.5 - w);
  }

  return w;
}
