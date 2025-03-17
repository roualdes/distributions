import { log_gamma } from "./log_gamma_fn";

export function xinbeta(p: number, q: number, alpha: number): number {
  const beta = log_gamma(p) + log_gamma(q) - log_gamma(p + q);

  let a;
  let acu;
  let adj;
  let fpu;
  let g;
  let h;
  let iex;
  let indx;
  let pp;
  let prev;
  let qq;
  let r;
  let s;
  let sae = -37.0;
  let sq;
  let t;
  let tx;
  let value;
  let w;
  let xin;
  let y;
  let yprev;

  fpu = Math.pow(10.0, sae);

  value = alpha;

  // Test for admissibility of parameters.
  if (p <= 0.0) {
    return Number.NaN;
  }

  if (q <= 0.0) {
    return Number.NaN;
  }

  if (alpha < 0.0 || 1.0 < alpha) {
    return Number.NaN;
  }

  // If answer is easy to determine, return immediately.
  if (alpha == 0.0) {
    return 0.0;
  }

  if (alpha == 1.0) {
    return 1.0;
  }

  // Change tail if necessary.
  if (0.5 < alpha) {
    a = 1.0 - alpha;
    pp = q;
    qq = p;
    indx = 1;
  } else {
    a = alpha;
    pp = p;
    qq = q;
    indx = 0;
  }

  // Calculate the initial approximation.
  r = Math.sqrt(-Math.log(a * a));

  y = r - (2.30753 + 0.27061 * r) / (1.0 + (0.99229 + 0.04481 * r) * r);

  if (1.0 < pp && 1.0 < qq) {
    r = (y * y - 3.0) / 6.0;
    s = 1.0 / (pp + pp - 1.0);
    t = 1.0 / (qq + qq - 1.0);
    h = 2.0 / (s + t);
    w =
      (y * Math.sqrt(h + r)) / h - (t - s) * (r + 5.0 / 6.0 - 2.0 / (3.0 * h));
    value = pp / (pp + qq * Math.exp(w + w));
  } else {
    r = qq + qq;
    t = 1.0 / (9.0 * qq);
    t = r * Math.pow(1.0 - t + y * Math.sqrt(t), 3);

    if (t <= 0.0) {
      value = 1.0 - Math.exp((Math.log((1.0 - a) * qq) + beta) / qq);
    } else {
      t = (4.0 * pp + r - 2.0) / t;

      if (t <= 1.0) {
        value = Math.exp((Math.log(a * pp) + beta) / pp);
      } else {
        value = 1.0 - 2.0 / (t + 1.0);
      }
    }
  }

  // Solve for X by a modified Newton-Raphson method, using the function BETAIN.
  r = 1.0 - pp;
  t = 1.0 - qq;
  yprev = 0.0;
  sq = 1.0;
  prev = 1.0;

  if (value < 0.0001) {
    value = 0.0001;
  }

  if (0.9999 < value) {
    value = 0.9999;
  }

  iex = Math.max(-5.0 / pp / pp - 1.0 / Math.pow(a, 0.2) - 13.0, sae);

  acu = Math.pow(10.0, iex);

  outer: while (true) {
    y = betain(value, pp, qq, beta);

    if (Number.isNaN(y)) {
      return Number.NaN;
    }

    xin = value;
    y = (y - a) * Math.exp(beta + r * Math.log(xin) + t * Math.log(1.0 - xin));

    if (y * yprev <= 0.0) {
      prev = Math.max(sq, fpu);
    }

    g = 1.0;

    while (true) {
      // Choose damping factor.
      while (true) {
        adj = g * y;
        sq = adj * adj;

        if (sq < prev) {
          tx = value - adj;

          if (tx >= 0.0 && tx <= 1.0) {
            break;
          }
        }
        g = g / 3.0;
      }

      // Check whether current estimate is acceptable.
      if (prev <= acu || y * y <= acu) {
        break outer;
      }

      if (tx != 0.0 && tx != 1.0) {
        break;
      }
      g = g / 3.0;
    }

    if (tx == value) {
      break;
    }

    value = tx;
    yprev = y;
  }

  if (indx) {
    value = 1.0 - value;
  }

  return value;
}

/**
 * Incomplete Beta function
 *
 * $I_x(p, q) = \frac{\Gamma(p + q)}{\Gamma(p)\Gamma(q)} \int_0^x t^{(p-1)}(1 - t)^{(q-1)}dt$
 *
 * @param x
 * @param p
 * @param q
 * @returns I_x(a, b)
 *
 * @category beta
 */
export function betain(x: number, p: number, q: number, beta: number): number {
  let indx;
  let ns;

  let acu = 0.1e-14;
  let ai;
  let cx;
  let pp;
  let psq;
  let qq;
  let rx;
  let temp;
  let term;
  let value;
  let xx;

  value = x;

  if (p <= 0.0 || q <= 0.0) {
    return value;
  }

  if (x < 0.0 || 1.0 < x) {
    return value;
  }

  if (x == 0.0 || x == 1.0) {
    return value;
  }

  psq = p + q;
  cx = 1.0 - x;

  if (p < psq * x) {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = 1;
  } else {
    xx = x;
    pp = p;
    qq = q;
    indx = 0;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = qq + cx * psq;

  rx = xx / cx;
  temp = qq - ai;
  if (ns == 0) {
    rx = xx;
  }

  while (true) {
    term = (term * temp * rx) / (pp + ai);
    value = value + term;
    temp = Math.abs(term);

    if (temp <= acu && temp <= acu * value) {
      value =
        (value *
          Math.exp(pp * Math.log(xx) + (qq - 1.0) * Math.log(cx) - beta)) /
        pp;

      if (indx) {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if (0 <= ns) {
      temp = qq - ai;
      if (ns == 0) {
        rx = xx;
      }
    } else {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}
