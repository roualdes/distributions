import { betaln, bratio, bratio_complement } from "./beta_inc_fn";

/** Invert the regularized incomplete beta ratio for finite, positive shapes. */
export function xinbeta(p: number, q: number, alpha: number): number {
  if (
    !Number.isFinite(p) ||
    p <= 0 ||
    !Number.isFinite(q) ||
    q <= 0 ||
    !Number.isFinite(alpha) ||
    alpha < 0 ||
    alpha > 1
  ) {
    return Number.NaN;
  }
  if (alpha === 0 || alpha === 1) return alpha;

  // Keep the root in [0, 0.5] so reflection cannot erase a small result.
  const reflect = alpha > bratio(p, q, 0.5);
  const a = reflect ? 1 - alpha : alpha;
  const pp = reflect ? q : p;
  const qq = reflect ? p : q;
  // Independently choose the smaller probability tail. The paired-coordinate
  // CDF preserves it even when 1 - value rounds to one.
  const upperTail = reflect ? alpha < 0.5 : alpha > 0.5;
  const target = upperTail ? (reflect ? alpha : 1 - alpha) : a;
  const residual = (x: number): number =>
    upperTail
      ? target - bratio_complement(pp, qq, x)
      : bratio(pp, qq, x) - target;
  // The specialized log-beta routine avoids cancellation of large log-gammas.
  const beta = betaln(pp, qq);
  const finish = (x: number): number => (reflect ? 1 - x : x);

  // Initial approximation from AS 109. Keep the logarithm of the probability
  // instead of squaring it, which would underflow in extreme tails.
  let r = Math.sqrt(-2 * Math.log(a));
  const y = r - (2.30753 + 0.27061 * r) / (1 + (0.99229 + 0.04481 * r) * r);
  let value: number;
  if (pp > 1 && qq > 1) {
    r = (y * y - 3) / 6;
    const s = 1 / (2 * pp - 1);
    const t = 1 / (2 * qq - 1);
    const h = 2 / (s + t);
    const w = (y * Math.sqrt(h + r)) / h - (t - s) * (r + 5 / 6 - 2 / (3 * h));
    value = pp / (pp + qq * Math.exp(2 * w));
  } else {
    r = 2 * qq;
    let t = 1 / (9 * qq);
    t = r * (1 - t + y * Math.sqrt(t)) ** 3;
    if (t <= 0) {
      value = -Math.expm1((Math.log1p(-a) + Math.log(qq) + beta) / qq);
    } else {
      t = (4 * pp + r - 2) / t;
      value =
        t <= 1
          ? Math.exp((Math.log(a) + Math.log(pp) + beta) / pp)
          : 1 - 2 / (t + 1);
    }
  }

  // Near zero, I_x(p,q) ~ x^p / (p B(p,q)). This estimate also preserves
  // representable roots far below the old fixed 0.0001 starting bound.
  const logTail = (Math.log(a) + Math.log(pp) + beta) / pp;
  if (logTail < Math.log(Number.MIN_VALUE) - Math.LN2) return finish(0);
  const tail = Math.exp(logTail);
  if (a < 1 && tail < 0.1 && qq * tail < 0.1) value = tail;
  if (a === 1) value = pp / (pp + qq);
  if (!Number.isFinite(value)) value = pp / (pp + qq);
  value = Math.max(Number.MIN_VALUE, Math.min(0.5, value));

  let lower = 0;
  let upper = 0.5;
  for (let iteration = 0; iteration < 2048; iteration++) {
    const error = residual(value);
    if (!Number.isFinite(error)) return Number.NaN;
    if (Math.abs(error) <= 8 * Number.EPSILON * target) return finish(value);

    if (error < 0) lower = value;
    else upper = value;

    // Compute the Newton correction in log space so that neither the density
    // nor its reciprocal needs to be representable on its own.
    const step =
      Math.sign(error) *
      Math.exp(
        Math.log(Math.abs(error)) +
          beta +
          (1 - pp) * Math.log(value) +
          (1 - qq) * Math.log1p(-value),
      );
    let next = value - step;
    if (!(next > lower && next < upper)) next = lower + (upper - lower) / 2;
    if (next === value || next === lower || next === upper) {
      // The remaining bracket has no representable interior point.
      const lowerError = Math.abs(residual(lower));
      const upperError = Math.abs(residual(upper));
      return finish(lowerError < upperError ? lower : upper);
    }
    value = next;
  }
  // Never return an unconverged estimate as a successful quantile.
  return Number.NaN;
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

  if (!Number.isFinite(p) || !Number.isFinite(q) || p <= 0.0 || q <= 0.0) {
    return Number.NaN;
  }

  if (!Number.isFinite(x) || !Number.isFinite(beta) || x < 0.0 || 1.0 < x) {
    return Number.NaN;
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
  // NS is an integer in AS 63; a fractional value never reaches zero.
  ns = Math.floor(qq + cx * psq);

  rx = xx / cx;
  temp = qq - ai;
  if (ns == 0) {
    rx = xx;
  }

  for (let iteration = 0; iteration < 10000; iteration++) {
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
      return value;
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

  return Number.NaN;
}
