import { expect, test, vi } from "vitest";
import {
  t_density,
  t_distribution,
  t_log_density,
  t_quantile,
  t_random,
} from "../src/index";
import {
  normal_density,
  normal_distribution,
  normal_log_density,
  normal_quantile,
} from "../src/normal";
import { Welford } from "../src/welford";

// R 4.4.1 dt(x, df), dt(x, df, log=TRUE), and pt(x, df).
test.each([
  [0, 1, 0.31830988618379069, -1.1447298858494002, 0.5],
  [1, 0.5, 0.11830465704039259, -2.1344921424125918, 0.69887838915867795],
  [-1, 2, 0.19245008972987526, -1.6479184330021646, 0.21132486540518713],
  [2.5, 5, 0.033326238887022824, -3.4014102377037108, 0.97275495032881187],
  [-3, 10, 0.011400549464542529, -4.4740937260991229, 0.006671827511284794],
  [0.1, 0.25, 0.20964921071430084, -1.5623195701504116, 0.52130954682675701],
  [10, 30, 5.3278145782741734e-11, -23.655494891676963, 0.99999999997712374],
  [
    -100, 100, 2.4076589692854326e-102, -233.9850245920266,
    2.4314972222696739e-102,
  ],
  [1e-10, 1e6, 0.39894218066587506, -0.91893878320467282, 0.5000000000398942],
  [
    -1e20, 1, 3.1830988618379062e-41, -93.24813360561123,
    3.1830988618379063e-21,
  ],
  [
    -1e200, 0.5, 1.6035048770710975e-301, -692.60592120954516,
    3.2070097541422692e-101,
  ],
])("Student t: density, log-density, and CDF x=%s df=%s", (x, df, d, ld, p) => {
  expect(Math.abs(t_density(x, df) / d - 1)).toBeLessThan(2e-12);
  expect(Math.abs(t_log_density(x, df) - ld)).toBeLessThan(2e-12);
  expect(Math.abs(t_distribution(x, df) / p - 1)).toBeLessThan(2e-12);
});

// R 4.4.1 qt reference values, including fractional and large df.
test.each([
  [0.025, 1, -12.706204736174703],
  [0.25, 1, -1],
  [0.975, 5, 2.5705818356363137],
  [0.99, 10, 2.7637694581126957],
  [1e-10, 30, -9.377489780407144],
  [1e-100, 0.5, -1.0284911563163961e199],
  [0.5 - 1e-12, 100, -2.5128469906434047e-12],
  [0.25, 0.01, -6.3641819284002167e28],
  [0.499999, 1e6, -2.5066289012237065e-6],
  [1e-200, 1, -3.1830988618379081e199],
  [1e-200, 1e8, -30.205663152517996],
  [0.025, 1e8, -1.9599640082627663],
  [1e-200, 1e12, -30.205594186476915],
  [0.75, 1e100, 0.67448975019608171],
])("Student t: quantile p=%s df=%s", (p, df, expected) => {
  expect(Math.abs(t_quantile(p, df) / expected - 1)).toBeLessThan(3e-12);
});

test("Student t: Cauchy and two-degree-of-freedom identities", () => {
  for (const x of [-100, -2, -0.01, 0, 0.01, 2, 100]) {
    expect(Math.abs(t_density(x, 1) * Math.PI * (1 + x * x) - 1)).toBeLessThan(
      2e-14,
    );
    expect(t_distribution(x, 1)).toBeCloseTo(0.5 + Math.atan(x) / Math.PI, 14);
    expect(t_distribution(x, 2)).toBeCloseTo(
      0.5 + x / (2 * Math.sqrt(2 + x * x)),
      14,
    );
  }
  // Unlike 0.5 + atan(x)/pi, this identity preserves the small Cauchy tail.
  expect(
    Math.abs(t_distribution(-1e200, 1) * 1e200 * Math.PI - 1),
  ).toBeLessThan(2e-13);
});

test("Student t: symmetry and CDF/quantile round trips", () => {
  for (const df of [0.01, 0.25, 0.5, 1, 2, 5, 30, 100, 1e6, 1e12]) {
    expect(t_distribution(0, df)).toBe(0.5);
    expect(t_quantile(0.5, df)).toBe(0);
    for (const p of [1e-10, 0.025, 0.25, 0.5 - 1e-12]) {
      const q = t_quantile(p, df);
      if (!Number.isFinite(q)) continue; // A heavy-tail quantile can overflow.
      expect(Math.abs(t_distribution(q, df) / p - 1)).toBeLessThan(3e-11);
    }
    for (const p of [0.025, 0.25]) {
      // Exact binary complements avoid introducing a tail rounding error.
      const lower = 1 - (1 - p);
      expect(
        Math.abs(t_quantile(lower, df) / -t_quantile(1 - p, df) - 1),
      ).toBeLessThan(2e-12);
    }
    const delta = 0.5 - 1e-12 - 0.5;
    expect(
      Math.abs((t_quantile(0.5 + delta, df) * t_density(0, df)) / delta - 1),
    ).toBeLessThan(2e-12);
  }
});

test("Student t: infinities, underflow, and invalid scalar arguments", () => {
  for (const df of [0.5, 1, 30, Infinity]) {
    expect(t_density(-Infinity, df)).toBe(0);
    expect(t_density(Infinity, df)).toBe(0);
    expect(t_log_density(-Infinity, df)).toBe(-Infinity);
    expect(t_log_density(Infinity, df)).toBe(-Infinity);
    expect(t_distribution(-Infinity, df)).toBe(0);
    expect(t_distribution(Infinity, df)).toBe(1);
    expect(t_quantile(0, df)).toBe(-Infinity);
    expect(t_quantile(1, df)).toBe(Infinity);
  }
  expect(t_density(1e200, 5)).toBe(0);
  expect(t_log_density(1e200, 5)).toBeCloseTo(-2759.2424174446073, 10);
  expect(t_quantile(1e-200, 0.01)).toBe(-Infinity);
  for (const df of [0, -1, NaN, -Infinity]) {
    expect(t_density(0, df)).toBeNaN();
    expect(t_log_density(0, df)).toBeNaN();
    expect(t_distribution(0, df)).toBeNaN();
    expect(t_quantile(0.5, df)).toBeNaN();
  }
  for (const f of [t_density, t_log_density, t_distribution])
    expect(f(NaN, 5)).toBeNaN();
  for (const p of [-1, 2, NaN, Infinity, -Infinity])
    expect(t_quantile(p, 5)).toBeNaN();
});

test("Student t: infinite df uses the existing standard normal functions", () => {
  for (const x of [-5, -0.5, 0, 0.5, 5]) {
    expect(t_density(x, Infinity)).toBe(normal_density(x, 0, 1));
    expect(t_log_density(x, Infinity)).toBe(normal_log_density(x, 0, 1));
    expect(t_distribution(x, Infinity)).toBe(normal_distribution(x, 0, 1));
  }
  for (const p of [1e-10, 0.025, 0.5, 0.975]) {
    expect(t_quantile(p, Infinity)).toBe(normal_quantile(p, 0, 1));
  }
});

test("Student t: random draws follow the CDF and finite moments", () => {
  // A fixed PRNG makes the statistical checks reproducible. Open-interval
  // uniforms also keep the underlying inverse-normal generator off endpoints.
  let state = 708;
  const random = vi.spyOn(Math, "random").mockImplementation(() => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return (state + 0.5) / 4294967296;
  });
  try {
    const N = 40000;
    for (const df of [0.5, 1, 5, 30, Infinity]) {
      const samples = t_random(N, df);
      expect(samples).toHaveLength(N);
      expect(samples.every(Number.isFinite)).toBe(true);
      for (const p of [0.1, 0.5, 0.9]) {
        const cutoff = t_quantile(p, df);
        const fraction = samples.filter(x => x <= cutoff).length / N;
        expect(Math.abs(fraction - p)).toBeLessThan(0.012);
      }
    }
    const moments = new Welford();
    for (const sample of t_random(N, 8)) moments.update(sample);
    expect(Math.abs(moments.mean())).toBeLessThan(0.04);
    expect(Math.abs(moments.std() / Math.sqrt(8 / 6) - 1)).toBeLessThan(0.04);
  } finally {
    random.mockRestore();
  }
});

test("Student t: random argument validation", () => {
  expect(t_random(0, 5)).toEqual([]);
  expect(t_random(0, Infinity)).toEqual([]);
  for (const df of [0, -1, NaN, -Infinity]) {
    expect(t_random(3, df)).toEqual([NaN, NaN, NaN]);
  }
  for (const N of [-1, 0.5, NaN, Infinity]) {
    expect(() => t_random(N, 5)).toThrow(RangeError);
  }
});
