import { expect, test, vi } from "vitest";
import {
  f_density,
  f_distribution,
  f_log_density,
  f_quantile,
  f_random,
  t_density,
  t_distribution,
  t_quantile,
} from "../src/index";
import { Welford } from "../src/welford";

// R 4.4.1 df(x, df1, df2), df(..., log=TRUE), and pf(x, df1, df2).
test.each([
  [
    0.1, 0.5, 3, 0.99843249576068793, -0.0015687340594157106,
    0.40872612733371255,
  ],
  [0.5, 1, 5, 0.40333954103270347, -0.90797653818772606, 0.48891591956971919],
  [1, 2, 7, 0.32273783965178487, -1.1309149272640773, 0.5850513490191338],
  [2.5, 5, 10, 0.093594754449685555, -2.3687809392732277, 0.89799772335573025],
  [10, 10, 20, 6.6529666038162549e-6, -11.920447696935106, 0.99999140588038204],
  [100, 30, 50, 6.8477265589960622e-36, -80.969146639796577, 1],
  [
    1e-10, 5, 10, 1.0368367495639181e-14, -32.200016810747684,
    4.1473469987000623e-25,
  ],
  [1e200, 1, 0.5, 1.6035048770711136e-251, -577.4766665598429, 1],
  [
    0.99, 100, 200, 2.3181390453976314, 0.84076472819084458,
    0.48459442033497774,
  ],
  [1, 1e6, 1e6, 199.47109033293728, 5.2956693152175172, 0.5],
  [
    1.01, 1000, 1000, 6.1670185901446306, 1.8192155104599173,
    0.56249141340369624,
  ],
])(
  "F: density, log-density, and CDF x=%s df1=%s df2=%s",
  (x, df1, df2, d, ld, p) => {
    expect(Math.abs(f_density(x, df1, df2) / d - 1)).toBeLessThan(2e-12);
    expect(Math.abs(f_log_density(x, df1, df2) - ld)).toBeLessThan(2e-12);
    expect(Math.abs(f_distribution(x, df1, df2) / p - 1)).toBeLessThan(2e-12);
  },
);

// R 4.4.1 qf reference values at well-conditioned probabilities.
test.each([
  [0.025, 1, 5, 0.001084777647615276],
  [0.25, 2, 7, 0.29983576605344997],
  [0.5, 5, 10, 0.93193316085104794],
  [0.975, 10, 20, 2.7736713751990809],
  [0.99, 30, 50, 2.0975934395718738],
  [0.975, 1000, 1000, 1.1320494987955811],
  [0.025, 100, 200, 0.70408485219391848],
])("F: quantile p=%s df1=%s df2=%s", (p, df1, df2, expected) => {
  expect(Math.abs(f_quantile(p, df1, df2) / expected - 1)).toBeLessThan(3e-12);
});

test("F: exact F(1,1) and F(2,2) identities and extreme tails", () => {
  for (const x of [1e-200, 1e-10, 0.1, 1, 10, 1e200]) {
    const cdf = (2 / Math.PI) * Math.atan(Math.sqrt(x));
    const ld = -Math.log(Math.PI) - 0.5 * Math.log(x) - Math.log1p(x);
    expect(Math.abs(f_distribution(x, 1, 1) / cdf - 1)).toBeLessThan(2e-13);
    expect(Math.abs(f_log_density(x, 1, 1) - ld)).toBeLessThan(3e-13);
    expect(f_distribution(x, 2, 2)).toBeCloseTo(x / (1 + x), 14);
    expect(
      Math.abs(f_log_density(x, 2, 2) / (-2 * Math.log1p(x)) - 1),
    ).toBeLessThan(2e-13);
  }
  for (const p of [1e-100, 1e-20, 1e-10, 0.025, 0.25, 0.5, 0.9, 1 - 1e-10]) {
    const tangent =
      p <= 0.5
        ? Math.tan((Math.PI * p) / 2)
        : 1 / Math.tan((Math.PI * (1 - p)) / 2);
    expect(
      Math.abs(f_quantile(p, 1, 1) / (tangent * tangent) - 1),
    ).toBeLessThan(3e-13);
    expect(Math.abs(f_quantile(p, 2, 2) / (p / (1 - p)) - 1)).toBeLessThan(
      2e-13,
    );
  }
  expect(f_distribution(1, Number.MAX_VALUE, Number.MAX_VALUE)).toBeCloseTo(
    0.5,
    14,
  );
  expect(
    Number.isFinite(f_log_density(1, Number.MAX_VALUE, Number.MAX_VALUE)),
  ).toBe(true);
});

test("F: reciprocal symmetry and CDF/quantile round trips", () => {
  for (const [df1, df2] of [
    [0.1, 0.5],
    [1, 5],
    [2, 7],
    [10, 20],
    [1000, 0.5],
    [1e6, 1e6],
  ]) {
    for (const x of [0.1, 1, 10]) {
      expect(
        f_distribution(x, df1, df2) + f_distribution(1 / x, df2, df1),
      ).toBeCloseTo(1, 13);
      const expectedLogDensity =
        f_log_density(1 / x, df2, df1) - 2 * Math.log(x);
      expect(
        Math.abs(f_log_density(x, df1, df2) - expectedLogDensity),
      ).toBeLessThan(5e-13 * Math.max(1, Math.abs(expectedLogDensity)));
    }
    for (const p of [1e-20, 1e-10, 0.025, 0.25, 0.5]) {
      const q = f_quantile(p, df1, df2);
      if (q === 0 || !Number.isFinite(q)) continue;
      expect(Math.abs(f_distribution(q, df1, df2) / p - 1)).toBeLessThan(2e-10);
    }
    for (const p of [0.25, 0.5]) {
      expect(
        Math.abs(f_quantile(p, df1, df2) * f_quantile(1 - p, df2, df1) - 1),
      ).toBeLessThan(2e-12);
    }
  }
});

test("F: two numerator degrees give an exact power survival", () => {
  for (const df2 of [0.5, 1, 10, 1000]) {
    const b = df2 / 2;
    for (const p of [1e-100, 1e-10, 0.025, 0.5, 0.975, 1 - 1e-10]) {
      const expected = b * Math.expm1(-Math.log1p(-p) / b);
      expect(Math.abs(f_quantile(p, 2, df2) / expected - 1)).toBeLessThan(
        2e-12,
      );
    }
    for (const x of [1e-100, 0.1, 1, 10]) {
      const expected = -Math.expm1(-b * Math.log1p(x / b));
      expect(Math.abs(f_distribution(x, 2, df2) / expected - 1)).toBeLessThan(
        2e-12,
      );
    }
  }
  // The transformed beta quantile underflows here, but the F quantile fits.
  const q = f_quantile(1e-100, 1, 1e300);
  expect(Math.abs(q / (Math.PI * 0.5e-200) - 1)).toBeLessThan(3e-13);
  expect(Math.abs(f_distribution(q, 1, 1e300) / 1e-100 - 1)).toBeLessThan(
    3e-13,
  );
});

test("F: square of a Student t has F(1,df) distribution", () => {
  for (const df of [0.5, 1, 2, 5, 30, 100]) {
    for (const x of [0.01, 0.5, 1, 3]) {
      expect(f_distribution(x * x, 1, df)).toBeCloseTo(
        2 * t_distribution(x, df) - 1,
        13,
      );
      expect(
        Math.abs(f_density(x * x, 1, df) / (t_density(x, df) / x) - 1),
      ).toBeLessThan(2e-12);
    }
    for (const p of [0.25, 0.5, 0.75]) {
      const expected = t_quantile((1 + p) / 2, df) ** 2;
      expect(Math.abs(f_quantile(p, 1, df) / expected - 1)).toBeLessThan(3e-12);
    }
  }
});

test("F: support and endpoint values", () => {
  for (const x of [-Infinity, -1]) {
    expect(f_density(x, 5, 10)).toBe(0);
    expect(f_log_density(x, 5, 10)).toBe(-Infinity);
    expect(f_distribution(x, 5, 10)).toBe(0);
  }
  expect(f_density(0, 1, 10)).toBe(Infinity);
  expect(f_density(0, 2, 10)).toBe(1);
  expect(f_density(0, 5, 10)).toBe(0);
  expect(f_log_density(0, 1, 10)).toBe(Infinity);
  expect(f_log_density(0, 2, 10)).toBe(0);
  expect(f_log_density(0, 5, 10)).toBe(-Infinity);
  expect(f_distribution(0, 5, 10)).toBe(0);
  expect(f_density(Infinity, 5, 10)).toBe(0);
  expect(f_log_density(Infinity, 5, 10)).toBe(-Infinity);
  expect(f_distribution(Infinity, 5, 10)).toBe(1);
  expect(f_quantile(0, 5, 10)).toBe(0);
  expect(f_quantile(1, 5, 10)).toBe(Infinity);
});

test("F: invalid and infinite degrees of freedom are rejected before endpoints", () => {
  for (const bad of [0, -1, NaN, Infinity, -Infinity]) {
    for (const [df1, df2] of [
      [bad, 5],
      [5, bad],
    ]) {
      for (const f of [f_density, f_log_density, f_distribution]) {
        expect(f(0, df1, df2)).toBeNaN();
        expect(f(Infinity, df1, df2)).toBeNaN();
      }
      for (const p of [0, 0.5, 1]) expect(f_quantile(p, df1, df2)).toBeNaN();
      expect(f_random(3, df1, df2)).toEqual([NaN, NaN, NaN]);
    }
  }
  for (const f of [f_density, f_log_density, f_distribution])
    expect(f(NaN, 5, 10)).toBeNaN();
  for (const p of [-1, 2, NaN, Infinity, -Infinity])
    expect(f_quantile(p, 5, 10)).toBeNaN();
  for (const N of [-1, 0.5, NaN, Infinity])
    expect(() => f_random(N, 5, 10)).toThrow(RangeError);
  expect(f_random(0, 5, 10)).toEqual([]);
});

test("F: random draws follow the CDF and finite moments", () => {
  let state = 709;
  const random = vi.spyOn(Math, "random").mockImplementation(() => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return (state + 0.5) / 4294967296;
  });
  try {
    const N = 30000;
    for (const [df1, df2] of [
      [0.5, 3],
      [1, 5],
      [5, 10],
      [30, 50],
    ]) {
      const samples = f_random(N, df1, df2);
      expect(samples).toHaveLength(N);
      expect(samples.every(x => Number.isFinite(x) && x >= 0)).toBe(true);
      for (const p of [0.1, 0.5, 0.9]) {
        const cutoff = f_quantile(p, df1, df2);
        const fraction = samples.filter(x => x <= cutoff).length / N;
        expect(Math.abs(fraction - p)).toBeLessThan(0.012);
      }
    }
    const moments = new Welford();
    for (const x of f_random(N, 10, 30)) moments.update(x);
    const mean = 30 / 28;
    const variance = (2 * 30 ** 2 * (10 + 30 - 2)) / (10 * 28 ** 2 * 26);
    expect(Math.abs(moments.mean() / mean - 1)).toBeLessThan(0.025);
    expect(Math.abs(moments.var() / variance - 1)).toBeLessThan(0.06);
  } finally {
    random.mockRestore();
  }
});
