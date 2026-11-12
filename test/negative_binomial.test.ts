import { expect, test, vi } from "vitest";
import {
  negative_binomial_density,
  negative_binomial_distribution,
  negative_binomial_log_density,
  negative_binomial_quantile,
  negative_binomial_random,
  geometric_density,
  geometric_distribution,
  geometric_log_density,
  geometric_quantile,
} from "../src/index";
import { Welford } from "../src/welford";

// R 4.4.1 dnbinom, log-dnbinom, and pnbinom; size is allowed to be fractional.
test.each([
  [0, 0.5, 0.2, 0.44721359549995793, -0.80471895621705014, 0.44721359549995793],
  [1, 2, 0.5, 0.25000000000000006, -1.3862943611198904, 0.5],
  [5, 5, 0.3, 0.051459672600000013, -2.9669568343718637, 0.15026833259999989],
  [
    20, 15, 0.7, 0.00023042379745892574, -8.3755903468941906,
    0.99977679681256781,
  ],
  [
    100, 100, 0.6, 0.0004753084181357931, -7.6515466632960329,
    0.99831521348008057,
  ],
  [
    3, 0.25, 0.4, 0.020130290321915965, -3.9055296170839182,
    0.97942560689717084,
  ],
  [
    7, 2.5, 0.8, 0.00013046115007924743, -8.9444350760652949,
    0.99995963561899814,
  ],
])(
  "Negative binomial: mass and CDF x=%s size=%s p=%s",
  (x, size, p, mass, logMass, cdf) => {
    expect(
      Math.abs(negative_binomial_density(x, size, p) / mass - 1),
    ).toBeLessThan(2e-12);
    expect(
      Math.abs(negative_binomial_log_density(x, size, p) - logMass),
    ).toBeLessThan(2e-12);
    expect(
      Math.abs(negative_binomial_distribution(x, size, p) / cdf - 1),
    ).toBeLessThan(2e-12);
    expect(negative_binomial_distribution(x + 0.9, size, p)).toBe(
      negative_binomial_distribution(x, size, p),
    );
  },
);

test.each([
  [0.025, 0.5, 0.2, 0],
  [0.25, 2, 0.5, 0],
  [0.5, 5, 0.3, 11],
  [0.9, 15, 0.7, 10],
  [0.975, 100, 0.6, 88],
  [1 - 1e-10, 0.25, 0.4, 38],
  [1e-20, 100, 0.5, 10],
])(
  "Negative binomial: R quantile u=%s size=%s p=%s",
  (u, size, p, expected) => {
    expect(negative_binomial_quantile(u, size, p)).toBe(expected);
  },
);

test("Negative binomial: integer size agrees with direct mass sums", () => {
  for (const size of [2, 3, 5, 10]) {
    for (const p of [1e-10, 0.01, 0.2, 0.5, 0.9]) {
      let choose = 1;
      let sum = 0;
      for (let k = 0; k <= 40; k++) {
        if (k > 0) choose *= (size + k - 1) / k;
        const expected = choose * p ** size * (1 - p) ** k;
        sum += expected;
        expect(
          Math.abs(negative_binomial_density(k, size, p) / expected - 1),
        ).toBeLessThan(2e-12);
        expect(
          Math.abs(negative_binomial_distribution(k, size, p) / sum - 1),
        ).toBeLessThan(2e-12);
      }
    }
  }
});

test("Negative binomial: geometric distribution is size one", () => {
  for (const p of [1e-10, 0.1, 0.5, 0.9, 1]) {
    for (const x of [-1, 0, 0.9, 1, 5, 100, Infinity]) {
      expect(negative_binomial_density(x, 1, p)).toBe(geometric_density(x, p));
      expect(negative_binomial_log_density(x, 1, p)).toBe(
        geometric_log_density(x, p),
      );
      expect(negative_binomial_distribution(x, 1, p)).toBe(
        geometric_distribution(x, p),
      );
    }
    for (const u of [0, 1e-10, 0.025, 0.5, 0.975, 1]) {
      expect(negative_binomial_quantile(u, 1, p)).toBe(
        geometric_quantile(u, p),
      );
    }
  }
});

test("Negative binomial: exact CDF jumps select the smallest integer", () => {
  for (const size of [2, 3, 5]) {
    let choose = 1;
    let sum = 0;
    for (let k = 0; k <= 12; k++) {
      if (k > 0) choose *= (size + k - 1) / k;
      sum += choose * 0.5 ** (size + k);
      expect(negative_binomial_quantile(sum, size, 0.5)).toBe(k);
      expect(negative_binomial_quantile(sum + 1e-12, size, 0.5)).toBe(k + 1);
    }
  }
  for (const size of [0.25, 0.5, 2, 10, 1000]) {
    for (const p of [0.01, 0.3, 0.9]) {
      let previous = -1;
      for (const u of [1e-20, 1e-10, 0.025, 0.25, 0.5, 0.9, 1 - 1e-10]) {
        const k = negative_binomial_quantile(u, size, p);
        expect(Number.isInteger(k)).toBe(true);
        expect(k).toBeGreaterThanOrEqual(previous);
        expect(
          negative_binomial_distribution(k, size, p) + 1e-13,
        ).toBeGreaterThanOrEqual(u);
        if (k > 0)
          expect(negative_binomial_distribution(k - 1, size, p)).toBeLessThan(
            u + 1e-13,
          );
        previous = k;
      }
    }
  }
});

test("Negative binomial: tiny probabilities and large counts remain finite", () => {
  expect(negative_binomial_quantile(0.5, 2, Number.MIN_VALUE)).toBe(Infinity);
  const q = negative_binomial_quantile(0.5, 2, 1e-20);
  expect(Number.isFinite(q)).toBe(true);
  expect(
    Math.abs(negative_binomial_distribution(q, 2, 1e-20) / 0.5 - 1),
  ).toBeLessThan(2e-12);
  // At equal large size and count with p=1/2, the central mass approaches
  // 1/sqrt(4*pi*size); the Stirling correction is below double precision here.
  const expected = -0.5 * (Math.log(4 * Math.PI) + Math.log(1e200));
  expect(
    Math.abs(negative_binomial_log_density(1e200, 1e200, 0.5) - expected),
  ).toBeLessThan(2e-13);
  expect(negative_binomial_log_density(20, 1000, 1e-10)).toBeCloseTo(
    -22929.8426639291,
    8,
  );
});

test("Negative binomial: support, degenerate case, and validation", () => {
  for (const x of [-Infinity, -1, 0.5, Infinity]) {
    expect(negative_binomial_density(x, 3, 0.4)).toBe(0);
    expect(negative_binomial_log_density(x, 3, 0.4)).toBe(-Infinity);
  }
  expect(negative_binomial_distribution(-0.1, 3, 0.4)).toBe(0);
  expect(negative_binomial_distribution(Infinity, 3, 0.4)).toBe(1);
  expect(negative_binomial_quantile(0, 3, 0.4)).toBe(0);
  expect(negative_binomial_quantile(1, 3, 0.4)).toBe(Infinity);
  expect(negative_binomial_density(0, 3, 1)).toBe(1);
  expect(negative_binomial_density(1, 3, 1)).toBe(0);
  expect(negative_binomial_distribution(0, 3, 1)).toBe(1);
  expect(negative_binomial_quantile(1, 3, 1)).toBe(0);
  expect(negative_binomial_random(3, 3, 1)).toEqual([0, 0, 0]);
  for (const [size, p] of [
    [0, 0.4],
    [-1, 0.4],
    [NaN, 0.4],
    [Infinity, 0.4],
    [3, 0],
    [3, -1],
    [3, 2],
    [3, NaN],
    [3, Infinity],
  ]) {
    for (const f of [
      negative_binomial_density,
      negative_binomial_log_density,
      negative_binomial_distribution,
    ])
      expect(f(0, size, p)).toBeNaN();
    for (const u of [0, 0.5, 1])
      expect(negative_binomial_quantile(u, size, p)).toBeNaN();
    expect(negative_binomial_random(2, size, p)).toEqual([NaN, NaN]);
  }
  for (const f of [
    negative_binomial_density,
    negative_binomial_log_density,
    negative_binomial_distribution,
  ])
    expect(f(NaN, 3, 0.4)).toBeNaN();
  for (const u of [-1, 2, NaN, Infinity])
    expect(negative_binomial_quantile(u, 3, 0.4)).toBeNaN();
  for (const N of [-1, 0.5, NaN, Infinity])
    expect(() => negative_binomial_random(N, 3, 0.4)).toThrow(RangeError);
  expect(negative_binomial_random(0, 3, 0.4)).toEqual([]);
});

test("Negative binomial: random gamma-Poisson mixture follows the CDF and moments", () => {
  let state = 712;
  const random = vi.spyOn(Math, "random").mockImplementation(() => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return (state + 0.5) / 4294967296;
  });
  try {
    const N = 30000;
    for (const [size, p] of [
      [0.5, 0.3],
      [2.5, 0.6],
      [5, 0.4],
      [100, 0.9],
    ]) {
      const samples = negative_binomial_random(N, size, p);
      expect(samples).toHaveLength(N);
      expect(samples.every(k => Number.isInteger(k) && k >= 0)).toBe(true);
      for (const k of [
        0,
        negative_binomial_quantile(0.5, size, p),
        negative_binomial_quantile(0.9, size, p),
      ]) {
        const fraction = samples.filter(x => x <= k).length / N;
        expect(
          Math.abs(fraction - negative_binomial_distribution(k, size, p)),
        ).toBeLessThan(0.012);
      }
    }
    const moments = new Welford();
    for (const k of negative_binomial_random(N, 5, 0.4)) moments.update(k);
    expect(Math.abs(moments.mean() / 7.5 - 1)).toBeLessThan(0.025);
    expect(Math.abs(moments.var() / 18.75 - 1)).toBeLessThan(0.06);
    expect(negative_binomial_random(3, 2, Number.MIN_VALUE)).toEqual([
      Infinity,
      Infinity,
      Infinity,
    ]);
  } finally {
    random.mockRestore();
  }
});
