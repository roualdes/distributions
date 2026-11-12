import { expect, test, vi } from "vitest";
import {
  geometric_density,
  geometric_distribution,
  geometric_log_density,
  geometric_quantile,
  geometric_random,
} from "../src/index";
import { Welford } from "../src/welford";

test("Geometric: probability masses and cumulative sums start at zero", () => {
  for (const p of [0.01, 0.2, 0.5, 0.9, 1]) {
    let cumulative = 0;
    for (let k = 0; k <= 50; k++) {
      const expected = p * (1 - p) ** k;
      cumulative += expected;
      expect(geometric_density(k, p)).toBeCloseTo(expected, 14);
      expect(geometric_distribution(k, p)).toBeCloseTo(cumulative, 13);
      expect(geometric_distribution(k + 0.9, p)).toBe(
        geometric_distribution(k, p),
      );
    }
  }
  expect(geometric_log_density(5, 0.5)).toBeCloseTo(-6 * Math.LN2, 14);
  expect(geometric_density(0, 0.3)).toBeCloseTo(0.3, 15);
});

test("Geometric: exact quantile jumps and smallest-integer inequalities", () => {
  for (let k = 0; k < 20; k++) {
    const jump = 1 - 2 ** (-k - 1);
    expect(geometric_quantile(jump, 0.5)).toBe(k);
    expect(geometric_quantile(jump + 1e-12, 0.5)).toBe(k + 1);
  }
  for (const p of [1e-10, 0.01, 0.3, 0.5, 0.9]) {
    let previous = -1;
    for (const u of [1e-20, 1e-10, 0.025, 0.25, 0.5, 0.9, 0.975, 1 - 1e-10]) {
      const k = geometric_quantile(u, p);
      expect(Number.isInteger(k)).toBe(true);
      expect(k).toBeGreaterThanOrEqual(previous);
      if (u > 0.5) {
        // The CDF itself can round identically at many adjacent counts when
        // the tail is small; compare log-survival to retain those differences.
        expect((k + 1) * Math.log1p(-p)).toBeLessThanOrEqual(Math.log1p(-u));
        if (k > 0) expect(k * Math.log1p(-p)).toBeGreaterThan(Math.log1p(-u));
      } else {
        expect(geometric_distribution(k, p)).toBeGreaterThanOrEqual(u);
        if (k > 0) expect(geometric_distribution(k - 1, p)).toBeLessThan(u);
      }
      previous = k;
    }
  }
});

test("Geometric: small probabilities and large counts retain precision", () => {
  expect(geometric_distribution(0, 1e-200)).toBe(1e-200);
  expect(geometric_log_density(10000, 0.5)).toBeCloseTo(-10001 * Math.LN2, 10);
  expect(geometric_density(10000, 0.5)).toBe(0);
  const k = geometric_quantile(0.5, 1e-20);
  expect(Math.abs(geometric_distribution(k, 1e-20) / 0.5 - 1)).toBeLessThan(
    1e-14,
  );
  expect(geometric_quantile(0.5, Number.MIN_VALUE)).toBe(Infinity);
});

test("Geometric: support, degenerate case, endpoints, and invalid arguments", () => {
  for (const x of [-Infinity, -1, 0.5, Infinity]) {
    expect(geometric_density(x, 0.3)).toBe(0);
    expect(geometric_log_density(x, 0.3)).toBe(-Infinity);
  }
  expect(geometric_distribution(-0.1, 0.3)).toBe(0);
  expect(geometric_distribution(Infinity, 0.3)).toBe(1);
  expect(geometric_quantile(0, 0.3)).toBe(0);
  expect(geometric_quantile(1, 0.3)).toBe(Infinity);
  expect(geometric_density(0, 1)).toBe(1);
  expect(geometric_density(1, 1)).toBe(0);
  expect(geometric_distribution(0, 1)).toBe(1);
  expect(geometric_quantile(1, 1)).toBe(0);
  expect(geometric_random(4, 1)).toEqual([0, 0, 0, 0]);
  for (const p of [0, -1, 2, NaN, Infinity, -Infinity]) {
    for (const f of [
      geometric_density,
      geometric_log_density,
      geometric_distribution,
    ])
      expect(f(0, p)).toBeNaN();
    for (const u of [0, 0.5, 1]) expect(geometric_quantile(u, p)).toBeNaN();
    expect(geometric_random(2, p)).toEqual([NaN, NaN]);
  }
  for (const f of [
    geometric_density,
    geometric_log_density,
    geometric_distribution,
  ])
    expect(f(NaN, 0.3)).toBeNaN();
  for (const u of [-1, 2, NaN, Infinity])
    expect(geometric_quantile(u, 0.3)).toBeNaN();
  for (const N of [-1, 0.5, NaN, Infinity])
    expect(() => geometric_random(N, 0.3)).toThrow(RangeError);
  expect(geometric_random(0, 0.3)).toEqual([]);
});

test("Geometric: random integer counts follow the CDF and moments", () => {
  let state = 711;
  const random = vi.spyOn(Math, "random").mockImplementation(() => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return (state + 0.5) / 4294967296;
  });
  try {
    const N = 40000;
    const p = 0.3;
    const samples = geometric_random(N, p);
    expect(samples).toHaveLength(N);
    expect(samples.every(k => Number.isInteger(k) && k >= 0)).toBe(true);
    const moments = new Welford();
    for (const k of samples) moments.update(k);
    expect(Math.abs(moments.mean() / ((1 - p) / p) - 1)).toBeLessThan(0.025);
    expect(Math.abs(moments.var() / ((1 - p) / (p * p)) - 1)).toBeLessThan(
      0.06,
    );
    for (const k of [0, 2, 5]) {
      const fraction = samples.filter(x => x <= k).length / N;
      expect(Math.abs(fraction - geometric_distribution(k, p))).toBeLessThan(
        0.012,
      );
    }
  } finally {
    random.mockRestore();
  }
});
