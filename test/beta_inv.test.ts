import { expect, test } from "vitest";
import { betaln, bratio, bratio_complement } from "../src/beta_inc_fn";
import { betain, xinbeta, xinbeta_complement } from "../src/beta_inv_fn";
import { beta_quantile } from "../src/beta";

// R 4.4.1 qbeta reference values; no runtime dependency on R.
test.each([
  [0.25, 0.5, 0.05, 1.8463996456027615e-5],
  [0.5, 0.5, 0.5, 0.49999999999999989],
  [1, 0.5, 0.05, 0.097500000000000017],
  [2.5, 0.5, 0.05, 0.43074146813810116],
  [5, 0.5, 0.5, 0.95331275466303356],
  [15, 0.5, 0.05, 0.87794059801520807],
  [50, 0.5, 2e-10, 0.66589330612500064],
  [500000, 0.5, 0.05, 0.99999615854663693],
  [1000, 1000, 1e-20, 0.3975320608971773],
  [10000, 20000, 0.975, 0.33867806945919066],
  [0.001, 0.1, 0.9, 3.1425485275408857e-42],
  [0.001, 20, 1 - 1e-10, 0.49786099662967664],
])("Inverse beta: reference a=%s b=%s p=%s", (a, b, p, expected) => {
  expect(Math.abs(xinbeta(a, b, p) / expected - 1)).toBeLessThan(3e-12);
});

test("Inverse beta: exact quantiles, including extreme tails", () => {
  for (const p of [1e-150, 1e-20, 1e-10, 0.025, 0.25, 0.5, 0.9]) {
    const expected = Math.sin((Math.PI * p) / 2) ** 2;
    expect(Math.abs(xinbeta(0.5, 0.5, p) / expected - 1)).toBeLessThan(3e-13);
    expect(Math.abs(xinbeta(3, 1, p) / p ** (1 / 3) - 1)).toBeLessThan(2e-13);
  }
  // An upper probability can still have a very small root in a skewed beta.
  const expected = Math.exp(Math.log(0.9) / 0.001);
  expect(Math.abs(xinbeta(0.001, 1, 0.9) / expected - 1)).toBeLessThan(3e-12);
  expect(xinbeta(0.001, 1, 0.1)).toBe(0); // True root rounds below MIN_VALUE.
});

test("Inverse beta: uses the CDF consistently across shapes and tails", () => {
  for (const a of [0.25, 0.5, 1, 5, 50, 1000]) {
    for (const b of [0.5, 1, 3, 100]) {
      for (const p of [1e-20, 1e-10, 0.025, 0.25, 0.5, 0.9]) {
        const x = xinbeta(a, b, p);
        expect(x).toBeGreaterThanOrEqual(0);
        expect(x).toBeLessThanOrEqual(1);
        // Near one, quantile rounding can dominate probability accuracy.
        if (x > 1 - 1e-6) continue;
        const tail = p > 0.5 ? bratio_complement(a, b, x) : bratio(a, b, x);
        const target = Math.min(p, 1 - p);
        expect(Math.abs(tail / target - 1)).toBeLessThan(1e-9);
      }
    }
  }
});

test("Inverse beta: fixes the Student t quantile regressions", () => {
  const cauchyBeta = beta_quantile(0.5, 0.5, 0.5);
  expect(Math.sqrt((1 - cauchyBeta) / cauchyBeta)).toBeCloseTo(1, 14);
  const t5Beta = beta_quantile(0.05, 2.5, 0.5);
  expect(Math.sqrt((5 * (1 - t5Beta)) / t5Beta)).toBeCloseTo(
    2.5705818356363155,
    13,
  );
});

test("Incomplete beta AS 63: integer recurrence counter", () => {
  expect(Math.abs(betain(0.5, 0.5, 0.5, betaln(0.5, 0.5)) - 0.5)).toBeLessThan(
    1e-14,
  );
  expect(
    Math.abs(betain(0.43, 2.5, 0.5, betaln(2.5, 0.5)) - bratio(2.5, 0.5, 0.43)),
  ).toBeLessThan(1e-14);
  expect(betain(NaN, 1, 1, 0)).toBeNaN();
});

test("Inverse beta: endpoints and invalid parameters terminate", () => {
  expect(xinbeta(0.5, 3, 0)).toBe(0);
  expect(xinbeta(0.5, 3, 1)).toBe(1);
  for (const bad of [0, -1, NaN, Infinity, -Infinity]) {
    expect(xinbeta(bad, 1, 0.5)).toBeNaN();
    expect(xinbeta(1, bad, 0.5)).toBeNaN();
  }
  for (const bad of [-1, 2, NaN, Infinity, -Infinity]) {
    expect(xinbeta(1, 1, bad)).toBeNaN();
  }
});

test("Inverse beta: complementary probabilities preserve small roots", () => {
  // Beta(a,1) has survival 1-x^a. Supplying this survival directly preserves
  // information that is lost in the subtraction 1 - survival.
  for (const a of [0.5, 1, 5, 1000]) {
    for (const p of [1e-100, 1e-10, 0.025, 0.5, 0.9]) {
      const x = xinbeta_complement(a, 1, p);
      const expected = Math.exp(Math.log1p(-p) / a);
      expect(Math.abs(x / expected - 1)).toBeLessThan(1e-12);
    }
  }
  const x = xinbeta_complement(0.5, 500000, 1e-100);
  expect(Math.abs(bratio_complement(0.5, 500000, x) / 1e-100 - 1)).toBeLessThan(
    1e-11,
  );
  expect(xinbeta_complement(2, 3, 0)).toBe(1);
  expect(xinbeta_complement(2, 3, 1)).toBe(0);
  expect(xinbeta_complement(2, 3, NaN)).toBeNaN();
  expect(xinbeta_complement(Infinity, 3, 0.5)).toBeNaN();
});
