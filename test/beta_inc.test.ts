import { expect, test } from "vitest";
import { bratio, bratio_complement } from "../src/beta_inc_fn";

test("Incomplete beta: power series includes its correction", () => {
  const expected = 0.99 ** 10 + 10 * 0.01 * 0.99 ** 9;
  expect(Math.abs(bratio(9, 2, 0.99) - expected)).toBeLessThan(1e-14);
});

test("Incomplete beta: recurrence uses integer term indices", () => {
  let expected = 0;
  let choose = 1;
  for (let k = 0; k <= 3; k++) {
    if (k > 0) choose *= (11 - k + 1) / k;
    expected += choose * 0.3 ** k * 0.7 ** (11 - k);
  }
  expect(Math.abs(bratio(8, 4, 0.7) - expected)).toBeLessThan(1e-14);
});

test("Incomplete beta: retains the asymptotic correction", () => {
  const expected = -Math.expm1(130 * Math.log(0.99));
  expect(Math.abs(bratio(1, 130, 0.01) - expected)).toBeLessThan(1e-14);
});

// R 4.4.1 pbeta reference values. These cover both power series and the
// gamma/error-function expansions, including the shapes used by Student's t.
test.each([
  [0.25, 0.5, 0.17, 0.49876020278406041],
  [1.5, 15, 0.1, 0.63996545556559103],
  [1.5, 20, 0.2, 0.97115898749431606],
  [0.01, 180, 0.01, 0.99934707745423834],
  [0.5, 1000, 0.1, 1],
  [180, 190, 0.5, 0.69864000975794682],
  [1000, 1000, 0.49, 0.18555265943151175],
  [10000, 20000, 0.33, 0.11020551912741867],
  [50, 0.5, 0.66589330612500064, 1.999999999999997e-10],
  [15, 0.5, 0.9, 0.077858670314667824],
])("Incomplete beta: reference a=%s b=%s x=%s", (a, b, x, expected) => {
  expect(Math.abs(bratio(a, b, x) / expected - 1)).toBeLessThan(2e-12);
});

test("Incomplete beta: exact distributions and double-precision tails", () => {
  for (const x of [1e-100, 1e-10, 0.01, 0.25, 0.5, 0.9]) {
    const arcsine = (2 / Math.PI) * Math.asin(Math.sqrt(x));
    expect(Math.abs(bratio(0.5, 0.5, x) / arcsine - 1)).toBeLessThan(1e-13);
    const power = x ** 3;
    expect(Math.abs(bratio(3, 1, x) / power - 1)).toBeLessThan(2e-13);
    const survival = -Math.expm1(0.001 * Math.log(x));
    expect(
      Math.abs(bratio_complement(0.001, 1, x) / survival - 1),
    ).toBeLessThan(1e-13);
  }
  // Exercises fpser's exponent cutoff: this tail fits in double precision.
  expect(
    Math.abs(bratio(15, 1e-16, 0.001) / 6.672922554580455e-63 - 1),
  ).toBeLessThan(2e-13);
});

// Independently computed with Python Decimal (80 digits) using the finite
// binomial sum I_x(a,b) = sum(k=a..a+b-1) choose(a+b-1,k) x^k (1-x)^(a+b-1-k).
test.each([
  [400, 30, 0.15, 2.3416586036110655e-287],
  [2000, 20, 0.7, 8.845087492090307e-275],
  [359, 20, 0.12, 1.153110369914359e-300],
])(
  "Incomplete beta: avoids recurrence underflow a=%s b=%s",
  (a, b, x, expected) => {
    expect(Math.abs(bratio(a, b, x) / expected - 1)).toBeLessThan(2e-12);
  },
);

test("Incomplete beta: boundaries and nonfinite arguments", () => {
  expect(bratio(2, 3, 0)).toBe(0);
  expect(bratio(2, 3, 1)).toBe(1);
  expect(bratio(0, 3, 0.5)).toBe(1);
  expect(bratio(2, 0, 0.5)).toBe(0);
  for (const bad of [NaN, Infinity, -Infinity]) {
    expect(bratio(bad, 1, 0.5)).toBeNaN();
    expect(bratio(1, bad, 0.5)).toBeNaN();
    expect(bratio(1, 1, bad)).toBeNaN();
  }
});

test("Incomplete beta: tiny-shape complement retains the digamma correction", () => {
  // R 4.4.1 pbeta(0.1, 1e-16, 0.25, lower.tail = FALSE).
  const expected = 5.8743259266780273e-16;
  expect(
    Math.abs(bratio_complement(1e-16, 0.25, 0.1) / expected - 1),
  ).toBeLessThan(1e-13);
});

test("Incomplete beta: huge shape times tiny x does not overflow internally", () => {
  const bx = Number.MAX_VALUE * Number.MIN_VALUE;
  const expected = 0.5 * bx * bx; // I_x(2,b) ~ (b*x)^2/2 in this regime.
  expect(
    Math.abs(bratio(2, Number.MAX_VALUE, Number.MIN_VALUE) / expected - 1),
  ).toBeLessThan(1e-12);
});
