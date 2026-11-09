import { expect, test } from "vitest";
import { bratio } from "../src/beta_inc_fn";

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
