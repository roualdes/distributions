import { expect, test } from "vitest";
import { Welford } from "../src/welford";

import {
  exponential_density,
  exponential_distribution,
  exponential_quantile,
  exponential_random,
} from "../src/exponential";

test("Exponential: density", () => {
  const d = exponential_density(0.453, 1.3);
  expect(d).toBeCloseTo(0.7214);
});

test("Exponential: distribution", () => {
  const p = exponential_distribution(2.5, 0.8);
  expect(p).toBeCloseTo(0.8646);
});

test("Exponential: quantile", () => {
  const q = exponential_quantile(0.1, 0.8);
  expect(q).toBeCloseTo(0.1317);
});

test("Exponential: random", () => {
  const N = 10_000;
  const lambda = 3.1;
  const x = exponential_random(N, lambda);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(1 / lambda, 1);
  expect(w.std()).toBeCloseTo(1 / lambda, 1);
});
