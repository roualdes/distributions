import { expect, test } from "vitest";
import { Welford } from "../src/welford";

import {
  uniform_density,
  uniform_distribution,
  uniform_quantile,
  uniform_random,
} from "../src/uniform";

test("Uniform: density", () => {
  const d = uniform_density(0.453, 0.0, 1.0);
  expect(d).toBeCloseTo(1.0);
});

test("Unfiorm: distribution", () => {
  const p = uniform_distribution(2.5, 2.0, 4.0);
  expect(p).toBeCloseTo(0.25);
});

test("Uniform: quantile", () => {
  const q = uniform_quantile(0.25, 2.0, 4.0);
  expect(q).toBeCloseTo(2.5);
});

test("Uniform: random", () => {
  const N = 10_000;
  const a = 2.0;
  const b = 7.7;

  const x = uniform_random(N, a, b);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(0.5 * (a + b), 1);

  const d = b - a;
  expect(w.std()).toBeCloseTo(Math.sqrt((d * d) / 12), 1);
});
