import { expect, test } from "vitest";
import { Welford } from "../src/welford";

import {
  normal_density,
  normal_distribution,
  normal_quantile,
  normal_random,
} from "../src/normal";

test("Normal: density", () => {
  const d = normal_density(-0.2345, 2.589, 20.34);
  expect(d).toBeCloseTo(0.01942561455842912);
});

test("Normal: distribution", () => {
  const p = normal_distribution(-0.2345, 6.589, 20.34);
  expect(p).toBeCloseTo(0.3686345403833694);
});

test("Normal: quantile", () => {
  const p = normal_distribution(-0.2345, 6.589, 20.34);
  const q = normal_quantile(p, 2.4, 4.2);
  expect(q).toBeCloseTo(0.9910176991150441);
});

test("Normal: random", () => {
  const N = 10_000;
  const m = -5.7;
  const s = 0.5;

  const x = normal_random(N, m, s);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(m, 1);
  expect(w.std()).toBeCloseTo(s, 1);
});
