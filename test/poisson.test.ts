import { relative_error } from "./tools";
import { Welford } from "../src/welford";

import { expect, test } from "vitest";

import {
  poisson_density,
  poisson_distribution,
  poisson_quantile,
  poisson_random,
} from "../src/poisson";

test("Poisson: density", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(poisson_density(1, 5), 0.03368973499542734),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_density(11, 17.54), 0.029217600644745594),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_density(13, 1.1), 1.845443105561357e-10),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_density(25, 50.1), 3.524876655343783e-5),
  ).toBeCloseTo(tolerance);
});

test("Poisson: distribution", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(poisson_distribution(1, 5), 0.04042768199451279),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(1, 3.7), 0.11620057441059513),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(2, 2), 0.6766764161830634),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(14, 7.25), 0.9922741119889988),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(5, 39), 9.91933214072047e-12),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(32, 39), 0.148157249294968),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(poisson_distribution(33, 39), 0.1907133228328563),
  ).toBeCloseTo(tolerance);
});

test("Poisson: random", () => {
  const N = 100_000;

  let lambda = 10;
  let x = poisson_random(N, lambda);

  const w = new Welford();
  for (let n = 0; n < N; ++n) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(lambda, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(lambda), 1);

  lambda = 50;
  x = poisson_random(N, lambda);

  w.reset();
  for (let n = 0; n < N; ++n) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(lambda, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(lambda), 1);

  lambda = 21.5;
  x = poisson_random(N, lambda);

  w.reset();
  for (let n = 0; n < N; ++n) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(lambda, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(lambda), 1);
});

test("Poisson: quantile", () => {
  const tolerance = 1e-12;

  expect(relative_error(poisson_quantile(0.75, 2), 3)).toBeCloseTo(tolerance);

  expect(relative_error(poisson_quantile(0.5, 3.7), 4)).toBeCloseTo(tolerance);

  expect(relative_error(poisson_quantile(0.77, 7.25), 9)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(poisson_quantile(0.15, 39), 33)).toBeCloseTo(tolerance);
});
