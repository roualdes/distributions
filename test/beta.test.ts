import { relative_error } from "./tools";
import { Welford } from "../src/welford";

import { expect, test } from "vitest";

import {
  beta_density,
  beta_distribution,
  beta_quantile,
  beta_random,
} from "../src/beta";

test("Beta: density", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(beta_density(0.25, 0.5, 0.5), 0.7351051938957226),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_density(0.25, 1.5, 3.5), 1.984784023518451),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_density(0.95, 1.5, 3.5), 0.004439938005136257),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_density(0.543, 31.5, 83.5), 2.860367620082414e-7),
  ).toBeCloseTo(tolerance);
});

test("Beta: distribution", () => {
  const tolerance = 1e-12;

  // https://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp
  // Regularized incomplete beta function
  // BetaRegularized

  expect(
    relative_error(beta_distribution(0.75, 0.5, 0.24), 0.43339455713545905),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(0.33, 1.5, 7.24), 0.8877525021108774),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(0.05, 2.4, 15.0), 0.1145884868689588),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(0.95, 15.0, 2.4), 0.885411513131041),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(0.3, 15.0, 1e-16), 1.331682743548004e-25),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(1e-16, 1e-16, 1), 0.9999999999999963),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_distribution(0.5, 180, 190), 0.6986400097579467),
  ).toBeCloseTo(tolerance);
});

test("Beta: quantile", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(beta_quantile(0.75, 0.5, 0.24), 0.9907914064651424),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_quantile(0.33, 1.5, 7.24), 0.09861696099455533),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_quantile(0.05, 2.4, 15.0), 0.03302727776083206),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(beta_quantile(0.95, 15.0, 2.4), 0.966972722239168),
  ).toBeCloseTo(tolerance);
});

test("Beta: random", () => {
  const N = 100_000;

  let a = 1.0;
  let b = 2.5;
  let x = beta_random(N, a, b);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / (a + b), 1);
  expect(w.std()).toBeCloseTo(
    Math.sqrt((a * b) / ((a + b) ** 2 * (a + b + 1))),
    1,
  );

  a = 0.7;
  b = 3.2;
  x = beta_random(N, a, b);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / (a + b), 1);
  expect(w.std()).toBeCloseTo(
    Math.sqrt((a * b) / ((a + b) ** 2 * (a + b + 1))),
    1,
  );

  a = 5.0;
  b = 5.2;
  x = beta_random(N, a, b);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / (a + b), 1);
  expect(w.std()).toBeCloseTo(
    Math.sqrt((a * b) / ((a + b) ** 2 * (a + b + 1))),
    1,
  );
});
