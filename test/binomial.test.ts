import { relative_error } from "./tools";
import { Welford } from "../src/welford";

import { expect, test } from "vitest";

import {
  binomial_density,
  binomial_distribution,
  binomial_quantile,
  binomial_random,
} from "../src/binomial";

test("Binomial: density 01", () => {
  const tolerance = 1e-12;

  expect(relative_error(binomial_density(1, 10, 0.5), 0.009765625)).toBeCloseTo(
    tolerance,
  );

  expect(
    relative_error(binomial_density(3, 30, 0.11), 0.23239900721727938),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(binomial_density(27, 30, 0.11), 3.752308635046128e-23),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(binomial_density(7, 18, 0.63), 0.022302706885899883),
  ).toBeCloseTo(tolerance);
});

test("Binomial: distribution", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(binomial_distribution(3.0, 10.0, 0.7), 0.01059207840000001),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(binomial_distribution(10, 11.0, 0.9), 0.68618940391),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(
      binomial_distribution(10, 24.0, 0.875),
      1.2944997659339275e-7,
    ),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(binomial_distribution(1, 11.0, 0.01), 0.9948202825096849),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(binomial_distribution(40, 130.0, 0.3), 0.6173916412888664),
  ).toBeCloseTo(tolerance);
});

test("Binomial: quantile", () => {
  const tolerance = 1e-12;

  expect(relative_error(binomial_quantile(0.01, 10, 0.7), 3)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(binomial_quantile(0.7, 11, 0.9), 11)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(binomial_quantile(1e-7, 24, 0.875), 10)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(binomial_quantile(0.62, 130, 0.3), 41)).toBeCloseTo(
    tolerance,
  );
});

test("Binomial: random", () => {
  const N = 50_000;

  let K = 10;
  let p = 0.75;
  let x = binomial_random(N, K, p);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(K * p, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(K * p * (1 - p)), 1);

  K = 50;
  p = 0.65;
  x = binomial_random(N, K, p);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(K * p, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(K * p * (1 - p)), 1);

  K = 20;
  p = 0.55;
  x = binomial_random(N, K, p);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(K * p, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(K * p * (1 - p)), 1);
});
