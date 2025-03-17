import { relative_error } from "./tools";
import { Welford } from "../src/welford";

import { expect, test } from "vitest";

import {
  gamma_density,
  gamma_distribution,
  gamma_quantile,
  gamma_random,
} from "../src/gamma";

test("Gamma: density 01", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(gamma_density(1, 1.5, 3.5), 0.22311380316036858),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_density(10, 10.5, 3.52), 7.892154824743942e-7),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_density(1, 10.5, 3.52), 0.014310169673460596),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_density(3.0, 15.0, 3.52), 0.22453261650890766),
  ).toBeCloseTo(tolerance);
});

test("Gamma: distribution", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(gamma_distribution(0.74, 1.65, 3.14), 0.7645943602808132),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(
      gamma_distribution(0.875, 24.0, 10.05),
      1.7068149854885828e-5,
    ),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_distribution(1.0, 11.0, 0.99), 9.078429600279503e-9),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(
      gamma_distribution(130.0, 13.0, 0.03),
      0.00021564345567396287,
    ),
  ).toBeCloseTo(tolerance);
});

test("Gamma: quantile", () => {
  const tolerance = 1e-12;

  expect(
    relative_error(gamma_quantile(0.25, 1.65, 3), 0.2366370744607247),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_quantile(0.05, 113.0, 0.15), 640.6776961186549),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_quantile(0.95, 3.51, 10.0), 0.7048131743374495),
  ).toBeCloseTo(tolerance);

  expect(
    relative_error(gamma_quantile(0.875, 24.0, 0.95), 31.26982740894499),
  ).toBeCloseTo(tolerance);
});

test("Gamma: random", () => {
  const N = 10_000;

  let a = 1.0;
  let b = 2.5;
  let x = gamma_random(N, a, b);

  const w = new Welford();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / b, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(a) / b, 1);

  a = 0.7;
  b = 3.2;
  x = gamma_random(N, a, b);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / b, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(a) / b, 1);

  a = 5.0;
  b = 5.2;
  x = gamma_random(N, a, b);

  w.reset();
  for (let n = 0; n < N; n++) {
    w.update(x[n]);
  }

  expect(w.mean()).toBeCloseTo(a / b, 1);
  expect(w.std()).toBeCloseTo(Math.sqrt(a) / b, 1);
});
