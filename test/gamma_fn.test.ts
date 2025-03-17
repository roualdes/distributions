import { gamma } from "../src/gamma_fn";
import { absolute_error, relative_error } from "./tools";

import { expect, test } from "vitest";

test("Gamma function", () => {
  const tolerance = 1e-12;

  // https://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp
  // gr = GammaRegularized
  // with arguments reversed and then 1.0 - gr
  expect(relative_error(gamma(2, 1), 1 - 0.13533528)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(2, 10), 1 - 0.9999535)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(3, 10), 1 - 0.99889751)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(4, 10), 1 - 0.99186776)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(5, 10), 1 - 0.96817194)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(1.5, 15), 1 - 0.9999999999176)).toBeCloseTo(
    tolerance,
  );
  expect(relative_error(gamma(2.5, 15), 1 - 0.99999993084686)).toBeCloseTo(
    tolerance,
  );
  expect(relative_error(gamma(3.5, 15), 1 - 0.99999574)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(4.5, 15), 1 - 0.99992634)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(11, 20), 1 - 0.99071054)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(12, 40), 1 - 0.999999999844129)).toBeCloseTo(
    tolerance,
  );
  expect(relative_error(gamma(30, 50), 1 - 0.99948111)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(39, 33), 1 - 0.14815725)).toBeCloseTo(tolerance);
  expect(relative_error(gamma(39, 32), 1 - 0.11214826)).toBeCloseTo(tolerance);
});
