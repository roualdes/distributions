import { log_gamma } from "../src/log_gamma_fn";
import { relative_error } from "./tools";

import { expect, test } from "vitest";

test("Log(Gamma) function", () => {
  const tolerance = 1e-12;

  expect(relative_error(log_gamma(1e-12), 27.6310211159)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(log_gamma(0.9999), 5.77297915613e-5)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(log_gamma(1.0001), -5.77133422205e-5)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(log_gamma(3.1), 0.787375083274)).toBeCloseTo(tolerance);

  expect(relative_error(log_gamma(6.3), 5.30734288962)).toBeCloseTo(tolerance);

  expect(relative_error(log_gamma(11.9999), 17.5020635801)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(log_gamma(12), 17.5023078459)).toBeCloseTo(tolerance);

  expect(relative_error(log_gamma(12.0001), 17.5025521125)).toBeCloseTo(
    tolerance,
  );

  expect(relative_error(log_gamma(27.4), 62.5755868211)).toBeCloseTo(tolerance);
});
