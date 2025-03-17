import { relative_error } from "./tools";
import { expect, test } from "vitest";

import { beta_inc } from "../src/beta_inc_fn";

test("beta_inc", () => {
  const tolerance = 1e-2;
  // TODO better testing will surely reveal my implementation is imperfect
  // expect(
  //   relative_error(beta_inc(10, 10, 0.5), 0.5))
  //   .toBeCloseTo(tolerance);
});
