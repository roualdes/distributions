import { expect, test, vi } from "vitest";
import {
  rayleigh_density,
  rayleigh_distribution,
  rayleigh_log_density,
  rayleigh_quantile,
  rayleigh_random,
} from "../src/index";
import { Welford } from "../src/welford";

// R 4.4.1 Weibull(shape=2, scale=sqrt(2)*s), equivalent to Rayleigh(s).
test.each([
  [0.453, 1.3, 0.25225762195505719, -1.3773044043275764, 0.058906443478925603],
  [0.1, 0.5, 0.39207946932270205, -0.93629073187415512, 0.019801326693244699],
  [1, 1, 0.60653065971263342, -0.5, 0.39346934028736652],
  [2.5, 3, 0.19629118829381004, -1.6281560676842868, 0.29335172214228367],
  [10, 2, 9.3166329301966926e-6, -11.583709268125844, 0.99999627334682795],
  [30, 1, 1.1081649205463026e-194, -446.59880261833774, 1],
])("Rayleigh: density, log-density, and CDF x=%s s=%s", (x, s, d, ld, p) => {
  expect(Math.abs(rayleigh_density(x, s) / d - 1)).toBeLessThan(3e-13);
  expect(Math.abs(rayleigh_log_density(x, s) - ld)).toBeLessThan(2e-13);
  expect(Math.abs(rayleigh_distribution(x, s) / p - 1)).toBeLessThan(1e-14);
});

test.each([
  [0.001, 0.5, 0.022366272974989974],
  [0.1, 1.3, 0.5967566865343471],
  [0.5, 1, 1.1774100225154747],
  [0.9, 3, 6.4378980788680433],
  [0.975, 2, 5.4324060629624773],
  [1 - 1e-10, 1, 6.7861404122225606],
])("Rayleigh: quantile p=%s s=%s", (p, s, expected) => {
  expect(Math.abs(rayleigh_quantile(p, s) / expected - 1)).toBeLessThan(1e-14);
});

test("Rayleigh: scale transformations and CDF/quantile round trips", () => {
  for (const s of [1e-300, 0.5, 1, 3, 1e300]) {
    for (const p of [1e-100, 1e-20, 0.025, 0.5, 0.975, 1 - 1e-10]) {
      const q = rayleigh_quantile(p, s);
      if (q === 0 || !Number.isFinite(q)) continue;
      expect(Math.abs(rayleigh_distribution(q, s) / p - 1)).toBeLessThan(2e-13);
    }
    for (const z of [0.1, 1, 5]) {
      expect(rayleigh_distribution(s * z, s)).toBeCloseTo(
        rayleigh_distribution(z, 1),
        14,
      );
      expect(
        Math.abs(
          rayleigh_log_density(s * z, s) +
            Math.log(s) -
            rayleigh_log_density(z, 1),
        ),
      ).toBeLessThan(3e-13);
    }
  }
});

test("Rayleigh: preserves small probabilities and finite log-densities", () => {
  expect(
    Math.abs(rayleigh_distribution(1e-150, 1) / 0.5e-300 - 1),
  ).toBeLessThan(2e-15);
  expect(rayleigh_distribution(rayleigh_quantile(Number.MIN_VALUE, 1), 1)).toBe(
    Number.MIN_VALUE,
  );
  expect(rayleigh_density(40, 1)).toBe(0);
  expect(rayleigh_log_density(40, 1)).toBe(Math.log(40) - 800);
  // Neither x^2 nor s^2 needs to fit when their ratio is ordinary.
  expect(rayleigh_log_density(1e200, 1e200)).toBeCloseTo(
    -Math.log(1e200) - 0.5,
    12,
  );
  expect(rayleigh_log_density(1e-200, 1e-200)).toBeCloseTo(
    -Math.log(1e-200) - 0.5,
    12,
  );
  expect(rayleigh_log_density(1.5e154, 1)).toBeCloseTo(-1.125e308, -295);
  expect(rayleigh_distribution(1e300, 1e-300)).toBe(1);
  expect(rayleigh_log_density(1e300, 1e-300)).toBe(-Infinity);
});

test("Rayleigh: support, boundaries, and invalid arguments", () => {
  for (const x of [-Infinity, -1, 0]) {
    expect(rayleigh_density(x, 2)).toBe(0);
    expect(rayleigh_log_density(x, 2)).toBe(-Infinity);
    expect(rayleigh_distribution(x, 2)).toBe(0);
  }
  expect(rayleigh_density(Infinity, 2)).toBe(0);
  expect(rayleigh_log_density(Infinity, 2)).toBe(-Infinity);
  expect(rayleigh_distribution(Infinity, 2)).toBe(1);
  expect(rayleigh_quantile(0, 2)).toBe(0);
  expect(rayleigh_quantile(1, 2)).toBe(Infinity);
  for (const s of [0, -1, NaN, Infinity, -Infinity]) {
    for (const f of [
      rayleigh_density,
      rayleigh_log_density,
      rayleigh_distribution,
    ]) {
      expect(f(0, s)).toBeNaN();
      expect(f(Infinity, s)).toBeNaN();
    }
    for (const p of [0, 0.5, 1]) expect(rayleigh_quantile(p, s)).toBeNaN();
    expect(rayleigh_random(3, s)).toEqual([NaN, NaN, NaN]);
  }
  for (const f of [
    rayleigh_density,
    rayleigh_log_density,
    rayleigh_distribution,
  ]) {
    expect(f(NaN, 2)).toBeNaN();
  }
  for (const p of [-1, 2, NaN, Infinity, -Infinity])
    expect(rayleigh_quantile(p, 2)).toBeNaN();
  for (const N of [-1, 0.5, NaN, Infinity])
    expect(() => rayleigh_random(N, 2)).toThrow(RangeError);
  expect(rayleigh_random(0, 2)).toEqual([]);
});

test("Rayleigh: random uniforms include zero without producing infinity", () => {
  const random = vi
    .spyOn(Math, "random")
    .mockReturnValueOnce(0)
    .mockReturnValueOnce(0.5)
    .mockReturnValueOnce(1 - Number.EPSILON / 2);
  try {
    const samples = rayleigh_random(3, 2);
    expect(samples[0]).toBe(0);
    expect(samples[1]).toBeCloseTo(2 * Math.sqrt(2 * Math.LN2), 14);
    expect(Number.isFinite(samples[2])).toBe(true);
  } finally {
    random.mockRestore();
  }
});

test("Rayleigh: random draws follow the CDF and theoretical moments", () => {
  let state = 710;
  const random = vi.spyOn(Math, "random").mockImplementation(() => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return (state + 0.5) / 4294967296;
  });
  try {
    const N = 40000;
    const s = 2.5;
    const samples = rayleigh_random(N, s);
    expect(samples).toHaveLength(N);
    expect(samples.every(x => Number.isFinite(x) && x >= 0)).toBe(true);
    const moments = new Welford();
    for (const x of samples) moments.update(x);
    expect(
      Math.abs(moments.mean() / (s * Math.sqrt(Math.PI / 2)) - 1),
    ).toBeLessThan(0.015);
    expect(
      Math.abs(moments.var() / ((s * s * (4 - Math.PI)) / 2) - 1),
    ).toBeLessThan(0.04);
    for (const p of [0.1, 0.5, 0.9]) {
      const cutoff = rayleigh_quantile(p, s);
      const fraction = samples.filter(x => x <= cutoff).length / N;
      expect(Math.abs(fraction - p)).toBeLessThan(0.01);
    }
  } finally {
    random.mockRestore();
  }
});
