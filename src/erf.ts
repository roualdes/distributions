export function erf(x: number): number {
  // John D. Cook
  // https://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;

  const sign = x === 0 ? 1 : Math.sign(x);
  const absx = Math.abs(x);

  if (absx < 1e-3) {
    return 1.1283791670955126 * x;
  }

  //  A & S 7.1.26
  const t = 1.0 / (1.0 + p * absx);
  const y =
    1.0 -
    ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-absx * absx);
  return sign * y;
}
