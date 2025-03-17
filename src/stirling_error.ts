import { log_gamma } from "./log_gamma_fn";
import { PI2 } from "./constants";

export function stirling_error(n: number): number {
  if (n < 16) {
    return log_gamma(n + 1.0) + n - n * Math.log(n) - 0.5 * Math.log(PI2 * n);
  }

  const nn = n * n;
  const S0 = 0.083333333333333333333; // 1/12
  const S1 = 0.00277777777777777777778; // 1/360
  const S2 = 0.00079365079365079365079365; // 1/1260
  const S3 = 0.000595238095238095238095238; // 1/1680
  const S4 = 0.0008417508417508417508417508; // 1/1188

  if (n > 500) {
    return (S0 - S1 / nn) / n;
  }

  if (n > 80) {
    return (S0 - (S1 - S2 / nn) / nn) / n;
  }

  if (n > 35) {
    return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  }

  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}
