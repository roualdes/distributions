export function absolute_error(actual: number, expect: number): number {
  return Math.abs(actual - expect);
}

export function relative_error(actual: number, expect: number): number {
  let abs_error = absolute_error(actual, expect);
  return abs_error / expect;
}
