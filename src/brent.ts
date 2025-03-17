// adapted from https://en.wikipedia.org/wiki/Brent's_method
// by Edward A. Roualdes on 2024-06-11

export function brent(
  a: number,
  b: number,
  u: number,
  tol: number,
  f: (x: number, u: number, ...theta: number[]) => number,
  ...theta: number[]
): number {
  let fa = f(a, u, ...theta);
  let fb = f(b, u, ...theta);

  if (fa * fb >= 0.0) {
    if (fa < fb) {
      return a;
    } else {
      return b;
    }
  }

  if (Math.abs(fa) < Math.abs(fb)) {
    let tmp = a;
    a = b;
    b = tmp;
    tmp = fa;
    fa = fb;
    fb = tmp;
  }

  let c = a;
  let fc = fa;
  let mflag = true;
  let s = 0.0;
  let fs = 0.0;
  let d = 0.0;
  // int tries = 0;

  while (Math.abs(fb) > tol && Math.abs(a - b) > tol) {
    if (Math.abs(fa - fc) > tol && Math.abs(fb - fc) > tol) {
      s = (a * fb * fc) / ((fa - fb) * (fa - fc));
      s += (b * fa * fc) / ((fb - fa) * (fb - fc));
      s += (c * fa * fb) / ((fc - fa) * (fc - fb));
    } else {
      s = b - (fb * (b - a)) / (fb - fa);
    }

    if (
      s < (3.0 * a + b) / 4.0 ||
      s > b ||
      (mflag && Math.abs(s - b) >= 0.5 * Math.abs(b - c)) ||
      (!mflag && Math.abs(s - b) >= 0.5 * Math.abs(c - d)) ||
      (mflag && Math.abs(b - c) < tol) ||
      (!mflag && Math.abs(c - d) < tol)
    ) {
      s = a + 0.5 * (b - a);
      mflag = true;
    } else {
      mflag = false;
    }

    fs = f(s, u, ...theta);
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (Math.abs(fa) < Math.abs(fb)) {
      let tmp = a;
      a = b;
      b = tmp;
      tmp = fa;
      fa = fb;
      fb = tmp;
    }

    // tries += 1;
    // if (tries > maxtries) {
    //   printf("Brent method: max tries\n") ;
    //   break;
    // }
  }
  return b;
}
