export class Welford {
  n: number = 0;
  m: number = 0;
  v: number = 0;

  update(x: number): void {
    this.n += 1;
    const a = 1 / this.n;
    const d = x - this.m;
    this.m += d * a;

    this.v += -this.v * a + d * d * a * (1.0 - a);
  }

  mean(): number {
    return this.m;
  }

  var(): number {
    return (this.v * this.n) / (this.n - 1.0);
  }

  std(): number {
    return Math.sqrt(this.var());
  }

  reset(): void {
    this.n = 0;
    this.m = 0;
    this.v = 0;
  }
}
