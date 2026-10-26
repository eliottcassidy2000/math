\\ Audit of HYP-9391: Fourier coefficients at primes p < N of 8.4.a.a, 8.6.a.a, 4.6.a.a and |b(p)|^2 for 32.3.c.a
\\ (b = q + 4i q^3 + 2 q^5 - 8i q^7 - 7 q^9 - ...).  Output lines: p a84 a86 a46 |b|^2.  N is set by the caller.
default(parisize, 1500000000);
f84 = mfeigenbasis(mfinit([8,4,1],0))[1];
f86 = mfeigenbasis(mfinit([8,6,1],0))[1];
f46 = mfeigenbasis(mfinit([4,6,1],0))[1];
fb  = mfeigenbasis(mfinit([32,3,-4],0))[1];
c84 = mfcoefs(f84, N); c86 = mfcoefs(f86, N); c46 = mfcoefs(f46, N); cb = mfcoefs(fb, N);
{
forprime(p = 3, N,
  bb = lift(cb[p + 1]);
  nb = if (type(bb) == "t_POL", polcoef(bb, 0)^2 + polcoef(bb, 1)^2, bb^2);
  print(p, " ", c84[p + 1], " ", c86[p + 1], " ", c46[p + 1], " ", nb));
}
quit;
