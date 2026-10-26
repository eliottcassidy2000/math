\\ exact Sturm spot checks in the range 1000 < q < 4000 (smallest-excess exceptions and nearby Ramanujan primes)
read("spectra_exact.gp");
{
  foreach ([1061, 1063, 1069, 2909, 2917], q,
    my(t0 = getabstime(), r = check(q));
    print(r, "   ms=", getabstime() - t0));
}
