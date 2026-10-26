\\ driver: exact Ramanujan census for primes 3 <= q < 1000 (A and the twisted A_chi), plus the two exceptional factors
read("spectra_exact.gp");
{
  my(badA = List(), badC = List(), n = 0, t0 = getabstime());
  forprime (q = 3, 999,
    my(r = check(q)); n++;
    if (r[3] != 1 || r[4], print("TRIVIAL-EIGENVALUE ANOMALY ", r));
    if (r[5][1] > 0 || r[5][2], listput(badA, [q, r[5]]));
    if (r[7], print("TWISTED +-4 ANOMALY ", r));
    if (r[8][1] > 0 || r[8][2], listput(badC, [q, r[8]]));
    if (q <= 211 || r[5][1] > 0 || r[8][1] > 0, print(r)));
  print("primes 3 <= q < 1000: ", n, " (", n - 1, " with q >= 5)");
  print("A (G_q): primes with an eigenvalue |mu| >= 2 sqrt 3 (exact Sturm count, [count, exactly-on-bound]): ", badA);
  print("A_chi (B_chi): primes with |mu| >= 2 sqrt 3: ", badC);
  print("elapsed ms: ", getabstime() - t0);
}
print("q=331 A exceptional factor [deg, mult, mu, lambda, 3/lambda]: ", exception(331, 1));
print("q=251 A_chi exceptional factor [deg, mult, mu, lambda, 3/lambda]: ", exception(251, 2));
