\\ Audit E (2026-10-14), THM-4656 part 3, computed in PARI/GP (independent of the Python scripts).
\\ pi(p) = order of phi in (Z[phi]/p)^*:  split p: lcm of znorder of the two roots of x^2-x-1 mod p;
\\                                        inert p: fforder of a root of x^2-x-1 in F_(p^2) (ffgen).
\\ (1) the primes p < 10^6 (p != 2, 5) with pi(p) a power of 2;  (2) p | L_(2^j) for some j <= 25 (L_2 = 3,
\\ L_(2n) = L_n^2 - 2), and then pi(p) = 2^(j+2);  (3) Mersenne primes M_q: phi^(2^q) = -1 in (Z/M_q)[x]/(x^2-x-1)
\\ for q = 3 mod 4, so pi(M_q) = 2^(q+1); for q = 1 mod 4 M_q is split and pi(M_q) is not a power of 2;
\\ (4) Lucas-Lehmer with seed 3 (s_0 = 3) for q = 3 mod 4, q <= 300: M_q | s_(q-2) iff M_q prime.
pis(p) = {
  my(k = kronecker(5, p), r, a);
  if (k == 1,
    r = polrootsmod(x^2 - x - 1, p);
    lcm(znorder(r[1]), znorder(r[2])),
    a = ffgen(Mod(1, p) * (x^2 - x - 1), 'a);   \\ irreducible over F_p since (5/p) = -1
    fforder(a));
}
ispow2(n) = (n > 0) && (n == 2^valuation(n, 2));
lucasdiv(p) = {
  my(L = Mod(3, p));
  for (j = 1, 25, if (L == 0, return(j)); L = L^2 - 2);
  0;
}
{
  my(G = List(), bad = 0, cnt = 0);
  forprime(p = 3, 10^6, if (p == 5, next);
    cnt++;
    my(pp = pis(p), j = lucasdiv(p));
    if (ispow2(pp), listput(G, [p, pp]));
    if (ispow2(pp) != (j > 0), bad++; print("MISMATCH ", p, " ", pp, " ", j));
    if (j > 0 && pp != 2^(j + 2), bad++; print("ORDER MISMATCH ", p, " ", pp, " ", j));
  );
  print("(1)-(2) primes checked: ", cnt, "; mismatches between [pi(p) = 2^k] and [p | L_(2^j), pi = 2^(j+2)]: ", bad);
  print("    golden-rigid primes p < 10^6 [p, pi(p)]: ", Vec(G));
  \\ (3) Mersenne primes
  foreach([2,3,5,7,13,17,19,31,61,89,107,127], q,
    my(M = 2^q - 1, f = Mod(Mod(1, M) * x, Mod(1, M) * (x^2 - x - 1)), e);
    e = f^(2^q);
    print("    q=", q, " (q mod 4 = ", q % 4, "), M_q mod 5 = ", M % 5, ", kronecker(5,M_q) = ", kronecker(5, M),
          ": phi^(2^q) == -1 ? ", e == -1, "; phi^(2^(q+1)) == 1 ? ", e^2 == 1,
          if (q <= 31, Str("; pi(M_q) = ", pis(M)), "")));
  \\ (4) Lucas-Lehmer seed 3
  my(okll = 1);
  forprime(q = 3, 300, if (q % 4 != 3, next);
    my(M = 2^q - 1, s = Mod(3, M));
    for (i = 1, q - 2, s = s^2 - 2);
    if ((s == 0) != isprime(M), okll = 0; print("LL seed 3 FAILS at q=", q)));
  print("(4) Lucas-Lehmer with seed 3 agrees with isprime(M_q) for all primes q = 3 mod 4, q < 300: ", okll);
  \\ seed 3 for q = 1 mod 4 (expected to fail as a test)
  my(f1 = List());
  forprime(q = 5, 130, if (q % 4 != 1, next);
    my(M = 2^q - 1, s = Mod(3, M));
    for (i = 1, q - 2, s = s^2 - 2);
    listput(f1, [q, isprime(M), s == 0]));
  print("    seed 3 with q = 1 mod 4 [q, M_q prime, M_q | s_(q-2)]: ", Vec(f1));
}
