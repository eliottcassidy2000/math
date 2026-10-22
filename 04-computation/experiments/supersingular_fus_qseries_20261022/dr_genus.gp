\\ Where supersingular counts enter the levels N = p^2, pq, p^2 q (opus-2026-10-22-ssj).
\\ Classical inputs (CITED): Deligne-Rapoport (X_0(lM) mod l, l not dividing M, is two copies
\\ of X_0(M) crossing transversally at the H(l, M) supersingular points), Katz-Mazur (X_0(l^2 M)
\\ mod l: three copies of X_0(M) indexed by the divisors 1, l, l^2, the middle one of multiplicity
\\ l - 1, all meeting at the supersingular points), Eichler's class number H(l, M).
\\ Tests (genera from PARI's mfdim, an independent implementation of the dimension formula):
\\   (G1) g(X_0(pq))  = 2 g(X_0(q)) + H(p, q) - 1 = 2 g(X_0(p)) + H(q, p) - 1
\\   (G2) g(X_0(p^2 q)) = 2 g(X_0(p^2)) + H(q, p^2) - 1
\\   (G3) H(l, 1) = ss(l) = deg ss_l (the census), so (G1) at q = 1 reads g(X_0(p)) = ss(p) - 1
\\ plus the cusp/divisor dictionary: #cusps(X_0(N)) = sum_{d | N} phi(gcd(d, N/d)), which equals
\\ tau(N) = F(N) + 2 exactly when N is squarefree.
\\ Run: gp -q -f dr_genus.gp > dr_genus.out
default(parisize, 64000000);
read("sspoly.gp");
PQ = 160;

H(l, M) = {
  my(fa = factor(M)[, 1], psi = M * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  nu2 = if (M % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (M % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (l - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, l)) + nu3 / 3 * (1 - kronecker(-3, l));
}
g0(N) = mfdim([N, 2], 1);
cusps(N) = sumdiv(N, d, eulerphi(gcd(d, N / d)));
FUS(N) = {
  my(D = select(d -> d > 1 && d < N, divisors(N)));
  [#D, #select(d -> isprime(d), D), #select(d -> issquarefree(d), D)];
}

{
  my(bad = 0, n1 = 0, n2 = 0, n3 = 0);
  forprime (p = 2, PQ,
    \\ (G3) for p >= 5 against the census polynomial
    if (p >= 5, n3++; if (H(p, 1) != poldegree(ssA(p)) || g0(p) != H(p, 1) - 1, bad++; print("G3 fails at ", p)));
    forprime (q = 2, PQ,
      if (q == p, next);
      n1++;
      if (g0(p * q) != 2 * g0(q) + H(p, q) - 1, bad++; print("G1 fails at ", [p, q]));
      if (p * p * q <= 40000, n2++;
        if (g0(p^2 * q) != 2 * g0(p^2) + H(q, p^2) - 1, bad++; print("G2 fails at ", [p, q])))));
  print("(G1) on ", n1, " ordered prime pairs < ", PQ, ";  (G2) on ", n2, " pairs with p^2 q <= 40000;  (G3) on ", n3, " primes");
  print("failures: ", bad);
  print("\nTable: the three asked shapes with their divisor counts (F,U,S), X_0(N) cusps, genus, and the");
  print("supersingular gluing numbers at each prime of bad reduction.");
  print("  shape  N      (F,U,S)    tau  #cusps  g(X_0(N))  gluing data");
  for (k = 1, 3,
    my(pp = [3, 5, 7][k], qq = [2, 11, 13][k]);
    my(N = pp^2); print("  p^2    ", N, "   ", FUS(N), "  ", numdiv(N), "   ", cusps(N), "      ", g0(N),
      "        at p: 3 components (divisors 1,p,p^2; middle mult. p-1) meeting at ss(p) = ", H(pp, 1), " points");
    N = pp * qq; print("  pq     ", N, "   ", FUS(N), "  ", numdiv(N), "   ", cusps(N), "      ", g0(N),
      "        at p: 2 x X_0(q) glued at H(p,q) = ", H(pp, qq), ";  at q: 2 x X_0(p) glued at H(q,p) = ", H(qq, pp));
    N = pp^2 * qq; print("  p^2q   ", N, "  ", FUS(N), "  ", numdiv(N), "   ", cusps(N), "      ", g0(N),
      "        at q: 2 x X_0(p^2) glued at H(q,p^2) = ", H(qq, pp^2), ";  at p: 3 x X_0(q) meeting at H(p,q) = ", H(pp, qq)));
  print("\nclosed forms (p, q distinct primes; (a/b) Kronecker symbols):");
  print("  H(p,q)   = (p-1)(q+1)/12 + (1/4)(1-(-4/p))(1+(-4/q)) + (1/3)(1-(-3/p))(1+(-3/q))");
  print("  H(q,p^2) = (q-1)p(p+1)/12 + (1/4)(1-(-4/q)) nu2(p^2) + (1/3)(1-(-3/q)) nu3(p^2),");
  print("             nu2(p^2) = 1+(-4/p) (0 at p = 2), nu3(p^2) = 1+(-3/p) (0 at p = 3)");
  print("  symmetric form of (G1): 2 ss(q) + H(p,q) = 2 ss(p) + H(q,p)  (ss(l) := H(l,1))");
  my(sym = 0);
  forprime (p = 2, PQ, forprime (q = 2, PQ, if (p != q, if (2 * H(q, 1) + H(p, q) != 2 * H(p, 1) + H(q, p), sym++))));
  print("  symmetric form fails on ", sym, " pairs");
  bad += sym;
  print("  cusp dictionary: #cusps(X_0(N)) = tau(N) = F+2 for squarefree N; at p^2: p+1 = F+2 + (p-2); at p^2 q: 2(p+1)");
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
