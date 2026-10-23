\\ audit_wave2 (independent auditor of monster_levels_square_sum_7_21_20261022.md, section 2)
\\ Genus of X_0(M)/W for squarefree M by Riemann-Hurwitz, with the fixed points of every
\\ Atkin-Lehner involution w_Q counted by CM points (class numbers), derived from first principles:
\\   w_Q fixes (E, C) iff E has an endomorphism a with ker a = C_Q, a^2 = Q*unit and a(C_{M/Q}) = C_{M/Q}.
\\   Q >= 5:  nu = h(-4Q) prod_{l | M/Q} (1 + (-4Q/l))
\\               + [Q = 3 mod 4] h(-Q) prod_{l | M/Q} b_l,  b_l = 1 + (-Q/l) (l odd),  b_2 = 3
\\            (for O = Z[(1+sqrt(-Q))/2], sqrt(-Q) = -1 mod 2O acts trivially on E[2], so all three
\\             lines are stable; for Z[sqrt(-Q)] one line is stable at 2).
\\   Q = 2:   nu = prod_{l | M/2} (1 + (-4/l)) + prod_{l | M/2} (1 + (-8/l))      (j = 1728 via 1+i; j = 8000 via sqrt(-2))
\\   Q = 3:   nu = prod_{l | M/3} (1 + (-12/l)) + prod_{l | M/3, l odd} (1 + (-3/l))
\\            (j = 54000 via sqrt(-3); j = 0 via 1-omega: at l = 2 the three lines are stable and are
\\             permuted transitively by omega, so they give one point).
\\ No modular symbols, no mf* functions are used for the lists.  The fixed-point formulas are then
\\ VALIDATED against the Lefschetz identity nu(w_Q) = 2 - 2 tr(W_Q | S_2(Gamma_0(M))) using PARI's
\\ Atkin-Lehner traces at small levels (M <= 150), a check of this file's formulas, not an input.
\\ Run: gp -q -f rh_genus_lists.gp > rh_genus_lists.out
default(parisize, 200000000);

g0(M) = {
  my(f = factor(M)[, 1], psi = prod(i = 1, #f, f[i] + 1),
     n2 = prod(i = 1, #f, 1 + kronecker(-4, f[i])), n3 = prod(i = 1, #f, 1 + kronecker(-3, f[i])));
  1 + psi / 12 - n2 / 4 - n3 / 3 - 2^(#f) / 2;
}

nufix(M, Q) = {
  my(Mp = M / Q, fl = if (Mp == 1, [], factor(Mp)[, 1]), s);
  if (Q == 2,
    return(prod(i = 1, #fl, 1 + kronecker(-4, fl[i])) + prod(i = 1, #fl, 1 + kronecker(-8, fl[i]))));
  if (Q == 3,
    return(prod(i = 1, #fl, 1 + kronecker(-12, fl[i]))
         + prod(i = 1, #fl, if (fl[i] == 2, 1, 1 + kronecker(-3, fl[i])))));
  s = qfbclassno(-4 * Q) * prod(i = 1, #fl, 1 + kronecker(-4 * Q, fl[i]));
  if (Q % 4 == 3,
    s += qfbclassno(-Q) * prod(i = 1, #fl, if (fl[i] == 2, 3, 1 + kronecker(-Q, fl[i]))));
  s;
}

\\ genus of X_0(M) / <w_Q : Q in Qs> where Qs is a subgroup (list of Hall divisors > 1, closed under *)
gquot(M, Qs) = {
  my(n = #Qs + 1, s = sum(i = 1, #Qs, nufix(M, Qs[i])), g = (2 * g0(M) - 2 - s) / (2 * n) + 1);
  if (type(g) != "t_INT", error("non-integral genus at M = ", M, " Qs = ", Qs, ": ", g));
  g;
}
hall(M) = select(Q -> Q > 1, divisors(M));                    \\ squarefree M: every divisor is Hall
gplusRH(M) = gquot(M, hall(M));

{
  my(bad = 0);
  print("== 1. validation of the fixed-point formulas: nu(w_Q) = 2 - 2 tr(W_Q | S_2), squarefree M <= 150 ==");
  my(cnt = 0, mism = 0);
  forstep (M = 2, 150, 1, if (!issquarefree(M), next);
    my(mf = mfinit([M, 2], 1), d = mfdim(mf), Qs = hall(M));
    if (d != g0(M), mism++; print("  genus mismatch at ", M));
    for (i = 1, #Qs, my(Q = Qs[i], tr = if (d == 0, 0, my(W = mfatkininit(mf, Q)); trace(W[2]) / W[3]));
      cnt++;
      if (nufix(M, Q) != 2 - 2 * tr, mism++; print("  nu mismatch M = ", M, " Q = ", Q, ": formula ", nufix(M, Q), " vs Lefschetz ", 2 - 2 * tr))));
  print("  ", cnt, " pairs (M, Q) checked, mismatches: ", mism);
  bad += mism;

  print("\n== 2. squarefree M <= 1000 with genus(X_0(M)^+) = 0 (Riemann-Hurwitz) ==");
  my(z = []);
  forstep (M = 1, 1000, 1, if (!issquarefree(M), next); if (M == 1 || gplusRH(M) == 0, z = concat(z, M)));
  print("  ", #z, " levels: ", z);
  \\ Conway-Norton's squarefree N+ list (auditor's recollection, 44 levels incl. 1), compared, not used as input
  my(cn = [1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,38,39,41,42,46,47,51,55,59,62,66,69,70,71,78,87,94,95,105,110,119]);
  print("  equals the recalled Conway-Norton squarefree N+ list: ", z == cn);
  if (z != cn, bad++);

  print("\n== 3. level-N rationality lists via Theorem 1.1: p prime, p not | N, p < 2000, genus(X_0(Np)^+) = genus(X_0(N)^+) ==");
  my(levels = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22],
     note = [[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71], [3,5,7,11,13,17,19,23,31,47], [2,5,7,11,13,17,23,29], [2,3,7,11,19],
             [5,7,11,13], [2,3,5,17], [3,7,11], [2,3,5], [2,3], [3,5], [2,7], [2,5], [3,5]]);
  for (k = 1, #levels, my(N = levels[k], gN = if (N == 1, 0, gplusRH(N)), L = []);
    forprime (p = 2, 2000, if (N % p == 0, next); if (gplusRH(N * p) == gN, L = concat(L, p)));
    print("  N = ", N, " (genus ", gN, "): ", L, if (L == note[k], "   = note's table", "   DIFFERS from note's table ", note[k]));
    if (L != note[k], bad++));

  print("\n== 4. row 2- (X = X_0(2), X^p = X_0(2p)/w_p) ==");
  my(L2m = []);
  forprime (p = 3, 2000, if (gquot(2 * p, [p]) == 0, L2m = concat(L2m, p)));
  print("  primes p < 2000 with genus(X_0(2p)/w_p) = 0: ", L2m, "  (Aricheta Table 1 / App. A: [3, 5, 7, 11, 23])");
  if (L2m != [3, 5, 7, 11, 23], bad++);
  print("  Aricheta Table 3 row 'gen X_0(2)^p' for p = 3..47: ", vector(14, i, gquot(2 * [3,5,7,11,13,17,19,23,29,31,37,41,43,47][i], [[3,5,7,11,13,17,19,23,29,31,37,41,43,47][i]])));
  print("  Aricheta Table 3 row 'gen X_0^+(2)^p' for p = 3..47: ", vector(14, i, gplusRH(2 * [3,5,7,11,13,17,19,23,29,31,37,41,43,47][i])));

  print("\n== 5. the note's drop-out genera ==");
  print("  genus X_0^+(2p), p = 29, 41, 59, 71: ", apply(p -> gplusRH(2 * p), [29, 41, 59, 71]), "   (note: 1, 1, 1, 1)");
  print("  genus X_0^+(6p), p = 17, 19, 23, 31, 47: ", apply(p -> gplusRH(6 * p), [17, 19, 23, 31, 47]), "   (note: 1, 1, 1, 2, 3)");
  print("  genus X_0^+(2p), p = 17, 19, 23, 31, 47: ", apply(p -> gplusRH(2 * p), [17, 19, 23, 31, 47]), "   (note: 0 for all five)");
  if (apply(p -> gplusRH(2 * p), [29, 41, 59, 71]) != [1,1,1,1] || apply(p -> gplusRH(6 * p), [17, 19, 23, 31, 47]) != [1,1,1,2,3], bad++);

  print("\n== 6. nesting is automatic: for squarefree N | N' (p not | N') genus(X_0^+(N'p)) >= genus(X_0^+(Np)) ==");
  my(viol = 0, tested = 0);
  forstep (Np = 2, 3000, 1, if (!issquarefree(Np), next);
    my(gp1 = gplusRH(Np));
    forprime (q = 2, 47, if (Np % q == 0, next); if (Np * q > 6000, break); tested++;
      if (gplusRH(Np * q) < gp1, viol++; print("  violation: ", Np, " -> ", Np * q))));
  print("  ", tested, " pairs (M, Mq) with M squarefree <= 3000, q <= 47, Mq <= 6000: violations ", viol);
  if (viol, bad++);

  print("\n== 7. growth of genus(X_0^+(Np)) (row A of section 5: 'grows like Np/24') ==");
  foreach ([1, 2, 6, 7, 21], N, my(p = nextprime(10^5 + N * 1000));
    if (N % p == 0, p = nextprime(p + 1));
    print("  N = ", N, ", p = ", p, ": genus = ", gplusRH(N * p), ";  Np/24 = ", round(N * p / 24),
          ";  psi(Np)/(12*2^omega(Np)) = ", round(prod(i = 1, #factor(N * p)[, 1], factor(N * p)[i, 1] + 1) / (12 * 2^omega(N * p)))));

  print("\nfailures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
