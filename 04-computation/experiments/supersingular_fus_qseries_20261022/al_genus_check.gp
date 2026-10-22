\\ Independent check of the F_p-rational supersingular counts through Atkin-Lehner traces
\\ (opus-2026-10-22-ssj).
\\ Deligne-Rapoport: X_0(Np) mod p (p not dividing N) is two copies of X_0(N) glued at the H(p, N)
\\ supersingular points; w_p swaps the copies and acts on the nodes as Frobenius. Hence
\\     g(X_0(Np)/w_p) = g(X_0(N)) + (H(p, N) - #X_0(N)^ss(F_p)) / 2,
\\ the level-N form of g(X_0(p)^+) = (#ss - #ss(F_p))/2.  The left side is computed here from PARI's
\\ Atkin-Lehner operator on S_2(Gamma_0(Np)) (g^+ = (dim + trace W_p)/2), an implementation
\\ independent of the Hauptmodul root counts; the right side uses the class-number formula
\\     #X_0(N)^ss(F_p) = (h(-4p) s_N(floor) + [p = 3 mod 4] h(-p) s_N(surface)) / 2
\\ with the local numbers s_N computed on (Z/N)^2 exactly as in level_fp_counts.gp.
\\ Run: gp -q -f al_genus_check.gp > al_genus_check.out
default(parisize, 512000000);
PMAX = 150;
LEVMAX = 300;   \\ cap on the level N p: mfatkininit cost and memory grow fast (a 512 MB stack overflowed below level 700)

H(l, M) = {
  if (M == 1, return((l - 1) / 12 + (1 - kronecker(-4, l)) / 4 + (1 - kronecker(-3, l)) / 3));
  my(fa = factor(M)[, 1], psi = M * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  nu2 = if (M % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (M % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (l - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, l)) + nu3 / 3 * (1 - kronecker(-3, l));
}
slocal(M, N) = {
  if (N == 1, return(1));
  my(cnt = 0, P, Q);
  for (a = 0, N - 1, for (b = 0, N - 1,
    if (N / gcd([a, b, N]) != N, next);
    P = [Mod(a, N), Mod(b, N)]~; Q = M * P;
    for (l = 0, N - 1, if (Q == l * P, cnt++; break))));
  cnt / eulerphi(N);
}
fpcount(p, N) = {
  my(h4 = qfbclassno(-4 * p), h1 = if (p % 4 == 3, qfbclassno(-p), 0), sf, ss = 0);
  sf = slocal(Mod([0, -p; 1, 0], max(N, 1)), N);
  if (p % 4 == 3, my(c = (1 + p) / 4); ss = slocal(2 * Mod([0, -c; 1, 1], max(N, 1)) - 1, N));
  (h4 * sf + h1 * ss) / 2;
}

{
  my(bad = 0, cnt = 0, rows = List());
  foreach ([1, 2, 4, 7, 8, 9], N,
    my(gN = if (N == 1, 0, mfdim([N, 2], 1)), okN = 0);
    forprime (p = 5, PMAX - 1,
      if (N % p == 0 || N * p > LEVMAX, next);
      my(mf = mfinit([N * p, 2], 1), d = mfdim(mf), tr, gplus, rhs, F);
      tr = if (d == 0, 0, trace(mfatkininit(mf, p)[2]));
      gplus = (d + tr) / 2;
      F = fpcount(p, N);
      rhs = gN + (H(p, N) - F) / 2;
      cnt++;
      if (gplus != rhs, bad++; print("  MISMATCH N=", N, " p=", p, ": g(X_0(Np)/w_p) = ", gplus, " vs ", rhs), okN++);
      if (p < 30, listput(rows, [N, p, d, gplus, H(p, N), F]))));
  print("checked ", cnt, " pairs (N in {1,2,4,7,8,9}, 5 <= p < ", PMAX, ", N p <= ", LEVMAX, "); mismatches: ", bad);
  print("sample rows [N, p, g(X_0(Np)), g(X_0(Np)/w_p), H(p,N), #X_0(N)^ss(F_p)]:");
  for (i = 1, #rows, print("   ", rows[i]));
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
