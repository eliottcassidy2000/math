\\ audit_wave2: the note's formula  genus(X_0^+(M)) = sum_{d | M} #{newforms of level d with all AL signs +1}
\\ (M squarefree) compared, as exact integers, with the independent Riemann-Hurwitz/class-number genus
\\ of rh_genus_lists.gp, at (a) every squarefree M <= 300 and (b) every level Np used for the note's
\\ thirteen lists (N in the table, p an Ogg prime not dividing N; M up to 22*71 = 1562).
\\ The newform side re-implements the note's counting with PARI's mfatkineigenvalues (newspace only,
\\ no mfatkininit); it is the method under audit, the R-H side is the independent comparison.
\\ Run: gp -q -f formula_vs_rh.gp > formula_vs_rh.out
default(parisize, 400000000);
g0(M) = { my(f = factor(M)[, 1], psi = prod(i = 1, #f, f[i] + 1),
     n2 = prod(i = 1, #f, 1 + kronecker(-4, f[i])), n3 = prod(i = 1, #f, 1 + kronecker(-3, f[i])));
  1 + psi / 12 - n2 / 4 - n3 / 3 - 2^(#f) / 2; }
nufix(M, Q) = { my(Mp = M / Q, fl = if (Mp == 1, [], factor(Mp)[, 1]), s);
  if (Q == 2, return(prod(i = 1, #fl, 1 + kronecker(-4, fl[i])) + prod(i = 1, #fl, 1 + kronecker(-8, fl[i]))));
  if (Q == 3, return(prod(i = 1, #fl, 1 + kronecker(-12, fl[i])) + prod(i = 1, #fl, if (fl[i] == 2, 1, 1 + kronecker(-3, fl[i])))));
  s = qfbclassno(-4 * Q) * prod(i = 1, #fl, 1 + kronecker(-4 * Q, fl[i]));
  if (Q % 4 == 3, s += qfbclassno(-Q) * prod(i = 1, #fl, if (fl[i] == 2, 3, 1 + kronecker(-Q, fl[i]))));
  s; }
gplusRH(M) = { if (M == 1, return(0)); my(Qs = select(Q -> Q > 1, divisors(M)));
  (2 * g0(M) - 2 - sum(i = 1, #Qs, nufix(M, Qs[i]))) / (2 * (#Qs + 1)) + 1; }

NP = Map();
newplus(d) = {
  my(v);
  if (d == 1, return(0));
  if (mapisdefined(NP, d, &v), return(v));
  my(mf = mfinit([d, 2], 0), cnt = 0);
  if (mfdim(mf) > 0,
    my(fa = factor(d)[, 1], ev = vector(#fa, i, mfatkineigenvalues(mf, fa[i])));
    for (o = 1, #ev[1], for (e = 1, #ev[1][o], if (prod(i = 1, #fa, ev[i][o][e] == 1), cnt++))));
  mapput(NP, d, cnt); cnt;
}
gplusNF(M) = sumdiv(M, d, newplus(d));

{
  my(bad = 0, n = 0, levs = List());
  forstep (M = 2, 300, 1, if (issquarefree(M), listput(levs, M)));
  my(Ns = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22], OGG = [2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]);
  foreach (Ns, N, foreach (OGG, p, if (N % p, listput(levs, N * p))));
  levs = Set(Vec(levs));
  my(maxg = 0, hist = Map());
  foreach (levs, M, my(a = gplusNF(M), b = gplusRH(M)); n++; maxg = max(maxg, a);
    if (a != b, bad++; print("  MISMATCH at M = ", M, ": newform count ", a, ", Riemann-Hurwitz ", b)));
  print("== newform-sign formula vs Riemann-Hurwitz, exact genus values ==");
  print("  ", n, " squarefree levels (all M <= 300 and every Np of the note's lists, max ", vecmax(levs), "); largest genus ", maxg, "; mismatches: ", bad);
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
