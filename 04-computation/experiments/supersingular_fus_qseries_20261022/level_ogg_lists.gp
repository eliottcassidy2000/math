\\ Level-N Ogg lists and the sporadic centralizer chain M > 2.B > 6.Fi22 (opus-2026-10-22-ssj, wave 2).
\\
\\ For squarefree N with X_0(N)^+ (quotient by all Atkin-Lehner involutions) of genus 0, a prime p not
\\ dividing N has the "rationality property" (every supersingular point of X_0(N)^+ mod p is F_p-rational)
\\ iff genus(X_0(Np)^+) = genus(X_0(N)^+) = 0 (Aricheta, SIGMA 15 (2019) 007, Theorem 1.1; for N = 1 this
\\ is Ogg's theorem). For squarefree level M the W-invariant part of S_2(Gamma_0(M)) is spanned, for each
\\ newform f of level d | M with all Atkin-Lehner signs +1, by exactly one vector of f's oldspace, so
\\     genus(X_0(M)^+) = sum_{d | M} #{newforms of level d (with Galois multiplicity) with eps_l = +1 for all l | d}.
\\ This is cross-checked against PARI's Atkin-Lehner traces on the full space for small M.
\\ Run: gp -q -f level_ogg_lists.gp > level_ogg_lists.out
default(parisize, 256000000);
OGG = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71];

\\ number of newforms of level d (counted with multiplicity) whose Atkin-Lehner signs are all +1
newplus(d) = {
  if (d == 1, return(0));
  my(mf = mfinit([d, 2], 0), fa, ev, cnt = 0, no);
  if (mfdim(mf) == 0, return(0));
  fa = factor(d)[, 1];
  ev = vector(#fa, i, mfatkineigenvalues(mf, fa[i]));   \\ ev[i][orbit][embedding]
  no = #ev[1];
  for (o = 1, no,
    for (e = 1, #ev[1][o],
      if (prod(i = 1, #fa, ev[i][o][e] == 1), cnt++)));
  cnt;
}
gplus(M) = { if (!issquarefree(M), error("squarefree levels only: ", M)); sumdiv(M, d, newplus(d)) };

\\ independent route: (1/|W|) sum over Hall divisors Q of trace(W_Q) on the full cusp space
gplus_trace(M) = {
  my(mf = mfinit([M, 2], 1), d = mfdim(mf), s = 0, Qs, W);
  if (d == 0, return(0));
  Qs = select(Q -> gcd(Q, M / Q) == 1, divisors(M));
  for (i = 1, #Qs, if (Qs[i] == 1, s += d, W = mfatkininit(mf, Qs[i]); s += trace(W[2])));
  s / #Qs;
}

{
  my(bad = 0, ctrl = 0);
  print("== cross-check of the newform count against Atkin-Lehner traces (squarefree M <= 150) ==");
  forstep (M = 2, 150, 1, if (!issquarefree(M), next);
    ctrl++;
    if (gplus(M) != gplus_trace(M), bad++; print("  MISMATCH at M = ", M, ": ", gplus(M), " vs ", gplus_trace(M))));
  print("  ", ctrl, " squarefree levels, mismatches: ", bad);

  print("\n== level-N Ogg lists: primes p in Ogg's list, p not dividing N, with genus(X_0(Np)^+) = 0 ==");
  my(levels = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22], lists = Map());
  for (k = 1, #levels,
    my(N = levels[k], L = List(), g0 = gplus(N));
    if (g0 != 0, print("  note: genus(X_0(", N, ")^+) = ", g0));
    for (i = 1, #OGG, my(p = OGG[i]); if (N % p == 0, next);
      if (gplus(N * p) == g0, listput(L, p)));
    mapput(lists, N, Vec(L));
    print("  N = ", N, ": ", Vec(L)));

  \\ Aricheta's Appendix A (CITED) for comparison
  my(table = [[2, [3,5,7,11,13,17,19,23,31,47]], [3, [2,5,7,11,13,17,23,29]], [5, [2,3,7,11,19]], [6, [5,7,11,13]],
              [7, [2,3,5,17]], [10, [3,7,11]], [11, [2,3,5]], [13, [2,3]], [14, [3,5]], [15, [2,7]], [21, [2,5]], [22, [3,5]]]);
  print("\n== comparison with Aricheta (2019), Appendix A ==");
  for (i = 1, #table, my(N = table[i][1], got = mapget(lists, N));
    if (got != table[i][2], bad++; print("  DIFFERS at N = ", N, ": ", got, " vs table ", table[i][2]), print("  N = ", N, ": agrees")));
  print("  N = 1: ", mapget(lists, 1), (if (mapget(lists, 1) == OGG, "  (Ogg's fifteen)", "  (MISMATCH with Ogg)")));
  if (mapget(lists, 1) != OGG, bad++);

  \\ sporadic groups (orders CITED from the ATLAS) and the centralizer chain
  my(orders = [["M", 808017424794512875886459904961710757005754368000000000],
               ["B", 4154781481226426191177580544000000],
               ["Fi24'", 1255205709190661721292800],
               ["HN", 273030912000000],
               ["Fi22", 64561751654400],
               ["He", 4030387200],
               ["HS", 44352000],
               ["M12", 95040],
               ["L3(3)", 5616]]);
  my(match = [[1, "M"], [2, "B"], [3, "Fi24'"], [5, "HN"], [6, "Fi22"], [7, "He"], [10, "HS"], [11, "M12"], [13, "L3(3)"]]);
  print("\n== level-N list together with the primes of N, against the prime set of the sporadic group ==");
  for (i = 1, #match,
    my(N = match[i][1], G = match[i][2], ord, ps, lst);
    for (k = 1, #orders, if (orders[k][1] == G, ord = orders[k][2]));
    ps = factor(ord)[, 1]~;
    lst = vecsort(concat(mapget(lists, N), factor(N)[, 1]~), , 8);
    print("  N = ", N, "  ", G, ": primes of |", G, "| = ", ps, "   level-N list + primes of N = ", lst, if (ps == lst, "   EQUAL", "   DIFFERENT")));
  print("\n== the owner's three sets as the filtration along levels 1 | 2 | 6 ==");
  my(L1 = mapget(lists, 1), L2 = vecsort(concat(mapget(lists, 2), [2])), L6 = vecsort(concat(mapget(lists, 6), [2, 3])));
  print("  level 1 (M):    ", L1);
  print("  level 2 (B):    ", L2);
  print("  level 6 (Fi22): ", L6);
  print("  M \\ B    = ", setminus(Set(L1), Set(L2)), "   (owner: {29,41,59,71})");
  print("  B \\ Fi22 = ", setminus(Set(L2), Set(L6)), "   (owner: {17,19,23,31,47})");
  print("  Fi22     = ", Set(L6), "   (owner: {2,3,5,7,11,13})");
  if (setminus(Set(L1), Set(L2)) != [29,41,59,71] || setminus(Set(L2), Set(L6)) != [17,19,23,31,47] || Set(L6) != [2,3,5,7,11,13], bad++);
  \\ the four levels where the filtration loses primes: genus of X_0(2p)^+ and X_0(6p)^+ for the dropped primes
  print("\n== genera where primes drop out ==");
  foreach ([29, 41, 59, 71], p, print("  genus(X_0(", 2 * p, ")^+) = ", gplus(2 * p)));
  foreach ([17, 19, 23, 31, 47], p, print("  genus(X_0(", 6 * p, ")^+) = ", gplus(6 * p), ",  genus(X_0(", 2 * p, ")^+) = ", gplus(2 * p)));
  \\ the owner's number and the three sets by valuations (orders CITED from the ATLAS)
  print("\n== |M|/|B|, the owner's number, and the valuation classes ==");
  my(oM = orders[1][2], oB = orders[2][2], r = oM / oB, num = 2^5 * 3^7 * 5^3 * 7^4 * 11 * 13^2);
  print("  |M| = ", factor(oM)~);
  print("  |B| = ", factor(oB)~);
  print("  |M|/|B| = ", factor(r)~);
  print("  owner's number 2^5*3^7*5^3*7^4*11*13^2 = ", num, ";  |M|/|B| / (29*41*59*71) = ", r / (29 * 41 * 59 * 71), if (r / (29 * 41 * 59 * 71) == num, "  EQUAL", "  DIFFERENT"));
  if (r / (29 * 41 * 59 * 71) != num, bad++);
  my(cl = vector(3, i, List()));
  for (i = 1, #OGG, my(q = OGG[i], a = valuation(oM, q), b = valuation(oB, q));
    if (a >= 2, listput(cl[1], q), if (b == 1, listput(cl[2], q), listput(cl[3], q))));
  print("  v_p|M| >= 2: ", Vec(cl[1]), ";  v_p|M| = v_p|B| = 1: ", Vec(cl[2]), ";  v_p|M| = 1, v_p|B| = 0: ", Vec(cl[3]));
  if (Vec(cl[1]) != [2,3,5,7,11,13] || Vec(cl[2]) != [17,19,23,31,47] || Vec(cl[3]) != [29,41,59,71], bad++);
  \\ hostile probe for the exponent "cycle": small-prime exponents of other indices
  my(extra = [["Th", 90745943887872000], ["Co1", 4157776806543360000]], sm = [2, 3, 5, 7, 11, 13]);
  print("  exponents at 2, 3, 5, 7, 11, 13 of the indices:");
  foreach (concat([[orders[2][1], oB], [orders[3][1], orders[3][2]], [orders[5][1], orders[5][2]], [orders[4][1], orders[4][2]], [orders[6][1], orders[6][2]]], extra), v,
    print("    [M : ", v[1], "] = ", vector(6, i, valuation(oM, sm[i]) - valuation(v[2], sm[i]))));
  print("    [B : Fi22] = ", vector(6, i, valuation(oB, sm[i]) - valuation(orders[5][2], sm[i])));
  \\ sizes of the successive differences along several level chains
  print("  difference sizes along level chains (list at each level, with the primes of the level):");
  my(full = (N) -> vecsort(concat(mapget(lists, N), factor(N)[, 1]~), , 8), chains = [[1, 2, 6], [1, 2, 10], [1, 3, 6], [1, 7, 21]]);
  for (i = 1, #chains, my(c = chains[i], A = if (c[1] == 1, OGG, full(c[1])), Bv = full(c[2]), Cv = full(c[3]));
    print("    ", c, ": ", [#A - #Bv, #Bv - #Cv, #Cv], "  nested: ", #setminus(Set(Bv), Set(A)) == 0 && #setminus(Set(Cv), Set(Bv)) == 0));
  \\ monotonicity (audit F6): genus(X_0^+(M)) = sum over d | M of newplus(d) is monotone under divisibility,
  \\ so for squarefree N | N' and p not dividing N', L(N') is inside L(N) minus the primes of N'
  print("\n== monotonicity lemma: L(N') inside L(N) for N | N' (squarefree) ==");
  my(mono = 0, pairs = 0);
  for (a = 1, #levels, for (b = 1, #levels,
    my(N = levels[a], Np = levels[b]);
    if (N == Np || Np % N, next);
    pairs++;
    my(LN = Set(mapget(lists, N)), LNp = Set(mapget(lists, Np)), ex = select(q -> Np % q != 0, LN));
    if (#setminus(LNp, Set(ex)) > 0, mono++; print("  violation: L(", Np, ") not inside L(", N, ")"))));
  print("  ", pairs, " divisibility pairs among the 13 levels; violations: ", mono);
  bad += mono;
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
