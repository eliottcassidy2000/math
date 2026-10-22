\\ F_p-rational supersingular points on X_0(2), X_0(4), X_0(8), X_1(8), X_1(9) against class
\\ numbers (opus-2026-10-22-ssj). Needs sspoly.gp and level_relations.gp.
\\ For each prime 5 <= p < PMAX: n_p = #{supersingular j in F_p}; for each such j the number of
\\ F_p-roots of the fibre polynomial A(x) - j B(x) on the coordinate line of the curve.
\\ Fit (exact, per residue class of p mod 8): count = gamma(p mod 8) n_p.
\\ X_1(9): expected 6 [p = 8 mod 9] n_p (eigenvalue argument, see the results note).
\\ Run: gp -q -f level_fp_counts.gp > level_fp_counts.out
default(parisize, 128000000);
read("sspoly.gp");
read("level_relations.gp");
PMAX = 2000;

\\ classical Hauptmoduln of X_0(2), X_0(4), checked on q-expansions first
{
  my(P = 300, E4 = 1 + 240 * sum(n = 1, P, sigma(n, 3) * 'q^n) + O('q^(P + 1)), Dl = 'q * eta('q + O('q^(P + 1)))^24, J = E4^3 / Dl,
     e1 = eta('q + O('q^P)), e2 = eta('q^2 + O('q^P)), e4 = eta('q^4 + O('q^P)), t2, t4);
  t2 = 'q^-1 * e1^24 / e2^24;
  t4 = 'q^-1 * e1^8 / e4^8;
  print("X_0(2): j = (t+256)^3/t^2 at t = (eta/eta_2)^24: ", valuation(J * t2^2 - (t2 + 256)^3 + O('q^(P - 10)), 'q) >= P - 10);
  print("X_0(4): j = (t^2+256t+4096)^3/(t^4 (t+16)) at t = (eta/eta_4)^8: ",
    valuation(J * t4^4 * (t4 + 16) - (t4^2 + 256 * t4 + 4096)^3 + O('q^(P - 10)), 'q) >= P - 10);
}

{
  curves = [["X_0(2)", ('x + 256)^3, 'x^2, 3], ["X_0(4)", ('x^2 + 256 * 'x + 4096)^3, 'x^4 * ('x + 16), 6],
            ["X_0(8)", P8, Q8, 12], ["X_1(8)", A8, B8, 24], ["X_1(9)", A9, B9, 36]];
}



eichler(p, N) = {
  my(fa = factor(N)[, 1], psi = N * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  nu2 = if (N % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (N % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (p - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, p)) + nu3 / 3 * (1 - kronecker(-3, p));
}
\\ local counts on (Z/N)^2 with Frobenius matrix M (entries mod N): typ 0 = pi-stable cyclic subgroups of
\\ order N (Gamma_0), typ 1 = +-classes of points of exact order N with M P = +-P (Gamma_1)
slocal(M, N, typ) = {
  my(cnt = 0, P, Q, ord);
  for (a = 0, N - 1, for (b = 0, N - 1,
    ord = N / gcd([a, b, N]);
    if (ord != N, next);
    P = [Mod(a, N), Mod(b, N)]~;
    Q = M * P;
    if (typ == 1, if (Q == P || Q == -P, cnt++),
      for (l = 0, N - 1, if (Q == l * P, cnt++; break)))));
  if (typ == 1, cnt / 2, cnt / eulerphi(N));
}
{
  my(data = vector(#curves, c, Map()), bad = 0, nprimes = 0);
  forprime (p = 5, PMAX - 1,
    my(s = ssA(p), roots = polrootsmod(s, p), np = #roots, h4 = qfbclassno(-4 * p), h1 = if (p % 4 == 3, qfbclassno(-p), 0));
    if (np != (h4 + h1) / 2, bad++);
    nprimes++;
    for (c = 1, #curves,
      if (c == 5 && p == 3, next);
      my(A = curves[c][2], B = curves[c][3], tot = 0, key = if (c == 5, p % 9, p % 8), lst);
      for (i = 1, #roots, my(j = lift(roots[i])); tot += #polrootsmod(A - j * B, p));
      if (mapisdefined(data[c], key), lst = mapget(data[c], key), lst = List());
      listput(lst, [p, np, h4, h1, tot]);
      mapput(data[c], key, lst)));
  print("\nprimes 5 <= p < ", PMAX, ": ", nprimes, ";  n_p = (h(-4p) + [p=3 mod 4] h(-p))/2 failures: ", bad);
  for (c = 1, 4,
    print("\n", curves[c][1], " (degree ", curves[c][4], " over the j-line): F_p-rational supersingular points");
    for (r = 1, 7, if (r % 2 == 0 || !mapisdefined(data[c], r), next);
      my(lst = Vec(mapget(data[c], r)), ratios, ok);
      \\ within a class mod 8, h(-4p) = 2 n_p (p = 1 mod 4), h(-4p) = h(-p) = n_p (p = 7 mod 8),
      \\ h(-4p) = 3 h(-p), n_p = 2 h(-p) (p = 3 mod 8): so test count = gamma n_p with gamma constant
      ratios = Set(vector(#lst, i, lst[i][5] / lst[i][2]));
      ok = (#ratios == 1);
      print("   p = ", r, " mod 8 (", #lst, " primes): count / n_p = ", if (ok, ratios[1], ratios), if (ok, "   constant on every prime", "   NOT constant"),
        ";  [p, n_p, h(-4p), h(-p), count]: ", lst[1..min(3, #lst)]);
      if (!ok, bad++)));
  \\ refinement for X_1(8): 2-adic eigenvalue argument (results note) predicts, for p = 7 mod 8,
  \\ 4 points per floor j and 8 [p = 15 mod 16] per surface j, i.e. 2 n_p (p = 7 mod 16), 6 n_p (p = 15 mod 16)
  print("\nX_1(8) refined by p mod 16 (prediction: 2 n_p for p = 7 mod 16, 6 n_p for p = 15 mod 16, 0 otherwise):");
  my(dev8 = 0, t8 = 0);
  for (r = 1, 7, if (r % 2 == 0 || !mapisdefined(data[4], r), next);
    my(lst = Vec(mapget(data[4], r)));
    for (i = 1, #lst, my(pp = lst[i][1], pred = if (pp % 16 == 7, 2, if (pp % 16 == 15, 6, 0)) * lst[i][2]);
      t8++; if (lst[i][5] != pred, dev8++; print("   deviation at p = ", pp, ": ", lst[i][5], " vs ", pred))));
  print("   deviations: ", dev8, " of ", t8, " primes");
  if (dev8 == 0, bad--);   \\ the constant-per-class-mod-8 ansatz for X_1(8) is superseded by the mod-16 rule
  bad += dev8;
  print("\nX_1(9): F_p-rational supersingular points against 6 [p = 8 mod 9] n_p:");
  my(dev = 0, tested = 0);
  for (r = 1, 8, if (!mapisdefined(data[5], r), next);
    my(lst = Vec(mapget(data[5], r)));
    for (i = 1, #lst, tested++; if (lst[i][5] != 6 * (r == 8) * lst[i][2], dev++));
    print("   p = ", r, " mod 9: ", #lst, " primes, counts ", vector(min(4, #lst), i, lst[i][5]), " vs predicted ", vector(min(4, #lst), i, 6 * (r == 8) * lst[i][2])));
  print("   deviations: ", dev, " of ", tested);
  bad += dev;
  \\ ---- local reduction: count = (h(-4p) s(floor) + [p = 3 mod 4] h(-p) s(surface)) / 2 ----
  \\ floor O = Z[pi], pi^2 = -p: pi acts on O/N by [0,-p;1,0]; surface O = Z[w], w = (1+pi)/2,
  \\ w^2 = w - (1+p)/4: w acts by [0,-c;1,1], c = (1+p)/4, and pi = 2w - 1 (p = 3 mod 4 only).
  \\ s0 = number of pi-stable cyclic subgroups of order N; s1 = number of +-classes of points of exact
  \\ order N with pi P = +-P.  The global counts above are compared prime by prime.
  print("\n== local reduction: count = (h(-4p) s(floor) + [p = 3 mod 4] h(-p) s(surface))/2 ==");
  my(levels = [[1, 2, 0], [2, 4, 0], [3, 8, 0], [4, 8, 1], [5, 9, 1]], ldev = 0, lt = 0, shown = Map());
  forprime (p = 5, PMAX - 1,
    my(h4 = qfbclassno(-4 * p), h1 = if (p % 4 == 3, qfbclassno(-p), 0));
    for (k = 1, #levels,
      my(c = levels[k][1], N = levels[k][2], typ = levels[k][3], key, lst, actual, sf, ss = 0, pred);
      if (c == 5 && p == 3, next);
      key = if (c == 5, p % 9, p % 8);
      lst = Vec(mapget(data[c], key));
      actual = select(v -> v[1] == p, lst)[1][5];
      sf = slocal(Mod([0, -p; 1, 0], N), N, typ);
      if (p % 4 == 3, my(cc = (1 + p) / 4, W = Mod([0, -cc; 1, 1], N)); ss = slocal(2 * W - 1, N, typ));
      pred = (h4 * sf + h1 * ss) / 2;
      lt++;
      if (pred != actual, ldev++; print("   local reduction FAILS: ", curves[c][1], " p = ", p, ": ", actual, " vs ", pred));
      my(tag = if (c == 5, Str(curves[c][1], " p=", p % 9, " mod 9"), Str(curves[c][1], " p=", p % 32, " mod 32")));
      if (!mapisdefined(shown, tag), mapput(shown, tag, [sf, if (p % 4 == 3, ss, "-")]))));
  print("   prime-by-prime comparisons: ", lt, ", failures: ", ldev);
  print("   local numbers [s(floor), s(surface)] by curve and p mod 32 (p mod 9 for X_1(9)):");
  my(K = Mat(shown));
  for (i = 1, #K[, 1], print("     ", K[i, 1], ": ", K[i, 2]));
  bad += ldev;
  \\ level-N analogue of Ogg's list: primes p where every supersingular point of the curve is F_p-rational
  print("\n== level-N Ogg lists: p < ", PMAX, " with every supersingular point F_p-rational ==");
  my(tot = [(p) -> eichler(p, 2), (p) -> eichler(p, 4), (p) -> eichler(p, 8), (p) -> 24 * (p - 1) / 12, (p) -> 36 * (p - 1) / 12]);
  for (c = 1, 5,
    my(L = List());
    forprime (p = 5, PMAX - 1, if (c == 5 && p == 3, next);
      my(key = if (c == 5, p % 9, p % 8), lst = Vec(mapget(data[c], key)), actual = select(v -> v[1] == p, lst)[1][5]);
      if (actual == tot[c](p), listput(L, p)));
    print("   ", curves[c][1], ": ", Vec(L)));
  \\ the remaining identity levels, by direct root counts on t7, t9 and the Rogers-Selberg coordinate M7
  my(more = [["X_0(7)", ('x^2 + 13 * 'x + 49) * ('x^2 + 245 * 'x + 2401)^3, 'x^7, 7, (p) -> eichler(p, 7)],
             ["X_0(9)", ('x + 9)^3 * ('x^3 + 243 * 'x^2 + 2187 * 'x + 6561)^3, 'x^9 * ('x^2 + 9 * 'x + 27), 3, (p) -> eichler(p, 9)],
             ["X_1(7)", A7, B7, 7, (p) -> 2 * (p - 1)]]);
  for (c = 1, #more,
    my(L = List(), A = more[c][2], B = more[c][3], bp = more[c][4]);
    forprime (p = 5, PMAX - 1, if (p == bp, next);
      my(roots = polrootsmod(ssA(p), p), tot = 0);
      for (i = 1, #roots, tot += #polrootsmod(A - lift(roots[i]) * B, p));
      if (tot == more[c][5](p), listput(L, p)));
    print("   ", more[c][1], ": ", Vec(L)));
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
