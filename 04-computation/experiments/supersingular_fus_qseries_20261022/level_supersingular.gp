\\ Supersingular points on the identity coordinates of X_1(N) and on X_0(N)
\\ (opus-2026-10-22-ssj). Needs level_relations.gp (from level_hauptmoduln.gp) and sspoly.gp.
\\
\\ Coordinates (all with j = A(x)/B(x), deg A > deg B, so x = infinity is a cusp):
\\   X_1(5): x = 1/M5 (Rogers-Ramanujan, control)   X_1(7): x = M7 (Rogers-Selberg)
\\   X_1(8): x = 1/M8 (Goellnitz-Gordon)             X_1(9): x = M9 (Bailey)
\\   X_0(7): t7 = (eta/eta_7)^4,  X_0(8): t8 = 1/M8 + M8,  X_0(9): t9 = (eta/eta_9)^3.
\\ For p >= 5, p not dividing N, the supersingular locus of the curve is the zero set of
\\   Q(x) = sum_k s_k A(x)^k B(x)^(d-k),   ss_p(J) = sum_k s_k J^k,  d = deg ss_p,
\\ (no zero at a cusp because gcd(A, B) = 1 mod p, checked).
\\ Tests:
\\  (T1) #ss(X_1(N)) = d_N (p-1)/12 with d_N = deg A = 12, 24, 24, 36   [no mod-12 correction]
\\  (T2) #ss(X_0(N)) = Eichler: (p-1) psi(N)/12 + nu2/4 (1-(-4/p)) + nu3/3 (1-(-3/p))
\\  (T3) over F_{p^2}: every irreducible factor of the squarefree Q over a generic j has
\\       degree k(p) = order of p in (Z/N)^x/{+-1}; exceptional fibres over j = 0, 1728 listed
\\  (T4) N = 7: Frobenius x -> x^(p^2) equals the diamond operator mu(x) = (x-1)/x or mu^2
\\       on every supersingular x, according to p mod 7
\\  (T5) F_p-rational supersingular points: X_0(7), X_0(9) per F_p-rational ss j against
\\       1 + (-p/l) (l = 7, 3); X_1(7) per F_p-rational ss j against 6 [p = 6 mod 7]; totals
\\       against the class numbers (h(-4p) + [p = 3 mod 4] h(-p))/2.
\\ Run: gp -q -f level_supersingular.gp > level_supersingular.out
default(parisize, 256000000);
read("sspoly.gp");
read("level_relations.gp");
PMAX = 400;

ordpm(p, N) = { my(k = 1, a = p % N); while (a != 1 && a != N - 1, a = (a * p) % N; k++); k };

\\ squarefree supersingular polynomial on a coordinate with j = A/B
ssQ(p, A, B) = {
  my(s = ssA(p), d = poldegree(s), Am = A * Mod(1, p), Bm = B * Mod(1, p), Q);
  if (poldegree(gcd(Am, Bm)) > 0, error("bad reduction of the j-map at p = ", p));
  Q = sum(k = 0, d, polcoeff(s, k) * Am^k * Bm^(d - k));
  Q / gcd(Q, Q');
}

\\ the j-fibre of an irreducible factor h over F_{p^2} (or F_p): 0, 1728, or -1 (other)
fibre(h, A, B, p) = {
  my(Am = A * Mod(1, p), Bm = B * Mod(1, p));
  if (h.mod == 0, 0);
  if (Mod(Am, h) == 0, return(0));
  if (Mod(Am - 1728 * Bm, h) == 0, return(1728));
  -1;
}

eichler(p, N) = {
  my(psi = N * prod(i = 1, #factor(N)[, 1], 1 + 1 / factor(N)[i, 1]), nu2, nu3, fa = factor(N)[, 1]);
  nu2 = if (N % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (N % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (p - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, p)) + nu3 / 3 * (1 - kronecker(-3, p));
}

nssFp(p) = { my(f = factor(ssA(p))); #select(h -> poldegree(h) == 1, f[, 1]) };

{
  my(curves, bad = 0, ntest = 0, nT3 = 0, nT4 = 0, exc = Map());
  curves = [[5, A5, B5, "X_1(5) x=1/M5"], [7, A7, B7, "X_1(7) x=M7"], [8, A8, B8, "X_1(8) x=1/M8"], [9, A9, B9, "X_1(9) x=M9"]];
  print("== (T1) supersingular points on X_1(N) and (T3) F_{p^2}-orbit lengths, primes 5 <= p < ", PMAX, " ==");
  for (c = 1, #curves,
    my(N = curves[c][1], A = curves[c][2], B = curves[c][3], dN = poldegree(curves[c][2]), cnt = 0, ok3 = 0, excs = List());
    forprime (p = 5, PMAX - 1,
      if (N % p == 0, next);
      my(Q = ssQ(p, A, B), k = ordpm(p, N), w = ffgen([p, 2], 'w), fac, gen, ex0, ex1728);
      cnt++;
      if (poldegree(Q) != dN * (p - 1) / 12, bad++; print("  T1 FAIL ", N, " p=", p, ": ", poldegree(Q)));
      \\ factor over F_{p^2}
      fac = factor(Q * w^0)[, 1];
      gen = List(); ex0 = List(); ex1728 = List();
      for (i = 1, #fac,
        my(h = fac[i], dg = poldegree(h), fb);
        \\ fibre test over F_{p^2}
        fb = if (subst(A, 'x, 'x) * w^0 % h == 0, 0, if ((A - 1728 * B) * w^0 % h == 0, 1728, -1));
        if (fb == -1, listput(gen, dg), if (fb == 0, listput(ex0, dg), listput(ex1728, dg))));
      if (#gen && (vecmin(Vec(gen)) != k || vecmax(Vec(gen)) != k), bad++; print("  T3 FAIL N=", N, " p=", p, " generic degrees ", Vec(gen), " k=", k), ok3++);
      if ((#ex0 && vecmax(Vec(ex0)) != vecmin(Vec(ex0))) || (#ex1728 && vecmax(Vec(ex1728)) != vecmin(Vec(ex1728))) || (#ex0 && Vec(ex0)[1] != k) || (#ex1728 && Vec(ex1728)[1] != k),
        listput(excs, [p, p % N, k, vecsort(Vec(ex0)), vecsort(Vec(ex1728))])));
    nT3 += ok3;
    print("  ", curves[c][4], ": d_N = ", dN, "; #ss = d_N (p-1)/12 on all ", cnt, " primes; generic F_{p^2}-degrees = k(p) on ", ok3);
    if (#excs, print("     exceptional fibres [p, p mod N, k(p), degrees over j=0, degrees over j=1728] (first 6 of ", #excs, "): ", Vec(excs)[1..min(6, #excs)]),
      print("     no exceptional fibre: degrees over j = 0, 1728 also equal k(p)"));
  );

  print("\n== (T2) supersingular points on X_0(N): Eichler's count ==");
  my(x0 = [[7, (t7 -> 0), "X_0(7)"], [8, 0, "X_0(8)"], [9, 0, "X_0(9)"]]);
  forprime (p = 5, PMAX - 1,
    \\ X_0(7): j = (t^2+13t+49)(t^2+245t+2401)^3 / t^7
    if (p != 7, my(Q = ssQ(p, ('x^2 + 13 * 'x + 49) * ('x^2 + 245 * 'x + 2401)^3, 'x^7));
      ntest++; if (poldegree(Q) != eichler(p, 7), bad++; print("  T2 FAIL X_0(7) p=", p)));
    my(Q8 = ssQ(p, P8, Q8));
    ntest++; if (poldegree(Q8) != eichler(p, 8), bad++; print("  T2 FAIL X_0(8) p=", p));
    my(Q9 = ssQ(p, ('x + 9)^3 * ('x^3 + 243 * 'x^2 + 2187 * 'x + 6561)^3, 'x^9 * ('x^2 + 9 * 'x + 27)));
    ntest++; if (poldegree(Q9) != eichler(p, 9), bad++; print("  T2 FAIL X_0(9) p=", p)));
  print("  X_0(7), X_0(8), X_0(9): Eichler count confirmed on ", ntest, " (curve, prime) pairs");
  print("  formulas: X_0(7): (2(p-1) + 4[p=2 mod 3])/3;  X_0(8): p-1;  X_0(9): p-1");

  print("\n== (T4) Frobenius x -> x^(p^2) on supersingular x equals a diamond operator ==");
  \\ diamond generators on the identity coordinates: N = 7: <2> = mu7(x) = (x-1)/x; N = 9: <2> = MU9;
  \\ N = 8: <3> = 1/x.  Record the exponent e with x^(p^2) = mu^e(x) for each residue class of p.
  my(gens = [[7, A7, B7, [-1, 1; -1, 0], 3, 2], [9, A9, B9, MU9, 3, 2], [8, A8, B8, [0, 1; 1, 0], 2, 3]]);
  for (c = 1, #gens,
    my(N = gens[c][1], A = gens[c][2], B = gens[c][3], mu = gens[c][4], ordmu = gens[c][5], gen = gens[c][6], tally = Map(), nT = 0);
    forprime (p = 5, PMAX - 1,
      if (N % p == 0, next);
      my(Q = ssQ(p, A, B), fac = factor(Q)[, 1], rule = -1, r = p % N, epred, gp = 1);
      \\ predicted exponent: <p> = <gen>^e in (Z/N)^x/(+-1)
      epred = -1;
      for (e = 0, ordmu - 1, if (gp == r || gp == N - r, epred = e; break); gp = (gp * gen) % N);
      for (i = 1, #fac,
        my(h = fac[i], X = Mod('x, h), F = X^(p^2), which = -1, Y = X);
        for (e = 0, ordmu - 1,
          if (F == Y, which = e; break);
          Y = (mu[1, 1] * Y + mu[1, 2]) / (mu[2, 1] * Y + mu[2, 2]));
        if (which == -1, bad++; print("  T4 FAIL N=", N, " p=", p));
        \\ fixed points of mu satisfy every rule; only record factors not fixed by mu
        if ((mu[1, 1] * X + mu[1, 2]) / (mu[2, 1] * X + mu[2, 2]) != X,
          if (rule == -1, rule = which, if (rule != which, bad++; print("  T4 inconsistent N=", N, " p=", p)))));
      if (rule != -1 && rule != epred, bad++; print("  T4 rule != <p> prediction N=", N, " p=", p, ": ", rule, " vs ", epred));
      nT++;
      if (mapisdefined(tally, r), if (mapget(tally, r) != rule, bad++; print("  T4 class rule changes N=", N, " p=", p)), mapput(tally, r, rule)));
    print("  X_1(", N, "): ", nT, " primes; x^(p^2) = <", gen, ">^e (x) with e by p mod ", N, " (predicted: <p> = <", gen, ">^e):");
    my(line = "     ");
    for (r = 1, N - 1, if (mapisdefined(tally, r), line = Str(line, "p=", r, ": e=", mapget(tally, r), "  ")));
    print(line));
  print("  the diamond <d> permutes the Klein forms f_a -> f_(da): on X_1(7) it cycles the three Rogers-Selberg");
  print("  products, on X_1(9) the three Bailey products, on X_1(8) it swaps the two Goellnitz-Gordon products.");

  print("\n== (T5) F_p-rational supersingular points ==");
  my(row = List(), n5 = 0, nspec = 0, dev = 0);
  forprime (p = 5, PMAX - 1,
    if (p == 7 || p == 3, next);
    my(s = ssA(p), roots = polrootsmod(s, p), tot7 = 0, tot9 = 0, tot17 = 0, gen7 = [], gen9 = [], gen17 = [], f7, f9, f17, nF = #roots, h4 = qfbclassno(-4 * p), h1 = if (p % 4 == 3, qfbclassno(-p), 0));
    for (i = 1, #roots,
      my(j = lift(roots[i]));
      \\ X_0(7) fibre: (t^2+13t+49)(t^2+245t+2401)^3 - j t^7
      f7 = #polrootsmod(('x^2 + 13 * 'x + 49) * ('x^2 + 245 * 'x + 2401)^3 - j * 'x^7, p);
      f9 = #polrootsmod(('x + 9)^3 * ('x^3 + 243 * 'x^2 + 2187 * 'x + 6561)^3 - j * 'x^9 * ('x^2 + 9 * 'x + 27), p);
      f17 = #polrootsmod(A7 - j * B7, p);
      tot7 += f7; tot9 += f9; tot17 += f17;
      if (j != 0 && j != 1728 % p,
        if (f7 != 1 + kronecker(-p, 7), bad++; print("  T5 FAIL X_0(7) p=", p, " j=", j));
        if (f9 != 1 + kronecker(-p, 3), bad++; print("  T5 FAIL X_0(9) p=", p, " j=", j));
        if (f17 != 6 * (p % 7 == 6), bad++; print("  T5 FAIL X_1(7) p=", p, " j=", j)),
        nspec++;
        if (f7 != 1 + kronecker(-p, 7) || f9 != 1 + kronecker(-p, 3) || f17 != 6 * (p % 7 == 6), dev++);
        listput(row, [p, j, f7, 1 + kronecker(-p, 7), f9, 1 + kronecker(-p, 3), f17, 6 * (p % 7 == 6)])));
    n5++;
    if (nF != (h4 + h1) / 2, bad++; print("  T5 FAIL class number count p=", p));
    if (tot7 != (1 + kronecker(-p, 7)) * nF, listput(row, [p, "X_0(7) total", tot7, (1 + kronecker(-p, 7)) * nF]));
    if (tot9 != (1 + kronecker(-p, 3)) * nF, listput(row, [p, "X_0(9) total", tot9, (1 + kronecker(-p, 3)) * nF])));
  print("  generic F_p-rational ss j (j != 0, 1728): fibre counts equal 1+(-p/7) on X_0(7), 1+(-p/3) on X_0(9),");
  print("  6[p = 6 mod 7] on X_1(7), on all ", n5, " primes (failures counted below)");
  print("  special fibres j in {0, 1728} that are F_p-rational and supersingular: ", nspec, "; deviations from the generic rule: ", dev);
  print("  totals: #X_0(7)^ss(F_p) = (1+(-p/7)) n_p and #X_0(9)^ss(F_p) = (1+(-p/3)) n_p, n_p = (h(-4p) + [p=3 mod 4] h(-p))/2,");
  print("          on all primes unless listed with the word total below");
  print("  fibres over j = 0, 1728 [p, j, X_0(7) got, generic rule, X_0(9) got, rule, X_1(7) got, rule]:");
  for (i = 1, min(#row, 40), print("     ", row[i]));
  if (#row > 40, print("     ... ", #row - 40, " more"));
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
