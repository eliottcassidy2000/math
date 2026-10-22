\\ Independent audit of sections 6.2 and 6.3 directly from elliptic curves (audit_ssj).
\\ Nothing here uses the Hauptmodul coordinates, the supersingular polynomial, or the local-number
\\ formula of the note, except where stated (Part 4 compares against them).
\\
\\ Part 1 (X_0(2), X_0(4), X_0(8), X_1(8), X_1(4)): for each supersingular j in F_p (ellap = 0) and one
\\   model E/F_p, the full group E[N] is computed inside E(F_{p^4}) = (Z/(p^2-1))^2 (ellgroup generators);
\\   Frobenius = coordinate-wise p-power; Aut(E) acts by explicit maps (x,y) -> (-x, i y) at j = 1728,
\\   (zeta x, y) at j = 0.  Level structures (cyclic subgroups of order N, or +-classes of points of exact
\\   order N) are enumerated as point sets; the F_p-points of the coarse curve over j are the
\\   Frobenius-stable Aut-orbits.  Totals and fibres are compared with the note's 6.3 table.
\\ Part 2 (X_1(N), N = 4,5,7,8,9): Tate normal form E_t (fine moduli for N >= 4): #X_1(N)^ss(F_p) =
\\   #{t in F_p : Delta(t) != 0, a_p(E_t) = 0}, all primes 5 <= p < 400 (p not dividing N).
\\ Part 3 (Proposition 6.2): on every supersingular point (E_t, P) of X_1(N), N = 5,7,8,9, with t in
\\   F_{p^2}[x]/(h) for each irreducible factor h of the fibre polynomial over F_{p^2}, the Tate-normal-form
\\   parameter of (E_t, [p]P) is t^(p^2); factor degrees vs k(p) and the exceptional fibres.
\\ Part 4: g(X_0(Np)/w_p) from Hecke traces on NEWFORMS (no Atkin-Lehner operator):
\\   g^+ = g(X_0(N)) + sum_{M | N} sigma_0(N/M) #{f new of level Mp : a_p(f) = -1}, compared with
\\   g(X_0(N)) + (H(p,N) - #X_0(N)^ss(F_p))/2 using the Part 1 counts (N = 2, 4, 8) and n_p (N = 1).
\\ Run: gp -q -f audit_fp_points.gp > audit_fp_points.out
default(parisize, 512000000);
PA = 104;     \\ Part 1 primes 5 <= p < PA
PB = 400;     \\ Part 2
PC = 100;     \\ Part 3

\\ ---------- generic Weierstrass arithmetic (a = [a1,a2,a3,a4,a6]), points [x,y] or [0] ----------
isO(P) = (#P == 1);
eneg(a, P) = if (isO(P), P, [P[1], -P[2] - a[1] * P[1] - a[3]]);
eadd(a, P, Q) = {
  my(l, n, x3, y3);
  if (isO(P), return(Q)); if (isO(Q), return(P));
  if (P[1] == Q[1],
    if (P[2] + Q[2] + a[1] * Q[1] + a[3] == 0, return([0]));
    l = (3 * P[1]^2 + 2 * a[2] * P[1] + a[4] - a[1] * P[2]) / (2 * P[2] + a[1] * P[1] + a[3]);
    n = (-P[1]^3 + a[4] * P[1] + 2 * a[5] - a[3] * P[2]) / (2 * P[2] + a[1] * P[1] + a[3]),
    l = (Q[2] - P[2]) / (Q[1] - P[1]);
    n = (P[2] * Q[1] - Q[2] * P[1]) / (Q[1] - P[1]));
  x3 = l^2 + a[1] * l - a[2] - P[1] - Q[1];
  y3 = -(l + a[1]) * x3 - n - a[3];
  [x3, y3];
}
emul(a, n, P) = { my(R = [0], B = P); if (n < 0, n = -n; B = eneg(a, P)); while (n, if (n % 2, R = eadd(a, R, B)); B = eadd(a, B, B); n \= 2); R };

\\ ---------- class numbers by reduced forms (own code) ----------
hred(D) = { my(h = 0, c); for (a = 1, sqrtint(-D \ 3), for (b = -a + 1, a, if ((b^2 - D) % (4 * a), next); c = (b^2 - D) / (4 * a);
  if (c < a || (c == a && b < 0) || gcd([a, b, c]) != 1, next); h++)); h };

\\ ---------- Part 1 ----------
\\ level structures on E[N] (N = 2^e) inside E(F_{p^4}); returns [#X_0(N) points over j, #X_1(N) points over j]
\\ integer key of a point over F_{p^4} (FFELT entries cannot be sorted directly)
fk(a) = { my(v = Vec(a.pol)); v = concat(vector(4 - #v), v); v };
ptkey(P) = if (#P == 1, [-1, 0, 0, 0, 0, 0, 0, 0], concat(fk(P[1]), fk(P[2])));
fibre_counts(p, j, N) = {
  my(w = ffgen([p, 4], 'w), A, B, E, G, g1, g2, m, ordN = List(), aut, Fr, res = vector(2));
  if (j == 0, A = 0; B = 1, if (j == 1728 % p, A = 1; B = 0, my(k = Mod(j, p) / (1728 - Mod(j, p))); A = lift(3 * k); B = lift(2 * k)));
  E = ellinit([A * w^0, B * w^0]);
  G = ellgroup(E, , 1);
  if (G[2] != [p^2 - 1, p^2 - 1], error("unexpected group ", G[2], " at p = ", p, " j = ", j));
  m = (p^2 - 1) / N;
  g1 = ellmul(E, G[3][1], m); g2 = ellmul(E, G[3][2], m);
  for (i = 0, N - 1, for (k = 0, N - 1, my(P = elladd(E, ellmul(E, g1, i), ellmul(E, g2, k)));
    if (ellmul(E, P, N / 2) != [0], listput(ordN, P))));          \\ exact order N = 2^e
  if (#ordN != N^2 - (N / 2)^2, error("order-N count"));
  if (j == 1728 % p, my(io = sqrt(-1 * w^0)); aut = (P -> if (#P == 1, P, [-P[1], io * P[2]])),
    if (j == 0, my(z = (-1 + sqrt(-3 * w^0)) / 2); aut = (P -> if (#P == 1, P, [z * P[1], P[2]])),
      aut = (P -> P)));
  Fr = (P -> if (#P == 1, P, [P[1]^p, P[2]^p]));
  for (i = 1, #ordN, if (ellisoncurve(E, aut(ordN[i])) == 0, error("aut not an automorphism")));
  \\ object keys: X_0(N): the cyclic subgroup <P>; X_1(N): the class {P, -P}
  my(key0 = (P -> vecsort(vector(N, i, ptkey(ellmul(E, P, i))))), key1 = (P -> vecsort([ptkey(P), ptkey(ellneg(E, P))])));
  for (typ = 1, 2,
    my(kf = if (typ == 1, key0, key1), done = Map(), nstab = 0);
    for (i = 1, #ordN,
      my(rep = ordN[i], k0 = kf(rep), orb = List(), cur = rep, fo, st = 0);
      if (mapisdefined(done, k0), next);
      for (r = 1, 12, my(kc = kf(cur)); if (!mapisdefined(done, kc), mapput(done, kc, 1); listput(orb, kc)); cur = aut(cur));
      fo = kf(Fr(rep));
      for (u = 1, #orb, if (orb[u] == fo, st = 1));
      nstab += st);
    res[typ] = nstab);
  res;    \\ [X_0(N) count, X_1(N) count] over j
}
surface(p, j) = {
  \\ for p = 3 mod 4 and j != 1728: End_{F_p} is the surface order iff E[2] is F_p-rational (both twists agree)
  if (p % 4 != 3, return(0));
  my(A, B);
  if (j == 0, A = 0; B = 1, my(k = Mod(j, p) / (1728 - Mod(j, p))); A = lift(3 * k); B = lift(2 * k));
  #polrootsmod('x^3 + A * 'x + B, p) == 3;
}
{
  my(bad = 0, nf = 0, fib = Map());
  print("== Part 1: F_p-rational supersingular points from E[N] in E(F_{p^4}), 5 <= p < ", PA, " ==");
  print("[p, p mod 16, n_p, h(-4p), h(-p), #X_0(2), #X_0(4), #X_0(8), #X_1(8), #X_1(4)]  vs note's 6.3 formulas");
  forprime (p = 5, PA - 1,
    my(ssj = List(), h4 = hred(-4 * p), h1 = if (p % 4 == 3, hred(-p), 0), tot = vector(5), pr);
    for (j = 0, p - 1, if (ellap(ellinit(ellfromj(Mod(j, p)))) == 0, listput(ssj, j)));
    foreach (ssj, j,
      my(c2 = fibre_counts(p, j, 2), c4 = fibre_counts(p, j, 4), c8 = fibre_counts(p, j, 8), typ);
      tot += [c2[1], c4[1], c8[1], c8[2], c4[2]];
      typ = if (j == 1728 % p && p % 4 == 3, "1728", if (j == 0, "j0", if (surface(p, j), "surf", "floor")));
      my(key = [typ, p % 32]);
      if (!mapisdefined(fib, key), mapput(fib, key, Set()));
      mapput(fib, key, setunion(mapget(fib, key), [[c2[1], c4[1], c8[1], c8[2]]])));
    pr = [(h4 + 3 * h1) / 2, 3 * h1 * (p % 4 == 3), 4 * h1 * (p % 8 == 7), (2 * (p % 16 == 7) + 6 * (p % 16 == 15)) * h1, 3 * h1 * (p % 4 == 3)];
    if (#ssj != (h4 + h1) / 2, bad++; print("  n_p mismatch p=", p));
    if (tot != pr, bad++; print("  MISMATCH p = ", p, ": got ", tot, " predicted ", pr));
    print("  ", [p, p % 16, #ssj, h4, h1, tot[1], tot[2], tot[3], tot[4], tot[5]]));
  print("  fibre counts [X_0(2), X_0(4), X_0(8), X_1(8)] by fibre type and p mod 32 (sets over all fibres):");
  my(K = Mat(fib)); for (i = 1, #K[,1], print("     ", K[i,1], ": ", K[i,2]));
  print("  failures in Part 1: ", bad);
  P1BAD = bad;
}

\\ ---------- Part 2: Tate normal forms ----------
tnf_bc(N, t) = {
  if (N == 4, return([t, 0]));
  if (N == 5, return([t, t]));
  if (N == 7, return([t^3 - t^2, t^2 - t]));
  if (N == 8, return([(2*t - 1) * (t - 1), (2*t - 1) * (t - 1) / t]));
  if (N == 9, my(c = t^2 * (t - 1)); return([c * (t^2 - t + 1), c]));
  error("N");
}
tnf_curve(N, t) = { my(bc = tnf_bc(N, t), b = bc[1], c = bc[2]); [1 - c, -b, -b, 0, 0] };
{
  my(bad = 0, rows = List());
  print("\n== Part 2: #X_1(N)^ss(F_p) by Tate normal form, 5 <= p < ", PB, " ==");
  \\ sanity: (0,0) has exact order N on E_t for a random rational t
  foreach ([4, 5, 7, 8, 9], N, my(a = tnf_curve(N, 3/7), P = [0, 0]);
    if (!isO(emul(a, N, P)) || sum(k = 1, N - 1, isO(emul(a, k, P))) > 0, bad++; print("  TNF order check fails N=", N)));
  forprime (p = 5, PB - 1,
    my(h4 = hred(-4 * p), h1 = if (p % 4 == 3, hred(-p), 0), np = (h4 + h1) / 2, got = vector(5), pred);
    my(Ns = [4, 5, 7, 8, 9]);
    for (i = 1, 5, my(N = Ns[i]); if (p == N, next);
      for (t = 0, p - 1, if (N == 8 && t == 0, next);
        my(E = ellinit(tnf_curve(N, Mod(t, p))));
        if (#E == 0, next);
        if (ellap(E) == 0, got[i]++)));
    pred = [3 * h1 * (p % 4 == 3), 4 * (p % 5 == 4) * np, 6 * (p % 7 == 6) * np, (2 * (p % 16 == 7) + 6 * (p % 16 == 15)) * h1, 6 * (p % 9 == 8) * np];
    if (p == 5, pred[2] = got[2]); if (p == 7, pred[3] = got[3]);
    if (got != pred, bad++; print("  MISMATCH p = ", p, ": ", got, " vs ", pred));
    if (p < 60 || p % 16 == 15 || p % 9 == 8, listput(rows, concat([p], got))));
  print("  [p, X_1(4) = X_0(4), X_1(5), X_1(7), X_1(8), X_1(9)] (sample):");
  for (i = 1, min(#rows, 30), print("     ", rows[i]));
  print("  predictions: X_0(4): 3[p=3(4)]h(-p); X_1(5): 4[p=4(5)]n_p; X_1(7): 6[p=6(7)]n_p; X_1(8): (2[p=7(16)]+6[p=15(16)])h(-p); X_1(9): 6[p=8(9)]n_p");
  print("  failures in Part 2: ", bad);
  P2BAD = bad;
}

\\ ---------- Part 3: Frobenius = diamond <p> in Tate-normal-form coordinates ----------
tnf_param2(N, a, Q) = {
  my(x0 = Q[1], y0 = Q[2], A1, A2, A3, A4, A6, s, B1, B2, B3, u, b, c, t);
  \\ translate (u,r,s,t) = (1, x0, 0, y0)
  A1 = a[1]; A2 = a[2] + 3 * x0; A3 = a[3] + x0 * a[1] + 2 * y0;
  A4 = a[4] + 2 * x0 * a[2] - y0 * a[1] + 3 * x0^2;
  A6 = a[5] + x0 * a[4] + x0^2 * a[2] + x0^3 - y0 * a[3] - y0^2 - x0 * y0 * a[1];
  if (A6 != 0, error("point not on curve"));
  \\ shear (1, 0, s, 0) with s = A4/A3:  a1' = A1 + 2s, a2' = A2 - s A1 - s^2, a3' = A3, a4' = A4 - s A3 = 0
  s = A4 / A3;
  B1 = A1 + 2 * s; B2 = A2 - s * A1 - s^2; B3 = A3;
  \\ scale by u = B3/B2: a_i -> a_i / u^i, making a2 = a3
  u = B3 / B2;
  b = -B2 / u^2; c = 1 - B1 / u;
  if (B3 / u^3 != -b, error("scaling"));
  if (N == 5, if (b != c, error("N=5 shape")); return(b));
  if (N == 7, t = b / c; if (c != t^2 - t, error("N=7 shape")); return(t));
  if (N == 8, t = b / c; if (b != (2*t - 1) * (t - 1), error("N=8 shape")); return(t));
  if (N == 9, t = c^2 / (b - c); if (c != t^2 * (t - 1), error("N=9 shape")); return(t));
}
jrat(N) = {
  my(a = tnf_curve(N, 't), b2, b4, b6, b8, c4, D);
  b2 = a[1]^2 + 4 * a[2]; b4 = 2 * a[4] + a[1] * a[3]; b6 = a[3]^2 + 4 * a[5];
  b8 = a[1]^2 * a[5] + 4 * a[2] * a[5] - a[1] * a[3] * a[4] + a[2] * a[3]^2 - a[4]^2;
  c4 = b2^2 - 24 * b4; D = -b2^2 * b8 - 8 * b4^3 - 27 * b6^2 + 9 * b2 * b4 * b6;
  c4^3 / D;
}
ordpm(p, N) = { my(k = 1, x = p % N); while (x != 1 && x != N - 1, x = (x * p) % N; k++); k };
{
  my(bad = 0, ntest = 0, exc = List(), Js = vector(9), shown = Map());
  print("\n== Part 3: Frobenius x -> x^(p^2) versus <p> on Tate normal forms, 5 <= p < ", PC, " ==");
  foreach ([5, 7, 8, 9], N, my(J = jrat(N), A = numerator(J), B = denominator(J), L);
    L = denominator(content([A, B])); A *= L; B *= L; Js[N] = [A, B];
    print("  N = ", N, ": deg j(t) = ", max(poldegree(A), poldegree(B)), " (deg num ", poldegree(A), ", deg den ", poldegree(B), ")"));
  forprime (p = 5, PC - 1,
    my(v = ffgen([p, 2], 'v), q = p^2, g = ffprimroot(v), z = 1 * v^0, ssl = List());
    \\ supersingular j in F_{p^2} by group order: trace = q + 1 - #E == 0 (mod p)
    for (i = 0, q - 1, my(jj = if (i == 0, 0 * v, z));
      if (i > 0, z *= g);
      if ((q + 1 - ellcard(ellinit(ellfromj(jj)))) % p == 0, listput(ssl, jj)));
    foreach ([5, 7, 8, 9], N, if (N % p == 0, next);
      my(Am = Js[N][1] * Mod(1, p), Bm = Js[N][2] * Mod(1, p), dN = max(poldegree(Js[N][1]), poldegree(Js[N][2])), k = ordpm(p, N), cnt = 0);
      if (poldegree(gcd(Am, Bm)) > 0, bad++; print("  bad reduction of j(t) at N=", N, " p=", p));
      foreach (ssl, j0,
        my(f = Am * v^0 - j0 * Bm * v^0, fs, fa, degs = List());
        fs = f / gcd(f, f');
        fa = factor(fs)[, 1];
        for (i = 1, #fa, my(h = fa[i], t0 = Mod('t, h), aa, Q, t1);
          listput(degs, poldegree(h));
          aa = tnf_curve(N, t0);
          Q = emul(aa, p % N, [0 * t0, 0 * t0]);
          t1 = tnf_param2(N, aa, Q);
          ntest++;
          if (t1 != t0^(p^2), bad++; print("  FROBENIUS != <p> at N=", N, " p=", p, " j0=", j0)));
        cnt += poldegree(fs);
        my(ds = vecsort(Vec(degs)));
        if (vecmin(ds) != k || vecmax(ds) != k,
          listput(exc, [N, p, p % N, k, if (j0 == 0, "j=0", if (j0 == 1728, "j=1728", "other")), ds])));
      if (cnt != dN * (p - 1) / 12, bad++; print("  #ss(X_1(", N, ")) at p=", p, ": ", cnt, " vs ", dN * (p - 1) / 12))));
  print("  tests of t(E_t, [p]P) = t^(p^2) (one per irreducible factor over F_{p^2}): ", ntest, "; failures: ", bad);
  print("  exceptional fibres (F_{p^2}-degrees not all k(p)), one line per [N, fibre, degrees] pattern [N, p, p mod N, k(p), fibre, degrees]:");
  for (i = 1, #exc, my(e = exc[i], kk = [e[1], e[5], e[6]]);
    if (!mapisdefined(shown, kk), mapput(shown, kk, 1); print("     ", e)));
  print("  number of exceptional fibres: N=5: ", #select(e -> e[1] == 5, Vec(exc)), ", N=7: ", #select(e -> e[1] == 7, Vec(exc)),
        ", N=8: ", #select(e -> e[1] == 8, Vec(exc)), ", N=9: ", #select(e -> e[1] == 9, Vec(exc)),
        "; fibres not over j = 0, 1728 among them: ", #select(e -> e[5] == "other", Vec(exc)));
  print("  failures in Part 3: ", bad);
  P3BAD = bad;
}

\\ ---------- Part 4: genus of X_0(Np)/w_p from newform Hecke traces ----------
H(l, M) = {
  my(fa = factor(M)[, 1], psi = M * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  if (M == 1, return((l - 1) / 12 + (1 - kronecker(-4, l)) / 4 + (1 - kronecker(-3, l)) / 3));
  nu2 = if (M % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (M % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (l - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, l)) + nu3 / 3 * (1 - kronecker(-3, l));
}
nminus(L, p) = { my(d = mfdim([L, 2], 0)); if (d == 0, 0, (d - mfcoef(mftraceform([L, 2], 0), p)) / 2) };
gplus_newforms(N, p) = { my(g = if (N == 1, 0, mfdim([N, 2], 1))); g + sumdiv(N, M, numdiv(N / M) * nminus(M * p, p)) };
countX0(p, N) = {
  \\ brute force: N = 1 -> n_p; N in {2,4,8} -> Part 1 method
  my(ssj = List(), tot = 0);
  for (j = 0, p - 1, if (ellap(ellinit(ellfromj(Mod(j, p)))) == 0, listput(ssj, j)));
  if (N == 1, return(#ssj));
  foreach (ssj, j, tot += fibre_counts(p, j, N)[1]);
  tot;
}
{
  my(bad = 0, n = 0, rows = List());
  print("\n== Part 4: g(X_0(Np)/w_p) from newform traces vs g(X_0(N)) + (H(p,N) - #X_0(N)^ss(F_p))/2 ==");
  foreach ([1, 2, 4, 8], N,
    forprime (p = 5, 103, if (N * p > 420, next);
      my(gp1 = gplus_newforms(N, p), gN = if (N == 1, 0, mfdim([N, 2], 1)), c = countX0(p, N), rhs);
      rhs = gN + (H(p, N) - c) / 2;
      n++;
      if (gp1 != rhs, bad++; print("  MISMATCH N=", N, " p=", p, ": ", gp1, " vs ", rhs));
      if (p < 30 || p == 23 || p == 47, listput(rows, [N, p, gp1, H(p, N), c]))));
  print("  pairs checked: ", n, "; sample [N, p, g(X_0(Np)/w_p), H(p,N), #X_0(N)^ss(F_p)]:");
  for (i = 1, #rows, print("     ", rows[i]));
  print("  failures in Part 4: ", bad);
  P4BAD = bad;
}
print("\nfailures: ", [P1BAD, P2BAD, P3BAD, P4BAD]);
if (P1BAD + P2BAD + P3BAD + P4BAD == 0, print("ALL AUDIT CHECKS PASSED"), print("AUDIT FOUND FAILURES"));
quit;
