\\ audit_wave2: DIRECT test of the rationality property (no genus formula, no modular symbols, no mf*).
\\ For each level N, an explicit Hauptmodul s of X_0^+(N) (quotient by ALL Atkin-Lehner involutions)
\\ is built from eta quotients / Eisenstein series; the curve X_0(N) is written as
\\     F(s, J) = prod_{Q | N} (J - j(Q tau)) = 0,     coefficients in Z[s]
\\ (j o w_Q = j(Q tau) as q-series, so the coefficients are W-invariant; that they ARE polynomials in s,
\\ with zero residual to high order, is the check that s is W-invariant with a single simple pole).
\\ The supersingular points of X_0^+(N) mod p are the roots of R_p(s) = Res_J(F(s, J), ss_p(J)).
\\ ss_p(J) is computed from Deuring's Legendre polynomial and cross-checked by point counting over F_{p^2}.
\\ Then, for every prime p < 200 with p not | N:
\\   (i)  p has the rationality property  <=>  every root of R_p lies in F_p;
\\   (ii) Theorem 1.1 is tested numerically:  #{non-F_p roots} / 2  =  genus(X_0^+(Np)) - genus(X_0^+(N)),
\\        with the genera from the independent Riemann-Hurwitz/class-number count of rh_genus_lists.gp.
\\ Run: gp -q -f direct_rationality.gp > direct_rationality.out
default(parisize, 400000000);
PREC = 140;

\\ ---------- Riemann-Hurwitz genus (copied from rh_genus_lists.gp, validated there) ----------
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

\\ ---------- supersingular polynomial: Deuring's Hasse polynomial in Legendre form ----------
ssp(p) = {
  if (p <= 3, return(Mod(1, p) * y));                       \\ only j = 0 (= 1728 when p = 3)
  my(m = (p - 1) / 2, H = sum(i = 0, m, binomial(m, i)^2 * z^i) * Mod(1, p),
     R = polresultant(H, (256 * (z^2 - z + 1)^3 - y * z^2 * (z - 1)^2) * Mod(1, p), z), fa = factor(R));
  my(P = prod(i = 1, #fa[, 1], fa[i, 1])); P / pollead(P);
}
\\ cross-check by point counting: the set of supersingular j in F_{p^2}
sscheck(p) = {
  my(g = ffgen([p, 2], 'w), P = lift(ssp(p)), Pff, fa, roots = [], cnt = []);
  Pff = sum(i = 0, poldegree(P), polcoef(P, i, y) * g^0 * y^i);
  fa = factor(Pff)[, 1];
  for (i = 1, #fa, if (poldegree(fa[i]) != 1, error("ss j not in F_p^2 at p = ", p)); roots = concat(roots, [(-polcoef(fa[i], 0, y) / polcoef(fa[i], 1, y)).pol]));
  for (k = 0, p^2 - 1, my(jj = subst(Pol(digits(k, p)), 'x, g) + 0 * g);
    if (ellissupersingular(ellinit(ellfromj(jj))), cnt = concat(cnt, [jj.pol])));
  [#roots, #cnt, Set(roots) == Set(cnt)];
}

\\ ---------- q-series ----------
etaq(d) = eta(q^d + O(q^PREC));
etaquo(dv, r) = { my(e = sum(i = 1, #dv, dv[i] * r[i])); if (e % 24, error("eta quotient not integral order"));
  q^(e / 24) * prod(i = 1, #dv, etaq(dv[i])^r[i]); }
D1 = q * etaq(1)^24;
E4 = 1 + 240 * sum(n = 1, PREC, sigma(n, 3) * q^n) + O(q^PREC);
J1 = E4^3 / D1;
jQ(Q) = subst(J1, q, q^Q);
E2(d) = 1 - 24 * sum(n = 1, PREC \ d, sigma(n) * q^(d * n)) + O(q^PREC);

\\ polynomial A in x with f - A(s) = O(q); returns [A, valuation of the residual]
topoly(f, s) = {
  my(A = 0, r = f, v);
  while ((v = valuation(r, q)) <= 0,
    my(c = polcoef(r, v, q), k = -v);                       \\ s = q^-1 + ..., so s^k = q^-k + ...
    A += c * x^k; r -= c * s^k);
  [A, r == 0, serprec(r, q)];                               \\ [polynomial, residual is a zero series?, its precision]
}

\\ Hauptmodul of X_0^+(N), normalised q^-1 + O(1)
hauptplus(N) = {
  my(t, c, k, f);
  if (N == 1, return(J1));
  if (isprime(N) && N != 11,
    k = 24 / (N - 1); t = etaquo([1, N], [k, -k]); c = N^(k / 2); return(t + c / t));
  if (N == 11,
    my(G = -(E2(1) - 11 * E2(11)) / 10, h = etaquo([1, 11], [2, 2]), s = G / h);
    return(s - polcoef(s, 0, q)));
  \\ composite squarefree N = e * (N/e): t = (eta(1) eta(e) / (eta(N/e) eta(N)))^k, invariant under w_e
  my(e = [0, 0, 0, 0, 0, 3, 0, 0, 0, 5, 0, 0, 0, 7, 5, 0, 0, 0, 0, 0, 3, 11][N],
     kk = [0, 0, 0, 0, 0, 6, 0, 0, 0, 4, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 1, 2][N],
     dv = [1, e, N / e, N], r = [kk, kk, -kk, -kk]);
  t = etaquo(dv, r);
  my(c2 = prod(i = 1, 4, (N / dv[i])^r[i]));
  if (!issquare(c2, &c) || type(c) != "t_INT", error("constant not integral ", c2));
  t + c / t;
}

curveF(N) = {
  my(s = hauptplus(N), Qs = divisors(N), js = apply(Q -> jQ(Q), Qs), n = #Qs, F = y^n, info = [], ok = 1);
  \\ integrality and shape of s
  my(cs = vector(PREC \ 2, i, polcoef(s, i - 2, q)));
  if (cs[1] != 1 || #select(c -> denominator(c) != 1, cs), ok = 0);
  for (k = 1, n,
    my(ek = esym(js, k), res = topoly(ek, s));
    if (!res[2] || res[3] < 20, ok = 0);
    info = concat(info, [[k, poldegree(res[1]), if (res[2], Str("0+O(q^", res[3], ")"), "NONZERO")]]);
    F += (-1)^k * res[1] * y^(n - k));
  [F, ok, info, s];
}
\\ elementary symmetric function e_k of a vector of series
esym(v, k) = { my(e = vector(#v + 1)); e[1] = 1; for (i = 1, #v, forstep (m = i + 1, 2, -1, e[m] += v[i] * e[m - 1])); e[k + 1]; }

\\ ---------- main ----------
{
  my(bad = 0);
  print("== A. supersingular polynomial (Deuring) against point counting over F_{p^2} ==");
  forprime (p = 5, 97, my(c = sscheck(p)); if (!c[3], bad++); print1("  p = ", p, ": ", c[1], " ss j, point count ", c[2], if (c[3], " ok", " MISMATCH"), if (p % 30 == 1, "\n", ";")));
  print();
  my(levels = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22],
     note = [[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71], [3,5,7,11,13,17,19,23,31,47], [2,5,7,11,13,17,23,29], [2,3,7,11,19],
             [5,7,11,13], [2,3,5,17], [3,7,11], [2,3,5], [2,3], [3,5], [2,7], [2,5], [3,5]]);
  print("\n== B. Hauptmodul s of X_0^+(N) and the curve F(s, J) = prod_{Q|N} (J - j(Q tau)) ==");
  my(curves = vector(#levels));
  for (k = 1, #levels, my(N = levels[k], C);
    if (N == 1, curves[k] = [y - x, 1, [], J1]; print("  N = 1: s = j"); next);
    C = curveF(N); curves[k] = C;
    if (!C[2], bad++);
    print("  N = ", N, ": s = ", truncate(C[4] + O(q^4)), " + ...;  [k, deg_s e_k, residual]: ", C[3], if (C[2], "  ok", "  FAILED")));
  print("\n== C. direct rationality test and Theorem 1.1, every prime p < 200 with p not | N ==");
  for (k = 1, #levels, my(N = levels[k], F = curves[k][1], L = [], th = 1, maxdeg = 0, gN = gplusRH(N));
    forprime (p = 2, 199, if (N % p == 0, next);
      my(R = polresultant(F * Mod(1, p), ssp(p), y), fa = factor(R), nrat = 0, nnon = 0);
      for (i = 1, #fa[, 1], my(d = poldegree(fa[i, 1])); maxdeg = max(maxdeg, d); if (d == 1, nrat++, nnon += d));
      if (nnon == 0, L = concat(L, p));
      if (nnon / 2 != gplusRH(N * p) - gN, th = 0; print("    Thm 1.1 identity FAILS at N = ", N, ", p = ", p, ": Q_p/2 = ", nnon / 2, ", genus = ", gplusRH(N * p))));
    print("  N = ", N, ": rationality list (p < 200) = ", L, if (L == note[k], "  = note", "  DIFFERS from note"),
          ";  Q_p/2 = genus(X_0^+(Np)) for all p < 200: ", if (th, "yes", "NO"), ";  max factor degree ", maxdeg);
    if (L != note[k] || !th || maxdeg > 2, bad++));
  print("\n== D. row 2- (X = X_0(2) itself): supersingular t-values, t = (eta/eta_2)^24, against Aricheta Table 2 ==");
  my(t2 = etaquo([1, 2], [24, -24]), P2 = topoly(J1 * t2^2, t2), L2 = [], th2 = 1);
  print("  j t^2 = ", P2[1], "  (residual zero: ", P2[2], ", precision ", P2[3], ")   equals (x + 256)^3: ", P2[1] == (x + 256)^3);
  forprime (p = 3, 199,
    my(S = ssp(p), R = 0, fa, nn = 0);
    R = sum(i = 0, poldegree(S), polcoef(S, i, y) * P2[1]^i * x^(2 * (poldegree(S) - i)));
    fa = factor(R * Mod(1, p));
    my(rad = prod(i = 1, #fa[, 1], fa[i, 1]));
    for (i = 1, #fa[, 1], if (poldegree(fa[i, 1]) > 1, nn += poldegree(fa[i, 1])));
    if (nn == 0, L2 = concat(L2, p));
    if (nn / 2 != (2 * g0(2 * p) - 2 - nufix(2 * p, p)) / 4 + 1, th2 = 0);
    if (p <= 29, print("    p = ", p, ": radical of ss poly on X_0(2) = ", lift(rad / pollead(rad)))));
  print("  X_0(2) rationality list p < 200: ", L2, if (L2 == [3, 5, 7, 11, 23], "  = Aricheta 2- [3,5,7,11,23]", "  DIFFERS"),
        ";  Q_p/2 = genus(X_0(2p)/w_p): ", if (th2, "yes", "NO"));
  if (L2 != [3, 5, 7, 11, 23] || !th2, bad++);
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
