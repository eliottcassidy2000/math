\\ Supersingular j-invariants in characteristic p: count, F_p / F_{p^2} split,
\\ class numbers, X_0(p)^+ genus, Ogg's list (opus-2026-10-22-ssj).
\\
\\ Route A (Hasse invariant = E_{p-1} mod p): write p - 1 = 12m + 4de + 6ep with
\\   de, ep in {0,1}. Since E_{p-1} == 1 (mod p) as a q-expansion and
\\   E_{p-1} = Delta^m E4^de E6^ep f(j) with deg f = m, the polynomial f mod p is the
\\   polynomial part (in j) of 1/(Delta^m E4^de E6^ep); then
\\   ss_p(j) = j^de (j - 1728)^ep f(j).
\\ Route B (Legendre line): ss_p is the squarefree part of
\\   Res_lambda( H_p(lambda), 256(lambda^2 - lambda + 1)^3 - j lambda^2 (lambda - 1)^2 ),
\\   H_p(lambda) = sum_k C(m2,k)^2 lambda^k, m2 = (p-1)/2  (Deuring/Igusa).
\\ Route C (direct): #{j in F_p : a_p(E_j) = 0} by point counting (ellap), p >= 5.
\\ Run: gp -q -f supersingular_census.gp > supersingular_census.out
default(parisize, 64000000);
PA = 3000;   \\ routes A and C, all primes 5 <= p < PA
PB = 700;    \\ route B cross-check, 5 <= p < PB

read("sspoly.gp");

ssB(p) = {
  my(m2 = (p - 1) / 2, H, R, g);
  H = sum(k = 0, m2, Mod(binomial(m2, k)^2, p) * 'x^k);
  R = polresultant(H, 256 * ('x^2 - 'x + 1)^3 - 'y * 'x^2 * ('x - 1)^2, 'x);
  R = subst(R, 'y, 'x);
  g = gcd(R, R');
  R = R / g;
  R / pollead(R);
}

countC(p) = {
  my(c = 0, E);
  for (j = 0, p - 1,
    E = ellinit(ellfromj(Mod(j, p)));
    if (ellap(E) == 0, c++));
  c;
}

\\ independent class number: reduced primitive forms of discriminant D < 0
hred(D) = {
  my(h = 0, c);
  for (a = 1, sqrtint(-D \ 3),
    for (b = -a + 1, a,
      if ((b^2 - D) % (4 * a), next);
      c = (b^2 - D) / (4 * a);
      if (c < a, next);
      if (c == a && b < 0, next);
      if (gcd([a, b, c]) != 1, next);
      h++));
  h;
}

kro(a, p) = kronecker(a, p);
g0(p) = 1 + (p + 1) / 12 - (1 + kro(-4, p)) / 4 - (1 + kro(-3, p)) / 3 - 1;

{
  my(bad = 0, nA = 0, nB = 0, nC = 0, ogg = [2, 3], tab = List(), mdeg, nlin, nquad, F, h1, h4, fx, pred, gp0, gplus, eps, ssfloor, okdeg);
  forprime (p = 5, PA - 1,
    my(sA = ssA(p), fac, degs);
    nA++;
    \\ degree law floor(p/12) + {0,1,1,2}
    eps = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2][p % 12];   \\ 0, 1, 1, 2 for p = 1, 5, 7, 11 mod 12
    ssfloor = p \ 12 + eps;
    if (poldegree(sA) != ssfloor, bad++; print("degree law fails at ", p));
    \\ factor over F_p
    fac = factor(sA);
    degs = vector(#fac[, 1], k, poldegree(fac[k, 1]));
    if (vecmax(fac[, 2]) > 1, bad++; print("ss_p not squarefree at ", p));
    mdeg = vecmax(degs);
    nlin = #select(d -> d == 1, degs);
    nquad = #select(d -> d == 2, degs);
    if (mdeg > 2, bad++; print("root outside F_{p^2} at ", p));
    \\ class numbers (PARI and independent count)
    h4 = qfbclassno(-4 * p);
    h1 = if (p % 4 == 3, qfbclassno(-p), 0);
    if (p < 3000, if (h4 != hred(-4 * p) || (p % 4 == 3 && h1 != hred(-p)), bad++; print("class number mismatch at ", p)));
    fx = h4 + h1;  \\ fixed points of the Fricke involution w_p on X_0(p)
    pred = fx / 2;
    if (nlin != pred, bad++; print("F_p count != class number prediction at ", p, ": ", nlin, " vs ", pred));
    gp0 = g0(p);
    if (gp0 != mfdim([p, 2], 1), bad++; print("genus formula vs mfdim at ", p));
    if (poldegree(sA) != gp0 + 1, bad++; print("#ss != g + 1 at ", p));
    gplus = (2 * gp0 + 2 - fx) / 4;
    if (nquad != gplus, bad++; print("#pairs != g(X_0(p)^+) at ", p, ": ", nquad, " vs ", gplus));
    if (nquad == 0, ogg = concat(ogg, p));
    if (p < PB, nB++; if (ssB(p) != sA, bad++; print("route A != route B at ", p)));
    if (p < 2000, nC++; if (countC(p) != nlin, bad++; print("route C (ellap) != linear factors at ", p)));
    if (p < 140 || p == 2999 || p == 2971 || setsearch([193, 239, 431, 599, 997, 1999], p),
      listput(tab, [p, p % 12, poldegree(sA), nlin, nquad, h1, h4, gplus]));
  );
  print("primes 5 <= p < ", PA, ": ", nA, " (route A); route B cross-check on ", nB, "; route C (ellap) on ", nC);
  print("columns: p, p mod 12, #ss = deg ss_p, #ss in F_p, #pairs in F_{p^2}\\F_p, h(-p) [p=3 mod 4], h(-4p), g(X_0(p)^+)");
  for (k = 1, #tab, print("  ", tab[k]));
  print("primes p < ", PA, " with every supersingular j in F_p: ", ogg);
  print("Ogg's list (Ogg 1975):                        ", [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]);
  if (ogg != [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71], bad++; print("Ogg list mismatch"));
  print("failures: ", bad);
  if (bad == 0, print("ALL CHECKS PASSED"));
}
quit;
