\\ Independent re-verification of the 17 DS15 identities used in qseries_identities.py, in PARI power
\\ series to O(q^K) with K = 1200 (the script used O(q^701) in pure Python), transcribed directly from
\\ the DS15 PDF (2nd ed. 2022) equations 2.7.1-3, 2.8.1, 2.8.4, 2.8.9-10, 2.9.1-3, 2.14.1-3, 2.27.1-4.
\\ Products are built as Pochhammer products (Jacobi triple product form), NOT as theta sums, so the
\\ product side is computed by a different route from the script (which summed f(a,b) directly).
\\ Also: the three theta-family relations and the sum-side corollaries, plus two hostile controls.
\\ Run: gp -q -f audit_qseries.gp > audit_qseries.out
default(parisize, 256000000);
K = 1200;
O1 = O('q^(K + 1));
\\ (a; q^s)_n with a = c q^e (c = +-1):  prod_{k<n} (1 - c q^(e + k s))
po(c, e, s, n) = { my(r = 1 + O1); for (k = 0, n - 1, if (e + k * s > K, break); r *= (1 - c * 'q^(e + k * s))); r };
poinf(c, e, s) = po(c, e, s, (K - e) \ s + 2);
\\ theta product f(-q^a, -q^(N-a)) = (q^a, q^(N-a), q^N; q^N)_inf ; f(q^a, q^(N-a)) = (-q^a, -q^(N-a), q^N; q^N)
thm(a, N) = poinf(1, a, N) * poinf(1, N - a, N) * poinf(1, N, N);
thp(a, N) = poinf(-1, a, N) * poinf(-1, N - a, N) * poinf(1, N, N);
fm(k) = poinf(1, k, k);                      \\ f(-q^k) = (q^k; q^k)_inf
psim = poinf(1, 1, 2) * poinf(1, 4, 4);      \\ psi(-q) = (q^2;q^2)/(-q;q^2) = (q;q^2)(q^4;q^4)  [DS15 (1.1.10)]
psim2 = fm(2) / poinf(-1, 1, 2);             \\ the same via (1.1.10) literally
psip = poinf(1, 2, 2) / poinf(1, 1, 2);      \\ psi(q) = (q^2;q^2)/(q;q^2)
\\ generic sum: sum_n q^E(n) * prod(num) / prod(den)
S(E, num, den) = { my(r = O1, n = 0, e); while ((e = E(n)) <= K, r += 'q^e * num(n) / den(n); n++); r };
chk(name, a, b) = { my(d = a - b, v = valuation(d + O('q^(K+1)), 'q)); if (d == 0 || v > K, print("  OK   ", name), print("  DIFF ", name, " first difference at q^", v)); (d == 0 || v > K) };
{
  my(ok = 1, RS = vector(3), GG = vector(3), JS = vector(3), B = vector(4), R14 = vector(6), D = vector(4));
  ok = ok && chk("psi(-q): (q;q^2)(q^4;q^4) = (q^2;q^2)/(-q;q^2)", psim, psim2);
  \\ 2.7.1-3 Rogers-Selberg
  RS[1] = S(n -> 2*n*(n+1), n -> 1, n -> po(1,2,2,n) * po(-1,1,1,2*n+1));
  RS[2] = S(n -> 2*n*(n+1), n -> 1, n -> po(1,2,2,n) * po(-1,1,1,2*n));
  RS[3] = S(n -> 2*n^2,     n -> 1, n -> po(1,2,2,n) * po(-1,1,1,2*n));
  for (a = 1, 3, ok = chk(Str("DS15 2.7.", a, "  f(-q^", a, ",-q^", 7-a, ")/f(-q^2)"), RS[a], thm(a, 7) / fm(2)) && ok);
  \\ 2.8.1, 2.8.4 Goellnitz-Gordon
  GG[1] = S(n -> n*(n+2), n -> po(-1,1,2,n), n -> po(1,2,2,n));
  GG[3] = S(n -> n^2,     n -> po(-1,1,2,n), n -> po(1,2,2,n));
  ok = chk("DS15 2.8.1  f(-q,-q^7)/psi(-q)", GG[1], thm(1, 8) / psim) && ok;
  ok = chk("DS15 2.8.4  f(-q^3,-q^5)/psi(-q)", GG[3], thm(3, 8) / psim) && ok;
  \\ 2.8.9, 2.8.10 Jackson-Slater
  JS[1] = S(n -> 2*n*(n+1), n -> 1, n -> po(1,1,1,2*n+1));
  JS[3] = S(n -> 2*n^2,     n -> 1, n -> po(1,1,1,2*n));
  ok = chk("DS15 2.8.9  f(q,q^7)/f(-q^2)", JS[1], thp(1, 8) / fm(2)) && ok;
  ok = chk("DS15 2.8.10 f(q^3,q^5)/f(-q^2)", JS[3], thp(3, 8) / fm(2)) && ok;
  \\ 2.9.1-3 Bailey mod 9
  B[1] = S(n -> 3*n*(n+1), n -> po(1,1,1,3*n+1), n -> po(1,3,3,n) * po(1,3,3,2*n+1));
  B[2] = S(n -> 3*n*(n+1), n -> po(1,1,1,3*n) * (1 - 'q^(3*n+2)), n -> po(1,3,3,n) * po(1,3,3,2*n+1));
  B[4] = S(n -> 3*n^2,     n -> po(1,1,1,3*n), n -> po(1,3,3,n) * po(1,3,3,2*n));
  foreach ([1, 2, 4], a, ok = chk(Str("DS15 2.9.", if (a == 4, 3, a), "  f(-q^", a, ",-q^", 9-a, ")/f(-q^3)"), B[a], thm(a, 9) / fm(3)) && ok);
  \\ 2.14.1-3 Rogers mod 14
  R14[2] = S(n -> n*(n+2), n -> 1, n -> po(1,1,2,n+1) * po(1,1,1,n));
  R14[4] = S(n -> n*(n+1), n -> 1, n -> po(1,1,2,n+1) * po(1,1,1,n));
  R14[6] = S(n -> n^2,     n -> 1, n -> po(1,1,2,n) * po(1,1,1,n));
  foreach ([2, 4, 6], a, ok = chk(Str("DS15 2.14.", a/2, " f(-q^", a, ",-q^", 14-a, ")/f(-q)"), R14[a], thm(a, 14) / fm(1)) && ok);
  \\ 2.27.1-4 Dyson mod 27
  D[1] = S(n -> n*(n+3), n -> po(1,3,3,n), n -> po(1,1,1,2*n+2) * po(1,1,1,n));
  D[2] = S(n -> n*(n+2), n -> po(1,3,3,n), n -> po(1,1,1,2*n+2) * po(1,1,1,n));
  D[3] = S(n -> n*(n+1), n -> po(1,3,3,n), n -> po(1,1,1,2*n+1) * po(1,1,1,n));
  D[4] = 1 + S(n -> (n+1)^2, n -> po(1,3,3,n), n -> po(1,1,1,2*n+1) * po(1,1,1,n+1));   \\ n -> n+1 shift of 2.27.4
  for (b = 1, 4, ok = chk(Str("DS15 2.27.", b, " f(-q^", 3*b, ",-q^", 27-3*b, ")/f(-q)"), D[b], thm(3*b, 27) / fm(1)) && ok);
  \\ hostile controls (must differ)
  print("  hostile controls (must DIFF):");
  chk("  RS_1 sum against f(-q^2,-q^5)/f(-q^2)", RS[1], thm(2, 7) / fm(2));
  chk("  Dyson_1 sum against f(-q^6,-q^21)/f(-q)", D[1], thm(6, 27) / fm(1));
  \\ theta families and sum-side corollaries
  print("  theta families / sum-side corollaries:");
  for (a = 1, 3, ok = chk(Str("  RS_", a, "(q^2) f(-q^4) = R14_", 2*a, "(q) f(-q)"), subst(RS[a], 'q, 'q^2) * fm(4) + O1, R14[2*a] * fm(1)) && ok);
  my(bb = [1, 2, 4]);
  for (i = 1, 3, my(a = bb[i]); ok = chk(Str("  B_", a, "(q^3) f(-q^9) = D_", a, "(q) f(-q)   [D_b has numerator theta_27(3b)]"), subst(B[a], 'q, 'q^3) * fm(9) + O1, D[a] * fm(1)) && ok);
  foreach ([1, 3], a, ok = chk(Str("  GG_", a, "(-q) psi(q) = JS_", a, "(q) f(-q^2)"), subst(GG[a], 'q, -'q) * psip, JS[a] * fm(2)) && ok);
  ok = chk("  Dyson 2.27.3 product = (q^9;q^9)/(q;q)", thm(9, 27) / fm(1), fm(9) / fm(1)) && ok;
  print(if (ok, "ALL AUDIT CHECKS PASSED (K = 1200)", "AUDIT FOUND FAILURES"));
}
quit;
