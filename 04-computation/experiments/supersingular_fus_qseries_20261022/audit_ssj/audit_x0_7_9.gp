\\ Independent audit of #X_0(N)^ss(F_p) for N = 5, 7, 9 (audit_ssj), as the number of p-Frobenius-stable
\\ orbits of the diamond group <d> (d = 2) on the supersingular points of X_1(N) in Tate normal form over
\\ ss j in F_p (X_0(N) = X_1(N)/<d>).  No class numbers, no X_0(N) Hauptmodul, no local-number formula.
\\ Then: (1 + (-p/l)) n_p (l = 7, 3), and for N = 5 the prediction (1 + (-p/5)) n_p of the note's
\\ theorem; and g(X_0(Np)/w_p) from newform Hecke traces against g(X_0(N)) + (H(p,N) - count)/2.
\\ Run: gp -q -f audit_x0_7_9.gp > audit_x0_7_9.out
default(parisize, 256000000);
PM = 400;
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
tnf_bc(N, t) = {
  if (N == 5, return([t, t]));
  if (N == 7, return([t^3 - t^2, t^2 - t]));
  if (N == 9, my(c = t^2 * (t - 1)); return([c * (t^2 - t + 1), c]));
}
tnf_curve(N, t) = { my(bc = tnf_bc(N, t), b = bc[1], c = bc[2]); [1 - c, -b, -b, 0, 0] };
tnf_param2(N, a, Q) = {
  my(x0 = Q[1], y0 = Q[2], A1, A2, A3, A4, A6, s, B1, B2, B3, u, b, c, t);
  A1 = a[1]; A2 = a[2] + 3 * x0; A3 = a[3] + x0 * a[1] + 2 * y0;
  A4 = a[4] + 2 * x0 * a[2] - y0 * a[1] + 3 * x0^2;
  A6 = a[5] + x0 * a[4] + x0^2 * a[2] + x0^3 - y0 * a[3] - y0^2 - x0 * y0 * a[1];
  if (A6 != 0, error("point not on curve"));
  s = A4 / A3; B1 = A1 + 2 * s; B2 = A2 - s * A1 - s^2; B3 = A3;
  u = B3 / B2; b = -B2 / u^2; c = 1 - B1 / u;
  if (N == 5, return(b));
  if (N == 7, return(b / c));
  if (N == 9, return(c^2 / (b - c)));
}
jrat(N) = {
  my(a = tnf_curve(N, 't), b2, b4, b6, b8, c4, D, J, A, B, L);
  b2 = a[1]^2 + 4 * a[2]; b4 = 2 * a[4] + a[1] * a[3]; b6 = a[3]^2 + 4 * a[5];
  b8 = a[1]^2 * a[5] + 4 * a[2] * a[5] - a[1] * a[3] * a[4] + a[2] * a[3]^2 - a[4]^2;
  c4 = b2^2 - 24 * b4; D = -b2^2 * b8 - 8 * b4^3 - 27 * b6^2 + 9 * b2 * b4 * b6;
  J = c4^3 / D; A = numerator(J); B = denominator(J); L = denominator(content([A, B])); [A * L, B * L];
}
hred(D) = { my(h = 0, c); for (a = 1, sqrtint(-D \ 3), for (b = -a + 1, a, if ((b^2 - D) % (4 * a), next); c = (b^2 - D) / (4 * a);
  if (c < a || (c == a && b < 0) || gcd([a, b, c]) != 1, next); h++)); h };
\\ count of F_p-points of X_0(N) = X_1(N)/<2> over the supersingular j in F_p
x0count(N, p, AB) = {
  my(cnt = 0, Am = AB[1] * Mod(1, p), Bm = AB[2] * Mod(1, p));
  for (j0 = 0, p - 1,
    if (ellap(ellinit(ellfromj(Mod(j0, p)))) != 0, next);
    my(f = Am - j0 * Bm, fs, fa);
    fs = f / gcd(f, f'); fa = factor(fs)[, 1];
    for (i = 1, #fa, my(h = fa[i], t0 = Mod('t, h), aa, orb, cur, Fr, st = 0);
      \\ the <2>-orbit of t0 inside F_p[t]/(h) = F_{p^deg h}
      orb = List(); cur = t0;
      for (r = 1, 6, my(isnew = 1); for (u = 1, #orb, if (orb[u] == cur, isnew = 0)); if (isnew, listput(orb, cur));
        aa = tnf_curve(N, cur); cur = tnf_param2(N, aa, emul(aa, 2, [0 * t0, 0 * t0])));
      Fr = t0^p;
      for (u = 1, #orb, if (orb[u] == Fr, st = 1));
      if (st, cnt += poldegree(h) / #orb)));
  cnt;
}
H(l, M) = {
  my(fa = factor(M)[, 1], psi = M * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  nu2 = if (M % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i])));
  nu3 = if (M % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (l - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, l)) + nu3 / 3 * (1 - kronecker(-3, l));
}
nminus(L, p) = { my(d = mfdim([L, 2], 0)); if (d == 0, 0, (d - mfcoef(mftraceform([L, 2], 0), p)) / 2) };
gplus_newforms(N, p) = mfdim([N, 2], 1) + sumdiv(N, M, numdiv(N / M) * nminus(M * p, p));
{
  my(bad = 0, rows = List(), ng = 0, AB = vector(9));
  foreach ([5, 7, 9], N, AB[N] = jrat(N));
  print("== #X_0(N)^ss(F_p) as Frobenius-stable <2>-orbits on Tate-normal-form X_1(N), 5 <= p < ", PM, " ==");
  forprime (p = 5, PM - 1,
    my(h4 = hred(-4 * p), h1 = if (p % 4 == 3, hred(-p), 0), np = (h4 + h1) / 2, got = vector(3), pred = vector(3), Ns = [5, 7, 9], ls = [5, 7, 3]);
    for (i = 1, 3, my(N = Ns[i]); if (N % p == 0, next);
      got[i] = x0count(N, p, AB[N]);
      pred[i] = (1 + kronecker(-p, ls[i])) * np);
    if (p == 5, got[1] = pred[1] = "-"); if (p == 7, got[2] = pred[2] = "-");
    if (got != pred, bad++; print("  MISMATCH p = ", p, ": ", got, " vs ", pred));
    if (p < 50, listput(rows, [p, np, got])));
  print("  [p, n_p, [#X_0(5), #X_0(7), #X_0(9)]] (sample p < 50):");
  for (i = 1, #rows, print("     ", rows[i]));
  print("  predictions (1 + (-p/l)) n_p with l = 5, 7, 3: failures ", bad);
  print("\n== genus relation for N = 7, 9 (and 5) via newform traces ==");
  foreach ([5, 7, 9], N,
    forprime (p = 5, 60, if (N % p == 0 || N * p > 450, next);
      my(c = x0count(N, p, AB[N]), lhs = gplus_newforms(N, p), rhs = mfdim([N, 2], 1) + (H(p, N) - c) / 2);
      ng++;
      if (lhs != rhs, bad++; print("  MISMATCH N=", N, " p=", p, ": ", lhs, " vs ", rhs))));
  print("  pairs checked: ", ng);
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
