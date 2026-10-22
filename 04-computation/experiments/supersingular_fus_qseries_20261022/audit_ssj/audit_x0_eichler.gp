\\ (copied helper block of audit_x0_7_9.gp) Eichler counts #X_0(N)^ss over F_pbar for N = 5, 7, 9 as the number of <2>-orbits on
\\ the supersingular locus of Tate-normal-form X_1(N) over ALL ss j in F_{p^2}; and the <2>-fixed points over j = 0 for N = 7
\\ (the note: 'the two linear factors lie over the two order-3 elliptic points of X_0(7)').
\\ Original header follows.
\\ Independent audit of #X_0(N)^ss(F_p) for N = 5, 7, 9 (audit_ssj), as the number of p-Frobenius-stable
\\ orbits of the diamond group <d> (d = 2) on the supersingular points of X_1(N) in Tate normal form over
\\ ss j in F_p (X_0(N) = X_1(N)/<d>).  No class numbers, no X_0(N) Hauptmodul, no local-number formula.
\\ Then: (1 + (-p/l)) n_p (l = 7, 3), and for N = 5 the prediction (1 + (-p/5)) n_p of the note's
\\ theorem; and g(X_0(Np)/w_p) from newform Hecke traces against g(X_0(N)) + (H(p,N) - count)/2.
\\ Run: gp -q -f audit_x0_eichler.gp > audit_x0_eichler.out
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

eich(p, N) = { my(fa = factor(N)[, 1], psi = N * prod(i = 1, #fa, 1 + 1 / fa[i]), nu2, nu3);
  nu2 = if (N % 4 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-4, fa[i]))); nu3 = if (N % 9 == 0, 0, prod(i = 1, #fa, 1 + kronecker(-3, fa[i])));
  (p - 1) * psi / 12 + nu2 / 4 * (1 - kronecker(-4, p)) + nu3 / 3 * (1 - kronecker(-3, p)) };
{
  my(bad = 0, fixed7 = List());
  my(AB = vector(9)); foreach ([5, 7, 9], N, AB[N] = jrat(N));
  forprime (p = 5, 60,
    my(v = ffgen([p, 2], 'v), q = p^2, g = ffprimroot(v), z = 1 * v^0, ssl = List());
    for (i = 0, q - 1, my(jj = if (i == 0, 0 * v, z)); if (i > 0, z *= g);
      if ((q + 1 - ellcard(ellinit(ellfromj(jj)))) % p == 0, listput(ssl, jj)));
    foreach ([5, 7, 9], N, if (N % p == 0, next);
      my(Am = AB[N][1] * Mod(1, p), Bm = AB[N][2] * Mod(1, p), orbits = 0);
      foreach (ssl, j0,
        my(f = Am * v^0 - j0 * Bm * v^0, fs = f / gcd(f, f'), fa = factor(fs)[, 1]);
        for (i = 1, #fa, my(h = fa[i], t0 = Mod('t, h), orb = List(), cur = t0, aa);
          for (r = 1, 6, my(isnew = 1); for (u = 1, #orb, if (orb[u] == cur, isnew = 0)); if (isnew, listput(orb, cur));
            aa = tnf_curve(N, cur); cur = tnf_param2(N, aa, emul(aa, 2, [0 * t0, 0 * t0])));
          orbits += poldegree(h) / #orb;
          if (N == 7 && j0 == 0 && poldegree(h) == 1, listput(fixed7, [p, #orb]))));
      if (orbits != eich(p, N), bad++; print("  Eichler count fails N=", N, " p=", p, ": ", orbits, " vs ", eich(p, N)))));
  print("Eichler counts #X_0(N)^ss = number of <2>-orbits on X_1(N)^ss (N = 5, 7, 9; 5 <= p <= 60): failures ", bad);
  print("N = 7, linear F_{p^2}-factors over j = 0 and their <2>-orbit sizes [p, orbit size] (1 = fixed = elliptic point of X_0(7)):");
  print("   ", Vec(fixed7));
  \\ the claim concerns k(p) = 3 (p != +-1 mod 7): then exactly two linear factors over j = 0, both <2>-fixed
  my(k3 = select(x -> x[1] % 7 != 1 && x[1] % 7 != 6, Vec(fixed7)));
  print("   for k(p) = 3 the linear factors over j = 0 are all <2>-fixed: ", #select(x -> x[2] != 1, k3) == 0, " (", #k3, " factors)");
  if (#select(x -> x[2] != 1, k3) == 0 && bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
