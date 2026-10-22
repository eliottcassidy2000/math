\\ Independent NUMERICAL audit of section 5 (Hauptmoduln) at points of the upper half plane (audit_ssj).
\\ The infinite products are evaluated directly at q = exp(2 pi i tau) and compared with ellj(tau);
\\ no q-expansion linear algebra is used.  Polynomials A_N, B_N, P8, Q8 are read from the committed
\\ level_relations.gp (data under audit).
\\ Also: Gamma_1(N)-invariance of the coordinates, the norm maps t7, t9, t8, the classical level-7/9
\\ j-maps, cusp factorizations / field discriminants, and the DIRECTION of the diamond action on
\\ M7, M9 under the standard complex convention (<d> induced by gamma in Gamma_0(N) whose LOWER-RIGHT
\\ entry is d, i.e. tau -> (C/(Z+Z tau), 1/N)).
\\ Run: gp -q -f audit_hauptmoduln_numeric.gp > audit_hauptmoduln_numeric.out
default(parisize, 256000000);
\p 80
read("../level_relations.gp");
TOL = 1e-50;
qq(t) = exp(2 * Pi * I * t);
\\ Klein-form numerator f_g = prod_{m>=1} (1 - q^(N(m-1)+g)) (1 - q^(Nm-g)), evaluated numerically
kf(N, g, q) = { my(r = 1., m = 1, a, b); while (1, a = q^(N*(m-1)+g); b = q^(N*m-g); r *= (1 - a) * (1 - b); if (abs(a) + abs(b) < 1e-90, break); m++); r };
M5(t) = { my(q = qq(t)); q * (kf(5,1,q) / kf(5,2,q))^5 };
M7(t) = { my(q = qq(t)); q^-1 * kf(7,2,q)^2 * kf(7,3,q) / kf(7,1,q)^3 };
M8(t) = { my(q = qq(t)); q * (kf(8,1,q) / kf(8,3,q))^2 };
M9(t) = { my(q = qq(t)); q^-1 * kf(9,2,q) * kf(9,4,q) / kf(9,1,q)^2 };
etaq(t) = eta(t, 1);                          \\ Dedekind eta with the q^(1/24)
rel(a, b) = abs(a - b) / max(1, abs(b));
act(g, t) = (g[1,1] * t + g[1,2]) / (g[2,1] * t + g[2,2]);
{
  my(bad = 0, taus = [0.13 + 0.6*I, 0.31 + 0.9*I, -0.2 + 0.7*I], e);
  print("== (1) j = A(x)/B(x) at three points tau (relative errors) ==");
  foreach (taus, t,
    my(j = ellj(t), x5 = 1 / M5(t), x7 = M7(t), x8 = 1 / M8(t), x9 = M9(t), t8, e5, e7, e8, e9, e8b);
    e5 = rel(subst(A5, 'x, x5) / subst(B5, 'x, x5), j);
    e7 = rel(subst(A7, 'x, x7) / subst(B7, 'x, x7), j);
    e8 = rel(subst(A8, 'x, x8) / subst(B8, 'x, x8), j);
    e9 = rel(subst(A9, 'x, x9) / subst(B9, 'x, x9), j);
    t8 = x8 + 1 / x8;
    e8b = rel(subst(P8, 'x, t8) / subst(Q8, 'x, t8), j);
    print("  tau = ", t, ":  N=5 ", e5, "  N=7 ", e7, "  N=8 ", e8, "  N=9 ", e9, "  X_0(8) via t8 ", e8b);
    if (vecmax([e5, e7, e8, e9, e8b]) > TOL, bad++));
  print("== (2) norm maps and classical level-7, level-9 j-maps ==");
  foreach (taus, t,
    my(x7 = M7(t), x9 = M9(t), t7 = (etaq(t) / etaq(7*t))^4, t9 = (etaq(t) / etaq(9*t))^3, ea, eb, ec, ed, ee, ef);
    ea = rel((x7^3 - 8*x7^2 + 5*x7 + 1) / (x7^2 - x7), t7);
    eb = rel((x9^3 - 6*x9^2 + 3*x9 + 1) / (x9^2 - x9), t9);
    ec = rel((t7^2 + 13*t7 + 49) * (t7^2 + 245*t7 + 2401)^3 / t7^7, ellj(t));
    ed = rel((t7^2 + 13*t7 + 49) * (t7^2 + 5*t7 + 1)^3 / t7, ellj(7*t));
    ee = rel((t9 + 9)^3 * (t9^3 + 243*t9^2 + 2187*t9 + 6561)^3 / (t9^9 * (t9^2 + 9*t9 + 27)), ellj(t));
    ef = rel((t9 + 3)^3 * (t9^3 + 9*t9^2 + 27*t9 + 3)^3 / (t9 * (t9^2 + 9*t9 + 27)), ellj(9*t));
    print("  tau = ", t, ": t7 ", ea, "  t9 ", eb, "  j(tau)[t7] ", ec, "  j(7tau)[t7] ", ed, "  j(tau)[t9] ", ee, "  j(9tau)[t9] ", ef);
    if (vecmax([ea, eb, ec, ed, ee, ef]) > TOL, bad++));
  print("== (3) Gamma_1(N)-invariance (gamma with a = d = 1 mod N) ==");
  my(t0 = 0.1 + 0.3*I);
  foreach ([[7, [1,0;7,1]], [7, [8,3;21,8]], [9, [1,0;9,1]], [9, [10,1;9,1]], [8, [1,0;8,1]], [8, [9,1;8,1]], [5, [6,1;5,1]]], gg,
    my(N = gg[1], g = gg[2], F = if (N == 7, M7, if (N == 9, M9, if (N == 8, M8, M5))), e1);
    e1 = rel(F(act(g, t0)), F(t0));
    print("  N = ", N, ", gamma = ", g, ": |M(gamma tau) - M(tau)|/|M| = ", e1);
    if (e1 > TOL, bad++));
  print("== (4) direction of the diamond action (standard complex convention: <d> <-> lower-right entry d) ==");
  \\ mu(x) = (x-1)/x is the note's '<2>' on M7 and M9; mu^2(x) = 1/(1-x)
  my(mu = (x -> (x - 1) / x), mu2 = (x -> 1 / (1 - x)));
  foreach ([[7, [4,1;7,2]], [7, [2,1;7,4]], [9, [5,1;9,2]], [9, [2,1;9,5]]], gg,
    my(N = gg[1], g = gg[2], F = if (N == 7, M7, M9), y = F(act(g, t0)), x = F(t0));
    print("  N = ", N, ", gamma = ", g, " (lower-right d = ", g[2,2] % N, ", upper-left a = ", g[1,1] % N, "):  |M(gamma tau) - mu(M)| = ",
      rel(y, mu(x)), ",  |M(gamma tau) - mu^2(M)| = ", rel(y, mu2(x))));
  my(g8 = [3,1;8,3], y8 = 1 / M8(act(g8, t0)), x8 = 1 / M8(t0));
  print("  N = 8, gamma = ", g8, ": |x(gamma tau) - 1/x| = ", rel(y8, 1 / x8), " (x = 1/M8; <3> is an involution, no direction issue)");
  if (rel(y8, 1 / x8) > TOL, bad++);
  print("== (5) factorizations, cusp fields ==");
  print("  A5 = (x^4+228x^3+494x^2-228x+1)^3: ", A5 == ('x^4 + 228*'x^3 + 494*'x^2 - 228*'x + 1)^3, ";  B5 = x (x^2-11x-1)^5: ", B5 == 'x * ('x^2 - 11*'x - 1)^5);
  print("  factor(A7): degrees/multiplicities ", vector(#factor(A7)[,1], i, [poldegree(factor(A7)[i,1]), factor(A7)[i,2]]));
  print("  factor(A8): ", vector(#factor(A8)[,1], i, [poldegree(factor(A8)[i,1]), factor(A8)[i,2]]), ";  A8 palindromic: ", Polrev(A8) == A8,
        ";  octic palindromic: ", my(o = factor(A8)[1,1]); Polrev(o) == o);
  print("  factor(A9): ", vector(#factor(A9)[,1], i, [poldegree(factor(A9)[i,1]), factor(A9)[i,2]]));
  print("  factor(B7) = ", factor(B7), ";  B8 = ", factor(B8), ";  B9 = ", factor(B9));
  print("  nfdisc: x^3-8x^2+5x+1: ", nfdisc('x^3 - 8*'x^2 + 5*'x + 1), " (poldisc ", poldisc('x^3 - 8*'x^2 + 5*'x + 1), ");  x^3-6x^2+3x+1: ",
        nfdisc('x^3 - 6*'x^2 + 3*'x + 1), " (poldisc ", poldisc('x^3 - 6*'x^2 + 3*'x + 1), ");  x^2-6x+1: ", nfdisc('x^2 - 6*'x + 1), " (poldisc ", poldisc('x^2 - 6*'x + 1), ")");
  print("  cyclic cubic fields are Q(zeta_7)^+ and Q(zeta_9)^+: ", nfisisom('x^3 - 8*'x^2 + 5*'x + 1, polsubcyclo(7, 3)) != 0, " ", nfisisom('x^3 - 6*'x^2 + 3*'x + 1, polsubcyclo(9, 3)) != 0);
  print("\nfailures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
