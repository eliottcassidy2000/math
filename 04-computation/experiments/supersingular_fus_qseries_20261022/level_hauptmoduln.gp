\\ Product sides of the identities as coordinates on X_1(N), N = 5 (control), 7, 8, 9
\\ (opus-2026-10-22-ssj).
\\ Klein-form numerators f_g = prod_{m>=1} (1 - q^{N(m-1)+g})(1 - q^{Nm-g}) = theta_N(g; q)/(q^N;q^N).
\\ Ratios of the identities' product sides are monomials in the f_g (common denominators cancel):
\\   N = 5  (Rogers-Ramanujan, control):  M5 = q   f_1^5 / f_2^5        = q (H/G)^5
\\   N = 7  (Rogers-Selberg S.31-33):     M7 = q^-1 f_2^2 f_3 / f_1^3   = q^-1 RS_2^2 RS_3 / RS_1^3
\\   N = 8  (Goellnitz-Gordon S.34/36):   M8 = q   (f_1 / f_3)^2         = q (GG_1/GG_3)^2
\\   N = 9  (Bailey S.40c-42c):           M9 = q^-1 f_2 f_4 / f_1^2     = q^-1 B_2 B_4 / B_1^2
\\ Test: find A, B in Q[x] with j = A(M)/B(M), deg A = d_N = [PSL2(Z) : Gamma_1(N)], solved on
\\ O(q^SOLVE) and verified on O(q^VERIFY); also the norm maps to the X_0(N) Hauptmoduln
\\ t_7 = (eta(tau)/eta(7 tau))^4 and t_9 = (eta(tau)/eta(9 tau))^3, and the diamond action on M7.
\\ Run: gp -q -f level_hauptmoduln.gp > level_hauptmoduln.out
default(parisize, 256000000);
SOLVE = 220; VERIFY = 520;
P = VERIFY + 60;

kf(N, g) = prod(m = 1, P \ N + 2, (1 - 'q^(N * (m - 1) + g)) * (1 - 'q^(N * m - g)) + O('q^P));
E4 = 1 + 240 * sum(n = 1, P, sigma(n, 3) * 'q^n) + O('q^(P + 1));
Dl = 'q * eta('q + O('q^(P + 1)))^24;
J = E4^3 / Dl;

\\ find j = A(M)/B(M) with deg A = d, deg B <= d - 1 (M has a simple pole at infinity)
rel(M, d, K) = {
  my(cols = List(), mat, ker, v, A, B, pw = vector(d + 1), lo, nrows);
  pw[1] = 1 + O('q^P);
  for (i = 1, d, pw[i + 1] = pw[i] * M);
  lo = -d - 1;
  nrows = K - lo;
  \\ columns: J*M^i (i = 0..d-1) for B, and -M^i (i = 0..d) for A
  for (i = 0, d - 1, listput(cols, J * pw[i + 1]));
  for (i = 0, d, listput(cols, -pw[i + 1]));
  mat = matrix(nrows, #cols, r, c, polcoeff(cols[c], lo + r - 1, 'q));
  ker = matker(mat);
  if (#ker != 1, return([#ker]));
  v = ker[, 1];
  B = sum(i = 0, d - 1, v[i + 1] * 'x^i);
  A = sum(i = 0, d, v[d + i + 1] * 'x^i);
  my(c = pollead(A)); A /= c; B /= c;
  [A, B];
}

verify(M, A, B, K) = {
  my(r = J * subst(B, 'x, M) - subst(A, 'x, M));
  r = r + O('q^K);
  valuation(r, 'q) >= K;
}

show(name, A, B) = {
  my(fa = factor(A), fb = factor(B));
  print("  ", name, ": deg A = ", poldegree(A), ", deg B = ", poldegree(B));
  print("    A = ", A);
  print("    B = ", fb);
}

{
  my(f, M, r, d);
  print("== Hauptmodul relations j = A(M)/B(M) (solved on O(q^", SOLVE, "), verified on O(q^", VERIFY, ")) ==");
  \\ N = 5 control
  f = vector(2, g, kf(5, g));
  M = 1 / ('q * f[1]^5 / f[2]^5);   \\ 1/M5 has a simple pole
  r = rel(M, 12, SOLVE); d = 12;
  if (#r != 2, print("N=5: no unique relation of degree 12: kernel dim ", r[1]),
    A5 = r[1]; B5 = r[2];
    print("  N = 5 (1/M5): relation found; verified to O(q^", VERIFY, "): ", verify(M, r[1], r[2], VERIFY));
    show("N = 5, x = 1/M5", r[1], r[2]));
  \\ N = 7
  f = vector(3, g, kf(7, g));
  M7 = 'q^-1 * f[2]^2 * f[3] / f[1]^3;
  r = rel(M7, 24, SOLVE);
  if (#r != 2, print("N=7: kernel dim ", r[1]),
    A7 = r[1]; B7 = r[2];
    print("  N = 7 (M7): relation found; verified to O(q^", VERIFY, "): ", verify(M7, A7, B7, VERIFY));
    show("N = 7, x = M7", A7, B7));
  \\ degree-23 attempt must fail (hostile control: M7 is not a coordinate of a smaller curve)
  r = rel(M7, 23, SOLVE);
  print("  hostile: degree 23 for M7 gives kernel dimension ", if (#r == 2, 1, r[1]), " (0 expected)");
  \\ N = 8
  f = vector(4, g, kf(8, g));
  M = 1 / ('q * (f[1] / f[3])^2);
  r = rel(M, 24, SOLVE);
  if (#r != 2, print("N=8 (1/M8, degree 24): kernel dim ", r[1]),
    A8 = r[1]; B8 = r[2];
    print("  N = 8 (1/M8): relation found; verified to O(q^", VERIFY, "): ", verify(M, A8, B8, VERIFY));
    show("N = 8, x = 1/M8", A8, B8));
  for (dd = 6, 23, r = rel(M, dd, SOLVE); if (#r == 2, print("  N = 8: ALSO a relation of degree ", dd)));
  \\ the diamond involution <3> of X_1(8) is x -> 1/x on x = 1/M8 (A8 palindromic, B8 symmetric);
  \\ t8 = x + 1/x is then a Hauptmodul of X_0(8): find j = P8(t8)/Q8(t8) of degree 12
  print("  A8 palindromic: ", Polrev(A8) == A8, ";  j-map invariant under x -> 1/x: ",
    subst(A8, 'x, 1/'x) / subst(B8, 'x, 1/'x) == A8 / B8);
  T8 = M + 1 / M;
  r = rel(T8, 12, SOLVE);
  if (#r != 2, print("N=8 X_0: kernel dim ", r[1]),
    P8 = r[1]; Q8 = r[2];
    print("  X_0(8): t8 = x + 1/x, j = P8(t8)/Q8(t8) verified to O(q^", VERIFY, "): ", verify(T8, P8, Q8, VERIFY));
    show("N = 8, X_0(8) coordinate t8", P8, Q8));
  \\ N = 9
  f = vector(4, g, kf(9, g));
  M9 = 'q^-1 * f[2] * f[4] / f[1]^2;
  r = rel(M9, 36, SOLVE);
  if (#r != 2, print("N=9: kernel dim ", r[1]),
    A9 = r[1]; B9 = r[2];
    print("  N = 9 (M9): relation found; verified to O(q^", VERIFY, "): ", verify(M9, A9, B9, VERIFY));
    show("N = 9, x = M9", A9, B9));
  r = rel(M9, 35, SOLVE);
  print("  hostile: degree 35 for M9 gives kernel dimension ", if (#r == 2, 1, r[1]), " (0 expected)");
}

\\ norm maps to X_0(N): t_N as rational functions of degree phi(N)/2 in M_N
norm_rel(T, M, d, K) = {
  \\ T = C(M)/D(M), deg C = d (T has a simple pole where M does), deg D <= d - 1
  my(cols = List(), mat, ker, v, C, D, pw = vector(d + 1), lo = -d - 1, nrows);
  pw[1] = 1 + O('q^P);
  for (i = 1, d, pw[i + 1] = pw[i] * M);
  nrows = K - lo;
  for (i = 0, d - 1, listput(cols, T * pw[i + 1]));
  for (i = 0, d, listput(cols, -pw[i + 1]));
  mat = matrix(nrows, #cols, r, c, polcoeff(cols[c], lo + r - 1, 'q));
  ker = matker(mat);
  if (#ker != 1, return([#ker]));
  v = ker[, 1];
  D = sum(i = 0, d - 1, v[i + 1] * 'x^i);
  C = sum(i = 0, d, v[d + i + 1] * 'x^i);
  my(c = pollead(C)); C /= c; D /= c;
  [C, D];
}

{
  my(e1, e7, e9, t7, t9, r, f, Mu, D2, num, den, mu, ker, cols, mat);
  print("\n== Norm maps X_1(N) -> X_0(N) ==");
  e1 = eta('q + O('q^P)); e7 = eta('q^7 + O('q^P)); e9 = eta('q^9 + O('q^P));
  t7 = 'q^-1 * e1^4 / e7^4;
  r = norm_rel(t7, M7, 3, SOLVE);
  if (#r == 2, print("  t7 = (eta/eta_7)^4 = C(M7)/D(M7): C = ", r[1], ",  D = ", r[2],
       "; verified: ", valuation(t7 * subst(r[2], 'x, M7) - subst(r[1], 'x, M7) + O('q^VERIFY), 'q) >= VERIFY),
    print("  t7: kernel dim ", r[1]));
  \\ check t7 against the classical level-7 j-maps (Klein/Fricke): with t = t7,
  \\   j(tau) = (t^2+13t+49)(t^2+245t+2401)^3/t^7   and   j(7 tau) = (t^2+13t+49)(t^2+5t+1)^3/t
  my(J7 = subst(J, 'q, 'q^7), J9 = subst(J, 'q, 'q^9), V7 = 7 * (VERIFY \ 7), V9 = 9 * (VERIFY \ 9));
  print("  classical j(tau)   = (t^2+13t+49)(t^2+245t+2401)^3/t^7 at t = t7 holds: ",
    valuation(J * t7^7 - (t7^2 + 13 * t7 + 49) * (t7^2 + 245 * t7 + 2401)^3 + O('q^(VERIFY - 7)), 'q) >= VERIFY - 7);
  print("  classical j(7 tau) = (t^2+13t+49)(t^2+5t+1)^3/t at t = t7 holds: ",
    valuation(J7 * t7 - (t7^2 + 13 * t7 + 49) * (t7^2 + 5 * t7 + 1)^3 + O('q^(V7 - 14)), 'q) >= V7 - 14);
  print("  (my first recollection j(tau) = (t^2+13t+49)(t^2+5t+1)^3/t is the j(7 tau) formula: ",
    valuation(J * t7 - (t7^2 + 13 * t7 + 49) * (t7^2 + 5 * t7 + 1)^3 + O('q^VERIFY), 'q) >= VERIFY, ")");
  t9 = 'q^-1 * e1^3 / e9^3;
  r = norm_rel(t9, M9, 3, SOLVE);
  if (#r == 2, print("  t9 = (eta/eta_9)^3 = C(M9)/D(M9): C = ", r[1], ",  D = ", r[2],
       "; verified: ", valuation(t9 * subst(r[2], 'x, M9) - subst(r[1], 'x, M9) + O('q^VERIFY), 'q) >= VERIFY),
    print("  t9: kernel dim ", r[1]));
  \\ level 9 with t = t9:  j(tau) = (t+9)^3 (t^3+243t^2+2187t+6561)^3 / (t^9 (t^2+9t+27)),
  \\                       j(9 tau) = (t+3)^3 (t^3+9t^2+27t+3)^3 / (t (t^2+9t+27))
  print("  classical j(tau)   = (t+9)^3 (t^3+243t^2+2187t+6561)^3/(t^9 (t^2+9t+27)) at t = t9 holds: ",
    valuation(J * t9^9 * (t9^2 + 9 * t9 + 27) - (t9 + 9)^3 * (t9^3 + 243 * t9^2 + 2187 * t9 + 6561)^3 + O('q^(VERIFY - 9)), 'q) >= VERIFY - 9);
  print("  classical j(9 tau) = (t+3)^3 (t^3+9t^2+27t+3)^3/(t(t^2+9t+27)) at t = t9 holds: ",
    valuation(J9 * t9 * (t9^2 + 9 * t9 + 27) - (t9 + 3)^3 * (t9^3 + 9 * t9^2 + 27 * t9 + 3)^3 + O('q^(V9 - 27)), 'q) >= V9 - 27);
  \\ diamond operator on M7: permuting the Klein forms 1 -> 2 -> 3 -> 1 sends M7 to D2 = f_1 f_3^2 / f_2^3
  f = vector(3, g, kf(7, g));
  D2 = f[1] * f[3]^2 / f[2]^3;
  \\ solve D2 * (c M7 + d) = a M7 + b
  cols = [D2 * M7, D2, -M7, -1 + O('q^P)];
  mat = matrix(SOLVE + 2, 4, r, c, polcoeff(cols[c], r - 2, 'q));
  ker = matker(mat);
  print("\n== Diamond action on the Rogers-Selberg coordinate ==");
  if (#ker != 1, print("  no Moebius relation: kernel dim ", #ker),
    my(v = ker[, 1], c0, d0, a0, b0);
    c0 = v[1]; d0 = v[2]; a0 = v[3]; b0 = v[4];
    mu = [a0, b0; c0, d0]; mu /= content(mu);
    print("  f_1 f_3^2/f_2^3 = (a M7 + b)/(c M7 + d) with [a, b; c, d] = ", mu);
    print("  verified to O(q^", VERIFY, "): ", valuation(D2 * (mu[2, 1] * M7 + mu[2, 2]) - (mu[1, 1] * M7 + mu[1, 2]) + O('q^VERIFY), 'q) >= VERIFY);
    print("  mu^3 is scalar (order 3 in PGL2): ", mu^3 == (mu^3)[1, 1] * matid(2));
    MU = mu);
  print("  (M7 has a simple pole at the cusp infinity; the permuted monomial is a unit there)");
  \\ N = 9: <2> permutes f_1 -> f_2 -> f_4 -> f_8 = f_1, sending M9 = q^-1 f_2 f_4/f_1^2 to f_1 f_4/f_2^2
  my(f9 = vector(4, g, kf(9, g)), D9 = f9[1] * f9[4] / f9[2]^2, cols9, mat9, ker9);
  cols9 = [D9 * M9, D9, -M9, -1 + O('q^P)];
  mat9 = matrix(SOLVE + 2, 4, r, c, polcoeff(cols9[c], r - 2, 'q));
  ker9 = matker(mat9);
  if (#ker9 != 1, print("  N = 9: no Moebius relation: kernel dim ", #ker9),
    my(v = ker9[, 1], mu9 = [v[3], v[4]; v[1], v[2]]);
    mu9 /= content(mu9);
    print("  N = 9: f_1 f_4/f_2^2 = (a M9 + b)/(c M9 + d) with [a, b; c, d] = ", mu9, "; verified to O(q^", VERIFY, "): ",
      valuation(D9 * (mu9[2, 1] * M9 + mu9[2, 2]) - (mu9[1, 1] * M9 + mu9[1, 2]) + O('q^VERIFY), 'q) >= VERIFY,
      "; order 3: ", mu9^3 == (mu9^3)[1, 1] * matid(2));
    MU9 = mu9);
  \\ N = 8: <3> swaps f_1 and f_3, sending x = 1/M8 = q^-1 (f_3/f_1)^2 to M8 = 1/x (A8 palindromic above)
  print("  N = 8: <3> acts as x -> 1/x on x = 1/M8 (equivalent to the palindromic j-map checked above)");
}

\\ save relations for the supersingular script (write appends, so remove any previous copy first)
system("rm -f level_relations.gp");
write("level_relations.gp", "A5=", A5, ";B5=", B5, ";P8=", P8, ";Q8=", Q8, ";A7=", A7, ";B7=", B7, ";A8=", A8, ";B8=", B8, ";A9=", A9, ";B9=", B9, ";MU=", MU, ";MU9=", MU9, ";");
print("\nrelations written to level_relations.gp");
quit;
