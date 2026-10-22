\\ Miscellaneous spot checks for the audit (audit_ssj).
\\ (M1) ledger C8: X_0(27) and the isogeny class 27a - j-invariants and CM discriminants.
\\ (M2) section 6.3 last paragraph: "At the identity levels the split is F_{p^2}-heavy for every p":
\\      fraction of F_p-rational supersingular points on X_0(7), X_0(8), X_0(9), X_1(7), X_1(8), X_1(9),
\\      using the formulas confirmed by audit_fp_points.gp / audit_x0_7_9.gp and own class numbers.
\\ (M3) X_0(14) and X_0(27): genus 1 (ledger C8), nu2 = nu3 = 0 (section 6.1).
\\ (M4) h(-36), h(-56) (ledger C9); Ogg list for p = 2, 3 is hard-coded in supersingular_census.gp.
\\ Run: gp -q -f audit_misc.gp > audit_misc.out
default(parisize, 64000000);
hred(D) = { my(h = 0, c); for (a = 1, sqrtint(-D \ 3), for (b = -a + 1, a, if ((b^2 - D) % (4 * a), next); c = (b^2 - D) / (4 * a);
  if (c < a || (c == a && b < 0) || gcd([a, b, c]) != 1, next); h++)); h };
{
  print("(M1) curves of conductor 27 (no elldata installed: explicit Weierstrass models; isogeny class of [0,0,1,0,-7]):");
  my(E0 = ellinit([0, 0, 1, 0, -7]), cls = ellisomat(E0), crv = cls[1]);
  for (i = 1, #crv, my(E = ellinit(crv[i][1]));
    print("   ", ellminimalmodel(E)[1..5], ":  conductor ", ellglobalred(E)[1], ",  j = ", E.j,
          ",  modular degree ", ellmoddegree(ellminimalmodel(E)), if (E.j == 0, ",  CM by disc -3", if (E.j == -12288000, ",  CM by disc -27", ""))));
  print("   the curve of modular degree 1 is X_0(27) itself (g(X_0(27)) = ", mfdim([27, 2], 1), ")");
  my(E14 = ellinit([1, 0, 1, 4, -6]));
  print("   [1,0,1,4,-6]: conductor ", ellglobalred(E14)[1], ", j = ", E14.j, ", modular degree ", ellmoddegree(E14), ";  g(X_0(14)) = ", mfdim([14, 2], 1));
  print("\n(M2) fraction of F_p-rational supersingular points, primes 5 <= p < 2000 where it is >= 1/2:");
  forprime (p = 5, 1999,
    my(h4 = hred(-4 * p), h1 = if (p % 4 == 3, hred(-p), 0), np = (h4 + h1) / 2, L = List());
    if (p != 7, listput(L, ["X_0(7)", (1 + kronecker(-p, 7)) * np, (2 * (p - 1) + 4 * (p % 3 == 2)) / 3]));
    listput(L, ["X_0(8)", 4 * (p % 8 == 7) * h1, p - 1]);
    listput(L, ["X_0(9)", (1 + kronecker(-p, 3)) * np, p - 1]);
    if (p != 7, listput(L, ["X_1(7)", 6 * (p % 7 == 6) * np, 2 * (p - 1)]));
    listput(L, ["X_1(8)", (2 * (p % 16 == 7) + 6 * (p % 16 == 15)) * h1, 2 * (p - 1)]);
    listput(L, ["X_1(9)", 6 * (p % 9 == 8) * np, 3 * (p - 1)]);
    for (i = 1, #L, if (2 * L[i][2] >= L[i][3], print("   p = ", p, ": ", L[i][1], "  #ss(F_p) = ", L[i][2], " of #ss = ", L[i][3],
       if (2 * L[i][2] > L[i][3], "   (F_p-heavy)", "   (exactly half)")))));
  print("\n(M3) nu2, nu3 of 14 and 27: ", [prod(i = 1, 2, 1 + kronecker(-4, [2, 7][i])), prod(i = 1, 2, 1 + kronecker(-3, [2, 7][i]))],
        " and [0 (4 does not divide 27 but 1 + (-4/3) = 0), 0 (9 | 27)]: ", [1 + kronecker(-4, 3), 0]);
  print("\n(M4) h(-36) = ", hred(-36), " (qfbclassno ", qfbclassno(-36), "), h(-56) = ", hred(-56), " (qfbclassno ", qfbclassno(-56), ")");
}
quit;
