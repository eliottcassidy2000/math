\\ Hostile probes for the integer coincidences around the moduli 7, 8, 9, 12, 14, 27, 60
\\ (opus-2026-10-22-ssj). Each probe states a tempting reading and its cheapest refutation.
\\ Run: gp -q -f numerology_checks.gp > numerology_checks.out
FUS(N) = {
  my(D = select(d -> d > 1 && d < N, divisors(N)));
  [#D, #select(d -> isprime(d), D), #select(d -> issquarefree(d), D)];
}
isdisc(D) = (D % 4 == 0 || D % 4 == 1) && !issquare(D);

{
  print("(N1) 'the moduli are |D| for class-number-one discriminants D':");
  foreach ([7, 8, 9, 12, 14, 27, 28], m,
    my(D = -m, D4 = -4 * m);
    print("   m = ", m, ":  D = -m ", if (isdisc(D), Str("is a discriminant, h = ", qfbclassno(D)), "is NOT a discriminant"),
      ";  h(-4m) = ", qfbclassno(D4)));
  print("   -> holds for 7, 8, 27 (and 12, 28), fails for 9 and 14 (the Bailey and Rogers moduli): no map, typed NUMEROLOGY.");

  print("\n(N2) (F, U, S) and the defect F - U - S of the moduli:");
  foreach ([7, 8, 9, 12, 14, 27, 60], m, my(v = FUS(m)); print("   ", m, " = ", factor(m), ":  (F,U,S) = ", v, ",  defect ", v[1] - v[2] - v[3]));
  print("   the asked shapes are realised as p^2 = 9 (Bailey), pq = 14 (Rogers), p^2q = 12 (the period of ss(p));");
  print("   the balanced shapes as p = 7 (Rogers-Selberg), p^3 = 8, 27 (Goellnitz-Gordon/Jackson-Slater, Dyson), p^2qr = 60.");
  print("   hostile: F(12) = phi(12) = 4, but F(N) != phi(N) at the other p^2 q values:");
  foreach ([12, 18, 20, 28, 44, 45, 50], m, print("     N = ", m, ": F = ", FUS(m)[1], ", phi = ", eulerphi(m)));

  print("\n(N3) 'p^2qr = 60 = |A_5| is the Rogers-Ramanujan icosahedral group, so F(60) counts its subgroups':");
  print("   subgroup orders of A_5: [1, 2, 3, 4, 5, 6, 10, 12, 60] (classical); proper nontrivial: 7 orders");
  print("   proper nontrivial divisors of 60: ", select(d -> d > 1 && d < 60, divisors(60)), " (F = 10): 15, 20, 30 are not subgroup orders");
  print("   |PSL_2(F_l)| = l(l^2-1)/gcd(2, l-1) for l = 2, 3, 5, 7, 11: ", vector(5, i, my(l = [2, 3, 5, 7, 11][i]); l * (l^2 - 1) / gcd(2, l - 1)),
    " with defects ", vector(5, i, my(l = [2, 3, 5, 7, 11][i], v = FUS(l * (l^2 - 1) / gcd(2, l - 1))); v[1] - v[2] - v[3]));
  print("   (pq, p^2q, p^2qr at l = 2, 3, 5 only; genuine content: #ss(X(l)) = (p-1)|PSL_2(F_l)|/12, a mass-formula fact)");

  print("\n(N4) 'the supersingular period 12 = 2^2 * 3 has shape p^2 q for a structural reason':");
  print("   the correction ss(p) - (p-1)/12 = (1/4)(1-(-4/p)) + (1/3)(1-(-3/p)) has conductor lcm(4, 3) = 12,");
  print("   and 4 = |disc Z[i]|, 3 = |disc Z[(1+sqrt(-3))/2]| are the two CM discriminants with extra units:");
  print("   disc of Z[i] = ", quaddisc(-1), ", disc of Z[omega] = ", quaddisc(-3), ";  so 12 = 2^2 * 3 is |D_-4| |D_-3|.");
  print("   this explains the factorisation of 12, not any role of F(12), U(12), S(12).");
}
quit;
