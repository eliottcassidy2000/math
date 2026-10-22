\\ Delfs-Galbraith input (a) of section 6.3, checked from curves (audit_ssj): for 5 <= p < 1000, enumerate the
\\ F_p-isomorphism classes of supersingular curves over F_p as (j, twist): j != 0, 1728 has the two quadratic
\\ twists (same End_{F_p}); j = 1728 (p = 3 mod 4) has y^2 = x^3 + x and y^2 = x^3 - x; j = 0 (p = 2 mod 3) has
\\ y^2 = x^3 + 1 and its quadratic twist (F_p^*/F_p^*6 has order gcd(6, p-1) = 2, F_p^*/F_p^*4 has order 2).
\\ End_{F_p} = surface Z[(1+sqrt(-p))/2] iff E[2] is F_p-rational (p = 3 mod 4), else floor Z[sqrt(-p)].
\\ Claim checked: #floor classes = h(-4p), #surface classes = [p = 3 mod 4] h(-p) (class numbers by reduced forms).
\\ Run: gp -q -f audit_dg_classes.gp > audit_dg_classes.out
hred(D) = { my(h = 0, c); for (a = 1, sqrtint(-D \ 3), for (b = -a + 1, a, if ((b^2 - D) % (4 * a), next); c = (b^2 - D) / (4 * a);
  if (c < a || (c == a && b < 0) || gcd([a, b, c]) != 1, next); h++)); h };
rat2(A, B, p) = #polrootsmod('x^3 + A * 'x + B, p) == 3;
{
  my(bad = 0, n = 0);
  forprime (p = 5, 999,
    my(fl = 0, su = 0);
    if (gcd(6, p - 1) != 2 && p % 3 == 2, bad++);
    for (j = 0, p - 1,
      my(E = ellinit(ellfromj(Mod(j, p))));
      if (ellap(E) != 0, next);
      if (j == 1728 % p,
        \\ two quartic twists y^2 = x^3 + x, y^2 = x^3 - x (gcd(4, p-1) = 2 for p = 3 mod 4)
        foreach ([1, -1], a, if (ellap(ellinit([a, 0], p)) != 0, bad++); if (p % 4 == 3 && rat2(a, 0, p), su++, fl++)),
        my(A, B);
        if (j == 0, A = 0; B = 1, my(k = Mod(j, p) / (1728 - Mod(j, p))); A = lift(3 * k); B = lift(2 * k));
        if (p % 4 == 3 && rat2(A, B, p), su += 2, fl += 2)));
    n++;
    if (fl != hred(-4 * p) || su != if (p % 4 == 3, hred(-p), 0), bad++; print("  mismatch p = ", p, ": floor ", fl, " vs ", hred(-4*p), ", surface ", su)));
  print("primes 5 <= p < 1000 checked: ", n, "; #floor classes = h(-4p), #surface classes = [p=3 mod 4] h(-p): failures ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
