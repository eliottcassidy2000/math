\\ Independent audit of section 2 for larger primes (audit_ssj).
\\ (1) For each listed prime p: count supersingular j in F_{p^2} \ F_p by computing the group order
\\     #E_j(F_{p^2}) with ellcard over F_{p^2} (t = p^2+1-#E == 0 mod p), j running over all of F_{p^2}\F_p;
\\     F_p part by ellap over F_p (only as a cross-check of the Python point count).
\\ (2) g(X_0(p)^+) = #{newforms f of level p with a_p(f) = -1} = (g - tr T_p | S_2(Gamma_0(p)))/2,
\\     with tr T_p from the trace form (mftraceform, Eichler-Selberg), NOT from the Fricke/Atkin-Lehner
\\     fixed-point formula used by supersingular_census.gp; compared with the number of conjugate pairs.
\\ Run: gp -q -f audit_census_ffp2.gp > audit_census_ffp2.out
default(parisize, 256000000);
L = [5, 7, 11, 13, 37, 43, 67, 97, 101, 103, 109, 113, 127, 131, 137, 139, 149, 163, 191, 193, 239, 331, 337, 383];
{
  my(bad = 0);
  print("[p, #ss in F_p, #pairs in F_{p^2}\\F_p, total, floor(p/12)+eps, g(X_0(p)), #{newforms a_p=-1}]");
  foreach (L, p,
    my(w = ffgen([p, 2], 'w), g = ffprimroot(w), z = g^0, n2 = 0, n1 = 0, q = p^2, deg, gg, T, trp, gplus, eps);
    \\ all j in F_{p^2}: enumerate as powers of a primitive root (plus 0)
    for (k = 0, q - 2,
      my(j = z, E, c);
      z *= g;
      if (j^p == j, next);                 \\ skip F_p (handled below)
      E = ellinit(ellfromj(j));
      c = ellcard(E);
      if ((q + 1 - c) % p == 0, n2++));
    for (j = 0, p - 1,
      if (ellap(ellinit(ellfromj(Mod(j, p)))) == 0, n1++));
    if (n2 % 2, bad++; print("odd number of non-F_p ss j at p = ", p));
    eps = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2][p % 12];
    deg = p \ 12 + eps;
    if (n1 + n2 != deg, bad++; print("total mismatch at p = ", p));
    gg = mfdim([p, 2], 1);
    T = mftraceform([p, 2], 0);
    trp = mfcoef(T, p);
    gplus = (mfdim([p, 2], 0) - trp) / 2;       \\ all of S_2(Gamma_0(p)) is new
    if (gplus != n2 / 2, bad++; print("pairs != #{a_p = -1} at p = ", p, ": ", n2 / 2, " vs ", gplus));
    print([p, n1, n2 / 2, n1 + n2, deg, gg, gplus]));
  print("failures: ", bad);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED"));
}
quit;
