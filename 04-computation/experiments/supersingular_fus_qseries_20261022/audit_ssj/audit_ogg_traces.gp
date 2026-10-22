\\ Ogg's list for 5 <= p < 3000 through Hecke traces (audit_ssj): g(X_0(p)^+) = #{newforms of level p with
\\ a_p = -1} = (g - tr T_p | S_2(Gamma_0(p)))/2 (trace form, Eichler-Selberg), and the census identity
\\ #pairs = (g + 1 - n_p)/2 with n_p = (h(-4p) + [p = 3 mod 4] h(-p))/2 (Delfs-Galbraith) - i.e. the two
\\ sides of supersingular_census.gp's g(X_0(p)^+) check, computed without the supersingular polynomial.
\\ Run: gp -q -f audit_ogg_traces.gp > audit_ogg_traces.out
default(parisize, 256000000);
{
  my(ogg = List(), bad = 0, n = 0);
  forprime (p = 5, 2999,
    my(g = mfdim([p, 2], 0), tr = if (g == 0, 0, mfcoef(mftraceform([p, 2], 0), p)), gplus = (g - tr) / 2,
       np = (qfbclassno(-4 * p) + if (p % 4 == 3, qfbclassno(-p), 0)) / 2);
    n++;
    if (gplus != (g + 1 - np) / 2, bad++; print("  mismatch at p = ", p, ": ", gplus, " vs ", (g + 1 - np) / 2));
    if (gplus == 0, listput(ogg, p)));
  print("primes checked: ", n, "; g(X_0(p)^+) from newform traces = (g + 1 - n_p)/2 failures: ", bad);
  print("primes 5 <= p < 3000 with g(X_0(p)^+) = 0: ", Vec(ogg));
  if (Vec(ogg) != [5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71], bad++);
  if (bad == 0, print("ALL AUDIT CHECKS PASSED (Ogg's list restricted to p >= 5)"));
}
quit;
