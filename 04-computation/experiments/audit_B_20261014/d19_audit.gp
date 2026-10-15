\\ Audit B (2026-10-14): independent checks of THM-4650.
\\ Jacobians are computed here from the classical invariants of the binary quartic (I, J; E: Y^2 = X^3 - 27 I X - 27 J),
\\ not with ellfromeqn, and the midline is re-derived with the OTHER conic parametrised.
default(parisize, 100000000);
quarticIJ(f, X) = {
  my(a = polcoef(f,4,X), b = polcoef(f,3,X), c = polcoef(f,2,X), d = polcoef(f,1,X), e = polcoef(f,0,X));
  [12*a*e - 3*b*d + c^2, 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3];
}
jac(f, X) = { my(v = quarticIJ(f, X)); ellinit([0,0,0,-27*v[1],-27*v[2]]); }
report(name, E) = {
  my(Em = ellminimalmodel(E), r = ellrank(Em), ar = ellanalyticrank(Em));
  print(name, ": minimal [a1..a6] = ", Em[1..5], ", conductor ", ellglobalred(Em)[1], ", j = ", Em.j,
        ", torsion ", elltors(Em)[2], ", ellrank ", r[1..2], ", analytic rank ", ar[1], " (leading L value ", ar[2], ")");
  Em;
}

\\ ---------------- midline x = 1/2, parametrising v = 2(1-y) (the second conic) by v = (s^2-1)/(2s)
v = (s^2-1)/(2*s);
f1 = simplify(4*s^2*(1 + (2 - v)^2));   \\ 1 + 4y^2 with 2y = 2 - v
print("midline quartic (other conic): ", f1);
Emid2 = report("midline Jacobian (other conic)", jac(f1, s));
Emid = ellinit([0,0,0,-11,14]);
print("  isomorphic to the session's y^2 = x^3 - 11x + 14: ", ellisisom(Emid2, Emid) != 0);
\\ session's quartic, our Jacobian formula
fm = t^4 - 8*t^3 + 18*t^2 + 8*t + 1;
{print("  session quartic check: 4t^2(4(1-y)^2+1) at y=(t^2-1)/(4t) = ", simplify(4*t^2*(4*(1-(t^2-1)/(4*t))^2+1)) == fm,
      "; its Jacobian iso to E_mid: ", ellisisom(ellminimalmodel(jac(fm, t)), Emid) != 0);}
\\ explicit 2-isogenies from the congruent-number curve
Ecn = ellinit([0,0,0,-1,0]);
{for(k = 1, 3, P = [[0,0],[1,0],[-1,0]][k]; iso = ellisogeny(Ecn, P); E2 = ellinit(iso[1]);
  print("  y^2 = x^3 - x, kernel <", P, ">: image ", ellminimalmodel(E2)[1..5], "; iso to E_mid: ", ellisisom(E2, Emid) != 0));}
print("  y^2 = x^3 - 16x iso to y^2 = x^3 - x: ", ellisisom(ellinit([0,0,0,-16,0]), Ecn));
report("congruent-number curve y^2 = x^3 - x", Ecn);
report("E_mid", Emid);
print("  j(E_mid) = 287496 = 66^3: ", Emid.j == 66^3);

\\ ---------------- side line y = 0, parametrising the OTHER conic: 1 - x = (s^2-1)/(2s)
w = (s^2-1)/(2*s);
f2 = simplify(4*s^2*(1 + (1 - w)^2));
print("side quartic (other conic): ", f2);
Es2 = report("side Jacobian (other conic)", jac(f2, s));
Eside = ellinit([0,0,0,-2,1]);
print("  isomorphic to the session's y^2 = x^3 - 2x + 1: ", ellisisom(ellminimalmodel(Es2), ellminimalmodel(Eside)) != 0);
report("E_side", Eside);
iferr(print("  ellidentify(E_side) = ", ellidentify(Eside)[1][1], "; ellidentify(E_mid) = ", ellidentify(Emid)[1][1]), err, print("  (elldata not installed)"));
print("  isogeny class of E_side (curves): ", apply(c -> c[1], ellisomat(Eside)[1]));

\\ ---------------- rectangle family: ratio r = 2a; triangle base r, height 1: 1+u^2 = sq, 1+(r-u)^2 = sq
\\ parametrise u = (s^2-1)/(2s)
fr(r) = simplify(4*s^2*(1 + (r - (s^2-1)/(2*s))^2));
{
  my(IJ = quarticIJ(fr(r), s), jj = 6912*IJ[1]^3/(4*IJ[1]^3 - IJ[2]^2));
  print("j of the ratio-r family from (I,J): ", factor(numerator(jj)), " / ", factor(denominator(jj)));
  print("  equals 16(16a^4+16a^2+1)^3/(a^2(a^2+1)) at a = r/2: ", jj == subst(16*(16*a^4+16*a^2+1)^3/(a^2*(a^2+1)), a, r/2));
  my(P = numerator(16*(16*a^4+16*a^2+1)^3 - 287496*a^2*(a^2+1)));
  print("  rational a with j = 287496 (CM by Z[2i]): factor = ", factor(P));
  my(P2 = numerator(16*(16*a^4+16*a^2+1)^3 - 1728*a^2*(a^2+1)));
  print("  rational a with j = 1728: factor = ", factor(P2));
}
print("rank table (ratio r = n/d, n <= 12, d <= 6, gcd 1): [r, conductor, torsion, ellrank, analytic rank]");
{
  my(R1 = List(), R0 = List(), T8 = List(), bad = List());
  for(n = 1, 12, for(d = 1, 6, if(gcd(n, d) == 1,
    r = n/d; E = ellminimalmodel(jac(fr(r), s));
    rk = ellrank(E); ar = ellanalyticrank(E)[1]; tor = elltors(E)[2];
    if(rk[1] != rk[2] || rk[1] != ar, listput(bad, [r, rk, ar]));
    if(rk[1] == 1, listput(R1, r), listput(R0, r));
    if(tor == [8], listput(T8, r));
    print("  ", r, "  N = ", ellglobalred(E)[1], "  tors ", tor, "  ellrank ", rk[1..2], "  an ", ar))));
  print("rank 1: ", Vec(R1));
  print("rank 0: ", Vec(R0));
  print("full elltors at 3/2, 5/6, 8/3: ", [elltors(ellminimalmodel(jac(fr(q), s)))[1..2] | q <- [3/2, 5/6, 8/3]]);
  print("torsion Z/8: ", Vec(T8));
  print("disagreements (ellrank bounds unequal or != analytic rank): ", Vec(bad));
  print("count: ", #R1 + #R0);
}
quit
