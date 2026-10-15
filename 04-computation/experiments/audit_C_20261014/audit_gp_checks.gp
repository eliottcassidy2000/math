\\ Audit C: PARI cross-checks (own script).
\\ (1) S_3(Gamma_0(32), (-4/.)): dimensions, newform orbits, CM flag, a_p for primes < 30000 (written for comparison
\\     with the independent eta/theta construction in audit_newform.py).
\\ (2) Elliptic curves: y^2 = x^3 - x, E_mid: y^2 = x^3 - 11x + 14 (THM-4650): conductor, isogeny class,
\\     a_p comparison with eta(4z)^2 eta(8z)^2.
default(parisizemax, 800000000);
mfnew = mfinit([32, 3, -4], 0);
mfcusp = mfinit([32, 3, -4], 1);
print("dim S_3^new(32,(-4/.)) = ", mfdim(mfnew), ";  dim S_3(32,(-4/.)) = ", mfdim(mfcusp));
L = mfeigenbasis(mfnew);
print("number of newform Galois orbits: ", #L, ";  fields: ", mffields(mfnew), ";  CM: ", mfisCM(L[1]));
co = mfcoefs(L[1], 30000);
P = primes([3, 30000]);
fn = "audit_gp_bp.txt";
{
  for (j = 1, #P,
    my(p = P[j], c = co[p + 1]);
    c = if (type(c) == "t_POLMOD", lift(c), c);
    write(fn, p, " ", polcoef(c, 0), " ", polcoef(c, 1)));
}
print("wrote a_p for ", #P, " primes 3 <= p <= 30000 as (constant, coefficient of y) with y^2 + 1 = 0");
E1 = ellinit([0, 0, 0, -1, 0]);
print("y^2 = x^3 - x: conductor ", ellglobalred(E1)[1], ", torsion ", elltors(E1)[2]);
Em = ellinit([0, 0, 0, -11, 14]);
print("E_mid: y^2 = x^3 - 11x + 14: conductor ", ellglobalred(Em)[1], ", j = ", Em.j, ", torsion ", elltors(Em)[2]);
{
  my(ok = 1);
  forprime (p = 3, 5000, if (ellap(Em, p) != ellap(E1, p), ok = 0));
  print("E_mid and y^2 = x^3 - x have equal a_p for all 3 <= p < 5000: ", ok);
}
iso = ellisomat(E1);
print("isogeny class of y^2 = x^3 - x, [a4, a6] of its members: ", vector(#iso[1], k, iso[1][k][1]));
f1 = mffrometaquo([4, 2; 8, 2]);
cf1 = mfcoefs(f1, 2000);
{
  my(ok = 1);
  forprime (p = 3, 2000, if (cf1[p + 1] != ellap(E1, p), ok = 0));
  print("eta(4z)^2 eta(8z)^2 has a_p = ellap(y^2 = x^3 - x) for all odd p < 2000: ", ok);
}
quit
