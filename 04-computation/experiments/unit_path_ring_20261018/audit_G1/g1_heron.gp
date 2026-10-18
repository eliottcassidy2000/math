\\ audit G1: Heron census of THM-4700(4) and the house numbers of 5(d)/(e), in PARI.
\\ run: gp -q -s 800000000 g1_heron.gp
read("g1_common.gp");
default(realprecision, 60);
emb = polroots(T2);
R = [one, 1/Nm(aa), 1/Nm(bb), 1/Nm(cc)]; nm = ["1", "d_a", "d_b", "d_c"];
\\ degree-one primes of F above p = 1 mod 15 (p >= 1000), with their residue maps
D1 = List();
{
forprime(p = 1000, 4000, if (p % 15 == 1,
  my(dec = idealprimedec(nf, p));
  for (k = 1, #dec, if (dec[k].f == 1, listput(D1, [p, nfmodprinit(nf, dec[k]), dec[k]]))));
  if (#D1 >= 120, break));
}
print("degree-one primes of F used: ", #D1, " (rational primes ", Set(vector(#D1, k, D1[k][1])), ")");
cert(y) = {
  my(yb = nfalgtobasis(nf, lift(y)));
  for (k = 1, #D1, my(p = D1[k][1], P = D1[k][3]);
    if (nfeltval(nf, yb, P) == 0, my(r = nfmodpr(nf, yb, D1[k][2])); if (!issquare(r), return([p, k]))));
  0;
}
print("side squared lengths (min, max over the 48 embeddings):");
for (i = 1, 4, my(vals = vector(#emb, k, real(ev(R[i], emb[k])))); print("  ", nm[i], ": ", vecmin(vals), " .. ", vecmax(vals)));
nreal = 0; ntri = 0; nexcl = 0; ntpfail = 0;
{
forvec(ix = [[1, 4], [1, 4], [1, 4]],
  my(R1 = R[ix[1]], R2 = R[ix[2]], R3 = R[ix[3]], H, Hv, tp, s, ct, name);
  H = 2*(R1*R2 + R2*R3 + R3*R1) - R1^2 - R2^2 - R3^2;
  Hv = vector(#emb, k, real(ev(H, emb[k])));
  tp = (vecmin(Hv) > 0);
  s = sqrtF(H/3);
  name = Str(nm[ix[1]], ", ", nm[ix[2]], ", ", nm[ix[3]]);
  if (s != 0,
    \\ H/3 = s^2 with s in F; realizable iff s can be taken in F+ (s real), which forces H totally positive
    print("  ", name, ": H/3 IS a square in F; s real: ", cj(s) == s, " ; H totally positive: ", tp, " ; min H = ", vecmin(Hv));
    if (cj(s) == s, ntri++),
    ct = cert(H/3); nexcl++; if (!tp, ntpfail++);
    print("  ", name, ": H/3 NOT a square in F (nfroots); quadratic non-residue certificate at a degree-one prime above p = ", ct[1], " ; H totally positive: ", tp, " ; min H = ", vecmin(Hv)));
  , 1);
}
print("realizable (H in 3 (F+)^*2): ", ntri, " ; excluded: ", nexcl, " ; excluded that also fail total positivity: ", ntpfail);
print("explicit constructions of statement 4:");
z10 = -z5^3;
print("  zeta10 = -zeta5^3: zeta10^2 = zeta5: ", z10^2 == z5, " ; zeta10 primitive 10th root: ", z10^5 == -1);
print("  N(1 - zeta10) = 1/N(a): ", Nm(1 - z10) == 1/Nm(aa), " ; (2 + z5 + z5^4)(2 + z5^2 + z5^3) = 1: ", (2 + z5 + z5^4)*(2 + z5^2 + z5^3) == 1);
print("  (1, d_a, d_a): |1/a| = |-zeta5/a| = d_a and |1/a + zeta5/a| = 1: ", Nm(1/aa) == R[2] && Nm(-z5/aa) == R[2] && Nm(1/aa + z5/aa) == 1);
print("  a-matching: units r with N(r - 1) = 1/N(a) are exactly zeta10^(+-1): roots of T^2 - (2 - 1/N(a)) T + 1 = ", [z10 + 1/z10 == 2 - 1/Nm(aa)]);
print("house numbers:");
print("  max |sigma(b)| = ", house(bb), " ; max |sigma(c)| = ", house(cc), " ; max |sigma(a)| = ", house(aa));
print("  min |sigma(b)| = ", vecmin(vector(#emb, k, abs(ev(bb, emb[k])))), " ; min |sigma(c)| = ", vecmin(vector(#emb, k, abs(ev(cc, emb[k])))));
quit;
