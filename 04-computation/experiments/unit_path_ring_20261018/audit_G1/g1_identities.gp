\\ audit G1: exact identities of THM-4700 (generators, walks, recipe, reflection links, mixed quadrilaterals) in PARI.
\\ run: gp -q -s 800000000 g1_identities.gp
read("g1_common.gp");
ok = 1;
chk(name, c) = ok = ok && c; print(if (c, "  ok   ", "  FAIL "), name);
print("(A) generators of O_F in Gamma");
chk("zeta15 unit", isU(zF));
v2 = (1 + s15)/4; chk("1/2 = v + conj v, v = (1 + sqrt-15)/4 unit", isU(v2) && v2 + cj(v2) == 1/2);
beta = vector(4, r, (z5^r - 3)/(z5^(5 - r) - 3));
chk("beta_r units and 11 sum beta_r - 38 = 1/11", prod(r = 1, 4, isU(beta[r])) && 11*vecsum(beta) - 38 == 1/11);
chk("rho = (5 + j)/6 unit and j = 6 rho - 5", isU(rho) && 6*rho - 5 == jF);
chk("u unit and w = (19 (u + conj u) + 6)/11", isU(uu) && (19*(uu + cj(uu)) + 6)/11 == wF);
print("(B) walks of statement 2");
walk2(n, y, sq) = my(vv = (1 + y*sq)/(2*n)); isU(vv) && vv + cj(vv) == 1/n;
wdelta(g, k, c, n) = my(D = cj(g)/g); isU(D) && k*(D + 1/D) + c == 1/n;
chk("1/2 two steps (D = 15, y = 1)", walk2(2, 1, s15));
chk("1/5 two steps (D = 11, y = 3)", walk2(5, 3, jF));
chk("1/13 two steps (D = 3, y = 15)", walk2(13, 15, s3));
chk("1/3 = delta + 1/delta + 2, g = (1 + j)/2 [4 steps]", wdelta((1 + jF)/2, 1, 2, 3));
chk("1/4 = delta + 1/delta + 2, g = (1 + sqrt-15)/2 [4 steps]", wdelta((1 + s15)/2, 1, 2, 4));
chk("1/7 = delta + 1/delta + 2, g = (1 + 3 sqrt-3)/2 [4 steps]", wdelta((1 + 3*s3)/2, 1, 2, 7));
chk("1/19 = delta + 1/delta + 2, g = (1 + 5 sqrt-3)/2 [4 steps]", wdelta((1 + 5*s3)/2, 1, 2, 19));
chk("1/23 = 2 delta + 2/delta - 3, g = (9 + j)/2 [7 steps]", wdelta((9 + jF)/2, 2, -3, 23));
chk("1/17 = 18 delta + 18/delta - 28, g = (11 + sqrt-15)/2 [64 steps]", wdelta((11 + s15)/2, 18, -28, 17));
print("(C) the cubic recipe");
chk("u^3 = m2 u^2 + m1 u + 1", uu^3 == m2*uu^2 + m1*uu + 1);
chk("Y unit and u^2 = a + b u + c Y", isU(YY) && uu^2 == aa + bb*uu + cc*YY);
print("  |M| = ", #Mset, " |C| = ", #Cset, " |A| = ", #Aset);
fac = [[aa, Aset], [bb, Cset], [cc, Cset]];
inS(p) = (p % 15 == 14) && (kronecker(p, 11) == -1);
{
for (i = 1, 3,
  my(k = fac[i][1], S = fac[i][2], pts, dens = 1, pdens = 1);
  pts = uniqv(concat(vector(#S, s, vector(#Mset, m, S[s] + k*Mset[m]))));
  for (t = 1, #pts,
    dens = lcm(dens, denominator(content(nfalgtobasis(nf, lift(pts[t])))));
    pdens = lcm(pdens, denominator(content(lift(pts[t])))));
  print("  factor ", ["A + aM", "C + bM", "C + cM"][i], ": ", #pts, " points; integral-basis denominator lcm ", factor(dens),
        " ; primes of S among them: ", select(p -> inS(p), factor(dens)[,1]~)));
}
print("  b, c: integral-basis denominators ", factor(denominator(content(nfalgtobasis(nf, lift(bb))))), " , ", factor(denominator(content(nfalgtobasis(nf, lift(cc))))));
print("(D) reflection links for k = b, c, every unit difference s of C, every unit q in M");
UC = List();
for (s = 1, #Cset, for (t = 1, #Cset, if (s != t && isU(Cset[t] - Cset[s]), listput(UC, Cset[t] - Cset[s]))));
UC = Vec(UC);
print("  unit differences in C (ordered pairs): ", #UC);
qs = select(m -> m != 0 && isU(m), Mset); print("  unit points q of M: ", #qs);
refl(k, s, q) = my(kq = k*q); (cj(k)/k)*cj(q)*(kq - s)/(cj(kq) - cj(s));
{
for (i = 2, 3,
  my(k = fac[i][1], good = 0, tot = 0, same = 0, qp);
  for (a = 1, #UC, for (b = 1, #qs,
     qp = refl(k, UC[a], qs[b]); tot++;
     if (isU(qp) && isU(UC[a] + k*(qp - qs[b])), good++);
     if (qp == qs[b], same++)));
  chk(Str("k = ", ["", "b", "c"][i], ": q' unit and |s + k(q' - q)| = 1 for all ", tot, " pairs (s, q); q' = q in ", same, " cases"), good == tot));
}
print("  N(b) is a square in F: ", sqrtF(Nm(bb)) != 0, " ; N(c) is a square in F: ", sqrtF(Nm(cc)) != 0, "  (if not, q' = q can never happen)");
print("(E) mixed quadrilaterals (THM-4700(5e)) for k = a, b, c");
tsplit(xx, s) = my(r = s3*s/Nm(xx)); [xx/2*(1 + r), xx/2*(1 - r)];
mixq(k) = {
  my(N = Nm(k), A = (k + cj(k))/2, B = (k - cj(k))/(2*s3), t, sg, z, kz, mu, nu, smu, snu, U, E, verts);
  if (B == 0, return(0));
  t = (1 + A)/(3*B); sg = (A + N)/(3*B);
  z = 2/(cj(k)*(1 + t*s3)); kz = k*z; mu = Nm(kz); nu = Nm(z);
  smu = 4*t/(1 + 3*t^2); snu = 4*sg/(N*(1 + 3*t^2));
  U = tsplit(z, snu); E = tsplit(-kz, smu);
  verts = [0, lift(U[1]), lift(z), lift(z + E[1]/k)];
  cj(A) == A && cj(B) == B && A + B*s3 == k && sg^2 - N*t^2 == (N - 1)/3 && smu^2 == mu*(4 - mu)/3 && snu^2 == nu*(4 - nu)/3
    && U[1] + U[2] + (E[1] + E[2])/k == 0 && isU(U[1]) && isU(U[2]) && isU(E[1]) && isU(E[2]) && z != 0 && #Set(verts) == 4 && mu == 4/(3*t^2 + 1)
    && cj(sg) == sg && cj(t) == t;
}
chk("k = a: conic point, closed walk u1 + u2 + (e1 + e2)/k = 0, units, 4 distinct vertices", mixq(aa));
chk("k = b: same", mixq(bb));
chk("k = c: same", mixq(cc));
print("  explicit Hilbert-symbol splittings: N((A + sqrtN)/B) = -3 and N(1 + sqrtN) = 1 - N hold identically since N = A^2 + 3B^2:");
chk("  N = A^2 + 3 B^2 for a, b, c", prod(i = 1, 3, my(k = fac[i][1], A = (k + cj(k))/2, B = (k - cj(k))/(2*s3)); Nm(k) == A^2 + 3*B^2));
print(if (ok, "ALL IDENTITY CHECKS PASSED", "SOME CHECK FAILED"));
quit;
