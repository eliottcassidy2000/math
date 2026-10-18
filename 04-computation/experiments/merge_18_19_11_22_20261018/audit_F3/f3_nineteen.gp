\\ audit F3 (independent referee), item (c): the 19-places of K' = Q(zeta5, j), F' = K'(w), and of F = Q(zeta15, j, w);
\\ valuations of b, c (golden_contact_circuit_20261017i.md section 4) and of u, u2; whether the spindle lies in F'.
\\ Usage: gp -q f3_nineteen.gp
default(parisize, "400M");
g = x^3 + 2*x^2 - 4*x - 4;
\\ ---------- K' = Q(zeta5, j)
C1 = polcompositum(polcyclo(5), x^2 + 11, 3);  \\ [R, zeta5, j, k], root = j + k*zeta5
R = C1[1]; z5 = C1[2]; jj = C1[3]; kk = C1[4];
K = nfinit(R);
print("K' = Q(zeta5, sqrt-11): degree ", poldegree(R), ", disc factorization ", factor(K.disc));
cK = lift(-jj + kk*z5^4);                       \\ c(root) = -j + k*zeta5^-1
if (subst(R, x, Mod(cK, R)) != 0 || subst(cK, x, Mod(cK, R)) != Mod(x, R), error("c on K'"));
cKv = nfalgtobasis(K, cK);
\\ coefficients of the u-square circuit
d = 1 + z5; a = 1 + z5^4;
m2 = -(20 + 11*jj)/19; m1 = (20 - 11*jj)/19;
cc = (a*m2 + 1 - a^2*m1 + a^3)/(d*a - 1);
bb = cc*d + a*m1 - a^2;
phi = 1 + z5 + z5^4;
if (d*a - 1 != phi, error("|d|^2-1 != phi"));
dec = idealprimedec(K, 19);
print("primes of K' above 19: ", #dec, " with [e, f] = ", vector(#dec, i, [dec[i].e, dec[i].f]));
{
for (i = 1, #dec,
  P = dec[i];
  cP = nfgaloisapply(K, cKv, P);
  modP = nfmodprinit(K, P);
  fr = (nfmodpr(K, cK, modP) == nfmodpr(K, x, modP)^19);
  printf("  P%d: c-fixed %d, c = Frobenius %d, v(b) = %d, v(c) = %d, v(b/c) = %d\n", i,
    idealhnf(K, cP) == idealhnf(K, P), fr, nfeltval(K, lift(bb), P), nfeltval(K, lift(cc), P), nfeltval(K, lift(bb/cc), P)));
}
print("norm of b: ", factor(norm(bb)), " ; norm of c: ", factor(norm(cc)));
\\ ---------- F' = K'(w)
C2 = polcompositum(R, g, 3);   \\ [R2, rootK, w, k2], theta = w + k2*rootK
R2 = C2[1]; rK = C2[2]; wF = Mod(lift(C2[3]), R2); k2 = C2[4];
toF(e) = Mod(subst(lift(e), x, lift(rK)), R2);
zF = toF(z5); jF = toF(jj);
th = Mod(x, R2);
if (th != wF + k2*toF(Mod(x, R)), error("tower F'"));
cth = wF + k2*(-jF + kk*zF^4);
cpol = lift(cth);
if (subst(R2, x, cth) != 0 || subst(cpol, x, cth) != th, error("c on F'"));
Fp = nfinit([R2, 1000]);
print("F' = Q(zeta5, sqrt-11, w): degree ", poldegree(R2), " signature ", Fp.sign);
cFv = nfalgtobasis(Fp, cpol);
hF = wF + 3;
uF = (hF - jF)/(hF + jF);
u2F = (hF + toF(phi) - jF)/(hF + toF(phi) + jF);
\\ unit vectors: x * c(x) = 1 with c(x) = x(theta -> c(theta))
cofF(e) = Mod(subst(lift(e), x, cth), R2);
print("u * c(u) = 1: ", uF*cofF(uF) == 1, " ; u2 * c(u2) = 1: ", u2F*cofF(u2F) == 1);
bF = toF(bb); cF = toF(cc);
decF = idealprimedec(Fp, 19);
print("primes of F' above 19: ", #decF, ", [e,f] = ", Set(vector(#decF, i, [decF[i].e, decF[i].f])));
nfix = 0;
{
for (i = 1, #decF,
  P = decF[i];
  cP = nfgaloisapply(Fp, cFv, P);
  isfix = (idealhnf(Fp, cP) == idealhnf(Fp, P));
  modP = nfmodprinit(Fp, P);
  wres = nfmodpr(Fp, lift(wF), modP);
  hres = nfmodpr(Fp, lift(hF), modP); jres = nfmodpr(Fp, lift(jF), modP);
  typ = if (hres == jres, "h=j", if (hres == -jres, "h=-j", if (hres == 7, "h=7", "?")));
  fr = (nfmodpr(Fp, cpol, modP) == nfmodpr(Fp, x, modP)^19);
  if (isfix, nfix++);
  printf("  P%2d: type %4s, w mod P = %s, c-fixed %d, c = Frob %d, v(b) %d, v(c) %d, v(u) %d, v(u2) %d\n", i, typ, wres,
    isfix, fr, nfeltval(Fp, lift(bF), P), nfeltval(Fp, lift(cF), P), nfeltval(Fp, lift(uF), P), nfeltval(Fp, lift(u2F), P)));
}
print("c-fixed primes of F' above 19: ", nfix);
\\ the spindle needs sqrt(-3): is x^2 + 3 solvable in F'?  (polcompositum / nfroots)
print("degrees of compositum F'(sqrt-3): ", apply(poldegree, polcompositum(R2, x^2 + 3)), " (48 means sqrt-3 is not in F')");
print("degrees of compositum K'(sqrt-3): ", apply(poldegree, polcompositum(R, x^2 + 3)), " (16 means sqrt-3 is not in K')");
print("done");
