\\ audit F3 (independent referee), item (a): conjugation-fixed primes of
\\ F = Q(zeta15)(sqrt-11)(w), w^3+2w^2-4w-4 = 0, by PARI idealprimedec on the ABSOLUTE degree-48 field
\\ and an explicit automorphism c (zeta -> zeta^-1, j -> -j, w -> w); conj-stability tested on ideals.
\\ Usage: gp -q f3_fixed_places.gp
default(parisize, "400M");
g = x^3 + 2*x^2 - 4*x - 4;
\\ tower: K = Q(zeta15, j), F = K(w)
C1 = polcompositum(polcyclo(15), x^2 + 11, 3);   \\ [R1, zeta, j, k1], phi = j + k1*zeta
R1 = C1[1]; zK = C1[2]; jK = C1[3]; k1 = C1[4];
C2 = polcompositum(R1, g, 3);                     \\ [R2, phi, w, k2], theta = w + k2*phi
R2 = C2[1]; phiF = C2[2]; wF = C2[3]; k2 = C2[4];
print("deg R1 = ", poldegree(R1), ", deg R2 = ", poldegree(R2));
zF = Mod(subst(lift(zK), x, lift(phiF)), R2);
jF = Mod(subst(lift(jK), x, lift(phiF)), R2);
wF = Mod(lift(wF), R2);
th = Mod(x, R2);
\\ sanity: generators satisfy their relations and theta = w + k2*(j + k1*zeta)
if (polcyclo(15, zF) != 0 || jF^2 + 11 != 0 || subst(g, x, wF) != 0 || th != wF + k2*(jF + k1*zF), error("tower"));
\\ complex conjugation as a polynomial in theta
cth = wF + k2*(-jF + k1*zF^14);
cpol = lift(cth);
if (subst(R2, x, cth) != 0, error("c is not an automorphism"));
if (subst(cpol, x, cth) != th, error("c is not an involution"));
print("c(theta) is an automorphism of order 2: OK");
\\ maximal order at primes <= B (the discriminant has only 2,3,5,11,37 as ramified primes)
B = 1300;
nf = nfinit([R2, B]);
print("nf degree ", nf.r1 + 2*nf.r2, " signature ", nf.sign);
cv = nfalgtobasis(nf, cpol);
\\ check the embedding claim: c agrees with complex conjugation at the first complex embedding
emb = nf.roots;
ok = 1; for (i = 1, #emb, if (abs(subst(cpol, x, emb[i]) - conj(emb[i])) > 1e-20 * (1 + abs(emb[i])), ok = 0));
print("c = complex conjugation at every complex embedding (", #emb, " roots listed): ", ok);
\\ residue-field action: for a fixed unramified prime P with f even, test c(theta) == theta^(p^(f/2)) mod P
frobtest(P) = {
  my(modP = nfmodprinit(nf, P), p = P.p, f = P.f, t1, t2);
  t1 = nfmodpr(nf, cpol, modP);
  t2 = nfmodpr(nf, x, modP)^(p^(f/2));
  t1 == t2;
}
fixedlist(p) = {
  my(dec = idealprimedec(nf, p), res = List(), P, cP);
  for (i = 1, #dec,
    P = dec[i];
    cP = nfgaloisapply(nf, cv, P);
    if (idealval(nf, idealhnf(nf, cP), P) > 0 && idealhnf(nf, cP) == idealhnf(nf, P),
      listput(res, [P.e, P.f, if (P.e == 1 && P.f % 2 == 0, frobtest(P), -1)])));
  [#dec, Vec(res)];
}
print("p : #primes of F above p ; c-fixed primes as [e, f, c = Frob^(f/2) on residue field?]");
{
forprime (p = 2, B,
  r = fixedlist(p);
  if (#r[2] > 0 || p <= 41,
    printf("%5d (p mod 15 = %2d, (p/11) = %2d, cubic mod p = %s): %3d primes; fixed: %s\n", p, p % 15, kronecker(p, 11), factormod(g, p)[,1]~, r[1], r[2])));
}
print("done");
