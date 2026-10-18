\\ audit G1: the primes above 19 in F' = Q(zeta5, sqrt-11, w) and F = F'(sqrt-3); delta_19 = (1 - 5 sqrt-3)/(1 + 5 sqrt-3).
\\ run: gp -q -s 800000000 g1_nineteen.gp
read("g1_common.gp");
\\ ---- F (degree 48)
cprid(P) = idealhnf(nf, nfgaloisapply(nf, cX, P));
dec = idealprimedec(nf, 19);
d19 = (1 - 5*s3)/(1 + 5*s3);
print("F: primes above 19: ", #dec, " ; residue degrees ", Set(apply(P -> P.f, dec)), " ; c-fixed: ", #select(P -> cprid(P) == idealhnf(nf, P), dec));
print("delta_19 unit vector: ", isU(d19), " ; (1 + 5 sqrt-3)/2 = 3 + 5 zeta3 integral of norm 19: ", (1 + 5*s3)/2 == 3 + 5*z3, " , ", Nm((1 + 5*s3)/2) == 19);
h7 = select(P -> nfeltval(nf, lift(wF - 4), P) > 0, dec);
print("F: h = 7 (w = 4) primes above 19: ", #h7);
{
for (k = 1, #h7, my(P = h7[k], cP = cprid(P), vk, vc, idx = 0);
  vk = nfeltval(nf, lift(d19), P);
  for (l = 1, #h7, if (idealhnf(nf, h7[l]) == cP, idx = l));
  vc = if (idx, nfeltval(nf, lift(d19), h7[idx]), "n/a");
  print("  P", k, ": c(P) = P", idx, " (fixed: ", idx == k, ") ; v_P(delta_19) = ", vk, " ; v_c(P)(delta_19) = ", vc,
        " ; v_P(b) = ", nfeltval(nf, lift(bb), P), " v_P(c) = ", nfeltval(nf, lift(cc), P), " v_P(u) = ", nfeltval(nf, lift(uu), P)));
}
print("F: valuations of delta_19 at all 24 primes above 19: ", vecsort(apply(P -> nfeltval(nf, lift(d19), P), dec)));
print("F: valuations of b at all 24 primes above 19: ", Set(apply(P -> nfeltval(nf, lift(bb), P), dec)), " ; of c: ", Set(apply(P -> nfeltval(nf, lift(cc), P), dec)));
\\ ---- F' (degree 24), built independently
f3 = x^3 + 2*x^2 - 4*x - 4;
r1 = polcompositum(polcyclo(5), x^2 + 11, 1)[1]; K1 = r1[1];
r2 = polcompositum(K1, f3, 1)[1]; T = r2[1];
[Tp, ap] = polredbest(T, 1);
thK = subst(lift(r2[2]), x, ap); wP = subst(lift(r2[3]), x, ap);
zP = subst(lift(r1[2]), x, thK); jP = subst(lift(r1[3]), x, thK);
if (polcyclo(5, zP) != 0 || jP^2 + 11 != 0 || subst(f3, x, wP) != 0, error("F' generators"));
Bz = vector(24, k, my(i = (k-1)%4, e = ((k-1)\4)%2, g = (k-1)\8); lift(zP^i * jP^e * wP^g));
Mz = matrix(24, 24, r, s, polcoef(Bz[s], r-1));
cx = Mz^(-1) * vectorv(24, r, r == 2);
cB = vector(24, k, my(i = (k-1)%4, e = ((k-1)\4)%2, g = (k-1)\8); zP^(5-i) * (-jP)^e * wP^g);
cXp = sum(k = 1, 24, cx[k] * cB[k]);
if (subst(Tp, x, cXp) != 0, error("c on F'"));
Bp = vector(24, k, my(i = (k-1)%4, e = ((k-1)\4)%2, g = (k-1)\8); lift(zP^i * ((1 + jP)/2)^e * [1, wP, wP^2/2][g+1]));
dBp = matdet(matrix(24, 24, r, s, trace(Mod(Bp[r]*Bp[s], Tp))));
dlocp = nfdisc([Tp, [2, 5, 11, 37]]);
print("F': disc of the basis ", factor(dBp), " ; p-maximal valuations at 2, 5, 11, 37: ", vector(4, k, valuation(dlocp, [2,5,11,37][k])));
nfp = nfinit([Tp, Bp]);
decp = idealprimedec(nfp, 19);
fixp = select(P -> idealhnf(nfp, nfgaloisapply(nfp, cXp, P)) == idealhnf(nfp, P), decp);
print("F': primes above 19: ", #decp, " ; c-fixed: ", #fixp, " ; residue degrees of fixed ", Set(apply(P -> P.f, fixp)), " ; w = 4 at each fixed prime: ", prod(k = 1, #fixp, nfeltval(nfp, lift(wP - 4), fixp[k]) > 0));
quit;
