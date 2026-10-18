\\ audit G1: census of conjugation-fixed primes of F = Q(zeta15, sqrt-11, w) with PARI idealprimedec in the degree-48 field.
default(parisize, 800000000);
V = read("g1_F.bin"); T2 = V[1]; zF = V[2]; jF = V[3]; wF = V[4]; cX = V[5]; nf = V[6];
f = x^3+2*x^2-4*x-4;
cprid(P) = idealhnf(nf, nfgaloisapply(nf, cX, P));
isfixed(P) = cprid(P) == idealhnf(nf, P);
inS(p) = (p % 15 == 14) && (kronecker(p, 11) == -1);
nroots(p) = #polrootsmod(f, p);
census(p) = {
  my(dec = idealprimedec(nf, p), fx = select(P -> isfixed(P), dec));
  [#dec, #fx, vecsort(apply(P -> P.f, fx)), vecsort(apply(P -> P.e, dec))];
}
print("== primes of S up to 3200 (plus the next split-completely ones): #primes of F above p, #c-fixed, residue degrees of fixed ones");
Slist = select(p -> inS(p), primes([2, 3200]));
print("S cap [2, 3200] = ", Slist);
bad = 0;
{
for (k = 1, #Slist, my(p = Slist[k], c = census(p), nr = nroots(p), pred = if (nr == 3, 24, 8),
     fpred = if (nr == 0, 6, 2));
  if (c[2] != pred || (c[2] > 0 && vecmax(c[3]) != fpred) || (c[2] > 0 && vecmin(c[3]) != fpred) || vecmax(c[4]) != 1, bad++);
  print("p = ", p, ": ", c[1], " primes, ", c[2], " fixed (predicted ", pred, "), fixed residue degrees ", Set(c[3]), " (predicted ", fpred, "), roots of f mod p: ", nr));
}
print("mismatches in S up to 3200: ", bad);
print("== control primes not in S (expect 0 fixed): ramified 2, 3, 5, 11, 37 and unramified others");
ctrl = [2, 3, 5, 11, 37, 7, 13, 17, 19, 23, 31, 41, 43, 47, 53, 59, 61, 89, 179, 269, 389, 419, 449, 509, 599, 719, 839, 929, 1049];
bad2 = 0;
for (k = 1, #ctrl, my(p = ctrl[k], c = census(p)); if (c[2] != 0, bad2++); print("p = ", p, " (p mod 15 = ", p % 15, ", (p/11) = ", kronecker(p, 11), "): ", c[1], " primes, ", c[2], " fixed, ramification indices ", Set(c[4])));
print("control mismatches: ", bad2);
print("== first primes of S with three roots of f mod p, p <= 20000");
three = select(p -> inS(p) && nroots(p) == 3, primes([2, 20000]));
print(three[1..min(#three, 12)], " ... (", #three, " in total up to 20000)");
print("roots mod 2549: ", polrootsmod(f, 2549));
quit;
