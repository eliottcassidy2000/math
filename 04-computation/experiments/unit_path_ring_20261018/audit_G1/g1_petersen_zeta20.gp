\\ audit G1: the golden unit-distance Petersen realization in Q(zeta20) (THM-4701(4)), exactly, and its reductions.
\\ run: gp -q g1_petersen_zeta20.gp
P20 = polcyclo(20); nf20 = nfinit(P20);
z = Mod(x, P20); ii = z^5; z5 = z^4;
cj20(y) = subst(lift(y), x, Mod(x^19, P20));
N20(y) = y * cj20(y);
s5 = 2*(z5 + z5^4) + 1; phi = 1 + z5 + z5^4;
print("I^2 = -1: ", ii^2 == -1, " ; s5^2 = 5: ", s5^2 == 5, " ; phi^2 = phi + 1: ", phi^2 == phi + 1, " ; phi = (1 + s5)/2: ", phi == (1 + s5)/2);
R = (z + 1/z)/s5; r = R/phi;
print("R^2 = (5 + s5)/10: ", R^2 == (5 + s5)/10, " ; r^2 = (5 - s5)/10: ", r^2 == (5 - s5)/10, " ; R^2 + r^2 = 1: ", R^2 + r^2 == 1, " ; R r = 1/s5: ", R*r == 1/s5, " ; R, r real: ", cj20(R) == R && cj20(r) == r);
pts = concat(vector(5, k, R*z5^(k-1)), vector(5, k, ii*r*z5^(k-1)));   \\ o_0..o_4, i_0..i_4
E = concat(concat(vector(5, k, [k, k % 5 + 1]), vector(5, k, [k, k + 5])), vector(5, k, [k + 5, (k + 1) % 5 + 6]));
print("edges: ", E);
print("all 15 edges are unit vectors: ", prod(e = 1, 15, N20(pts[E[e][1]] - pts[E[e][2]]) == 1));
\\ graph check: the 15 pairs form the Petersen graph (3-regular, girth 5, 10 vertices)
deg = vector(10); for (e = 1, 15, deg[E[e][1]]++; deg[E[e][2]]++); print("degrees: ", deg);
\\ non-edges: no other pair is at unit distance (so the realization is also induced)
nonedge = 0; for (a = 1, 10, for (b = a + 1, 10, if (N20(pts[a] - pts[b]) == 1, nonedge++))); print("pairs at unit distance in total: ", nonedge, " (15 = exactly the edges)");
\\ norms of all 45 differences: primes that can make two points collide under reduction
pr = Set(); for (a = 1, 10, for (b = a + 1, 10, my(nn = norm(pts[a] - pts[b])); pr = setunion(pr, Set(factor(numerator(nn))[,1]~)); pr = setunion(pr, Set(factor(denominator(nn))[,1]~))));
print("primes dividing numerators or denominators of N_(Q(zeta20)/Q)(p_a - p_b), all 45 pairs: ", pr);
print("primes dividing norms of the 10 points: ", Set(concat(vector(10, a, factor(numerator(norm(pts[a])))[,1]~))), " / denominators ", Set(concat(vector(10, a, factor(denominator(norm(pts[a])))[,1]~))));
\\ reduction at the prime above 2
dec2 = idealprimedec(nf20, 2); print("primes above 2: ", #dec2, " with e = ", dec2[1].e, ", f = ", dec2[1].f);
mp = nfmodprinit(nf20, dec2[1]);
red2(y) = nfmodpr(nf20, lift(y), mp);
rp = vector(10, a, red2(pts[a]));
w = red2(R);
print("R -> w with w^3 = 1, w != 1: ", w^3 == 1 && w != 1, " ; r -> w^2: ", red2(r) == w^2, " ; i -> 1: ", red2(ii) == 1, " ; phi -> w^2: ", red2(phi) == w^2);
print("images distinct: ", #Set(rp) == 10, " ; none in mu_5 or 0: ", prod(a = 1, 10, rp[a] != 0 && rp[a]^5 != 1), " ; edges map to mu_5 (x^5 = 1): ", prod(e = 1, 15, (rp[E[e][1]] - rp[E[e][2]])^5 == 1));
print("complex conjugation reduces to x -> x^4 at the prime above 2: ", red2(cj20(z)) == red2(z)^4);
\\ reductions at q = 19 (mod 20)
{
foreach([19, 59, 79, 139, 179, 199, 239, 359], q,
  my(dec = idealprimedec(nf20, q), okall = 1);
  for (k = 1, #dec, my(P = dec[k], m = nfmodprinit(nf20, P), im, fixed);
    fixed = idealhnf(nf20, nfgaloisapply(nf20, x^19, P)) == idealhnf(nf20, P);
    im = vector(10, a, nfmodpr(nf20, lift(pts[a]), m));
    okall = okall && fixed && P.f == 2 && #Set(im) == 10 && prod(e = 1, 15, my(d = im[E[e][1]] - im[E[e][2]]); d^(q + 1) == 1)
            && nfmodpr(nf20, lift(cj20(z)), m) == nfmodpr(nf20, lift(z), m)^q);
  print("q = ", q, ": ", #dec, " primes, all c-fixed with residue F_(q^2), conj -> Frob_q, reduction injective and edges -> mu_(q+1): ", okall));
}
quit;
