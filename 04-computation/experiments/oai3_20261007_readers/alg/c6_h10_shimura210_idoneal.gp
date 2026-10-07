\\ C6: the five contact points of openai/math #004's height estimate.
\\ B = quaternion algebra over Q of discriminant D = 210, maximal order; X = Shimura curve (genus 5);
\\ W = Atkin-Lehner group (Z/2)^4; X/W = P^1 with five branch values labelled m = 30, 42, 70, 105, 210.
\\ Fixed points of w_m (Ogg / Eichler): sum over orders O containing an element of norm m generating Q(sqrt(-m)),
\\   #Fix(w_m) = sum_O h(O) * prod_{p | D/m} (1 - kronecker(disc O, p)),  O embeddable in B (non-split at p | D).
D = 210; P = factor(D)[,1];
fixcount(m) = {
  \\ orders O = Z[f w_K] of discriminant d = f^2 dK; optimal embedding into the maximal order of B needs
  \\ K non-split at every p | D, and at p | D the local order must be maximal (p not dividing f);
  \\ local factor at p | D/m is 1 - (Eichler symbol) = 1 - kronecker(dK, p) when p does not divide f.
  my(discs = if (m == 2, [-4, -8], if (m % 4 == 3, [-m, -4*m], [-4*m])), tot = 0, detail = List());
  for (i = 1, #discs, my(d = discs[i], dK = coredisc(d), f = sqrtint(d/dK), ok = 1, fac = 1);
    for (j = 1, #P, my(p = P[j]); if (kronecker(dK, p) == 1 || f % p == 0, ok = 0));
    if (ok, for (j = 1, #P, my(p = P[j]); if ((D/m) % p == 0, fac *= (1 - kronecker(dK, p))));
            my(h = qfbclassno(d)); tot += h*fac; listput(detail, [d, h, fac, quadclassunit(d).cyc])));
  [tot, Vec(detail)]
};
total = 0;
fordiv(D, m, if (m > 1, my(r = fixcount(m)); total += r[1]; if (r[1] > 0, print("m = ", m, ": #Fix(w_m) = ", r[1], "   [disc, h, local factor, class group] = ", r[2]))));
print("sum of fixed points = ", total, "   (Riemann-Hurwitz: 5 branch values x |W|/2 = 5 x 8 = 40)");
\\ idoneal check: n idoneal <=> class group of disc -4n is elementary 2-abelian (one class per genus)
idoneal(n) = {my(c = quadclassunit(-4*n).cyc); for (i = 1, #c, if (c[i] > 2, return(0))); 1};
L = select(n -> idoneal(n), [1..2000]);
print("idoneal numbers <= 2000 (class group of -4n is 2-torsion): ", #L, " numbers: ", L);
print("labels 30, 42, 70, 105, 210 idoneal? ", vector(5, i, idoneal([30,42,70,105,210][i])));
\\ D = 210 itself and the other 4-prime discriminants: which squarefree D = product of 4 primes have all labels idoneal?
