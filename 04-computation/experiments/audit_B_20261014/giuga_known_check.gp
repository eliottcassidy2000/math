\\ Audit B: the 13 known Giuga numbers (12 terms of OEIS A007850 + the 97-digit 10-factor one in its comments).
\\ Factorisations are recomputed with factor() for the 12 terms; for the 13th the listed factors are proved prime (isprime,
\\ a primality proof) and multiplied.  Checks: Giuga condition, index sum 1/p - 1/g, arithmetic derivative, Korselt score,
\\ and PROVED primality (isprime) of g - 1 and g + 1.
G = [30, 858, 1722, 66198, 2214408306, 24423128562, 432749205173838, 14737133470010574, 550843391309130318, 244197000982499715087866346, 554079914617070801288578559178, 1910667181420507984555759916338506];
P13 = [2,3,11,23,31,47059,2217342227,1729101023519,8491659218261819498490029296021,58254480569119734123541298976556403];
{
  for(i = 1, 13,
    my(ps, g);
    if(i <= 12, g = G[i]; ps = factor(g)[,1]~; if(vecmax(factor(g)[,2]) > 1, error("not squarefree")),
                ps = P13; g = prod(j = 1, #ps, ps[j]); if(!vecmin(apply(isprime, ps)), error("factor not prime")));
    my(giuga = vecmin(apply(p -> (g/p - 1) % p == 0, ps)));
    my(idx = sum(j = 1, #ps, 1/ps[j]) - 1/g);
    my(der = sum(j = 1, #ps, g/ps[j]));
    my(kor = sum(j = 1, #ps, (g - 1) % (ps[j] - 1) == 0));
    printf("G%2d: %3d digits, omega %2d, Giuga %d, index %s, g' = g+1: %d, Korselt score %d, g-1 prime %d, g+1 prime %d\n",
      i, #Str(g), #ps, giuga, Str(idx), der == g + 1, kor, isprime(g - 1), isprime(g + 1)));
}
quit
