\\ D19 midline curve vs the congruent-number curve y^2 = x^3 - x (n = 1: Fermat's right-triangle theorem).
Emid = ellinit([0,0,0,-11,14]);
Ecn = ellinit([0,0,0,-1,0]);
Eside = ellinit([0,0,0,-2,1]);
R = ellisomat(Emid);
L = R[1];
print("isogeny class of y^2 = x^3 - 11x + 14 (as [a4, a6] or curves): ", L);
print("isogeny degree matrix: ", R[2]);
{
found = 0;
for(i = 1, #L,
  c = L[i];
  E = if(type(c) == "t_VEC" && #c == 2, ellinit([0,0,0,c[1],c[2]]), ellinit(c));
  if(E.j == Ecn.j && ellisomorphism(E, Ecn) != 0, found = i));
print("member isomorphic over Q to y^2 = x^3 - x: index ", found);
}
print("a_p(Emid) == a_p(y^2 = x^3 - x) for all primes 3 <= p < 1000: ", prod(k = 2, primepi(1000), ellap(Emid, prime(k)) == ellap(Ecn, prime(k))));
print("side curve y^2 = x^3 - 2x + 1: conductor ", ellglobalred(Eside)[1], ", torsion ", elltors(Eside), ", class size ", #ellisomat(Eside)[1]);
quit
