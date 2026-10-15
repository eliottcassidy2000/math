\\ D19 midline curve vs the congruent-number curve y^2 = x^3 - x (n = 1: Fermat's right-triangle theorem).
\\ (Fixed after audit B: the isogeny-class entries are [[a4, a6], map, map]; gp 2.17 tests isomorphism with ellisisom.)
Emid = ellinit([0,0,0,-11,14]);
Ecn = ellinit([0,0,0,-1,0]);
Eside = ellinit([0,0,0,-2,1]);
R = ellisomat(Emid);
L = R[1];
print("isogeny class of y^2 = x^3 - 11x + 14, as [a4, a6]: ", vector(#L, i, L[i][1]));
print("isogeny degree matrix: ", R[2]);
{
found = 0;
for(i = 1, #L,
  E = ellinit([0, 0, 0, L[i][1][1], L[i][1][2]]);
  if(ellisisom(E, Ecn) != 0, found = i));
print("member isomorphic over Q to y^2 = x^3 - x: index ", found, ", degree ", if(found, R[2][1, found], "-"),
      "; isomorphism ", if(found, ellisisom(ellinit([0, 0, 0, L[found][1][1], L[found][1][2]]), Ecn), "-"));
}
\\ explicit 2-isogeny from y^2 = x^3 - x with kernel <(-1, 0)> (Velu), and its image vs Emid
{
[Eimg, phi] = ellisogeny(Ecn, [-1, 0]);
Eimg = ellinit(Eimg);
print("Velu image of y^2 = x^3 - x by <(-1,0)>: [a1,a2,a3,a4,a6] = ", Eimg[1..5], "; isomorphic to Emid: ", ellisisom(Eimg, Emid) != 0, "; map ", phi);
}
print("a_p(Emid) == a_p(y^2 = x^3 - x) for all primes 3 <= p < 1000: ", prod(k = 2, primepi(1000), ellap(Emid, prime(k)) == ellap(Ecn, prime(k))));
print("side curve y^2 = x^3 - 2x + 1: conductor ", ellglobalred(Eside)[1], ", torsion ", elltors(Eside), ", class size ", #ellisomat(Eside)[1]);
quit
