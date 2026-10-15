\\ Midline x = 1/2 of the 1 x a rectangle: 4y^2 + 1 = square (y = (t^2-1)/(4t)) and 4(a-y)^2 + 1 = square.
\\ Equivalently: a rational triangle with base 2a and height 1 (MacLeod's base/height ratio 2a).
Q = 4*t^2*(4*(a-(t^2-1)/(4*t))^2+1);
print("Q_a(t) = ", Q);
E = ellinit(ellfromeqn(y^2 - subst(Q, t, x)));
print("j(E_a) = ", factor(numerator(E.j)), " / ", factor(denominator(E.j)));
print("j at a = 1: ", subst(E.j, a, 1), "   (j(2i) = 287496)");
\\ ranks for small ratios r = 2a = num/den
{
for(n = 1, 12, for(d = 1, 6, if(gcd(n, d) == 1,
  aa = n/(2*d); Ea = ellinit(ellfromeqn(y^2 - subst(subst(Q, a, aa), t, x)));
  if(Ea == [], next);
  Em = ellminimalmodel(Ea);
  r = ellrank(Em);
  print("ratio 2a = ", n, "/", d, ": conductor ", ellglobalred(Em)[1], ", torsion ", elltors(Em)[1], " = Z^", elltors(Em)[2], ", rank in [", r[1], ",", r[2], "]")
)));
}
quit
