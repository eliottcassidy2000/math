\\ Guy D19 special curves.  Midline x = 1/2: y = (t^2-1)/(4t) gives 4y^2+1 = ((t^2+1)/(2t))^2, and the second
\\ condition 4(1-y)^2+1 = square becomes w^2 = Qmid(t).  Side y = 0: x = (s^2-1)/(2s) gives x^2+1 = square, and
\\ (1-x)^2+1 = square becomes w^2 = Qside(s).  Each is a genus-1 quartic; rank and torsion of its Jacobian decide
\\ its rational points.
default(parisize, 200000000);
Qmid = 4*t^2*(4*(1-(t^2-1)/(4*t))^2+1);
Qside = 4*t^2*((1-(t^2-1)/(2*t))^2+1);
print("Qmid(t)  = ", Qmid);
print("Qside(t) = ", Qside);
{
for(i = 1, 2,
  Q = if(i == 1, Qmid, Qside);
  nm = if(i == 1, "midline x = 1/2", "side y = 0");
  E = ellinit(ellfromeqn(y^2 - subst(Q, t, x)));
  Em = ellminimalmodel(E);
  print("---- ", nm);
  print("  minimal model [a1,a2,a3,a4,a6] = ", Em[1..5], ", conductor ", ellglobalred(Em)[1], ", j = ", Em.j);
  print("  torsion: ", elltors(Em));
  r = ellrank(Em);
  print("  ellrank: rank in [", r[1], ", ", r[2], "]");
  print("  analytic rank: ", ellanalyticrank(Em)[1]);
);
}
quit
