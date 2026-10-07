\\ C3: 2-adic integral preimage counts N(F(P)) for Haar-random P in Z_2^3, via the fibre cubic
\\ phi(Y) = -2Y^3 + 3bY^2 - 18aY + (18ab - b^3 - 27a^2c), x = (b-y)/(y^2-by+3a), z = (2x-3x^2y-c)/x^3.
default(parisize, 200000000);
F(P) = {my(X=P[1],Y=P[2],Z=P[3],U=1+X*Y); [U^3*Z + Y^2*U*(4+3*X*Y), Y + 3*X*U^2*Z + 3*X*Y^2*(4+3*X*Y), 2*X - 3*X^2*Y - X^3*Z]};
\\ sanity: the collision
print("F(0,0,-1/4) = ", F([0,0,-1/4]), "  F(1,-3/2,13/2) = ", F([1,-3/2,13/2]), "  F(-1,3/2,13/2) = ", F([-1,3/2,13/2]));
count(P, prec) = {
  my(w = F(P), a = w[1], b = w[2], c = w[3], phi, R, N = 0, pts = List());
  phi = -2*'Y^3 + 3*b*'Y^2 - 18*a*'Y + (18*a*b - b^3 - 27*a^2*c);
  R = polrootspadic(phi, 2, prec);
  for (i = 1, #R,
    my(yy = R[i], den = yy^2 - b*yy + 3*a, xx, zz);
    if (valuation(yy - b, 2) >= prec - 5, xx = 0; zz = a - 4*b^2,
        xx = (b - yy)/den; zz = (2*xx - 3*xx^2*yy - c)/xx^3);
    if (valuation(xx,2) >= 0 && valuation(yy,2) >= 0 && valuation(zz,2) >= 0, N++; listput(pts, [xx,yy,zz])));
  [N, #R, pts]
};
setrand(20261007);
T = 4000; H = vector(4); HR = vector(4);
for (t = 1, T, my(P = vector(3, i, random(2^40)), r = count(P, 80)); H[r[1]+1]++; HR[r[2]+1]++);
print("source-Haar sample of ", T, " points P: #(integral preimages of F(P)) histogram N=0..3: ", H);
print("   (number of Q_2-roots of the fibre cubic, 0..3): ", HR);
print("   predicted from mod-2^k census: N=1: 7/16 = ", 7/16.*T, ", N=2: 3/8 = ", 3/8.*T, ", N=3: 3/16 = ", 3/16.*T);
\\ explicit example of a 2-adic merge
setrand(7);
{for (t = 1, 200, my(P = vector(3, i, random(2^12)), r = count(P, 60));
  if (r[1] >= 2, print("example P = ", P, ": N = ", r[1]);
     for (j = 1, #r[3], my(Q = r[3][j]); print("   preimage mod 2^20: ", vector(3, i, lift(Q[i] + O(2^20)))));
     break));}
\\ rational collisions inside Z_(2)^3: P integer, other fibre points rational with odd denominators
{found = 0;
forvec(P = [[-6,6],[-6,6],[-6,6]],
  my(w = F(P), a = w[1], b = w[2], c = w[3], phi, q, D);
  phi = -2*'Y^3 + 3*b*'Y^2 - 18*a*'Y + (18*a*b - b^3 - 27*a^2*c);
  q = phi \ ('Y - P[2]);
  if (poldegree(q) == 2 && polcoeff(phi % ('Y - P[2]), 0) == 0,
    D = poldisc(q);
    if (D != 0 && issquare(D),
      my(rts = nfroots(, q), pts = List());
      for (j = 1, #rts, my(yy = rts[j], den = yy^2 - b*yy + 3*a, xx, zz);
        if (yy == b, xx = 0; zz = a - 4*b^2, xx = (b - yy)/den; zz = (2*xx - 3*xx^2*yy - c)/xx^3);
        listput(pts, [xx, yy, zz]));
      my(ok = 1); for (j = 1, #pts, for (i = 1, 3, if (valuation(pts[j][i], 2) < 0, ok = 0)));
      if (ok && found < 8, found++; print("Z_(2)-collision: F", P, " = F", pts[1], " = F", pts[2], " = ", w,
          "   check: ", F(pts[1]) == w && F(pts[2]) == w)))));
print("rational Z_(2)^3 triple collisions found (|coords| <= 6 for P): ", found);}
