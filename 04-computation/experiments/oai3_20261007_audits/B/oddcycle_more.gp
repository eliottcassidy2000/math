default(parisizemax, 1000000000);
read("oddcycle_lib.gp");
\\ quadratic field Q(sqrt d): test generated subgroup of many unit vectors
testquad(d) = {
  my(r = quadgen(4*d), V = List(), W);
  for(a=-3,3, for(b=-2,2, for(c=-3,3, for(e=-2,2,
     my(A=a+b*r, C=c+e*r, den=A^2+C^2);
     if(den != 0, my(x=(A^2-C^2)/den, y=2*A*C/den); listput(V, [real(x), imag(x), real(y), imag(y)]~))))));
  W = Vec(Set(Vec(V)));
  my(t = oddtest(W));
  printf("Q(sqrt%d)^2: %d unit vectors; odd closed walk forced: %d\n", d, #W, t[1]);
};
testquad(5); testquad(13); testquad(37); testquad(17); testquad(7); testquad(11); testquad(23); testquad(3); testquad(2);
\\ explicit 11-cycle in Q(sqrt11)^2: 6 alternating steps (+-sqrt11/6, 5/6), then 5 steps (0,-1)
r = quadgen(44); P = [[0,0]]; s = [[r/6,5/6],[-r/6,5/6],[r/6,5/6],[-r/6,5/6],[r/6,5/6],[-r/6,5/6],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1]];
for(i=1,#s, P = concat(P, [P[#P] + s[i]]));
print("11-cycle closes: ", P[#P] == [0,0], "; all steps unit: ", vector(#s,i, s[i][1]^2+s[i][2]^2)==vector(#s,i,1), "; distinct vertices: ", #Set(P[1..#s]) == #s);
