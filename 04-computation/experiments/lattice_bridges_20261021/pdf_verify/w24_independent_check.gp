read("Lraw_cand2.gp");
print("dims ", matsize(G), " symmetric ", G == G~);
print("all diagonal even ", vecmin(vector(#G, i, (G[i,i] % 2) == 0)));
d = matdet(G); print("det = 7^", valuation(d, 7), " * ", d / 7^valuation(d,7));
print("min diag ", vecmin(vector(#G, i, G[i,i])));
print("positive definite (qfsign) ", qfsign(G));
v = [0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]~;
print("witness norm v~*G*v = ", v~*G*v);
quit
