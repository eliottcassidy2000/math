default(parisizemax, 2000000000);
read("Lraw_cand2.gp");
d = matdet(G); print("cand 2: det = 7^", valuation(d,7), " cofactor ", d/7^valuation(d,7));
U = qflllgram(G); R = U~*G*U;
dg = vector(96,i,R[i,i]); print("cand 2: LLL diag min = ", vecmin(dg), "  sorted first 10: ", vecsort(dg)[1..10]);
\ dual: 7*G^-1 should be integral and similar to G (7-modular check)
Gd = 7*G^-1; print("cand 2: 7*G^-1 integral: ", denominator(Gd)==1);
quit
