\\ independent cross-check in PARI/GP using its own basis of M_48(SL2(Z))
mf = mfinit([1,48,1],4);
B = mfbasis(mf);
print("dim M_48 = ", #B);
M = matrix(13, #B, r, c, mfcoefs(B[c],12)[r]);
A = M[1..5,];
c = matsolve(A, [1,0,0,0,0]~);
f = M*c;
print("gp extremal dim96 q^0..q^8: ", vector(9,i,f[i]));
print("q^5 = ", f[6], "  q^6 = ", f[7]);
quit
