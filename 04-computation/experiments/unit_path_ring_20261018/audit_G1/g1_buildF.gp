\\ audit G1: build F = Q(zeta15, sqrt-11, w) as an absolute field with an explicit integral basis, independently of the session code.
default(parisize, 800000000);
f = x^3+2*x^2-4*x-4;
r1 = polcompositum(polcyclo(15), x^2+11, 1)[1];
K1 = r1[1]; zK = r1[2]; jK = r1[3];
r2 = polcompositum(K1, f, 1)[1];
T = r2[1]; thK = r2[2]; wT = r2[3];
[T2, a] = polredbest(T, 1);
\\ express generators mod T2
thK2 = subst(lift(thK), x, a); wF = subst(lift(wT), x, a);
zF = subst(lift(zK), x, thK2); jF = subst(lift(jK), x, thK2);
if (polcyclo(15, zF) != 0 || jF^2 + 11 != 0 || subst(f, x, wF) != 0, error("generator check failed"));
print("generators: Phi15(z) = 0, j^2 = -11, f(w) = 0 : ok");
\\ complex conjugation: z -> z^-1, j -> -j, w -> w ; image of the primitive root x of T2
\\ x = polynomial in (z, j, w): get it by solving: we know a = image of theta (root of T); invert: x as polmod is Mod(x, T2)
\\ c is a field automorphism; determine c(x) by expressing x in terms of z, j, w via linear algebra over Q
Bzjw = vector(48, k, my(i = (k-1)%8, e = ((k-1)\8)%2, g = (k-1)\16); lift(zF^i * jF^e * wF^g));
Mz = matrix(48, 48, r, s, polcoef(Bzjw[s], r-1));
cx = Mz^(-1) * vectorv(48, r, r == 2);   \\ coordinates of x in basis z^i j^e w^g
cBzjw = vector(48, k, my(i = (k-1)%8, e = ((k-1)\8)%2, g = (k-1)\16); zF^(15-i) * (-jF)^e * wF^g);
cX = sum(k = 1, 48, cx[k] * cBzjw[k]);
if (subst(T2, x, cX) != 0, error("c(x) not a root"));
cmap(y) = subst(lift(y), x, cX);
if (cmap(cmap(Mod(x, T2))) != Mod(x, T2), error("c not an involution"));
if (cmap(zF) * zF != 1 || cmap(jF) != -jF || cmap(wF) != wF, error("c action wrong"));
print("complex conjugation c: involution, z -> 1/z, j -> -j, w -> w : ok");
\\ integral basis O_K (x) Z[1, w, w^2/2]
B = vector(48, k, my(i = (k-1)%8, e = ((k-1)\8)%2, g = (k-1)\16); lift(zF^i * ((1 + jF)/2)^e * [1, wF, wF^2/2][g+1]));
dB = matdet(matrix(48, 48, r, s, trace(Mod(B[r]*B[s], T2))));
dexp = (nfdisc(polcyclo(15))^2 * 11^8)^3 * 148^16;
print("disc(B) = ", factor(dB), " ; expected disc(O_K)^3 disc(O_w)^16 = ", factor(dexp), " ; equal: ", dB == dexp);
t0 = getabstime();
dloc = nfdisc([T2, [2, 3, 5, 11, 37]]);
print("p-maximal discriminant at 2,3,5,11,37 (PARI round 4): valuations ", vector(5, k, valuation(dloc, [2,3,5,11,37][k])), " vs disc(B): ", vector(5, k, valuation(dB, [2,3,5,11,37][k])), "  [", getabstime()-t0, " ms]");
nf = nfinit([T2, B]);
print("nfinit done: degree ", poldegree(nf.pol), " disc = ", factor(nf.disc), " [", getabstime()-t0, " ms]");
writebin("g1_F.bin", [T2, zF, jF, wF, cX, nf]);
print("saved g1_F.bin");
quit;
