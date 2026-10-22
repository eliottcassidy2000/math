\\ Audit of THM-4755: rho(T) = det(dI - S)/|det(i d I - S)|, written from scratch.
\\ Tournaments are given by adjacency matrices A (A[i,j]=1 iff i->j).
default(realprecision, 60);
circA(d, D) = matrix(d, d, i, j, if(i != j && setsearch(Set(D), (j - i) % d), 1, 0));
skew(A) = A - A~;
isTour(A) = { my(d = #A); for(i = 1, d, if(A[i,i], return(0)); for(j = i + 1, d, if(A[i,j] + A[j,i] != 1, return(0)))); 1 };
isReg(A) = { my(d = #A, r = vecsum(A[1,])); for(i = 1, d, if(vecsum(A[i,]) != r, return(0))); 2*r == d - 1 };
rho(A) = { my(d = #A, S = skew(A), cp = charpoly(S, 'x)); subst(cp, 'x, d) / abs(subst(cp, 'x, I*d)) };
lexprod(AQ, AB) = { my(q = #AQ, b = #AB, M = matrix(q*b, q*b), i1, j1, k1, l1); for(r = 1, q*b, for(c = 1, q*b, i1 = (r-1)\b + 1; j1 = (r-1)%b + 1; k1 = (c-1)\b + 1; l1 = (c-1)%b + 1; M[r,c] = if(i1 != k1, AQ[i1,k1], AB[j1,l1]))); M };
Rd(d) = circA(d, vector((d-1)/2, r, r));
print("--- Wormald table, exact ---");
W = [[3, 2], [5, 44/19], [7, 3032/1231], [9, 437584/171769], [11, 108859552/41833969], [13, 41637757888/15769379963]];
for(t = 1, #W, my(d = W[t][1], r = rho(Rd(d))); print("R_", d, ": rho = ", r, "   Wormald ", W[t][2], "   equal: ", r == W[t][2]));
C3 = Rd(3); R9p = lexprod(C3, C3);
print("R'_9 = C3[C3]: tournament ", isTour(R9p), " regular ", isReg(R9p), " rho = ", rho(R9p), "  Wormald 5488/2197 equal: ", rho(R9p) == 5488/2197);
\\ Z_9{1,3,4,7} equals C3[C3]? (compare rho, and an explicit check of the 3-adic rule)
Z9 = circA(9, [1,3,4,7]); print("Z_9{1,3,4,7}: tournament ", isTour(Z9), " rho = ", rho(Z9));
\\ Wormald's 7-vertex example: vertices a,b,c,a',b',c',d = 1..7
A7 = matrix(7,7); arcs = [[1,2],[2,3],[3,1],[4,5],[5,6],[6,4],[1,4],[2,5],[3,6],[1,7],[2,7],[3,7]];
for(t = 1, #arcs, A7[arcs[t][1], arcs[t][2]] = 1);
\\ complete to a regular tournament: forced as derived by hand; verify uniqueness by brute force over remaining pairs
rem = List(); for(i = 1, 7, for(j = i+1, 7, if(A7[i,j] + A7[j,i] == 0, listput(rem, [i,j]))));
sols = List(); for(m = 0, 2^#rem - 1, my(B = A7); for(t = 1, #rem, my(p = rem[t]); if(bittest(m, t-1), B[p[1],p[2]] = 1, B[p[2],p[1]] = 1)); if(isReg(B), listput(sols, B)));
print("7-vertex example: completions to a regular tournament = ", #sols, "; rho = ", vector(#sols, s, rho(sols[s])), "  Wormald 3088/1287");
print("--- Wormald p.6 nine-vertex circulants ---");
Ds = [[1,3,4,5], [1,2,4,5], [1,2,3,5], [2,3,4,5], [1,3,5,6], [1,4,5,6], [2,3,5,6]]; vals = [1.4582, 1.6811, 2.4557, 1.5618, 1.5529, 1.9325, 1.8117];
for(t = 1, #Ds, my(D = Ds[t], A = circA(9, D), both = select(r -> setsearch(Set(D), 9 - r), D)); print("D = ", D, ": tournament ", isTour(A), ", pairs r & 9-r inside D: ", both, ", Wormald value ", vals[t], if(isTour(A), Str(", rho = ", rho(A), " = ", rho(A)*1.0), "")));
print("Jensen bound d=9: (10/8)^4 = ", (10/8)^4*1.0);
print("--- Large rotational, decimal ---");
Wl = [[21, "2.7215188", "2.8542501"], [41, "2.7867038", "2.8555079"], [81, "2.820801899", "2.855842467"], [101, "2.827747088", "2.855883644"]];
for(t = 1, #Wl, my(d = Wl[t][1], r = rho(Rd(d)) * 1.0, cf = ((d+1)^d + (d-1)^d) / real((d+I)^d + (d-I)^d)); printf("R_%d: rho_T = %.15f (Wormald %s...), closed form %.15f, rho_R = rho*e^(1/d) = %.15f (Wormald %s...)\n", d, r, Wl[t][2], cf*1.0, r*exp(1/d), Wl[t][3]));
print("cosh(1)/cos(1) = ", cosh(1)/cos(1));
print("--- Ternary tower S_k = C3^[k] ---");
Sk = C3; for(k = 1, 4, if(k > 1, Sk = lexprod(C3, Sk)); my(r = rho(Sk), pr = prod(i = 0, k-1, ((3^(2*i+1)+1)/(3^(2*i+1)-1))^(3^i))); print("k=", k, " n=", #Sk, ": tournament ", isTour(Sk), " regular ", isReg(Sk), "; rho(det) == level product: ", r == pr, "; rho = ", r*1.0));
\\ also the 3-adic Cayley model Cay(Z/3^k, unit part = 1 mod 3) has the same rho
sgn(z, m) = { my(v = 0); z = z % m; if(z == 0, return(0)); while(z % 3 == 0, z /= 3; v++); z % 3 };
for(k = 1, 4, my(m = 3^k, A = matrix(m, m, i, j, if(i != j && sgn(j - i, m) == 1, 1, 0))); print("Cay(Z/3^", k, ", sign 1): tournament ", isTour(A), "; rho equals level product: ", rho(A) == prod(i = 0, k-1, ((3^(2*i+1)+1)/(3^(2*i+1)-1))^(3^i))));
\\ prodinf with exact factors can be slow; do a log-sum instead
lk = suminf(i = 0, 3^i * log((3^(2*i+1)+1)/(3^(2*i+1)-1)));
print("kappa_S = ", exp(lk));
print("series 1 + 2 sum_{m>=3 odd} 3^-m/(m(1-3^(1-2m))) -> exp = ", exp(1 + 2*suminf(t = 1, my(m = 2*t+1); 3.^(-m)/(m*(1 - 3.^(1-2*m))))));
r9 = 5488/2197 * exp(1/9); print("(5488/2197) e^(1/9) = ", r9, "  relative excess of kappa_S: ", (exp(lk) - r9)/r9);
print("--- Lexicographic law spot checks ---");
P7 = circA(7, [1,2,4]); R5 = Rd(5);
lexlaw(AQ, AB) = { my(q = #AQ, b = #AB, d = q*b, SB = skew(AB), cp = charpoly(SB, 'x), ev = polroots(cp)); my(f = prod(t = 1, #ev, my(s = imag(ev[t])); if(abs(s) > 1e-20, ((d^2 + s^2)/(d^2 - s^2))^(1/2), 1))); rho(AQ) * f^q };
pairsL = [[R5, C3, "R5[C3]"], [C3, R5, "C3[R5]"], [P7, C3, "P7[C3]"], [C3, P7, "C3[P7]"], [R5, P7, "R5[P7]"]];
for(t = 1, #pairsL, my(AQ = pairsL[t][1], AB = pairsL[t][2], L = lexprod(AQ, AB)); printf("%s: rho(det) = %.20f  law = %.20f\n", pairsL[t][3], rho(L)*1.0, real(lexlaw(AQ, AB))));
