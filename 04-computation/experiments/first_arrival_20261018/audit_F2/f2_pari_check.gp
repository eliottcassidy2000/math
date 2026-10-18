\\ Audit F2 -- independent PARI/GP cross-check of the 19-adic valuations of THM-4690(5).
\\ b, c, a = conj(d), the spindle and the sheet moves live in K = Q(zeta15, sqrt(-11)) (degree 16);
\\ u lives in L = Q(w, sqrt(-11)) (degree 6); u2 in Q(w, sqrt(-11), sqrt5) (degree 12).
\\ Valuations are PARI's nfeltval at idealprimedec(., 19): exact, no reduction ambiguity.
\\ F/K and F/L are unramified above 19, so these are the valuations at the primes of F above them.
default(parisize, 200000000);

\\ ---------------- K = Q(zeta15, j)
pc = polcompositum(polcyclo(15), x^2 + 11, 1)[1];
RK = pc[1]; zt = pc[2]; jj = pc[3];
K = nfinit(RK);
print("K: degree ", poldegree(RK), ", zeta of order 15: ", zt^15 == 1 && zt^5 != 1 && zt^3 != 1, ", j^2 = -11: ", jj^2 == -11);
alpha = zt^3; phi = 1 + alpha + alpha^4; d = 1 + alpha; a = 1 + alpha^4;
m2 = -(20 + 11*jj)/19; m1 = (20 - 11*jj)/19;
c = (1 + a*m2 - a^2*m1 + a^3)/(d*a - 1);
b = c*d + a*m1 - a^2;
PK = idealprimedec(K, 19);
print("primes of K above 19: ", #PK, " ; (e, f) = ", vector(#PK, i, [PK[i].e, PK[i].f]));
print("v(b) at the K-primes: ", vector(#PK, i, nfeltval(K, b, PK[i])));
print("v(c) at the K-primes: ", vector(#PK, i, nfeltval(K, c, PK[i])));
print("v(b/c): ", vector(#PK, i, nfeltval(K, b/c, PK[i])), " ; v(a): ", vector(#PK, i, nfeltval(K, a, PK[i])));
eta6 = 1 + zt^5; rho = (5 + jj)/6;
Mset = [0, 1, eta6, 1 + eta6, rho, rho*eta6, rho*(1 + eta6)];
Dl = List();
{
for (i = 1, 7, for (k = 1, 7,
  my(t = Mset[i] - Mset[k], new = 1);
  for (l = 1, #Dl, if (Dl[l] == t, new = 0; break));
  if (new, listput(Dl, t))));
}
print("|M - M| (with 0) = ", #Dl);
print("nonzero spindle differences: v = 0 at all K-primes: ", prod(l = 1, #Dl, if (Dl[l] == 0, 1, prod(i = 1, #PK, nfeltval(K, Dl[l], PK[i]) == 0))));
cnt = vector(9); allint = 0; npairs = 0;
{
for (i2 = 1, #Dl, for (i3 = 1, #Dl,
  if (Dl[i2] == 0 && Dl[i3] == 0, next);
  npairs++;
  my(yy = b*Dl[i2] + c*Dl[i3], s = 0);
  if (yy == 0, s = #PK, s = sum(i = 1, #PK, nfeltval(K, yy, PK[i]) >= 0));
  cnt[s + 1]++;
  if (s == #PK, allint++)));
}
print("sheet moves: ", npairs, " nonzero pairs; [#K-primes where b dm2 + c dm3 is integral, count]: ", vector(9, s, [s - 1, cnt[s]]));
print("pairs 19-integral: ", allint);
\\ complex conjugation and the splitting of the primes of K above 19
G = nfgaloisconj(K);
cc = 0;
{
for (i = 1, #G,
  my(s = Mod(G[i], RK));
  if (subst(lift(zt), x, s) == zt^14 && subst(lift(jj), x, s) == -jj, cc = G[i]));
}
print("complex conjugation found: ", cc != 0);
fix = sum(i = 1, #PK, idealhnf(K, nfgaloisapply(K, cc, PK[i])) == idealhnf(K, PK[i]));
print("primes of K above 19 fixed by complex conjugation: ", fix, " of ", #PK);
\\ explicit delta = conj(gamma)/gamma with a pole at PK[1]
q = PK[1]; qc = nfgaloisapply(K, cc, q);
gam = nfbasistoalg(K, q.gen[2]); kk = 0;
while (nfeltval(K, gam, qc) > 0 || nfeltval(K, gam, q) < 1, kk++; gam = gam + 19*kk);
del = nfbasistoalg(K, nfgaloisapply(K, cc, gam)) / gam;
delc = nfbasistoalg(K, nfgaloisapply(K, cc, del));
print("delta = conj(gamma)/gamma: delta*conj(delta) = 1: ", del*delc == 1, " ; v at the K-primes above 19: ", vector(#PK, i, nfeltval(K, del, PK[i])));

\\ ---------------- L = Q(w, j)
pl = polcompositum(x^3 + 2*x^2 - 4*x - 4, x^2 + 11, 1)[1];
RL = pl[1]; ww = pl[2]; jl = pl[3];
L = nfinit(RL);
h = ww + 3; u = (h - jl)/(h + jl);
print("L: degree ", poldegree(RL), " ; N_{L/Q}(h + j) = ", nfeltnorm(L, h + jl), " ; N_{L/Q}(h - j) = ", nfeltnorm(L, h - jl), " ; 76^2 = ", 76^2);
PL = idealprimedec(L, 19);
print("primes of L above 19: ", #PL, " ; (e, f) = ", vector(#PL, i, [PL[i].e, PL[i].f]));
{
for (i = 1, #PL,
  my(t = if (nfeltval(L, h - 7, PL[i]) > 0, "h=7", if (nfeltval(L, h - jl, PL[i]) > 0, "h=+j", if (nfeltval(L, h + jl, PL[i]) > 0, "h=-j", "?"))));
  print("  L-prime ", i, " type ", t, " : v(h - j) = ", nfeltval(L, h - jl, PL[i]), ", v(h + j) = ", nfeltval(L, h + jl, PL[i]), ", v(u) = ", nfeltval(L, u, PL[i])));
}
\\ ---------------- Q(w, j, sqrt5) for u2
p3 = polcompositum(RL, x^2 - x - 1, 1)[1];
R3 = p3[1]; th = p3[2]; ph = p3[3];
L3 = nfinit(R3);
w3 = subst(lift(ww), x, th); j3 = subst(lift(jl), x, th);
h3 = w3 + 3; u2 = (h3 + ph - j3)/(h3 + ph + j3);
P3 = idealprimedec(L3, 19);
print("Q(w,j,sqrt5): degree ", poldegree(R3), " ; primes above 19: ", #P3, " ; v(u2) = ", vector(#P3, i, nfeltval(L3, u2, P3[i])));
print("   and v(u) there = ", vector(#P3, i, nfeltval(L3, (h3 - j3)/(h3 + j3), P3[i])), " ; phi^2 = phi + 1: ", ph^2 == ph + 1);
quit;
