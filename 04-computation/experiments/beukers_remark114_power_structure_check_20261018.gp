\\ Beukers' Remark 1.14(2) families (Cooper s7, s10, s18, Apery2*central, Franel*central) and inert CM primes.
\\ Checks for 05-knowledge/results/beukers_remark114_power_structure_and_inert_cm_20261018.md
\\ (session cloud-2026-10-18-M1).
\\
\\ Each family: F(t) = sum g_n t^n solves a third-order equation, t(q) is a Hauptmodul of an Atkin-Lehner extension
\\ of Gamma0(N), Theta(q) = F(t(q)) has weight 2, and t^a Theta^b is an eta product (a = b - 1).
\\ T_p(x) = (1 - Theta(q^p)/Theta x) prod_k (p^2 - Theta((tau+k)/p)/Theta x) is Beukers' polynomial (his (6)),
\\ T_{p,m} its coefficient of x^m, t^sigma = t(q^p) the excellent lift.
\\
\\ Sections: I (identities), T (T_{p,p+1} = +-Q^a and the degree count), F (F == F_p F(t^sigma) mod p^2),
\\ N (Newton polygons of the full T_p(x, t0) at CM points), CM (rational CM values of all thirteen families of
\\ Beukers' Appendix C that have them: inert primes, split primes, ramified primes), S (Z.-W. Sun's Apery-sum conjecture).
\\
\\ Run: gp -q 04-computation/experiments/beukers_remark114_power_structure_check_20261018.gp

default(parisizemax, 8*10^9);
default(realprecision, 80);

LT = 58; PMAXT = 61;    \\ section T: primes p <= PMAXT, t-adic precision LT
K = 3; PMAXF = 110;     \\ section F: primes p < PMAXF, coefficients of t^0 .. t^(K p)
PMAXCM = 400;           \\ section CM: primes 5 <= p < PMAXCM

\\ n^3 g_n = A(n-1) g_(n-1) - B(n-2) g_(n-2) for the operator theta^3 - t A(theta) + t^2 B(theta)
gseq(A, B, M) = {
  my(P1 = m -> A[1]*m^3 + A[2]*m^2 + A[3]*m + A[4], P2 = m -> B[1]*m^3 + B[2]*m^2 + B[3]*m + B[4]);
  my(g = vector(M+1)); g[1] = 1;
  for(n = 1, M, g[n+1] = (P1(n-1)*g[n] - if(n >= 2, P2(n-2)*g[n-1], 0))/n^3);
  g
};
\\ [g, q(t)]: q(t) = t exp(F1/F) from the Frobenius solution F log t + F1
mirror(A, B, M) = {
  my(P1 = m -> A[1]*m^3 + A[2]*m^2 + A[3]*m + A[4], P2 = m -> B[1]*m^3 + B[2]*m^2 + B[3]*m + B[4]);
  my(D1 = m -> 3*A[1]*m^2 + 2*A[2]*m + A[3], D2 = m -> 3*B[1]*m^2 + 2*B[2]*m + B[3]);
  my(g = vector(M+1), f = vector(M+1));
  g[1] = 1; f[1] = 0;
  for(n = 1, M,
    g[n+1] = (P1(n-1)*g[n] - if(n >= 2, P2(n-2)*g[n-1], 0))/n^3;
    my(rhs = 3*n^2*g[n+1] - D1(n-1)*g[n] + if(n >= 2, D2(n-2)*g[n-1], 0));
    f[n+1] = (P1(n-1)*f[n] - if(n >= 2, P2(n-2)*f[n-1], 0) - rhs)/n^3);
  [g, x*exp(Ser(f, x)/Ser(g, x))]
};
bn(m, k) = if(m < 0, 0, binomial(m, k));
etq(m, M) = subst(eta(x + O(x^(M\m + 1))), x, x^m);
E(m, tau) = eta(m*tau, 1);
\\ q-series [t(q), Theta(q)] with Theta = (eta product / t^a)^(1/b), and t(tau) numerically
tTh7(M)  = { my(h = x*(etq(7,M)/etq(1,M))^4, t = h/(1+13*h+49*h^2), e = x^2*(etq(1,M)*etq(7,M))^6); [t, (e/t^2)^(1/3)] };
tTh10(M) = { my(h = x*etq(2,M)*etq(10,M)^3/(etq(1,M)^3*etq(5,M)), t = h*(1+4*h)*(1+5*h)/(1+8*h+20*h^2)^2,
                e = x^3*(etq(1,M)*etq(2,M)*etq(5,M)*etq(10,M))^4); [t, (e/t^3)^(1/4)] };
tTh18(M) = { my(h = x*etq(2,M)*etq(3,M)*etq(18,M)^2/(etq(1,M)^2*etq(6,M)*etq(9,M)),
                t = h*(1+2*h)*(1+3*h)*(1+3*h+3*h^2)*(1+6*h+12*h^2)/(1+6*h+6*h^2)^4,
                e = x^3*(etq(3,M)*etq(6,M))^8); [t, (e/t^3)^(1/4)] };
tThA5(M) = { my(h = x*(etq(5,M)/etq(1,M))^6, t = h/(1+22*h+125*h^2), e = x*(etq(1,M)*etq(5,M))^4); [t, (e/t)^(1/2)] };
tThF6(M) = { my(a = etq(2,M)*etq(3,M), b = etq(1,M)*etq(6,M), Th = a^7/b^5 - x*b^7/a^5);
             [x*(a*b)^2/Th^2, Th] };
tThAp(M) = { my(a = etq(2,M)*etq(3,M), b = etq(1,M)*etq(6,M)); [x*(b/a)^12, a^7/b^5] };
t7n(tau)  = { my(h = (E(7,tau)/E(1,tau))^4); h/(1+13*h+49*h^2) };
t10n(tau) = { my(h = E(2,tau)*E(10,tau)^3/(E(1,tau)^3*E(5,tau))); h*(1+4*h)*(1+5*h)/(1+8*h+20*h^2)^2 };
t18n(tau) = { my(h = E(2,tau)*E(3,tau)*E(18,tau)^2/(E(1,tau)^2*E(6,tau)*E(9,tau)));
              h*(1+2*h)*(1+3*h)*(1+3*h+3*h^2)*(1+6*h+12*h^2)/(1+6*h+6*h^2)^4 };
tA5n(tau) = { my(h = (E(5,tau)/E(1,tau))^6); h/(1+22*h+125*h^2) };
tF6n(tau) = { my(a = E(2,tau)*E(3,tau), b = E(1,tau)*E(6,tau), Th = a^7/b^5 - b^7/a^5); (a*b)^2/Th^2 };

\\ [name, A, B, tTh, a, b, N, binomial formula]
{FAM = [
 ["s7", [26,39,21,4], [-27,-81,-78,-24], tTh7, 2, 3, 7,
  n -> sum(k=0, n, binomial(n,k)^2*binomial(n+k,k)*binomial(2*k,n))],
 ["s10", [12,18,10,2], [-64,-192,-188,-60], tTh10, 3, 4, 10,
  n -> sum(k=0, n, binomial(n,k)^4)],
 ["s18", [28,42,26,6], [192,576,564,180], tTh18, 3, 4, 18,
  n -> sum(k=0, n\2, (-1)^k*binomial(n,k)*binomial(2*k,k)*binomial(2*n-2*k,n-k)*(bn(2*n-3*k-1,n)+bn(2*n-3*k,n)))],
 ["Apery2*central", [44,66,34,6], [-16,-48,-44,-12], tThA5, 1, 2, 5,
  n -> binomial(2*n,n)*sum(k=0, n, binomial(n,k)^2*binomial(n+k,k))],
 ["Franel*central", [28,42,22,4], [-128,-384,-352,-96], tThF6, 1, 2, 6,
  n -> binomial(2*n,n)*sum(k=0, n, binomial(n,k)^3)]
];}

\\ ---------------------------------------------------------------------------------------------- I
MI = 120;
print("== I: identities (n <= ", MI, ")");
okI = 1;
{for(i = 1, #FAM, my(f = FAM[i], S = f[4](MI), g = gseq(f[2], f[3], MI), R = mirror(f[2], f[3], MI), r);
  r = [vector(MI+1, n, f[8](n-1)) == g, denominator(content(truncate(R[2]))) == 1,
       denominator(content(truncate(S[2]))) == 1, valuation(subst(truncate(Ser(g, x)), x, S[1]) - S[2], x) > MI,
       valuation(subst(truncate(S[1]), x, R[2]) - x, x) > MI];
  okI = okI && vecmin(r);
  printf("%-15s binomial formula = ODE: %d;  q(t) in Z[[t]]: %d;  Theta in Z[[q]]: %d;  F(t(q)) = Theta: %d;  t(q(t)) = t: %d\n",
    f[1], r[1], r[2], r[3], r[4], r[5]));}
print("I all: ", okI);

\\ ---------------------------------------------------------------------------------------------- T
minval(d, p, L) = { my(v = 99, cc); for(k = 0, L, cc = polcoeff(d, k); if(cc != 0, v = min(v, valuation(cc, p)))); v };
print("\n== T: T_{p,p+1} = (-1)^(p+1) Q^a, Q in 1 + p Z[t], deg Q <= floor((p+1)/b), deg T_{p,p} <= floor(p a/b),");
print("       and T_{p,p} + T_{p,p+1} F_p == 0 mod p^2 as polynomials; primes p <= ", PMAXT, ", p !| N");
okT = 1;
{for(i = 1, #FAM, my(f = FAM[i], MQ = PMAXT*LT, S = f[4](MQ), Th = S[2], a = f[5], b = f[6], N = f[7], qt, logTh, invTh, g, ok = 1, res = []);
  qt = mirror(f[2], f[3], LT)[2]; g = gseq(f[2], f[3], PMAXT);
  logTh = log(Th); invTh = 1/Th;
  forprime(p = 2, PMAXT, if(N % p == 0, next);
    my(L = min(MQ\p, LT), ThL, Np, Up, Rr, ep, T1, T0, Q, d1, d0, dQ, ok1, Fp, rel, b0 = (p*a)\b, b1 = ((p+1)*a)\b, bQ = (p+1)\b);
    ThL = Th + O(x^(L+1));
    Np = exp(p*sum(m=1, L, polcoeff(logTh, p*m)*x^m) + O(x^(L+1)));    \\ prod_k Theta((tau+k)/p) = exp(p U_p log Theta)
    Up = sum(m=0, L, polcoeff(invTh, p*m)*x^m) + O(x^(L+1));           \\ U_p(1/Theta)
    Rr = subst(ThL, x, x^p)/ThL; ep = Np/ThL^p;
    T1 = subst(truncate((-1)^(p+1)*Rr*ep), x, qt + O(x^(L+1)));
    T0 = subst(truncate((-1)^p*ep*(1 + p^3*Rr*ThL*Up)), x, qt + O(x^(L+1)));
    d1 = poldegree(truncate(T1)); d0 = poldegree(truncate(T0));
    Q = ((-1)^(p+1)*T1)^(1/a); dQ = poldegree(truncate(Q));
    Fp = sum(n = 0, p-1, g[n+1]*x^n);
    rel = truncate(T0) + truncate(T1)*Fp;
    ok1 = (L >= b1 + 4) && denominator(content(truncate(T1))) == 1 && denominator(content(truncate(T0))) == 1
          && denominator(content(truncate(Q))) == 1 && d0 <= b0 && d1 <= b1 && dQ <= bQ
          && (content(truncate(Q) - 1) % p == 0) && (content(rel) % p^2 == 0);
    if(!ok1, printf("  FAIL %s p=%d\n", f[1], p));
    ok = ok && ok1; res = concat(res, [[p, d0, b0, dQ, bQ]]));
  okT = okT && ok;
  printf("%-15s (a,b) = (%d,%d): %d;  [p, deg T_pp, floor(pa/b), deg Q, floor((p+1)/b)] for the last three p: %s\n",
    f[1], a, b, ok, Str(res[#res-2..#res])));}
print("T all: ", okT);

\\ ---------------------------------------------------------------------------------------------- F
print("\n== F: F(t) == F_p(t) F(t^sigma) mod p^2, t^sigma = t(q(t)^p), coefficients up to t^(", K, "p), primes p < ", PMAXF);
okF = 1;
{for(i = 1, #FAM, my(f = FAM[i], N = f[7], g = gseq(f[2], f[3], K*PMAXF + 2), tq = f[4](K*PMAXF + 2)[1], qtall,
      okc = 1, okg = 1, okb = 1, notexact = [], badp = []);
  qtall = serreverse(truncate(tq) + O(x^(K*PMAXF + 2)));
  forprime(p = 2, PMAXF,
    my(Lf = K*p, P3 = p^3, qtm, qp, acc, ts, Fs, Ft, Fpp, vF);
    qtm = Mod(1, P3)*(truncate(qtall) + O(x^(Lf+1)));
    qp = qtm^p; acc = qp; ts = 0;
    for(j = 1, K, ts += polcoeff(tq, j)*acc; acc *= qp);
    Fs = sum(j = 0, K, g[j+1]*ts^j);
    Ft = Mod(1, P3)*(sum(n = 0, Lf, g[n+1]*x^n) + O(x^(Lf+1)));
    Fpp = Mod(1, P3)*sum(n = 0, p-1, g[n+1]*x^n);
    vF = minval(lift(Ft - Fpp*Fs), p, Lf);
    if(N % p == 0, badp = concat(badp, [[p, vF]]); next);
    okc = okc && (vF >= 2); if(vF != 2, notexact = concat(notexact, [[p, vF]]));
    okg = okg && ((g[p+1] - g[2]) % p^2 == 0);
    okb = okb && (f[8](p) == g[p+1]));
  okF = okF && okc && okg && okb;
  printf("%-15s p !| N: v_p >= 2 for all p: %d (v_p != 2 only at [p, v]: %s);  g_p == g_1 mod p^2: %d;  primes p | N, [p, v_p]: %s\n",
    f[1], okc, Str(notexact), okg, Str(badp)));}
print("F all: ", okF);

\\ ---------------------------------------------------------------------------------------------- N
\\ Full T_p(x, t) from power sums: sum_k Theta((tau+k)/p)^j = p U_p(Theta^j), Newton's identities for e_j.
Tpoly(tTh, p, r) = {
  my(L = p + 8, M = p*L + 5, S = tTh(M), tq = S[1], Th = S[2] + O(x^(M+1)), qt, ThL, R, Sp = vector(p), e = vector(p+1), T = vector(p+2), Thj = 1, ok = 1);
  qt = serreverse(tq + O(x^(L+1)));
  ThL = Th + O(x^(L+1));
  R = subst(Th + O(x^(L*p+1)), x, x^p)/ThL;
  for(j = 1, p, Thj = Thj*Th; Sp[j] = p*sum(m = 0, L, polcoeff(Thj, p*m)*x^m)/ThL^j + O(x^(L+1)));
  e[1] = 1 + O(x^(L+1));
  for(j = 1, p, e[j+1] = sum(i = 1, j, (-1)^(i-1)*e[j-i+1]*Sp[i])/j);
  for(m = 0, p+1, my(c = if(m <= p, e[m+1]*p^(2*(p-m)), 0) + if(m >= 1, R*e[m]*p^(2*(p-m+1)), 0), P);
    P = truncate((-1)^m*subst(c, x, qt + O(x^(L+1))));
    ok = ok && denominator(content(P)) == 1 && poldegree(P) <= floor(m*r) && L >= floor(m*r) + 3;
    T[m+1] = P);
  [T, ok]
};
newton(v) = {  \\ v[m+1] = valuation of the coefficient of x^m (99999 if zero); returns [[root valuation, multiplicity], ...]
  my(n = #v - 1, pts = [], hull = List(), out = []);
  for(m = 0, n, if(v[m+1] < 99999, pts = concat(pts, [[m, v[m+1]]])));
  foreach(pts, P, while(#hull >= 2 && (hull[#hull][2] - hull[#hull-1][2])*(P[1] - hull[#hull-1][1])
                                       >= (P[2] - hull[#hull-1][2])*(hull[#hull][1] - hull[#hull-1][1]), listpop(hull));
    listput(hull, P));
  for(i = 1, #hull - 1, out = concat(out, [[-(hull[i+1][2] - hull[i][2])/(hull[i+1][1] - hull[i][1]), hull[i+1][1] - hull[i][1]]]));
  out
};
print("\n== N: Newton polygon of T_p(x, t0) at CM values t0 (root valuations [v, multiplicity]); predicted:");
print("       inert [[2p/(p+1), p+1]],  split [[2, p], [0, 1]],  ramified [[(2p-1)/p, p], [1, 1]]");
okN = 1;
{NT = [["s7", tTh7, 2/3, 7, [5, 11, 13], [[-1, -7], [1/2, -12], [1/27, -28]]],
       ["s10", tTh10, 3/4, 10, [7, 11, 13], [[1, -15], [-1/2, -4], [1/16, -40]]],
       ["s18", tTh18, 3/4, 18, [5, 7, 11], [[1/16, -72], [1/12, -36], [-1/4, -180]]],
       ["Apery2*central", tThA5, 1/2, 5, [7, 11, 13], [[1/4, -11], [1/18, -16]]],
       ["Franel*central", tThF6, 1/2, 6, [5, 7, 11], [[-1/4, -3], [1/32, -24]]],
       ["Apery", tThAp, 1, 6, [5, 7, 11, 13], [[1, -8], [-1, -3]]]];}
{foreach(NT, F, foreach(F[5], p, if(F[4] % p == 0, next);
   my(TT = Tpoly(F[2], p, F[3]), line = "");
   okN = okN && TT[2];
   foreach(F[6], c, my(t0 = c[1], D = c[2], k = kronecker(D, p), v, nw, pred, P, fa);
     if(denominator(t0) % p == 0, next);
     v = vector(p+2, i, my(z = subst(TT[1][i], 'x, t0)); if(z == 0, 99999, valuation(z, p)));
     nw = newton(v);
     pred = if(k == -1, [[2*p/(p+1), p+1]], k == 1, [[2, p], [0, 1]], [[(2*p-1)/p, p], [1, 1]]);
     okN = okN && (nw == pred);
     P = sum(m = 0, p+1, subst(TT[1][m+1], 'x, t0)*'y^m); fa = factorpadic(P, p, 20);
     line = concat(line, Str("  t0=", t0, " (D/p)=", k, ": ", nw, " p-adic factor degrees ", apply(poldegree, fa[,1]~), "^", fa[,2]~)));
   printf("%-15s p=%2d  T_{p,m} in Z[t], deg <= floor(m r): %d %s\n", F[1], p, TT[2], line)));}
print("N all (every Newton polygon as predicted): ", okN);

\\ ---------------------------------------------------------------------------------------------- CM
\\ tau = (u + sqrt(v))/w with v < 0 rational; D = discriminant of the primitive form of tau.
disc(u, v, w) = { my(f = [w^2, -2*u*w, u^2 - v]); f = f*denominator(f); f = f/content(f); f[2]^2 - 4*f[1]*f[3] };
classify(g, t0, D, N, AL) = {
  my(nsplit = 0, hit = 0, nin = 0, zin = 0, nram = 0, vram1 = 0, bad = [], badin = []);
  forprime(p = 5, PMAXCM, if(N % p == 0 || denominator(t0) % p == 0, next);
    my(P3 = p^3, v = lift(Mod(sum(n = 0, p-1, g[n+1]*Mod(t0, P3)^n), P3)), k = kronecker(D, p));
    if(k == -1, nin++; if(v % p^2 == 0, zin++, badin = concat(badin, p)); next);
    if(k == 0, nram++; if(v % p == 0 && v % p^2 != 0, vram1++); next);
    nsplit++; my(ok = 0, P2 = p^2);
    foreach(AL, e, my(m = 4*e*p);
      for(Y = 0, sqrtint(m\abs(D)), my(r = m - abs(D)*Y^2);
        if(issquare(r), my(X = sqrtint(r));
          foreach([1, -1], s, if(Mod(v, P2) == Mod(s*(X^2 - 2*e*p), P2)/e, ok = 1)))));
    if(ok, hit++, bad = concat(bad, p)));
  [nsplit, hit, nin, zin, nram, vram1, bad, badin]
};
\\ [name, ODE [A, B], binomial formula, t(tau), N, Atkin-Lehner e, CM list [t0, u, v, w], twisted]
\\ Corrected entries (Beukers arXiv:2403.03301v3 prints): 9.1 t((1+i)/2) = -1/8 (-1/2); 9.2 t((1+3i)/2) = -1/12288
\\ (-1/17088), t(sqrt(-11/2)) = 1/2509056 (1/9801); 9.3 t((1+sqrt(-11))/6) = 1/64 (1/1024); 9.14 t(sqrt(-3/5)) = 1/147
\\ (1/588); 9.15 t((1+sqrt(-7/3))/2) = -1/112 (-1/116).  Twisted 9.7, 9.8: the printed t(alpha) = t0 refers to the
\\ untwisted sequence; for the eta quotient t(q) used here it reads t(alpha + 1/2) = -t0, which is what is listed below.
{CMF = [
 ["9.1 C(2n,n)^3", [[64,96,48,8],[0,0,0,0]], n -> binomial(2*n,n)^3, tau -> (E(1,tau)*E(4,tau)/E(2,tau)^2)^24, 4, [1,4],
  [[1/64,0,-1,2], [-1/64,1,-2,2], [1/256,0,-3,2], [1/16,1,-3,4], [-1/8,1,-1,2], [-1/512,1,-4,2], [1,3,-7,8], [1/4096,0,-7,2]]],
 ["9.2 (4n)!/n!^4", [[256,384,176,24],[0,0,0,0]], n -> (4*n)!/n!^4,
  tau -> my(h = (E(2,tau)/E(1,tau))^24); h/(1+64*h)^2, 2, [1,2],
  [[1/256,0,-1/2,1], [-1/1024,1,-5,2], [1/2304,0,-3/2,1], [-1/144,1,-3,2], [-1/12288,1,-9,2], [1/648,0,-1,1],
   [1/20736,0,-5/2,1], [1/81,1,-7,4], [-1/82944,1,-13,2], [-1/3969,1,-7,2], [1/614656,0,-9/2,1], [1/2509056,0,-11/2,1],
   [-1/6635520,1,-25,2], [-1/199148544,1,-37,2], [1/24591257856,0,-29/2,1]]],
 ["9.3 (3n)!(2n)!/n!^5", [[108,162,78,12],[0,0,0,0]], n -> (3*n)!*(2*n)!/n!^5,
  tau -> my(h = (E(3,tau)/E(1,tau))^12); h/(1+27*h)^2, 3, [1,3],
  [[1/108,0,-1/3,1], [-1/27,3,-15,6], [1/216,0,-2/3,1], [1/8,1,-2,3], [-1/1728,1,-17/3,2], [-1/192,1,-3,2],
   [1/64,1,-11,6], [-1/8640,1,-25/3,2], [1/3375,0,-5/3,1], [-1/110592,1,-41/3,2], [-1/326592,1,-49/3,2]]],
 ["9.5 Apery", [[34,51,27,5],[1,3,3,1]], n -> sum(k=0, n, binomial(n,k)^2*binomial(n+k,k)^2),
  tau -> (E(1,tau)*E(6,tau)/(E(2,tau)*E(3,tau)))^12, 6, [1,6], [[1,2,-2,6], [-1,3,-3,6]]],
 ["9.7 Domb", [[-20,-30,-18,-4],[64,192,192,64]], n -> (-1)^n*sum(k=0, n, binomial(n,k)^2*binomial(2*k,k)*binomial(2*n-2*k,n-k)),
  tau -> (E(2,tau)*E(6,tau)/(E(1,tau)*E(3,tau)))^6, 6, [1,3],
  [[-1,42,-60,48], [1/2,12,-12,12], [-1/4,18,-12,24], [-1/8,8,-8,12], [1/8,4,-8/3,4], [-1/16,6,-12,12], [1/32,4,-16/3,4], [-1/64,30,-60,60]]],
 ["9.8 Almkvist-Zudilin", [[-14,-21,-13,-3],[81,243,243,81]], n -> sum(k=0, n\3, (-3)^(n-3*k)*(n+k)!/(k!^4*(n-3*k)!)),
  tau -> (E(3,tau)*E(6,tau)/(E(1,tau)*E(2,tau)))^4, 6, [1,2],
  [[1,12,-8,12], [-1/3,10,-4,12], [-1/9,18,-12,24], [1/9,4,-8/3,4], [-1/27,2,-4,4], [1/81,4,-8,4]]],
 ["9.9 Almkvist (eps)", [[24,36,20,4],[16,48,48,16]], n -> sum(k=ceil(n/2), n, binomial(n,k)^2*binomial(2*k,n)^2),
  tau -> (E(1,tau)*E(8,tau)/(E(2,tau)*E(4,tau)))^8, 8, [1,8], [[1,5,-7,16], [-1/4,2,-2,4], [1/4,1,-1,4], [1/16,1,-7,8]]],
 ["9.10 Almkvist (zeta)", [[18,27,15,3],[-27,-81,-81,-27]],
  n -> sum(k=0, n, sum(l=0, k, binomial(n,k)^2*binomial(n,l)*binomial(k,l)*binomial(k+l,n))),
  tau -> (E(1,tau)*E(9,tau)/E(3,tau)^2)^6, 9, [1,9],
  [[1,5,-11,18], [-1,4,-2,9], [-1/3,1,-1/3,2], [1/9,1,-3,6], [1/27,0,-2,3], [-1/27,3,-11,6]]],
 ["9.11 Cooper s7", [[26,39,21,4],[-27,-81,-78,-24]], FAM[1][8], t7n, 7, [1,7],
  [[-1,1,-1/7,2], [1/2,2,-3,7], [1/8,3,-19,14], [-1/8,1,-5/7,2], [1/24,1,-27,14], [1/27,0,-1/7,1], [-1/64,1,-13/7,2],
   [1/125,0,-4/7,1], [-1/10648,1,-61/7,2]]],
 ["9.12 Cooper s10", [[12,18,10,2],[-64,-192,-188,-60]], FAM[2][8], t10n, 10, [1,2,5,10],
  [[1,5,-15,20], [-1/2,2,-1,5], [1/4,2,-6,10], [-1/4,1,-1/5,2], [-1/9,1,-3/5,2], [1/12,1,-9,10], [1/16,0,-1/10,1],
   [-1/20,1,-1,2], [1/36,0,-3/10,1], [-1/64,1,-9/5,2], [1/196,0,-7/10,1], [-1/324,1,-17/5,2], [1/1296,0,-13/10,1],
   [1/5776,0,-19/10,1]]],
 ["9.13 Cooper s18", [[28,42,26,6],[192,576,564,180]], FAM[3][8], t18n, 18, [1,2,9,18],
  [[-1/4,3,-5,6], [1/9,3,-7,12], [-1/9,3,-7,6], [1/12,1,-1,6], [1/16,0,-2,6], [1/18,0,-1,3], [1/36,0,-10,6],
   [-1/36,3,-13,6], [1/144,0,-22,6], [-1/180,3,-25,6], [-1/576,3,-37,6]]],
 ["9.14 Apery2*central", [[44,66,34,6],[-16,-48,-44,-12]], FAM[4][8], tA5n, 5, [1,5],
  [[1/4,3,-11,10], [-1/3,1,-3/5,2], [-1/28,1,-7/5,2], [1/36,1,-19,10], [1/18,-1,-4,5], [1/72,0,-2/5,1],
   [1/147,0,-3/5,1], [-1/828,1,-23/5,2], [-1/15228,1,-47/5,2]]],
 ["9.15 Franel*central", [[28,42,22,4],[-128,-384,-352,-96]], FAM[5][8], tF6n, 6, [1,2,3,6],
  [[-1/4,1,-1/3,2], [1/16,1,-5,6], [-1/16,1,-1,2], [1/5,3,-15,12], [1/32,0,-1/6,1], [1/96,0,-1/2,1], [1/50,0,-1/3,1],
   [-1/112,1,-7/3,2], [-1/49,1,-5/3,2], [1/320,0,-5/6,1], [-1/400,1,-11/3,2], [1/896,0,-7/6,1], [-1/2704,1,-19/3,2],
   [1/10400,0,-13/6,1], [1/39200,0,-17/6,1]]]
];}
print("\n== CM: rational CM values t0 = t(alpha); primes 5 <= p < ", PMAXCM, ", p !| N den(t0)");
print("   [inert primes, F_p(t0) == 0 mod p^2] [split primes, F_p(t0) == s(X^2 - 2ep)/e mod p^2 with 4ep = X^2 + |D|Y^2, e | N Atkin-Lehner]");
print("   [ramified primes, v_p(F_p(t0)) = 1]");
okCMin = 1; okCMsp = 1; okCMnum = 1; nvals = 0; ninert = 0; nsplit = 0; nram = 0; nram1 = 0;
{foreach(CMF, F, my(g = gseq(F[2][1], F[2][2], PMAXCM), okf = (vector(41, n, F[3](n-1)) == g[1..41]));
  printf("%s (level %d, e in %s; binomial formula = ODE for n <= 40: %d)\n", F[1], F[5], Str(F[6]), okf);
  okCMnum = okCMnum && okf;
  foreach(F[7], c, my(t0 = c[1], tau = (c[2] + sqrt(c[3]))/c[4], D = disc(c[2], c[3], c[4]), tn = F[4](tau), r, numok);
    numok = abs(tn - t0) < 10^-40*max(1, abs(t0));
    r = classify(g, t0, D, F[5], F[6]);
    okCMnum = okCMnum && numok; okCMin = okCMin && (r[4] == r[3]); okCMsp = okCMsp && (r[2] == r[1]);
    nvals++; ninert += r[3]; nsplit += r[1]; nram += r[5]; nram1 += r[6];
    printf("  t0 = %-14s D = %-6d |t(alpha) - t0| < 1e-40: %d   inert %2d zero %2d%s   split %2d matched %2d%s   ramified %d v=1 %d\n",
      Str(t0), D, numok, r[3], r[4], if(#r[8], Str(" FAIL ", r[8]), ""), r[1], r[2], if(#r[7], Str(" FAIL ", r[7]), ""), r[5], r[6])));}
{printf("CM: %d values; t(alpha) = t0 and formulas: %d;  inert (%d prime-value pairs) all zero mod p^2: %d;  split (%d pairs) all of the predicted shape: %d;  ramified: v_p = 1 at %d of %d\n",
  nvals, okCMnum, ninert, okCMin, nsplit, okCMsp, nram1, nram);}
\\ Z.-W. Sun's conjecture (Guo-Zeng, JNT 132 (2012), Conjecture 1.2) for the Apery numbers A_k at t0 = 1
print("\n== S: sum_{k<p} A_k == 4x^2 - 2p (mod p^2) if p = x^2 + 2y^2, == 0 (mod p^2) if p == 5, 7 (mod 8); 3 <= p < ", PMAXCM);
okS = 1; nS = [0, 0];
{my(av = gseq([34,51,27,5], [1,3,3,1], PMAXCM));
 forprime(p = 3, PMAXCM, my(S = sum(k = 0, p-1, Mod(av[k+1], p^2)));
   if(p % 8 == 1 || p % 8 == 3, my(x0 = qfbsolve(Qfb(1, 0, 2), p)[1]); okS = okS && (S == 4*x0^2 - 2*p); nS[1]++,
      okS = okS && (S == 0); nS[2]++));}
print("S: ", okS, "  (", nS[1], " primes p = 1, 3 mod 8 and ", nS[2], " primes p = 5, 7 mod 8)");
print("\nALL (I, T, F, N, CM, S): ", okI && okT && okF && okN && okCMnum && okCMin && okCMsp && nvals == 118 && okS);
quit;
