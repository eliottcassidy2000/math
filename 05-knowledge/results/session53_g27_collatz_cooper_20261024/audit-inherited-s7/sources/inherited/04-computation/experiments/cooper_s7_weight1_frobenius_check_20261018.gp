\\ Cooper's sporadic sequence s7 and its weight-one square root c7 (OEIS A279619).
\\ Checks for 05-knowledge/results/cooper_s7_weight_one_root_frobenius_20261018.md
\\ (session cloud-2026-10-18-M1).
\\
\\   s7(n) = sum_k C(n,k)^2 C(n+k,k) C(2k,n),  F(t) = sum s7(n) t^n,  G(t) = F(t)^(1/2) = sum c7(n) t^n,
\\   h = (eta7/eta1)^4,  t = h/(1+13h+49h^2),  s = (1-49h^2)/(1+13h+49h^2),
\\   theta = 1 + 2 sum (k/7) q^k/(1-q^k)  (theta series of x^2+xy+2y^2).
\\
\\ Sections: S (exact q-series identities), T (Beukers' weight-one polynomial T_p for theta),
\\ A/B (Frobenius congruences with the excellent lift t^sigma(q) = t(q^p)), C (Lucas congruences for c7),
\\ CM (mod p^2 values at the nine rational CM points), O (the Chan-Cooper-Sica congruence
\\ c7(np) == c7(n) mod p^2, numerically), D (the ingredients of its proof, including p = 2 and p = 7),
\\ E (inert primes: the twisted congruence d7(pn) == c7(n) mod p^2 for H = G/s = sum d7(n) t^n, and its proof).
\\
\\ Run: gp -q 04-computation/experiments/cooper_s7_weight1_frobenius_check_20261018.gp

default(parisizemax, 8*10^9);

PMAXA = 110;     \\ sections A/B: primes p < PMAXA
K = 4;           \\ sections A/B: compare coefficients of t^0 .. t^(K p)
M = K*PMAXA;     \\ q-adic precision for sections S, A, B
MT = 1500;       \\ q-adic precision for section T
PMAXT = 43;      \\ section T: primes p <= PMAXT
PMAXC = 60;      \\ sections C, O: primes p < PMAXC
NC = 3000;       \\ sections C, O: indices n < NC
PMAXCM = 400;    \\ section CM: primes p < PMAXCM

s7(n) = sum(k=0, n, binomial(n,k)^2*binomial(n+k,k)*binomial(2*k,n));

c7list(N) = {
  my(c = vector(N+1), num);
  c[1] = 1; c[2] = 2;
  for(n = 1, N-1,
    num = (26*n^2 + 13*n + 2)*c[n+1] + 3*(3*n-1)*(3*n-2)*c[n];
    if(num % (n+1)^2, error("c7 not integral at n = ", n+1));
    c[n+2] = num/(n+1)^2);
  c
};

chi(p) = kronecker(-7, p);

\\ ------------------------------------------------------------------ S
qseries(P) = {
  my(E, h, tq, sq, th);
  E = (m -> subst(eta(x + O(x^(P\m + 1))), x, x^m));           \\ prod_n (1 - q^(m n))
  h = x*(E(7)/E(1))^4;
  tq = h/(1 + 13*h + 49*h^2);
  sq = (1 - 49*h^2)/(1 + 13*h + 49*h^2);
  th = 1 + 2*sum(n=1, P, sumdiv(n, d, kronecker(d, 7))*x^n) + O(x^(P+1));
  [tq, sq, th, x*(E(1)*E(7))^3]
};

print("== S: exact identities (q-precision ", M, ")");
QS = qseries(M); tq = QS[1]; sq = QS[2]; th = QS[3]; e17 = QS[4];
S = vector(M+1, n, s7(n-1));
C = c7list(M);
{print("S1 c7 * c7 = s7 (Cauchy square) for n <= ", M, ": ",
      vector(M+1, n, sum(j=0, n-1, C[j+1]*C[n-j])) == S);}
Fy = sum(n=0, M, S[n+1]*x^n) + O(x^(M+1));      \\ F as a series in t (variable x)
Gy = sum(n=0, M, C[n+1]*x^n) + O(x^(M+1));      \\ G as a series in t
sy = sqrt(1 - 26*x - 27*x^2 + O(x^(M+1)));      \\ s as a series in t
print("S2 F(t(q)) = theta^2: ", valuation(subst(truncate(Fy), x, tq) - th^2, x) > M);
print("S3 G(t(q)) = theta: ", valuation(subst(truncate(Gy), x, tq) - th, x) > M);
print("S4 t theta^3 = (eta1 eta7)^3: ", valuation(tq*th^3 - e17, x) > M);
print("S5 s^2 = 1 - 26 t - 27 t^2: ", valuation(sq^2 - (1 - 26*tq - 27*tq^2), x) > M);
print("S6 q dt/dq / t = theta^2 s: ", valuation(x*deriv(tq, x)/tq - th^2*sq, x) >= M);
print("S7 s(t) = sqrt(1-26t-27t^2) has integral t-coefficients: ", denominator(content(truncate(sy))) == 1);
print("S8 s(t(q)) = s(q): ", valuation(subst(truncate(sy), x, tq) - sq, x) > M);
{print("S9 theta^6 (1+t)^2 = theta^6 + 2 (eta1 eta7)^3 theta^3 + (eta1 eta7)^6 (holomorphic, weight 6): ",
      valuation(th^6*(1 + tq)^2 - (th^6 + 2*e17*th^3 + e17^2), x) > M);}

\\ ------------------------------------------------------------------ T
\\ T_p(x) = (1 - R x) prod_k (p - chi(p) r_k x),  R = theta(q^p)/theta(q),  r_k = theta(zeta^k q^(1/p))/theta(q).
\\ a = coefficient of x^p, b = coefficient of x^(p+1), as q-series to precision L = MT\p.
print("\n== T: Beukers' weight-one polynomial for theta (q-precision ", MT, ")");
QT = qseries(MT); thT = QT[3]; tqT = QT[1];
logthT = log(thT);
invthT = 1/thT;
Tstruct(p) = {
  my(L = MT\p, c = chi(p), thL, N, Tr, R, a, b, rel, qt, at, bt, st, A0, okA, okB, degA, degB, boundA, boundB);
  thL = thT + O(x^(L+1));
  N = exp(p*sum(m=1, L, polcoeff(logthT, p*m)*x^m) + O(x^(L+1)));    \\ prod_k theta(zeta^k q^(1/p))
  Tr = p*sum(m=0, L, polcoeff(invthT, p*m)*x^m) + O(x^(L+1));        \\ sum_k 1/theta(zeta^k q^(1/p))
  R = subst(thL, x, x^p)/thL;
  b = -(-c)^p*R*N/thL^p;
  a = (-c)^p*N/thL^p - R*p*(-c)^(p-1)*(N/thL^p)*thL*Tr;
  rel = lift(Mod(1, p^2)*truncate(a + b/R));
  qt = serreverse(tqT + O(x^(L+1)));
  at = subst(truncate(a), x, qt);
  bt = subst(truncate(b), x, qt);
  st = sqrt(1 - 26*x - 27*x^2 + O(x^(L+1)));
  if(c == 1, A0 = at; boundA = p\3, A0 = at/st; boundA = p\3 - 1);
  degA = poldegree(truncate(A0)); degB = poldegree(truncate(bt)); boundB = (p+1)\3;
  okA = (denominator(content(truncate(A0))) == 1) && (degA <= boundA);
  okB = (denominator(content(truncate(bt))) == 1) && (degB <= boundB) && ((truncate(bt) - c) % p == 0);
  printf("p=%3d chi=%2d L=%3d  a%s: integral polynomial of degree %2d <= %2d: %d   b: integral, degree %2d <= %2d, b == chi mod p: %d   a + b theta/theta(q^p) == 0 mod p^2: %d\n",
         p, c, L, if(c == 1, "  ", "/s"), degA, boundA, okA, degB, boundB, okB, rel == 0);
  okA && okB && (rel == 0)
};
okT = 1;
forprime(p = 2, PMAXT, if(p != 7, okT = okT && Tstruct(p)));
print("T all: ", okT);

\\ ------------------------------------------------------------------ A, B
\\ t^sigma(t) = t(q(t)^p) modulo p^3.  Theorem A: F == F_p(t) F(t^sigma) mod p^2.
\\ Theorem B: split p: G == G_p(t) G(t^sigma);  inert p: H := G/s satisfies H == H_p(t) G(t^sigma) mod p^2
\\ and H_p mod p^2 has degree <= floor(p/3) - 1 + floor((p+1)/3).
print("\n== A/B: Frobenius congruences, primes p < ", PMAXA, ", p != 7, coefficients up to t^(", K, " p)");
qtA = serreverse(tq);
Hy = Gy/sy;
minval(d, p, L) = {
  my(v = 99, cc);
  for(k = 0, L, cc = polcoeff(d, k); if(cc != 0, v = min(v, valuation(cc, p))));
  v
};
frob(p) = {
  my(L = K*p, P3 = p^3, qtm, qp, acc, ts, Fs, Gs, Ft, Fp, Gt, Gp, Ht, Hp, vF, vG, vN, deghp, boundH, res);
  qtm = Mod(1, P3)*(truncate(qtA) + O(x^(L+1)));
  qp = qtm^p; acc = qp; ts = 0;
  for(j = 1, L\p, ts += polcoeff(tq, j)*acc; acc *= qp);
  Fs = sum(j=0, L\p, S[j+1]*ts^j);
  Gs = sum(j=0, L\p, C[j+1]*ts^j);
  Ft = Mod(1, P3)*(truncate(Fy) + O(x^(L+1)));
  Gt = Mod(1, P3)*(truncate(Gy) + O(x^(L+1)));
  Fp = Mod(1, P3)*sum(n=0, p-1, S[n+1]*x^n);
  vF = minval(lift(Ft - Fp*Fs), p, L);
  \\ naive lift t -> t^p, for comparison
  vN = minval(lift(Ft - Fp*subst(truncate(Ft), x, x^p)), p, L);
  if(chi(p) == 1,
    Gp = Mod(1, P3)*sum(n=0, p-1, C[n+1]*x^n);
    vG = minval(lift(Gt - Gp*Gs), p, L);
    deghp = -1; boundH = -1,
    Ht = Mod(1, P3)*(truncate(Hy) + O(x^(L+1)));
    Hp = Mod(1, P3)*sum(n=0, p-1, polcoeff(Hy, n)*x^n);
    vG = minval(lift(Ht - Hp*Gs), p, L);
    boundH = p\3 - 1 + (p+1)\3;
    deghp = poldegree(lift(Mod(1, p^2)*lift(Hp))));
  res = [p, chi(p), vF, vN, vG, deghp, boundH];
  res
};
resA = [];
forprime(p = 2, PMAXA - 1, if(p != 7, resA = concat(resA, [frob(p)])));
print("[p, chi(p), v_p(F - F_p F^sigma), v_p(F - F_p F(t^p)), v_p(G-side), deg(H_p mod p^2), bound] (99 = exact zero):");
for(i = 1, #resA, print("  ", resA[i]));
{print("A: v(F - F_p F^sigma) >= 2 for all p: ", #select((r -> r[3] < 2), resA) == 0,
      ";  = 2 exactly except at: ", apply((r -> r[1]), select((r -> r[3] != 2), resA)));}
print("A': naive lift t^p gives only v = 1 at: ", apply((r -> r[1]), select((r -> r[4] == 1), resA)));
{print("B: v(G-side) >= 2 for all p: ", #select((r -> r[5] < 2), resA) == 0,
      ";  inert degree bound holds: ", #select((r -> r[2] == -1 && r[6] > r[7]), resA) == 0);}

\\ ------------------------------------------------------------------ C
print("\n== C: Lucas congruences for c7 modulo p, n < ", NC);
CC = c7list(NC);
lucas(p) = {
  my(m, pr);
  for(n = 0, NC - 1,
    m = n; pr = 1;
    while(m, pr = pr*CC[m % p + 1] % p; m \= p);
    if((CC[n+1] - pr) % p, return(n)));
  -1
};
lucasall() = {
  my(f, expect, ok = 1);
  forprime(p = 2, PMAXC - 1,
    f = lucas(p); expect = (p == 7) || (chi(p) == 1);
    printf("  p=%2d chi=%2d  Lucas holds: %d (first failure n = %d)  expected: %d\n", p, chi(p), f == -1, f, expect);
    ok = ok && ((f == -1) == expect));
  ok
};
okC = lucasall();
print("C: Lucas holds exactly for p = 7 and chi(p) = 1: ", okC);
w6 = th^6*(1 + tq)^2;
{print("C7 theta^6 (1+t)^2 == 1 mod 7 up to q^", M, ": ", lift(Mod(1,7)*truncate(w6 - 1)) == 0,
      ";  first coefficients (Sturm bound 4): ", vector(5, k, polcoeff(w6, k-1)));}
{print("C7 G(t) == (1+t)^2 G(t^7) mod 7 up to t^", M, ": ",
      lift(Mod(1,7)*truncate(Gy - (1 + x)^2*subst(truncate(Gy), x, x^7))) == 0);}

\\ ------------------------------------------------------------------ CM
\\ Rule: if 4p = X^2 - D Y^2 is solvable then F_p(t0) == X^2 - 2p, else 28p = X^2 - D Y^2 and F_p(t0) == 2p - X^2/7 (mod p^2).
print("\n== CM: F_p(t0) mod p^2 at the nine rational CM values, split primes 5 <= p < ", PMAXCM);
SCM = vector(PMAXCM, n, s7(n-1));
repX(D, m) = {
  my(y = 0, r, X);
  while(-D*y^2 <= 4*m,
    r = 4*m + D*y^2;
    if(issquare(r, &X) && (X - D*y) % 2 == 0, return(X));
    y++);
  -1
};
cmcheck(t0, D) = {
  my(good = 0, bad = [], P2, X1, X7, val, pred, n1 = 0, n7 = 0);
  forprime(p = 5, PMAXCM - 1,
    if(p == 7 || denominator(t0) % p == 0 || D % p == 0 || kronecker(D, p) != 1, next);
    P2 = p^2;
    val = lift(sum(n = 0, p-1, SCM[n+1]*Mod(t0, P2)^n));
    X1 = repX(D, p);
    if(X1 >= 0, pred = lift(Mod(X1^2 - 2*p, P2)); n1++,
      X7 = repX(D, 7*p);
      if(X7 < 0, error("no representation"));
      pred = lift(Mod(2*p, P2) - Mod(X7^2, P2)/7); n7++);
    if(val == pred, good++, bad = concat(bad, p)));
  printf("  t0 = %9s  D = %5d: match %3d / %3d   (N(gamma) = p: %3d, N(gamma) = 7p: %3d)  failures: %s\n",
         Str(t0), D, good, good + #bad, n1, n7, Str(bad));
  #bad == 0
};
cmlist = [[-1, -7], [1/2, -12], [1/8, -19], [-1/8, -35], [1/24, -27], [1/27, -28], [-1/64, -91], [1/125, -112], [-1/10648, -427]];
okCM = 1;
for(i = 1, #cmlist, okCM = okCM && cmcheck(cmlist[i][1], cmlist[i][2]));
print("CM all: ", okCM);

\\ ------------------------------------------------------------------ O
print("\n== O: the Chan-Cooper-Sica congruence c7(np) == c7(n) mod p^2, numerically (n < ", NC, "/p)");
ccsall() = {
  my(ok1, ok2, ok3);
  forprime(p = 2, PMAXC - 1,
    ok1 = 1; ok2 = 1; ok3 = 1;
    for(n = 1, (NC - 1)\p,
      if((CC[n*p+1] - CC[n+1]) % p, ok1 = 0);
      if((CC[n*p+1] - CC[n+1]) % p^2, ok2 = 0);
      if((CC[n*p+1] - CC[n+1]) % p^3, ok3 = 0));
    printf("  p=%2d chi=%2d  mod p: %d  mod p^2: %d  mod p^3: %d\n", p, chi(p), ok1, ok2, ok3));
};
ccsall();

\\ ------------------------------------------------------------------ D
\\ Proof ingredients of the Chan-Cooper-Sica congruence c7(pn) == c7(n) mod p^2.
\\ Split p >= 11:
\\   (E)  prod(1-q^(pn)) / prod(1-q^n)^p == E_{p-1}^(-1/2) mod p^2,
\\   (F)  p eps_p == 3 ( theta^(p-1) G_p(t) - (E_{p-1}(tau) + E_{p-1}(7 tau))/2 ) mod p^2   (q-series),
\\   (W)  Psi_p = (E_{p-1}(tau) + E_{p-1}(7 tau))/(2 theta^(p-1)) is == A_p(t) mod p^2, A_p a polynomial of degree <= (p-1)\3,
\\   (P)  eps_p G_p == 3 (G_p - A_p)/p mod p, a polynomial of degree <= p\3 + (p+1)\3 with constant term 0.
\\ p = 2: t(q^2) == t^2 (1 + 2t) mod 4.  p = 7: f = theta^3 s = 1 - (7/8) A1 - (49/8) A2,  U_7 f == f mod 49,
\\ t(q^7) == t^7 (1 + 14 t + 28 t^2) mod 49  (both t-identities via holomorphic forms and Sturm bounds).
print("\n== D: proof ingredients for c7(pn) == c7(n) mod p^2");
Pq = eta(x + O(x^(M+1)));
EisD(k) = 1 - (2*k/bernfrac(k))*sum(n=1, M, sigma(n, k-1)*x^n) + O(x^(M+1));
iszeromod(f, m) = { my(g = truncate(f)); for(i = 0, poldegree(g), if(polcoeff(g, i) % m, return(0))); 1 };
ccsproof(p) = {
  my(P2 = p^2, k = p-1, E1, E7, okE, eps, Gpt, okF, Psi, At, degA, okW, epsT, GpT, lhs, okP, degP, c0, prec = M - 10);
  E1 = EisD(k); E7 = subst(E1, x, x^7);
  okE = iszeromod((subst(Pq, x, x^p)/Pq^p)^2*E1 - 1, P2);
  eps = (subst(tq, x, x^p)/tq^p - 1)/p;
  Gpt = sum(n=0, p-1, C[n+1]*tq^n);
  okF = iszeromod(p*eps - 3*(th^k*Gpt - (E1 + E7)/2), P2);
  Psi = (E1 + E7)/(2*th^k);
  At = lift(Mod(1, P2)*subst(truncate(Psi), x, qtA)) + O(x^prec);
  degA = poldegree(truncate(At));
  okW = (degA <= (p-1)\3);
  epsT = Mod(1, p)*subst(truncate(eps), x, qtA) + O(x^prec);
  GpT = sum(n=0, p-1, C[n+1]*x^n);
  lhs = epsT*GpT;
  okP = iszeromod(GpT - truncate(At), p) && (lift(truncate(lhs - Mod(1,p)*3*((GpT - truncate(At))/p))) == 0);
  degP = poldegree(lift(truncate(lhs))); c0 = lift(polcoeff(lhs, 0));
  printf("  p=%2d: (E) %d  (F) %d  (W) deg A_p = %2d <= %2d: %d  (P) %d, deg(eps_p G_p mod p) = %2d <= %2d, constant term %d\n",
         p, okE, okF, degA, (p-1)\3, okW, okP, degP, p\3 + (p+1)\3, c0);
  okE && okF && okW && okP && (degP <= p\3 + (p+1)\3) && (c0 == 0)
};
okD = 1;
forprime(p = 11, 43, if(chi(p) == 1, okD = okD && ccsproof(p)));
print("D split primes 11..43 all ingredients: ", okD);
A1 = sum(n=1, M, sumdiv(n, d, kronecker(-7, d)*d^2)*x^n) + O(x^(M+1));
A2 = sum(n=1, M, sumdiv(n, d, kronecker(-7, n/d)*d^2)*x^n) + O(x^(M+1));
fw3 = th^3*sq;
print("D7 f = theta^3 s = 1 - (7/8) A1 - (49/8) A2 (weight 3, Gamma0(7), chi_-7) to q^", M, ": ", valuation(fw3 - (1 - 7/8*A1 - 49/8*A2), x) > M);
{print("D7 U_7 f == f mod 49 to q^", M\7, ": ", #select((n -> (polcoeff(fw3, 7*n) - polcoeff(fw3, n)) % 49 != 0), vector(M\7 + 1, n, n-1)) == 0);}
th7 = subst(th, x, x^7); e17_7 = subst(e17, x, x^7);
F30 = th^27*e17_7 - e17^7*th^6*th7^3 - 14*e17^8*th^3*th7^3 - 28*e17^9*th7^3;
{print("D7 theta^27 theta(q^7)^3 (t(q^7) - t^7 - 14 t^8 - 28 t^9) == 0 mod 49 to q^", M, " (weight 30 on Gamma0(49), Sturm bound 140): ",
      iszeromod(F30, 49));}
th2 = subst(th, x, x^2); e17_2 = subst(e17, x, x^2);
F12 = th^9*e17_2 - e17^2*th^3*th2^3 - 2*e17^3*th2^3;
{print("D2 theta^9 theta(q^2)^3 (t(q^2) - t^2 - 2 t^3) == 0 mod 4 to q^", M, " (weight 12 on Gamma0(14), Sturm bound 24): ",
      iszeromod(F12, 4));}
eps7T = lift(Mod(1, 7)*subst(truncate((subst(tq, x, x^7)/tq^7 - 1)/7), x, qtA)) + O(x^(M-10));
print("D7 eps_7 mod 7 = ", truncate(eps7T), ";  eps_7 (1+t)^2 mod 7 = ", lift(Mod(1,7)*truncate(eps7T*(1+x)^2)));
eps2T = lift(Mod(1, 2)*subst(truncate((subst(tq, x, x^2)/tq^2 - 1)/2), x, qtA)) + O(x^(M-10));
print("D2 eps_2 mod 2 = ", truncate(eps2T), ";  G_2 = ", C[1] + C[2]*x);

\\ ------------------------------------------------------------------ E
\\ Inert primes.  H = G/s = sum d7(n) t^n, 1/s = sum u(n) t^n with u(n) = [x^n] (x^2 + 13x + 49)^n.
\\ Theorem E: d7(pn) == c7(n) mod p^2 for every inert p.  Ingredients for inert p >= 5:
\\   (F') p eps_p == 3 ( theta^(p-1) s H_p(t) - (E_{p-1}(tau) + E_{p-1}(7 tau))/2 ) mod p^2,
\\   (W') Psi_p == s B_p(t) mod p^2, B_p a polynomial of degree <= (p-1)\3 - 1,
\\   (P') eps_p H_p == 3 (H_p - B_p)/p mod p, a polynomial of degree <= p\3 - 1 + (p+1)\3 with constant term 0.
\\ p = 3: t(q^3) == t^3 mod 9 (weight 12 on Gamma0(21), Sturm bound 32), so eps_3 == 0 mod 3.
print("\n== E: inert primes, twisted congruence d7(pn) == c7(n) mod p^2, where H = G/s = sum d7(n) t^n");
ulist(N) = {
  my(u = vector(N+1), num);
  u[1] = 1; u[2] = 13;
  for(n = 1, N-1,
    num = (26*n + 13)*u[n+1] + 27*n*u[n];
    if(num % (n+1), error("u not integral at n = ", n+1));
    u[n+2] = num/(n+1));
  u
};
UU = ulist(NC);
{print("E0 1/s = sum u(n) t^n with u(n) = [x^n] (x^2+13x+49)^n (n <= 60) and the recurrence (n <= ", M, "): ",
      vector(61, n, polcoeff((x^2 + 13*x + 49)^(n-1), n-1)) == UU[1..61],
      " ", vector(M+1, n, polcoeff(1/sy, n-1)) == UU[1..M+1]);}
twistnum(p) = {
  my(P3 = p^3, cm, um, d, diff, ok1 = 1, ok2 = 1, ok3 = 1);
  cm = vector(NC, i, Mod(CC[i], P3)); um = vector(NC, i, Mod(UU[i], P3));
  for(n = 1, (NC - 1)\p,
    d = sum(k = 0, n*p, cm[k+1]*um[n*p-k+1]);
    diff = lift(d - cm[n+1]);
    if(diff % p, ok1 = 0); if(diff % p^2, ok2 = 0); if(diff, ok3 = 0));
  [ok1, ok2, ok3]
};
twistall() = {
  my(r, ok = 1);
  forprime(p = 2, PMAXC - 1,
    r = twistnum(p);
    printf("  p=%2d chi=%2d  d7(pn) == c7(n):  mod p: %d  mod p^2: %d  mod p^3: %d\n", p, chi(p), r[1], r[2], r[3]);
    if(chi(p) == -1, ok = ok && r[2]));
  ok
};
okEn = twistall();
print("E numerics: mod p^2 holds for every inert p < ", PMAXC, ": ", okEn);
Dl = vector(M+1, n, polcoeff(Hy, n-1));
twistproof(p) = {
  my(P2 = p^2, k = p-1, E1, E7, okE, eps, Hpt, okF, Psi, Bt, degB, okW, epsT, HpT, lhs, okP, degP, c0, prec = M - 10);
  E1 = EisD(k); E7 = subst(E1, x, x^7);
  okE = iszeromod((subst(Pq, x, x^p)/Pq^p)^2*E1 - 1, P2);
  eps = (subst(tq, x, x^p)/tq^p - 1)/p;
  Hpt = sum(n=0, p-1, Dl[n+1]*tq^n);
  okF = iszeromod(p*eps - 3*(th^k*sq*Hpt - (E1 + E7)/2), P2);
  Psi = (E1 + E7)/(2*th^k);
  Bt = lift(Mod(1, P2)*(subst(truncate(Psi), x, qtA)/sy)) + O(x^prec);
  degB = poldegree(truncate(Bt));
  okW = (degB <= (p-1)\3 - 1);
  epsT = Mod(1, p)*subst(truncate(eps), x, qtA) + O(x^prec);
  HpT = sum(n=0, p-1, Dl[n+1]*x^n);
  lhs = epsT*HpT;
  okP = iszeromod(HpT - truncate(Bt), p) && (lift(truncate(lhs - Mod(1,p)*3*((HpT - truncate(Bt))/p))) == 0);
  degP = poldegree(lift(truncate(lhs))); c0 = lift(polcoeff(lhs, 0));
  printf("  p=%2d: (E) %d  (F') %d  (W') deg B_p = %2d <= %2d: %d  (P') %d, deg(eps_p H_p mod p) = %2d <= %2d, constant term %d\n",
         p, okE, okF, degB, (p-1)\3 - 1, okW, okP, degP, p\3 - 1 + (p+1)\3, c0);
  okE && okF && okW && okP && (degP <= p\3 - 1 + (p+1)\3) && (c0 == 0)
};
okEp = 1;
forprime(p = 5, 43, if(chi(p) == -1, okEp = okEp && twistproof(p)));
print("E inert primes 5..41 all ingredients: ", okEp);
th3 = subst(th, x, x^3); e17_3 = subst(e17, x, x^3);
{print("E3 theta^9 theta(q^3)^3 (t(q^3) - t^3) == 0 mod 9 to q^", M, " (weight 12 on Gamma0(21), Sturm bound 32): ",
      iszeromod(th^9*e17_3 - e17^3*th3^3, 9));}
quit;
