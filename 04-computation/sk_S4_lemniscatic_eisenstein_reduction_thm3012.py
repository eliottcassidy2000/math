"""THM-3012 addendum 5 referee: the EXACT lemniscatic reduction of S(4).

Claim chain (each step checked here):

 (A) M_{2m} = int_0^1 k^{2m} K(k) dk obeys, EXACTLY,
        4 m^2 M_{2m} = (2m-1)^2 M_{2m-2} + 1,   M_0 = 2G,
     from the Legendre ODE (k(1-k^2)K')' = kK plus two integrations by parts.
     Hence M_{2m} = 2 c_m G + b_m,  c_m = [(1/2)_m/m!]^2,  b_m rational,
        4 m^2 b_m = (2m-1)^2 b_{m-1} + 1,  b_0 = 0.

 (B) 1/(2-k^2) = (1/2) sum_m (k^2/2)^m  gives
        I := int_0^1 K(k)/(2-k^2) dk = (2/pi) K(1/sqrt2) G + beta,
        beta = (1/2) sum_m b_m 2^{-m}.

 (C) B(x) = sum_m b_m x^m solves the INHOMOGENEOUS hypergeometric equation
        x(1-x)B'' + (1-2x)B' - B/4 = 1/(4(1-x)),
     whose homogeneous solutions are (2/pi)K(sqrt x), (2/pi)K(sqrt(1-x)) and
     whose Wronskian is exactly the LEGENDRE RELATION, W = -1/(pi x(1-x)).
     Variation of parameters (f/W = -pi/(4(1-x)) is elementary) gives
        B(x) = (pi/4)[ y1 int_0^x y2 dt/(1-t) - y2 int_0^x y1 dt/(1-t) ],
     and at x = 1/2 (where y1 = y2) the k -> k' involution collapses it to
        beta = (K(1/sqrt2)/pi) * D,   D = int_{1/sqrt2}^1 [K(k) - K'(k)] dk/k.

 (D) Therefore, with K = K(1/sqrt2) = Gamma(1/4)^2/(4 sqrt(pi)),
        S(4) = (2 sqrt2/pi) I = (2 sqrt2 K/pi^2) * Lambda,  Lambda = 2G + D,
        S(4) = sqrt2 Gamma(1/4)^2 Lambda / (2 pi^{5/2}).

 (E) In the modular coordinate x = lambda(tau), tau = i K'/K, the point x=1/2 is
     tau = i and  K dk/k = (i pi^2/4) theta3^2 theta4^4 dtau, so
        D = (pi^2/4) int_1^infty (s-1) theta2^4 theta3^2 (i s) ds.
     theta2^4 theta3^2 = 16 sum_n B_n q^n with B_n = sum_{d|n} chi_{-4}(n/d) d^2
     (weight-3 EISENSTEIN, no cusp-form part), so
        D = 4 sum_{n odd} chi_{-4}(n) / (n^2 (e^{pi n} - 1)) = 4 sum_{m>=1} Ti_2(e^{-pi m}),
        Lambda = 2 sum_{n odd} (-1)^{(n-1)/2} coth(n pi/2)/n^2.

 (F) Residues of  pi cot(pi z) sech(pi z)/z^2  give the single Fricke relation
        Lambda = 5 pi^2/24 - (1/2) sum_{m>=1} sech(pi m)/m^2.

Run:  python3 sk_S4_lemniscatic_eisenstein_reduction_thm3012.py
"""
import mpmath as mp
import itertools
import time
from fractions import Fraction

FAIL = []


def check(name, cond, extra=""):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name} {extra}")
    if not cond:
        FAIL.append(name)


# --------------------------------------------------------------- C1: S(4) values
print("C1  S(4) by three independent routes")
mp.mp.dps = 60
t0 = time.time()
# the defining 3F2 at unit argument is the expensive one -- 40 digits is plenty to
# pin the two integral routes to it; C5 then carries the pair to 200 digits.
with mp.workdps(40):
    S4_hyper = +mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/4], [1, mp.mpf(5)/4], 1)
S4_ell = mp.quad(lambda k: mp.ellipk(k**2)/(2-k**2), [0, 1])*2*mp.sqrt(2)/mp.pi
S4_smooth = mp.quad(lambda p: (p + mp.asinh(mp.sin(p)))/mp.sqrt(1+mp.sin(p)**2),
                    [0, mp.pi/2])*2/mp.pi
print("     3F2 (40 digits) =", mp.nstr(S4_hyper, 38))
check("elliptic route  == 3F2 (35 digits)", abs(S4_ell - S4_hyper) < mp.mpf(10)**-35)
check("smooth  route  == 3F2 (35 digits)", abs(S4_smooth - S4_hyper) < mp.mpf(10)**-35,
      f"[{time.time()-t0:.0f}s]")
check("elliptic route == smooth route (55 digits)",
      abs(S4_ell - S4_smooth) < mp.mpf(10)**-55)
S4_hyper = S4_smooth   # high-precision reference for the later checks

# ------------------------------------------------------- C2: the moment recurrence
print("C2  exact K-moment recurrence  4m^2 M_2m = (2m-1)^2 M_2m-2 + 1")
G = mp.catalan
b = Fraction(0); c = Fraction(1); okc = True
for m in range(1, 8):
    b = (Fraction((2*m-1)**2)*b + 1)/Fraction(4*m*m)
    c = c*Fraction((2*m-1)**2, 4*m*m)
    pred = 2*mp.mpf(c.numerator)/c.denominator*G + mp.mpf(b.numerator)/b.denominator
    num = mp.quad(lambda k: mp.ellipk(k**2)*k**(2*m), [0, 1])
    okc &= abs(pred - num) < mp.mpf(10)**-50
check("M_2m = 2 c_m G + b_m for m<=7", okc)


def beta_series(dps):
    with mp.workdps(dps+30):
        bb = Fraction(0); tot = mp.mpf(0); m = 0
        while True:
            m += 1
            bb = (Fraction((2*m-1)**2)*bb + 1)/Fraction(4*m*m)
            t = mp.mpf(bb.numerator)/mp.mpf(bb.denominator)/mp.mpf(2)**m
            tot += t
            if m > 60 and t < mp.mpf(10)**(-(dps+25)):
                break
        return +(tot/2)


# ------------------------------------------------- C3: the Legendre-relation split
print("C3  I = (2/pi) K(1/r2) G + beta   and   beta = (K/pi) D")
mp.mp.dps = 90
K = mp.ellipk(mp.mpf(1)/2)
G = mp.catalan
I = mp.quad(lambda k: mp.ellipk(k**2)/(2-k**2), [0, 1])
beta = beta_series(90)
check("I - (2/pi)K G = beta", abs(I - 2*K*G/mp.pi - beta) < mp.mpf(10)**-80)
D1 = mp.quad(lambda k: mp.ellipk(k**2)/k, [1/mp.sqrt(2), 1])
D2 = mp.quad(lambda k: mp.ellipk(1-k**2)/k, [1/mp.sqrt(2), 1])
D2b = mp.quad(lambda k: k*mp.ellipk(k**2)/(1-k**2), [0, 1/mp.sqrt(2)])
check("k->k' involution on the second piece", abs(D2 - D2b) < mp.mpf(10)**-80)
D = D1 - D2
check("beta = (K/pi) D", abs(beta - K*D/mp.pi) < mp.mpf(10)**-80)
Lam = 2*G + D
check("S(4) = (2 r2 K/pi^2) Lambda",
      abs(2*mp.sqrt(2)*K*Lam/mp.pi**2 - S4_hyper) < mp.mpf(10)**-55)
g14 = mp.gamma(mp.mpf(1)/4)
check("S(4) = r2 Gamma(1/4)^2 Lambda/(2 pi^{5/2})",
      abs(mp.sqrt(2)*g14**2*Lam/(2*mp.pi**mp.mpf(2.5)) - S4_hyper) < mp.mpf(10)**-55)

# ------------------------------------------ C4: the two weight-3 Eisenstein identities
print("C4  theta identities (exact integer power series to q^159)")
N = 161


def _mul(a, bb):
    out = [0]*N
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(bb):
                if i+j >= N:
                    break
                if bj:
                    out[i+j] += ai*bj
    return out


def chi(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


t2 = [0]*N; t3 = [0]*N; t4 = [0]*N
n = 0
while n*(n+1) < N:
    t2[n*(n+1)] += 1; n += 1
t3[0] = t4[0] = 1
n = 1
while n*n < N:
    t3[n*n] += 2; t4[n*n] += 2*(-1)**n; n += 1
t2_4 = _mul(_mul(t2, t2), _mul(t2, t2))
t3_2 = _mul(t3, t3)
t4_4 = _mul(_mul(t4, t4), _mul(t4, t4))
h_over_q = [16*x for x in _mul(t2_4, t3_2)]
fser = _mul(t3_2, t4_4)
okh = okf = True
for n in range(1, N-1):
    divs = [d for d in range(1, n+1) if n % d == 0]
    B = sum(chi(n//d)*d*d for d in divs)
    A = sum(chi(d)*d*d for d in divs)
    okh &= (h_over_q[n-1] == 16*B)
    okf &= (fser[n] == -4*A)
check("theta2^4 theta3^2 = 16 sum B_n q^n  (Eisenstein, n<=159)", okh)
check("theta3^2 theta4^4 = 1 - 4 sum A_n q^n (Eisenstein, n<=159)", okf)

# ------------------------------------------------- C5: Lambert / coth / sech forms
print("C5  D and Lambda normal forms at dps 200")
mp.mp.dps = 200
DPS = 200
D_ref = mp.pi*beta_series(DPS)/mp.ellipk(mp.mpf(1)/2)
with mp.workdps(DPS+30):
    Dq = mp.mpf(0); n = -1
    while True:
        n += 2
        t = 4*chi(n)/(n**2*(mp.e**(mp.pi*n)-1))
        Dq += t
        if abs(t) < mp.mpf(10)**(-(DPS+20)):
            break
    Dq = +Dq
    U = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.sech(mp.pi*m)/m**2
        U += t
        if t < mp.mpf(10)**(-(DPS+20)):
            break
    U = +U
check("D = 4 sum chi(n)/(n^2(e^{pi n}-1))", abs(D_ref - Dq) < mp.mpf(10)**-190)
Lam200 = 2*mp.catalan + D_ref
check("Lambda = 5pi^2/24 - (1/2) sum sech(pi m)/m^2",
      abs(Lam200 - (5*mp.pi**2/24 - U/2)) < mp.mpf(10)**-190)
g14 = mp.gamma(mp.mpf(1)/4)
S4_closed = mp.sqrt(2)*g14**2*Lam200/(2*mp.pi**mp.mpf(2.5))
S4_sm = mp.quad(lambda p: (p + mp.asinh(mp.sin(p)))/mp.sqrt(1+mp.sin(p)**2),
                [0, mp.pi/2])*2/mp.pi
check("S(4) closed-form-modulo-Lambda vs quadrature at 200 digits",
      abs(S4_closed - S4_sm) < mp.mpf(10)**-190,
      f"resid {mp.nstr(S4_closed-S4_sm, 3)}")
print("     Lambda =", mp.nstr(Lam200, 60))
print("     S(4)   =", mp.nstr(S4_closed, 60))

# --------------------------------------------------------------- C6: POSITIVE CONTROLS
print("C6  positive controls (a null result is meaningless without these)")
mp.mp.dps = 130
r = mp.sqrt; pi = mp.pi
g13 = mp.gamma(mp.mpf(1)/3); g18 = mp.gamma(mp.mpf(1)/8); g38 = mp.gamma(mp.mpf(3)/8)
z3 = mp.zeta(3)

# C6a  blind rediscovery of pi*S(3) = sqrt3 log(5+2 sqrt6) - 2 atan(sqrt2/5)
MULT = {'1': mp.mpf(1), 'r2': r(2), 'r3': r(3), 'r5': r(5), 'r6': r(6), 'r7': r(7),
        'r10': r(10), 'r15': r(15), 'r2/2': r(2)/2, 'r3/3': r(3)/3,
        '2^(1/4)': mp.mpf(2)**(mp.mpf(1)/4)}
LOG = {'log2': mp.log(2), 'log3': mp.log(3), 'log5': mp.log(5),
       'log(1+r2)': mp.log(1+r(2)), 'log(2+r3)': mp.log(2+r(3)),
       'log(5+2r6)': mp.log(5+2*r(6)), 'log(3+2r2)': mp.log(3+2*r(2)),
       'log(phi)': mp.log((1+r(5))/2), 'atan(r2)': mp.atan(r(2)),
       'atan(r2/5)': mp.atan(r(2)/5), 'atan(1/2)': mp.atan(mp.mpf(1)/2),
       'atan(1/3)': mp.atan(mp.mpf(1)/3), 'pi': pi, 'one': mp.mpf(1)}
prods = [(f"{m}*{l}", MULT[m]*LOG[l]) for m in MULT for l in LOG]
piS3 = pi*mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/3], [1, mp.mpf(4)/3], 1)
hits3 = []
for cc in itertools.combinations(prods, 2):
    rel = mp.pslq([piS3, cc[0][1], cc[1][1]], maxcoeff=10**5, maxsteps=40000,
                  tol=mp.mpf(10)**-120)
    if rel and rel[0] != 0:
        hits3.append((rel, [cc[0][0], cc[1][0]]))
check("C6a blind rediscovery of pi*S(3) (11781 pairs, 154 products)", bool(hits3),
      str(hits3[:1]))

# C6b  the weight-2 Eisenstein Eichler value at tau = i (one weight BELOW the target)
P1 = mp.nsum(lambda n: 1/(n*(mp.e**(2*pi*n)-1)), [1, mp.inf])
LP = {'pi': pi, 'log2': mp.log(2), 'logpi': mp.log(pi), 'log(g14)': mp.log(g14),
      'one': mp.mpf(1), 'G': mp.catalan, 'log3': mp.log(3), 'zeta3': z3}
h1 = [(mp.pslq([P1]+[x[1] for x in cc], maxcoeff=10**8, maxsteps=80000,
               tol=mp.mpf(10)**-100), [x[0] for x in cc])
      for cc in itertools.combinations(list(LP.items()), 4)]
h1 = [x for x in h1 if x[0] and x[0][0] != 0]
check("C6b weight-2 Eisenstein Eichler value at tau=i is FOUND (Gamma(1/4))",
      bool(h1), str(h1[:1]))

# C6c  Ramanujan coth sums + Gamma(1/4)^8 K-moments
P2 = mp.nsum(lambda n: mp.coth(pi*n)/n**3, [1, mp.inf])
check("C6c1 sum coth(pi n)/n^3 = 7 pi^3/180",
      bool(mp.pslq([P2, pi**3], maxcoeff=10**6, tol=mp.mpf(10)**-95)))


def Ksafe(m):
    m = mp.mpf(m)
    return mp.ellipk(1 - mp.mpf(10)**(-mp.mp.dps) if m >= 1 else m)


IK3 = mp.quad(lambda t: Ksafe(mp.sin(t)**2)**3*mp.sin(t), [0, pi/2])
check("C6c2 int_0^1 K'^3 dk = Gamma(1/4)^8/(128 pi^2)",
      abs(IK3 - g14**8/(128*pi**2)) < mp.mpf(10)**-40,
      f"resid {mp.nstr(IK3 - g14**8/(128*pi**2), 3)}")
IK3b = mp.quad(lambda k: Ksafe(k**2)**3, [0, 1])
check("C6c3 int_0^1 K^3 dk = 3 Gamma(1/4)^8/(1280 pi^2)",
      abs(IK3b - 3*g14**8/(1280*pi**2)) < mp.mpf(10)**-40)
check("C6c4 K(sqrt2-1) = sqrt(sqrt2+1) G(1/8)G(3/8)/(2^{13/4} sqrt(pi))",
      abs(mp.ellipk((r(2)-1)**2) - r(r(2)+1)*g18*g38/(mp.mpf(2)**(mp.mpf(13)/4)*r(pi)))
      < mp.mpf(10)**-40)
check("C6c5 K((r3-1)/(2r2)) = 3^{1/4} Gamma(1/3)^3/(2^{7/3} pi)",
      abs(mp.ellipk(((r(3)-1)/(2*r(2)))**2)
          - mp.mpf(3)**mp.mpf(0.25)*g13**3/(mp.mpf(2)**(mp.mpf(7)/3)*pi))
      < mp.mpf(10)**-40)

# ------------------------------- C7: the general-x by-product (addendum 5, R6)
print("C7  int_0^1 K/(1-x k^2) dk = (4/pi)K(sqrt x) G + B(x)  for general x")
mp.mp.dps = 60


def _y1(t):
    return (2/mp.pi)*mp.ellipk(t)


def _y2(t):
    return (2/mp.pi)*mp.ellipk(1-t)


def _B(x):
    a = mp.quad(lambda t: _y2(t)/(1-t), [0, x])
    bq = mp.quad(lambda t: _y1(t)/(1-t), [0, x])
    return (mp.pi/4)*(_y1(x)*a - _y2(x)*bq)


ok7 = True
for xs in ('0.05', '0.3', '0.5', '0.9', '-0.5'):
    x = mp.mpf(xs)
    lhs = mp.quad(lambda k: mp.ellipk(k**2)/(1 - x*k**2), [0, 1])
    ok7 &= abs(lhs - ((4/mp.pi)*mp.ellipk(x)*mp.catalan + _B(x))) < mp.mpf(10)**-50
check("C7 general-x moment formula at x = 0.05, 0.3, 0.5, 0.9, -0.5", ok7)

print()
if FAIL:
    print("*** FAILURES:", FAIL)
else:
    print("ALL THM-3012 ADDENDUM-5 REFEREE CHECKS PASSED")
