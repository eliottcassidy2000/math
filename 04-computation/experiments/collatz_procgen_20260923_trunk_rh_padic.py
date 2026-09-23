#!/usr/bin/env python3
"""
collatz_procgen_20260923_trunk_rh_padic.py

Exact 3-adic facts about the Collatz trunk  T(i) = (4^i - 1)/3  and its relation
to the Kubota-Leopoldt 3-adic zeta function (session collatz-procgen-20260922,
trunk/RH lane, 2026-09-23).

Everything here is exact modular arithmetic on 3-adic integers (Python ints
mod 3^P) or exact rationals (Fractions).  No floating point.

Sections (printed in this order):
  P1  3-adic log/exp self-tests; lambda = log_3(4), v_3(lambda) = 1
  P2  T is an isometry of Z_3 onto Z_3; the trunk chart i(m) = log(3m+1)/log(4)
  P3  fixed points {0, 1, -1/2} (Strassmann), special values T(1/2) = -1,
      T(-1) = -1/4, T(i*) = 1/2 with i* = log(5/2)/log(4)
  P4  minus sheet: exit set -T(n), positive-integer trunk S(i) = 2T(i)+1,
      Jacobsthal interleaving, conjugate chart -T(-i) fixing {0,-1,1/2}
  P5  backward E-game in the logarithmic charts: translation law, kappa formula
  P6  the functional-equation involution s -> 1-s in the trunk chart
  P7  Mazur measure E_{1,c} on Z_3: moments (c = 2 gives Genocchi numbers)
  P8  Iwasawa power series G(X) of zeta_3 in the variable X = 4^s - 1 = 3 T(s):
      coefficients, interpolation check, G(0) a unit, zero-freeness,
      zeta_3 at s = 1/2 and at the E-game point s = i*
  P9  |(s-1) zeta_3(s)|_3 = 3 on s = 1-k (von Staudt-Clausen), k <= KMAX
  P10 archimedean special values: eta(1-2i) = -3 T(i) zeta(1-2i),
      Genocchi G_{2i} = -6 T(i) B_{2i}
  P11 DRIFT/SHEET control: the qx+1 trunks are Iwasawa coordinates too
  P12 half-integer values T(n/2) = ((-2)^n - 1)/3; the exits of -1

Runtime: well under a minute; memory < 100 MB.
"""
import math
import random
import sys
import time
from fractions import Fraction

T0 = time.time()
random.seed(20260923)

P = 60                  # working precision: integers mod 3^P
MOD = 3 ** P


# ----------------------------------------------------------------------------
# basic 3-adic helpers
# ----------------------------------------------------------------------------
def v3(n):
    """3-adic valuation of a nonzero int or Fraction (inf for 0)."""
    if isinstance(n, Fraction):
        if n == 0:
            return math.inf
        return v3(n.numerator) - v3(n.denominator)
    n = int(n)
    if n == 0:
        return math.inf
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def vmod(x, prec):
    """valuation of an element of Z/3^prec given by an int representative"""
    x %= 3 ** prec
    if x == 0:
        return prec          # means ">= prec"
    return v3(x)


def to_Z3(q, prec=P):
    """Fraction / int with 3-free denominator -> int mod 3^prec"""
    q = Fraction(q)
    assert q.denominator % 3 != 0, q
    m = 3 ** prec
    return (q.numerator % m) * pow(q.denominator % m, -1, m) % m


def log1p(t, prec):
    """log(1+t) mod 3^prec for an int t with 3 | t."""
    assert t % 3 == 0
    m = 3 ** prec
    res = 0
    j = 1
    vt = vmod(t, prec + 5) if t % (3 ** (prec + 5)) else prec + 5
    while True:
        vj = v3(j)
        if j * vt - vj >= prec and j - math.log(j, 3) > prec + 2:
            break
        num = pow(t, j, 3 ** (prec + vj))
        term = (num // 3 ** vj) * pow((j // 3 ** vj) % m, -1, m)
        res = (res + term) % m if j % 2 == 1 else (res - term) % m
        j += 1
    return res % m


def exp3(y, prec):
    """exp(y) mod 3^prec for an int y with 3 | y."""
    assert y % 3 == 0
    m = 3 ** prec
    res = 0
    vfact = 0          # v_3(n!)
    ufact = 1          # unit part of n!
    for n in range(0, 2 * prec + 12):
        if n > 0:
            vn = v3(n)
            vfact += vn
            ufact = (ufact * (n // 3 ** vn)) % (3 ** (prec + 2))
        num = pow(y, n, 3 ** (prec + vfact))
        term = (num // 3 ** vfact) * pow(ufact % m, -1, m)
        res = (res + term) % m
    return res


PP = P + 4                              # a little extra precision internally
LAM = log1p(3, PP + 1)                  # lambda = log_3(4)  (mod 3^(P+5))
LAM3 = (LAM // 3) % 3 ** PP             # lambda/3, a 3-adic unit
LAM3INV = pow(LAM3, -1, 3 ** PP)
LOGM2 = log1p(-3, PP + 1)               # log_3(-2)  (-2 = 1 - 3)


def pow4(s, prec=P):
    """4^s for s in Z_3 (int rep mod 3^PP) -> int mod 3^(prec+1)"""
    y = (s * LAM) % 3 ** (PP + 1)
    return exp3(y, prec + 1)


def trunk_exact(s, prec=P):
    x = (pow4(s, prec) - 1) % 3 ** (prec + 1)
    assert x % 3 == 0
    return (x // 3) % 3 ** prec


def ichart(m, prec=P):
    """i(m) = log(3m+1)/log(4) mod 3^prec (m in Z_3)"""
    L = log1p((3 * m) % 3 ** (prec + 2), prec + 2)        # v >= 1
    return ((L // 3) * LAM3INV) % 3 ** prec


def teich_sign(x):
    """omega(x) in {+1,-1} for a 3-adic unit x (int rep)"""
    return 1 if x % 3 == 1 else -1


def logunit(x, prec=P):
    """log<x> for a 3-adic unit x, where <x> = omega(x)^(-1) x in 1+3Z_3"""
    y = x if x % 3 == 1 else (-x)
    return log1p((y - 1) % 3 ** (prec + 2), prec + 2)


def base3(x, n):
    ds = []
    for _ in range(n):
        ds.append(x % 3)
        x //= 3
    return ''.join(str(d) for d in reversed(ds))


def hdr(t):
    print()
    print('=' * 78)
    print(t)
    print('=' * 78)


# ----------------------------------------------------------------------------
hdr('P1  3-adic log/exp self-tests; lambda = log_3(4)')
print(f'working precision P = {P} ternary digits (all checks below are mod 3^(P-3) or finer)')
print(f'v_3(lambda) = {vmod(LAM, P)}   (expected 1)')
print(f'v_3(lambda - 3) = {vmod(LAM - 3, P)}   (expected 2, used by Strassmann)')
print(f'lambda/3 mod 3^12 = {LAM3 % 3**12}  (base 3, low digits last: {base3(LAM3, 12)})')
ok = True
for _ in range(200):
    a = 3 * random.randrange(3 ** (P - 1))
    b = 3 * random.randrange(3 ** (P - 1))
    # exp(a+b) = exp(a) exp(b);  log(exp(a)) = a
    lhs = exp3((a + b) % 3 ** P, P)
    rhs = exp3(a, P) * exp3(b, P) % 3 ** P
    if lhs != rhs:
        ok = False
    if log1p((exp3(a, P) - 1) % 3 ** P, P) % 3 ** (P - 2) != a % 3 ** (P - 2):
        ok = False
print(f'exp(a+b)=exp(a)exp(b) and log(exp(a))=a on 200 random a,b in 3Z_3: {"OK" if ok else "FAIL"}')
print(f'4^(1/2) = exp(lambda/2) == -2 mod 3^{P-2}: '
      f'{(pow4(to_Z3(Fraction(1, 2), PP)) + 2) % 3**(P-2) == 0}')
print(f'log(-2) = lambda/2 mod 3^{P-2}: {(LOGM2 - LAM * to_Z3(Fraction(1,2), PP+1)) % 3**(P-2) == 0}')

# ----------------------------------------------------------------------------
hdr('P2  T(i) = (4^i-1)/3 is an isometry of Z_3 onto Z_3; the trunk chart')
PC = P - 4                   # comparison precision
samples = [random.randrange(3 ** PP) for _ in range(300)]
samples += [to_Z3(Fraction(a, b), PP) for a, b in
            [(1, 2), (-1, 2), (1, 4), (-1, 4), (1, 5), (5, 7), (2, 1), (-1, 1), (0, 1), (1, 1)]]
bad = 0
npairs = 0
TV = [trunk_exact(x) for x in samples]
for ix in range(120):
    for iy in range(120, len(samples)):
        d = vmod(samples[ix] - samples[iy], PC)
        if d >= PC - 2:
            continue
        npairs += 1
        if vmod(TV[ix] - TV[iy], PC) != d:
            bad += 1
print(f'v_3(T(x)-T(y)) == v_3(x-y) on {npairs} pairs (incl. 1/2, -1/2, 1/4, 1/5, 5/7): '
      f'{"all OK" if bad == 0 else str(bad)+" FAILURES"}')
bad = 0
for m in samples:
    if (trunk_exact(ichart(m % 3 ** PP, PP)) - m) % 3 ** PC:
        bad += 1
    if (ichart(trunk_exact(m, PP), PP) - m) % 3 ** PC:
        bad += 1
print(f'i(T(s)) = s and T(i(m)) = m on {len(samples)} samples: {"OK" if bad == 0 else "FAIL"}')
bad = 0
for i in range(0, 200):
    Ti = (4 ** i - 1) // 3
    if (trunk_exact(i) - Ti) % 3 ** PC:
        bad += 1
    if vmod(Ti, 40) != min(v3(i), 40) if i else False:
        bad += 1
print(f'T(i) agrees with the integers (4^i-1)/3 for i<200 and v_3(T(i)) = v_3(i): '
      f'{"OK" if bad == 0 else "FAIL"}')
print('first trunk values T(1..8) =', [(4 ** i - 1) // 3 for i in range(1, 9)],
      '; illegal (3 | T(i) iff 3 | i):', [(4 ** i - 1) // 3 for i in (3, 6, 9)])

# ----------------------------------------------------------------------------
hdr('P3  fixed points and special values')
half = to_Z3(Fraction(1, 2), PP)
mhalf = to_Z3(Fraction(-1, 2), PP)
for name, s, expect in [('T(0)', 0, Fraction(0)), ('T(1)', 1, Fraction(1)),
                        ('T(-1/2)', mhalf, Fraction(-1, 2)),
                        ('T(1/2)', half, Fraction(-1)),
                        ('T(-1)', to_Z3(-1, PP), Fraction(-1, 4)),
                        ('T(2)', 2, Fraction(5)), ('T(-2)', to_Z3(-2, PP), Fraction(-5, 16))]:
    val = trunk_exact(s)
    print(f'  {name:8s} == {str(expect):6s} mod 3^{PC}: {(val - to_Z3(expect)) % 3**PC == 0}')
istar = ichart(to_Z3(Fraction(1, 2), PP), PP)
print(f'  i* := i(1/2) = log(5/2)/log(4);  T(i*) == 1/2 mod 3^{PC}: '
      f'{(trunk_exact(istar) - half) % 3**PC == 0}')
print(f'  i* mod 3^30 (base 3, most significant first): {base3(istar, 30)}')
print(f'  v_3(i* - 1/2) = {vmod(istar - half, PC)},  v_3(i* - 1) = {vmod(istar - 1, PC)},'
      f'  v_3(i*) = {vmod(istar, PC)}')
print('  i* is irrational: if i* = a/b then (5/2)^b = 4^a, i.e. 5^b = 2^(2a+b), impossible for b != 0.')
# Strassmann for f(i) = T(i) - i = (exp(lambda i) - 1)/3 - i
print('  Strassmann for f(i) = T(i) - i = sum_n a_n i^n:  a_1 = lambda/3 - 1,  a_n = lambda^n/(3 n!)')
vals = []
for n in range(1, 31):
    if n == 1:
        vn = vmod(LAM3 - 1, PC)
    else:
        vn = n * 1 - 1 - v3(math.factorial(n))
    vals.append(vn)
print('   v_3(a_n), n=1..30:', vals)
mn = min(vals)
last = max(n for n, v in enumerate(vals, 1) if v == mn)
print(f'   min valuation {mn}, attained last at n = {last}; tail bound v_3(a_n) >= (n-1)/2 > 1 for n >= 4')
print(f'   => f has at most {last} zeros in Z_3; 0, 1, -1/2 are zeros  => Fix(T) = {{0, 1, -1/2}} exactly.')
# finite Hensel count as a sanity check
counts = []
for r in range(1, 10):
    c = 0
    for i in range(3 ** r):
        if (trunk_exact(i, r + 1) - i) % 3 ** r == 0:
            c += 1
    counts.append(c)
print('   #{i mod 3^r : T(i) = i mod 3^r}, r = 1..9:', counts)
dv = []
for z, name in [(0, '0'), (1, '1'), (mhalf, '-1/2')]:
    der = (pow4(z) * LAM3 - 1) % 3 ** PC      # T'(z) - 1 = 4^z lambda/3 - 1
    dv.append(vmod(der, PC))
    print(f"   v_3(T'({name}) - 1) = {dv[-1]}")
print(f'   consistency: a simple zero z with v_3(f\'(z)) = d accounts for 3^d residues mod 3^r (r large);'
      f' 3^{dv[0]} + 3^{dv[1]} + 3^{dv[2]} = {sum(3**d for d in dv)} = the stable count {counts[-1]}.')

# ----------------------------------------------------------------------------
hdr('P4  minus sheet: exit set, positive-integer trunk, Jacobsthal interleaving')
ok = all((2 ** (2 * n + 1) * Fraction(-1, 2) + 1) / 3 == Fraction(1 - 4 ** n, 3) for n in range(40))
print(f'  E_- reverse move (2^k x + 1)/3 from the minus hostile point -1/2 with k = 2n+1 gives -T(n): {ok}')
ok = all((2 ** (2 * n) * (-1) + 1) / Fraction(3) == Fraction(1 - 4 ** n, 3) for n in range(40))
print(f'  ... and from -1 with k = 2n gives -T(n) as well: {ok}')
S = [(2 * 4 ** n + 1) // 3 for n in range(10)]
print('  positive-integer minus trunk S(n) = (2^(2n+1)+1)/3 = 2T(n)+1 (3S-1 a power of 2):', S)
J = [(2 ** t - (-1) ** t) // 3 for t in range(16)]
print('  Jacobsthal J_t = (2^t - (-1)^t)/3:', J)
ok = all(((-2) ** t - 1) // 3 * 1 == (-1) ** t * J[t] for t in range(16)) and \
     all(J[2 * n] == (4 ** n - 1) // 3 and J[2 * n + 1] == S[n] for n in range(8))
print(f'  ((-2)^t - 1)/3 = (-1)^t J_t;  J_(2n) = T(n) (plus trunk), J_(2n+1) = S(n) (minus trunk): {ok}')
print('  so both trunks are the lattice points of ONE chart, t(m) = log(3m+1)/log(-2) = 2 i(m).')


def Sfun(s):
    return (2 * trunk_exact(s) + 1) % 3 ** PC


for name, s, expect in [('S(0)', 0, 1), ('S(1)', 1, 3), ('S(-1)', to_Z3(-1, PP), Fraction(1, 2)),
                        ('S(1/2)', half, -1), ('S(-1/2)', mhalf, 0)]:
    print(f'  {name:8s} == {str(expect):5s}: {(Sfun(s) - to_Z3(expect, PC)) % 3**PC == 0}')
# unique fixed point of S by Newton
i0 = 2
for _ in range(12):
    f = (Sfun(i0) - i0) % 3 ** PC
    fp = (2 * pow4(i0) * LAM3 - 1) % 3 ** PC
    i0 = (i0 - f * pow(fp, -1, 3 ** PC)) % 3 ** PC
vS = [0, vmod(2 * LAM3 - 1, PC)] + [n - 1 - v3(math.factorial(n)) for n in range(2, 12)]
print(f'  S(i) - i = sum_n b_n i^n with b_0 = 1, b_1 = 2 lambda/3 - 1, b_n = 2 lambda^n/(3 n!):'
      f' v_3(b_n), n=0..11 = {vS}')
print(f'  => (Strassmann, N = 1) S has exactly one fixed point i_S, i_S = 2 mod 3; '
      f'S(i_S) = i_S: {(Sfun(i0) - i0) % 3**PC == 0}; i_S mod 3^20 = {base3(i0, 20)}')
for z, name in [(0, '0'), (to_Z3(-1, PP), '-1'), (half, '1/2')]:
    val = (-trunk_exact((-z) % 3 ** PP)) % 3 ** PC
    print(f'  conjugate chart Tt(i) = -T(-i): Tt({name}) == {name}: {(val - z) % 3**PC == 0}')
print(f'  Tt(-1/2) = -T(1/2) = 1: {((-trunk_exact(half)) - 1) % 3**PC == 0}')

# ----------------------------------------------------------------------------
hdr('P5  the backward E-game in logarithmic charts (translation law = kappa formula)')
print('  departure chart delta(x) = log(2^eps x)/log 4 (eps = 0 if x=1 mod 3, 1 if x=2 mod 3),')
print('  arrival chart i(y) = log(3y+1)/log 4.  Claim: a legal move k = 2n+eps(x) lands at')
print('  y = (2^k x - 1)/3 with i(y) = n + delta(x), and v_3(y - h) = v_3(i(y) - i(h)) for every h.')
bad = 0
ntest = 0
HIS = [(h, ichart(h, PP - 2)) for h in (1, half, to_Z3(Fraction(43, 32), PP), to_Z3(Fraction(59, 64), PP))]
for _ in range(300):
    x = random.randrange(3 ** PP)
    if x % 3 == 0:
        continue
    eps = 0 if x % 3 == 1 else 1
    delta = ((logunit((2 ** eps * x) % 3 ** (PP + 2), PP + 2) // 3) * LAM3INV) % 3 ** PP
    # note: 2^eps x = 1 mod 3, so <2^eps x> = 2^eps x
    for n in range(0, 25):
        k = 2 * n + eps
        num = (pow(2, k, 3 ** (PP + 1)) * x - 1) % 3 ** (PP + 1)
        assert num % 3 == 0
        y = num // 3
        iy = ichart(y, PP - 2)
        ntest += 1
        if (iy - (n + delta)) % 3 ** (PC - 2):
            bad += 1
        for h, ih in HIS:
            d1 = vmod(y - h, PC - 2)
            d2 = vmod(iy - ih, PC - 2)
            if d1 < PC - 4 and d1 != d2:
                bad += 1
print(f'  checked {ntest} moves (and 4 targets h each): {"all OK" if bad == 0 else str(bad)+" FAILURES"}')
print(f'  delta(1) = delta(1/2) = 0 (both primary hostile points sit at the chart origin): '
      f'{logunit(1) == 0 and logunit((2 * half) % 3**PP) == 0}')
print('  precision k(x) = v_3(x - h(x)) = 1 + v_3(delta(x)):', end=' ')
bad = 0
for _ in range(500):
    x = random.randrange(3 ** PP)
    if x % 3 == 0:
        continue
    eps = 0 if x % 3 == 1 else 1
    h = 1 if eps == 0 else half
    delta = ((logunit((2 ** eps * x) % 3 ** (PP + 2), PP + 2) // 3) * LAM3INV) % 3 ** PP
    k1 = vmod(x - h, PC)
    if k1 < PC - 3 and k1 != 1 + vmod(delta, PC):
        bad += 1
print('OK' if bad == 0 else f'{bad} FAILURES')
# the foundry's kappa formula, and kappa(3m+1) = 2 i(m)
bad = 0
cnt = 0
for _ in range(60):
    r = random.randrange(1, 3 ** 12)
    if r % 3 == 0:
        continue
    e = 0 if r % 3 == 1 else 1
    r1 = r if e == 0 else (-r) % 3 ** PP
    kap = ((log1p((r1 - 1) % 3 ** (PP + 2), PP + 2) // 3) *
           pow((LOGM2 // 3) % 3 ** PP, -1, 3 ** PP)) % 3 ** PP
    for K in range(0, 400):
        cnt += 1
        lhs = v3(2 ** K - r) if 2 ** K != r else math.inf
        rhs = 0 if (K - e) % 2 else 1 + vmod(K - kap, PC)
        if lhs != rhs and not (lhs == math.inf):
            bad += 1
print(f'  kappa formula v_3(2^K - r) = [K=e mod 2](1 + v_3(K - kappa(r))), '
      f'kappa(r) = log(r_1)/log(-2): {cnt} cases, {"OK" if bad == 0 else str(bad)+" FAILURES"}')
bad = 0
for _ in range(200):
    m = random.randrange(3 ** PP)
    kap = ((log1p((3 * m) % 3 ** (PP + 2), PP + 2) // 3) *
           pow((LOGM2 // 3) % 3 ** PP, -1, 3 ** PP)) % 3 ** PP
    if (kap - 2 * ichart(m, PP)) % 3 ** PC:
        bad += 1
print(f'  kappa(3m+1) = 2 i(m) (kappa IS the doubled trunk coordinate): {"OK" if bad == 0 else "FAIL"}')

# ----------------------------------------------------------------------------
hdr('P6  the functional-equation involution s -> 1-s in the trunk chart')
print('  M(m) = (1-m)/(3m+1).  Claim: T(1-s) = M(T(s)) for all s in Z_3.')
bad = 0
for s in samples[:200]:
    t = trunk_exact(s, PP)
    lhs = trunk_exact((1 - s) % 3 ** PP, PP)
    rhs = ((1 - t) * pow((3 * t + 1) % 3 ** PP, -1, 3 ** PP)) % 3 ** PP
    if (lhs - rhs) % 3 ** PC:
        bad += 1
print(f'  checked on 200 samples: {"OK" if bad == 0 else "FAIL"}')


def M(q):
    q = Fraction(q)
    return (1 - q) / (3 * q + 1)


print(f'  M is an involution: M(M(m)) = m (exact, 20 rationals): '
      f'{all(M(M(Fraction(a, b))) == Fraction(a, b) for a, b in [(1,2),(2,5),(7,11),(-4,5),(3,8)]*4)}')
print(f'  fixed points of M in Q: m = -1 (= T(1/2)) and m = 1/3 (not in Z_3); M(-1) = {M(-1)}')
print(f'  M(1) = {M(1)}, M(1/2) = {M(Fraction(1,2))}, M(0) = {M(0)}, M(-1/2) = {M(Fraction(-1,2))}')
print('  => the hostile set {1, 1/2} of the backward E-game is NOT M-invariant (1 -> 0 illegal, 1/2 -> 1/5');
print('     generic): the E-game has no s <-> 1-s symmetry in the trunk chart.')

# ----------------------------------------------------------------------------
hdr('P7  Mazur measure E_{1,c} on Z_3 and its moments')


def bernoulli_list(n):
    B = [Fraction(1)]
    for m in range(1, n + 1):
        s = Fraction(0)
        for j in range(m):
            s += math.comb(m + 1, j) * B[j]
        B.append(-s / (m + 1))
    return B        # B_1 = -1/2


KMAX = 300
B = bernoulli_list(KMAX + 2)
assert B[1] == Fraction(-1, 2) and B[2] == Fraction(1, 6) and B[12] == Fraction(-691, 2730)


def mazur_twice(c, N):
    """return list e2[a] = 2*E_{1,c}(a + 3^N Z_3) (an integer) for a in [0, 3^N)"""
    Mm = 3 ** N
    cinv = pow(c, -1, Mm)
    e2 = []
    for a in range(Mm):
        b = (cinv * a) % Mm
        # E = (a - c b)/M + (c-1)/2, and (a - c b) is divisible by M
        q, r = divmod(a - c * b, Mm)
        assert r == 0
        e2.append(2 * q + (c - 1))
    return e2


for c in (2, 4):
    N = 9
    e2 = mazur_twice(c, N)
    rows = []
    for k in range(1, 13):
        tot = sum(pow(a, k - 1) * e2[a] for a in range(3 ** N))
        riem = Fraction(tot, 2)
        true = (1 - Fraction(c) ** k) * B[k] / k
        rows.append((k, v3(riem - true) if riem != true else math.inf))
        # restricted to units: extra Euler factor (1 - 3^(k-1))
    print(f'  c = {c}, N = {N}: v_3( sum_a a^(k-1) E(a+3^N Z_3) - (1-c^k)B_k/k ), k=1..12:',
          [r[1] for r in rows])
    rows = []
    for k in range(1, 13):
        tot = sum(pow(a, k - 1) * e2[a] for a in range(3 ** N) if a % 3)
        riem = Fraction(tot, 2)
        true = (1 - Fraction(c) ** k) * (1 - Fraction(3) ** (k - 1)) * B[k] / k
        rows.append(v3(riem - true) if riem != true else math.inf)
    print(f'        restricted to Z_3^x vs (1-c^k)(1-3^(k-1))B_k/k:', rows)
print('  (valuations ~ N - O(log k): the Riemann sums converge 3-adically to the stated moments;')
print('   c = 2 gives the Genocchi moments (1-2^k)B_k/k = G_k/(2k) on Z_3.)')

# ----------------------------------------------------------------------------
hdr('P8  Iwasawa power series of zeta_3 in the trunk variable X = 4^s - 1 = 3 T(s)')
print('  zeta_3(1-s) = G(X)/X,  G(X) = int_{Z_3^x} (1+X)^{l(x)} x^{-1} dE_{1,4}(x),  l(x) = log<x>/log 4')
print('  (Washington normalization L_3(s,1); RJW Thm 4.1: int x^k zeta_p = (1-p^(k-1)) zeta(1-k)).')
NI = 11
M_ = 3 ** NI
e2 = mazur_twice(4, NI)
# discrete log base 4 on 1 + 3Z mod 3^NI:  <a> = 4^(l(a)),  l(a) mod 3^(NI-1)
dlog = {}
x = 1
for j in range(3 ** (NI - 1)):
    dlog[x] = j
    x = (x * 4) % M_
NCOEF = 10
R = 3 ** (NI - 1)
coef = []
inv2 = pow(2, -1, 3 ** (NI + 2))
modn = 3 ** (NI + 2)
UNITS = []
for a in range(1, M_):
    if a % 3 == 0:
        continue
    br = a if a % 3 == 1 else (M_ - a)
    UNITS.append((dlog[br % M_], pow(a, -1, modn) * e2[a] % modn))
for n in range(NCOEF):
    tot = 0
    for L, w in UNITS:
        tot += math.comb(L, n) * w
    cn = (tot * inv2) % modn
    prec_n = NI - 1 - v3(math.factorial(n))
    coef.append((cn % 3 ** prec_n, prec_n))
for n, (cn, pr) in enumerate(coef):
    print(f'   c_{n} mod 3^{pr} = {cn:>8d}   (base 3: {base3(cn, pr)})   v_3(c_{n}) = '
          f'{vmod(cn, pr) if cn % 3**pr else ">=" + str(pr)}')
c0 = coef[0][0]
print(f'  G(0) = c_0 is a 3-adic UNIT (c_0 mod 3 = {c0 % 3}) => G is a unit of Z_3[[X]]:')
print('   mu = lambda = 0, and G(X) != 0 for every X in the open unit disc of C_3.')
# residue check: G(0) = -(2/3) log 4 = -2 lambda/3
chk = (-2 * LAM3) % 3 ** coef[0][1]
print(f'  residue check G(0) = -(2/3) log 4  (pole of zeta_3 at s=1 with residue 1-1/3): '
      f'{(c0 - chk) % 3**coef[0][1] == 0}')


def Gval(Xv, prec):
    """evaluate sum c_n X^n mod 3^prec (X int with 3 | X)"""
    s = 0
    for n, (cn, pr) in enumerate(coef):
        s += cn * pow(Xv, n, 3 ** (prec + 2))
    return s % 3 ** prec


def Lomega(k):
    """L(omega, 1-k) = -B_{k,omega}/k for the quadratic character mod 3 (k odd)"""
    Bk_poly = lambda xx: sum(math.comb(k, j) * B[j] * Fraction(xx) ** (k - j) for j in range(k + 1))
    Bko = Fraction(3) ** (k - 1) * (Bk_poly(Fraction(1, 3)) - Bk_poly(Fraction(2, 3)))
    return -Bko / k


print('  interpolation check: G(4^k-1)/(4^k-1) vs zeta_3(1-k) = (1-3^(k-1))zeta(1-k) (k even),')
print('                                        vs L(omega,1-k) (k odd; RJW Thm 5.1, chi = omega):')
for k in range(1, 13):
    X = 4 ** k - 1
    if k % 2 == 0:
        target = (1 - Fraction(3) ** (k - 1)) * (-B[k] / k)
    else:
        target = Lomega(k)
    num_target = target * X            # = G(X), a 3-adic unit
    ok = (Gval(X, 6) - to_Z3(num_target, 6)) % 3 ** 6 == 0
    print(f'   k={k:2d}: zeta_3(1-k) = {str(target):>28s}   v_3 = {v3(target):3d}   '
          f'G(4^k-1) = X*zeta_3(1-k) mod 3^6: {"OK" if ok else "MISMATCH"}')
Xhalf = (pow4(half, 20) - 1) % 3 ** 21
print(f'  X at s = 1/2 is 4^(1/2)-1 = -3 (= 3*T(1/2) = 3*(-1)): {(Xhalf + 3) % 3**20 == 0}')
g = Gval(-3, 7)
print(f'  zeta_3(1/2) = G(-3)/(-3):  G(-3) mod 3^7 = {g} (unit: {g % 3 != 0}) => |zeta_3(1/2)|_3 = 3, nonzero')
g2 = Gval(to_Z3(Fraction(3, 2), 12), 7)
print(f'  E-game point s = i*: X = 3 T(i*) = 3/2;  G(3/2) mod 3^7 = {g2} (unit: {g2 % 3 != 0})'
      f' => zeta_3(1-i*) = (2/3) G(3/2) != 0, |.|_3 = 3')

# ----------------------------------------------------------------------------
hdr('P9  |(s-1) zeta_3(s)|_3 = 3 at s = 1-k (von Staudt-Clausen), exact')
bad_even = [k for k in range(2, KMAX + 1, 2) if v3((1 - Fraction(3) ** (k - 1)) * B[k]) != -1]
print(f'  even k <= {KMAX}: v_3((s-1) zeta_3(s)) = v_3((1-3^(k-1)) B_k) = -1 for all: {not bad_even}')
bad_odd = []
for k in range(1, 121, 2):
    if v3(Lomega(k) * k) != -1:
        bad_odd.append(k)
print(f'  odd k <= 119:  v_3(k L(omega,1-k)) = v_3(B_(k,omega)) = -1 for all: {not bad_odd}')
print('  Proof (all k): for even k, (3-1) | k so 3 B_k = -1 mod 3; odd k similarly via B_(k,omega).')
print('  Since the points 1-k are dense in Z_3 and (s-1)zeta_3(s) is continuous, |(s-1)zeta_3(s)|_3 = 3')
print('  on all of Z_3: zeta_3 has NO zeros (in particular none at s = 1/2, none at s = i*).')

# ----------------------------------------------------------------------------
hdr('P10 archimedean special values: eta and Genocchi numbers carry the trunk')
ok = True
for i in range(1, 41):
    Ti = Fraction(4 ** i - 1, 3)
    zeta_neg = -B[2 * i] / (2 * i)                 # zeta(1-2i)
    eta = (1 - Fraction(2) ** (2 * i)) * zeta_neg  # eta(s) = (1-2^(1-s)) zeta(s) at s = 1-2i
    if eta != -3 * Ti * zeta_neg:
        ok = False
print(f'  eta(1-2i) = (1 - 4^i) zeta(1-2i) = -3 T(i) zeta(1-2i), i = 1..40: {ok}')
# Genocchi numbers from the generating function 2t/(e^t+1)
NG = 60
ex = [Fraction(1, math.factorial(n)) for n in range(NG + 2)]
den = ex[:]                 # e^t + 1
den[0] += 1
# solve den * g = 2t  (power series)
g = [Fraction(0)] * (NG + 1)
rhs = [Fraction(0)] * (NG + 1)
rhs[1] = Fraction(2)
for n in range(NG + 1):
    s = rhs[n] - sum(den[j] * g[n - j] for j in range(1, n + 1))
    g[n] = s / den[0]
Gen = [g[n] * math.factorial(n) for n in range(NG + 1)]
ok1 = all(Gen[n] == 2 * (1 - Fraction(2) ** n) * B[n] for n in range(2, NG + 1))
ok2 = all(Gen[2 * i] == -6 * Fraction(4 ** i - 1, 3) * B[2 * i] for i in range(1, NG // 2 + 1))
print(f'  Genocchi numbers (2t/(e^t+1)): G_2..G_12 = {[int(Gen[n]) for n in range(2, 13, 2)]}')
print(f'  G_n = 2(1-2^n)B_n (n<={NG}): {ok1};   G_(2i) = -6 T(i) B_(2i) (i<={NG//2}): {ok2}')
print('  These are values at NEGATIVE integers; zeta(1-2i) != 0; the factor 1-4^i vanishes only at')
print('  i = pi*k*sqrt(-1)/log 2, i.e. s = 1-2i on the line Re s = 1 (zeros of 1-2^(1-s)), never on Re s = 1/2.')

# ----------------------------------------------------------------------------
hdr('P11 DRIFT/SHEET control: the trunk of every qx+1 map is an Iwasawa coordinate (q non-Wieferich)')
print('  trunk of qx+1 = {x : qx+1 is a power of 2} = {(u^i - 1)/q}, u = 2^ord_q(2).  LTE: v_q(u^i - 1) =')
print('  v_q(u - 1) + v_q(i); so i -> (u^i-1)/q is an isometry of Z_q onto Z_q iff v_q(u - 1) = 1.')
for q in (3, 5, 7, 11, 13, 17, 31, 127, 1093):
    o = 1
    while pow(2, o, q) != 1:
        o += 1
    u = 2 ** o
    vu = v3(u - 1) if q == 3 else None
    vq = 0
    t = u - 1
    while t % q == 0:
        t //= q
        vq += 1
    ok = True
    for i in range(1, 400):
        t = pow(u, i) - 1
        vv = 0
        while t % q == 0:
            t //= q
            vv += 1
        vi = 0
        ii = i
        while ii % q == 0:
            ii //= q
            vi += 1
        if vv != vq + vi:
            ok = False
    kind = ('u = 1+q (Mersenne prime)' if u == 1 + q else f'u = {u}' if u < 10**6 else f'u = 2^{o}')
    print(f'   q = {q:5d}: ord_q(2) = {o:4d}, {kind:26s} v_q(u-1) = {vq}  LTE check i<400: {ok}  '
          f'=> trunk {"IS" if vq == 1 else "is NOT"} an isometry onto Z_q'
          + ('' if vq == 1 else ' (Wieferich prime: image q^' + str(vq - 1) + 'Z_q)'))
print('  The 5x+1 and 7x+1 maps (expected to have divergent orbits) carry the same structure; 7x+1 even has')
print('  u = 1+q exactly, like 3x+1.  The minus sheet uses the negated trunk.  So the identity is DRIFT- and')
print('  SHEET-blind: it cannot carry information that distinguishes 3x+1 from its controls.')

# ----------------------------------------------------------------------------
hdr('P12 half-integer values: T(n/2) = ((-2)^n - 1)/3 and the exits of -1')
ok = True
for n in range(-20, 41):
    val = trunk_exact(to_Z3(Fraction(n, 2), PP))
    exp_ = (Fraction(-2) ** n - 1) / 3
    if (val - to_Z3(exp_)) % 3 ** PC:
        ok = False
print(f'  T(n/2) = ((-2)^n - 1)/3 for -20 <= n <= 40 (since 4^(1/2) = -2 in 1+3Z_3): {ok}')
ok = all(Fraction((-2) ** (2 * n + 1) - 1, 3) == -Fraction(2 ** (2 * n + 1) + 1, 3) for n in range(30))
print(f'  hence T(n + 1/2) = -S(n) = -(2^(2n+1)+1)/3: the plus trunk passes through the negated minus trunk: {ok}')
dm1 = ((logunit((2 * (MOD - 1)) % 3 ** (PP + 2), PP + 2) // 3) * LAM3INV) % 3 ** PC
dmh = ((logunit(mhalf % 3 ** (PP + 2), PP + 2) // 3) * LAM3INV) % 3 ** PC
print(f'  departure chart: delta(-1) = 1/2: {(dm1 - to_Z3(Fraction(1, 2), PC)) % 3**PC == 0};  '
      f'delta(-1/2) = -1/2: {(dmh - to_Z3(Fraction(-1, 2), PC)) % 3**PC == 0};  delta(1) = delta(1/2) = 0')
ok = all((2 ** (2 * n + 1) * (-1) - 1) // 3 == -(2 ** (2 * n + 1) + 1) // 3 for n in range(30))
print(f'  plus-game reverse moves from -1 with k = 2n+1 land at -S(n) = T(n + 1/2): {ok}')
print('  So "s = 1/2" in the trunk chart is the first exit of the point -1 (a hostile point of the forward game')
print('  and of the minus sheet), because log(-2)/log(4) = 1/2.  The number 1/2 here is the exponent relating the')
print('  two generators -2 and 4 of 1+3Z_3.')

print()
print(f'[padic] done in {time.time() - T0:.1f} s')
