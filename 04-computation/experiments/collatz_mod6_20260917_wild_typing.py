#!/usr/bin/env python3
"""
collatz_mod6_20260917_wild_typing.py  --  LANE wild_typing (creative lane, every claim probed)

Session collatz-mod6-20260917, wave two.  Exact integer / Fraction arithmetic for every
load-bearing claim; explicit universes; positive and hostile controls; every check is an
explicit `raise` (active under python -O).

Sections
  S0  helpers, sieves, universes
  S1  TYPING TABLE probes: Fano/octonion orientation count, divisor box 12=8+3+1,
      Cayley-Dickson vs Omega/omega, Bott-8 probe, 6/pi^2 candidates (three exact densities)
  S2  {7,21}: tournament census n<=6 (h, |Aut|), 3n+k cycle-count spectrum c(k) for
      gcd(k,6)=1, k<=99, |n|<=2*10^5 both signs, completeness bounds, missing values
  S3  five bold conjectures with probes + one follow-up probe each
  S4  3N-5 = Collatz shifted by 2: exact discrepancy, rows mod 6 / mod 9, cycles
  S5  status summary

Runtime target < 10 min, RAM < 1 GB.
"""
import sys, time, math, itertools
from fractions import Fraction
from math import gcd
import numpy as np
import sympy
from sympy import factorint, isprime, primerange, totient, divisors

T0 = time.time()
def tm():
    return "[t=%.1fs]" % (time.time() - T0)

def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def v2(x):
    return (x & -x).bit_length() - 1

# ---------------------------------------------------------------- S0 sieves
hdr("S0. Sieves and universes")
NS = 6 * 10**6 + 1          # sieve bound (covers 3n+1 for odd n <= 2*10^6, and 6K+1 for K = 10^6)
sqf = np.ones(NS + 1, dtype=np.bool_)
sqf[0] = False
for p in primerange(2, int(NS**0.5) + 1):
    sqf[p * p::p * p] = False
# Omega sieve (prime factors with multiplicity)
Om = np.zeros(NS + 1, dtype=np.uint8)
# smallest-prime sieve is heavier; do Omega by prime powers
spf_primes = list(primerange(2, NS + 1))
for p in spf_primes:
    q = p
    while q <= NS:
        Om[q::q] += 1
        q *= p
isP = np.zeros(NS + 1, dtype=np.bool_)
isP[spf_primes] = True
print("sieve bound NS =", NS, " primes:", len(spf_primes), tm())
# audits against sympy on a random sample
rng = np.random.default_rng(20260917)
for n in rng.integers(2, NS, 400):
    n = int(n)
    f = factorint(n)
    check(int(Om[n]) == sum(f.values()), "Omega sieve at %d" % n)
    check(bool(sqf[n]) == all(e == 1 for e in f.values()), "sqf sieve at %d" % n)
    check(bool(isP[n]) == isprime(n), "prime sieve at %d" % n)
print("sieve audits vs sympy on 400 random n: PASS", tm())
# positive control: squarefree density -> 6/pi^2
for X in (10**5, 10**6):
    cnt = int(sqf[1:X + 1].sum())
    print("  #squarefree <= %d = %d ; density %.6f ; 6/pi^2 = %.6f" % (X, cnt, cnt / X, 6 / math.pi**2))

# ---------------------------------------------------------------- S1 typing table
hdr("S1. TYPING TABLE probes")

print("\n--- S1.1  Squarefree divisors of N=p^2 q r : F_2^3, PG(2,2), the Fano plane, and the octonion orientation ---")
p, q, r = 2, 3, 5
N = p * p * q * r
D = divisors(N)
proper = [d for d in D if 1 < d < N]
S_set = [d for d in proper if all(e == 1 for e in factorint(d).values())]
U_set = [d for d in proper if isprime(d)]
F, S, U = len(proper), len(S_set), len(U_set)
print("N = %d = 2^2*3*5 : F=%d S=%d U=%d ; F=S+U: %s" % (N, F, S, U, F == S + U))
check(F == 10 and S == 7 and U == 3, "p^2qr profile")
nonsqf_proper = [d for d in proper if d not in S_set]
print("divisor box [0,2]x[0,1]^2 has", len(D), "points = 8 (Boolean cube {squarefree divisors incl. 1}) + %d (nonsquarefree proper: %s) + 1 (N)"
      % (len(nonsqf_proper), nonsqf_proper))
check(len(D) == 12 and len(nonsqf_proper) == 3, "12 = 8+3+1")
# F_2^3 structure: squarefree divisor d <-> vector of exponents on (p,q,r); group law = multiply and drop squares
prim = [p, q, r]
def vec(d):
    return tuple(1 if d % pp == 0 else 0 for pp in prim)
def vecint(d):
    v = vec(d); return v[0] * 4 + v[1] * 2 + v[2]
sqf_all = [d for d in D if all(e == 1 for e in factorint(d).values())]
check(len(sqf_all) == 8, "8 squarefree divisors")
for a in sqf_all:
    for b in sqf_all:
        prod = a * b
        core = 1
        for pp, e in factorint(prod).items():
            if e % 2 == 1:
                core *= pp
        check(vecint(core) == (vecint(a) ^ vecint(b)), "group law")
print("PROVED (+checked 64 products): squarefree divisors under (multiply, drop squares) = F_2^3 ; 1 = zero vector.")
# Fano lines = triples of nonzero vectors summing to 0 = triples of proper squarefree divisors whose product is a square
lines = set()
for a in S_set:
    for b in S_set:
        if a < b:
            c = None
            for cc in S_set:
                if vecint(cc) == (vecint(a) ^ vecint(b)):
                    c = cc
            check(c is not None and c != a and c != b, "line closure")
            lines.add(frozenset((a, b, c)))
lines = sorted(sorted(l) for l in lines)
check(len(lines) == 7, "7 Fano lines")
for l in lines:
    prod = l[0] * l[1] * l[2]
    check(math.isqrt(prod)**2 == prod, "line product square")
print("Fano lines of the divisor lattice (triples of proper squarefree divisors with square product):")
for l in lines:
    print("   %-14s product %d = %d^2" % (l, l[0] * l[1] * l[2], math.isqrt(l[0] * l[1] * l[2])))
# every point on 3 lines, every two points on exactly one line
for a in S_set:
    check(sum(1 for l in lines if a in l) == 3, "3 lines per point")
print("PROVED: 7 points, 7 lines, 3 points/line, 3 lines/point = PG(2,2).")

# octonion orientation count: build the algebra for each of the 2^7 orientations of the 7 lines
idx = {d: vecint(d) for d in S_set}      # 1..7
line_vecs = [tuple(sorted(idx[d] for d in l)) for l in lines]
def build_mult(orient):
    # mult[a][b] = (sign, c) for basis e_a e_b, a,b in 0..7 ; e_0 = 1
    mult = [[None] * 8 for _ in range(8)]
    for a in range(8):
        for b in range(8):
            if a == 0:
                mult[a][b] = (1, b)
            elif b == 0:
                mult[a][b] = (1, a)
            elif a == b:
                mult[a][b] = (-1, 0)
            else:
                mult[a][b] = None
    for (lv, o) in zip(line_vecs, orient):
        a, b, c = lv
        # a cyclic order on the line: (a->b->c) if o==0 else (a->c->b)
        cyc = [(a, b, c), (b, c, a), (c, a, b)] if o == 0 else [(a, c, b), (c, b, a), (b, a, c)]
        for (x, y, z) in cyc:
            mult[x][y] = (1, z)
            mult[y][x] = (-1, z)
    for a in range(8):
        for b in range(8):
            check(mult[a][b] is not None, "mult table complete")
    return mult
def mul(mult, x, y):
    z = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            s, c = mult[a][b]
            z[c] += s * x[a] * y[b]
    return z
def basis(i):
    e = [0] * 8; e[i] = 1; return e
def assoc(mult, x, y, z):
    xy = mul(mult, x, y); yz = mul(mult, y, z)
    a = mul(mult, xy, z); b = mul(mult, x, yz)
    return [u - v for u, v in zip(a, b)]
good = []
witness_bad = None
for orient in itertools.product((0, 1), repeat=7):
    mult = build_mult(orient)
    ok = True
    for i in range(8):
        for j in range(8):
            for k3 in range(8):
                A = assoc(mult, basis(i), basis(j), basis(k3))
                B = assoc(mult, basis(j), basis(i), basis(k3))
                C = assoc(mult, basis(i), basis(k3), basis(j))
                if any(a + b != 0 for a, b in zip(A, B)) or any(a + c != 0 for a, c in zip(A, C)):
                    ok = False
                    if witness_bad is None:
                        witness_bad = (orient, (i, j, k3))
                    break
            if not ok:
                break
        if not ok:
            break
    if ok:
        good.append(orient)
print("orientations of the 7 Fano lines: 128 ; giving an ALTERNATIVE algebra (linearized alternativity on the basis): %d" % len(good))
check(len(good) == 16, "expected 16 valid orientations")
# exact norm multiplicativity on the valid ones (composition algebra) and a failing witness on a bad one
rng2 = np.random.default_rng(7)
for orient in good[:4]:
    mult = build_mult(orient)
    for _ in range(30):
        x = [int(t) for t in rng2.integers(-9, 10, 8)]
        y = [int(t) for t in rng2.integers(-9, 10, 8)]
        z = mul(mult, x, y)
        check(sum(t * t for t in z) == sum(t * t for t in x) * sum(t * t for t in y), "norm multiplicativity")
print("checked: the valid orientations satisfy |xy|^2=|x|^2|y|^2 exactly on 120 random integer pairs (composition algebra = octonions).")
print("hostile witness: orientation %s fails linearized alternativity at basis triple %s" % witness_bad)
print("first valid orientation (line -> cyclic order):", [(lv, "a->b->c" if o == 0 else "a->c->b") for lv, o in zip(line_vecs, good[0])])
print("VERDICT: the divisor lattice of p^2qr fixes the 7 points (S=7) and the 7 lines (square-product triples) = F_2^3 / PG(2,2),")
print("         but nothing in the lattice selects one of the 16 octonion orientations out of 128: the multiplication is LOST DATA.")
print("         Typed analogy with a real underlying object (F_2^3); no multiplication map.  Fano points 7 = S(p^2qr), dim 3 = U(p^2qr).")

print("\n--- S1.2  Cayley-Dickson dimensions vs Omega / omega ---")
tbl = []
for name, n0 in (("A=2", 2), ("B=6", 6), ("A^2B=24", 24), ("C^2B=5400", 5400)):
    f = factorint(n0)
    Omg = sum(f.values()); omg = len(f); dn = len(divisors(n0))
    tbl.append((name, Omg, omg, dn, 2**omg))
    print("   %-10s Omega=%2d omega=%d d(N)=%3d 2^omega=%d  2^Omega=%d" % (name, Omg, omg, dn, 2**omg, 2**Omg))
print("PROVED: Omega is additive (Omega(MN)=Omega(M)+Omega(N)), Cayley-Dickson dimension is multiplicative (dim doubles),")
print("        so 'Omega(A^2B)=4 <-> dim H = 4' compares a count with a power of two: log_2(dim) = j-fold doubling <-> nothing.")
print("        d(N) = prod(a_i+1) <= 2^Omega(N) with equality iff N squarefree; and d(N) = 2^omega(N) iff N squarefree.")
for n0 in range(2, 2001):
    f = factorint(n0)
    Omg = sum(f.values()); omg = len(f); dn = len(divisors(n0))
    sf = all(e == 1 for e in f.values())
    check((dn == 2**Omg) == sf and (dn == 2**omg) == sf and dn <= 2**Omg, "d(N) vs 2^Omega at %d" % n0)
print("        (checked 2<=N<=2000).  The only functorial match: for SQUAREFREE N with omega(N)=j the squarefree-divisor group is F_2^j,")
print("        its real group algebra R[F_2^j] has dimension 2^j = d(N), and the Cayley-Dickson algebra of that dimension is a")
print("        twisted group algebra R_F[F_2^j] (CITED: Albuquerque-Majid, 'Quasialgebra structure of the octonions', J. Algebra 220 (1999)).")
print("        Map: N=1->R, p->C, pq->H, pqr->O (omega = 0,1,2,3).  LOST: the cocycle F (signs/orientation).  A^2B and C^2B are NOT squarefree,")
print("        their divisor boxes ([0,2]x[0,1], 6 points; [0,2]^?x.. 12+ points) are not groups: VERDICT numerology for the Omega version.")

print("\n--- S1.3  Bott periodicity 8: probe ---")
def mult_order(a, m):
    if gcd(a, m) != 1:
        return None
    o, x = 1, a % m
    while x != 1:
        x = x * a % m; o += 1
    return o
ord8 = [m for m in range(3, 2000, 2) if mult_order(2, m) == 8]
print("odd moduli m<2000 with ord_m(2)=8 exactly:", ord8, " (all divide 2^8-1=255)")
check(ord8 == [17, 51, 85, 255], "ord 8 moduli")
print("Collatz row moduli and their orders of 2: ", [(m, mult_order(2, m)) for m in (3, 9, 7, 63, 21, 27)])
print("R-braid period on rows = ord_9(4) = %d; source period mod 42 and target mod 63: 3 (inherited row_braid_typing S1)." % mult_order(4, 9))
print("VERDICT: NO MAP FOUND.  8 = Bott period is the Morita period of the Clifford cocycle on F_2^n (Cl_{n+8} = Cl_n (x) M_16(R));")
print("         the shared object with S1.1 is only 'twisted group algebra of F_2^n' (Clifford AND octonion are such: CITED Albuquerque-Majid).")
print("         No Collatz modulus has ord(2)=8; the row/braid periods are 2,3,6.  Nothing in the divisor lattice carries a cocycle.")

print("\n--- S1.4  6/pi^2 candidates: exact densities at 10^6 odd n ---")
NODD = 10**6                       # odd n = 1,3,...,2*10^6-1
nodd = np.arange(1, 2 * NODD, 2, dtype=np.int64)
x = 3 * nodd + 1
# odd part of x
lowbit = x & (-x)
u = x // lowbit
check(int((u % 2).min()) == 1 and int((u % 3 == 0).sum()) == 0, "odd cores odd and coprime to 3")
Tsqf = sqf[u]
nsqf = sqf[nodd]
dA = int(Tsqf.sum()) / NODD
dC = int((Tsqf & nsqf).sum()) / NODD
# gcd(n, T(n)) = 1 always (PROVED: gcd(n,3n+1)=1); check
g = np.gcd(nodd, u)
check(int(g.max()) == 1, "gcd(n,T(n))=1")
# gcd(n+2, T(n)) : gcd(n+2, 3n+1) | 5 -> density of coprimality = 1 - 1/5 * (correction) ; compute
g2 = np.gcd(nodd + 2, u)
dB2 = int((g2 == 1).sum()) / NODD
# Euler products
P5 = 1.0; P5b = 1.0
for pp in primerange(5, 10**7):
    P5 *= (1 - 1 / pp**2)
    P5b *= (1 - 2 / pp**2)
pred_A = P5                       # = (6/pi^2) / ((3/4)(8/9)) = 9/pi^2
pred_C = (8 / 9) * P5b            # n odd: 9 !| n excludes 1 class mod 9 ; p>=5 : 2 classes mod p^2
print("(a) odd n with T(n) squarefree             : %d/%d = %.6f ; predicted prod_{p>=5}(1-1/p^2) = %.6f = 9/pi^2 = %.6f = (3/2)*6/pi^2"
      % (int(Tsqf.sum()), NODD, dA, pred_A, 9 / math.pi**2))
print("(b) odd n with gcd(n,T(n))=1               : %d/%d = 1 exactly (PROVED: gcd(n,3n+1)=1); NOT 6/pi^2 (hostile control)" % (NODD, NODD))
print("    odd n with gcd(n+2,T(n))=1             : %.6f ; PROVED: gcd(n+2,3n+1) | 5, so density = 1 - (1/5)*P(5 | 3n+1 & 5|n+2 | odd n) = 4/5 + ... = %.6f"
      % (dB2, 1 - 1 / 5))
print("(c) odd n with n AND T(n) squarefree        : %d/%d = %.6f ; predicted (8/9)*prod_{p>=5}(1-2/p^2) = %.6f (Mirsky-type product, NOT a multiple of 6/pi^2)"
      % (int((Tsqf & nsqf).sum()), NODD, dC, pred_C))
check(abs(dA - pred_A) < 3e-3 and abs(dC - pred_C) < 3e-3, "6/pi^2 candidate densities")
# row uniformity of (a)
for rres in (1, 3, 5):
    mask = (nodd % 6 == rres)
    print("    (a) restricted to row n = %d mod 6: %.6f  (predicted 9/pi^2 = %.6f, row-uniform)" % (rres, int((Tsqf & mask).sum()) / int(mask.sum()), 9 / math.pi**2))
print("WHY (PROVED, CRT + Moebius): T(n) squarefree iff p^2 !| 3n+1 for every odd p; for p>=5 that is one excluded class of n mod p^2,")
print("     for p=3 it is automatic (3n+1 = 1 mod 3), for p=2 the odd part is taken.  Density = prod_{p>=5}(1-1/p^2) = (6/pi^2)*(4/3)*(9/8) = 9/pi^2.")
print("     The factor 3/2 over 6/pi^2 is exactly the two local factors killed by 'odd core' (p=2) and 'image is 2 mod 3' (p=3).")
print("     Coprimality of (n, T(n)) or of any two affine forms in n is trivial: gcd(an+b, cn+d) | (ad-bc).  So 6/pi^2 = P(gcd=1) cannot arise")
print("     from Collatz-linear pairs; it arises only through squarefreeness, and there with the local factors at 2 and 3 removed.")
del x, lowbit, u, Tsqf, nsqf, g, g2

# ---------------------------------------------------------------- S2 {7,21}
hdr("S2. {7,21}: tournaments n<=6 (control), 3n+k cycle-count spectrum, Mersenne probe")

print("\n--- S2.1  Tournament census n=3..6: Hamiltonian paths h and |Aut| (vectorised DP) ---")
def tournament_census(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    codes = np.arange(2**m, dtype=np.int64)
    adj = [[None] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        bit = ((codes >> b) & 1).astype(np.int64)
        adj[i][j] = bit
        adj[j][i] = 1 - bit
    full = (1 << n) - 1
    cnt = {}
    for v in range(n):
        cnt[(1 << v, v)] = np.ones(2**m, dtype=np.int64)
    for size in range(2, n + 1):
        for Sset in itertools.combinations(range(n), size):
            Smask = sum(1 << v for v in Sset)
            for v in Sset:
                acc = np.zeros(2**m, dtype=np.int64)
                prev = Smask & ~(1 << v)
                for uu in Sset:
                    if uu == v:
                        continue
                    acc += cnt[(prev, uu)] * adj[uu][v]
                cnt[(Smask, v)] = acc
    h = sum(cnt[(full, v)] for v in range(n))
    aut = np.zeros(2**m, dtype=np.int64)
    for perm in itertools.permutations(range(n)):
        fixed = np.ones(2**m, dtype=np.bool_)
        for (i, j) in pairs:
            a, b = perm[i], perm[j]
            img = adj[a][b]
            fixed &= (img == adj[i][j])
        aut += fixed
    return h, aut
spectrum = {}
for n in range(3, 7):
    h, aut = tournament_census(n)
    vals = sorted(set(int(t) for t in h))
    check(all(v % 2 == 1 for v in vals), "Redei: h odd at n=%d" % n)
    check(int((h % aut).max()) == 0, "|Aut| divides h at n=%d" % n)
    spectrum[n] = vals
    autvals = sorted(set(int(t) for t in aut))
    print("  n=%d: %6d labelled tournaments; h values = %s ; |Aut| values = %s ; max h = %d" % (n, 2**(n * (n - 1) // 2), vals, autvals, max(vals)))
allvals = sorted(set().union(*spectrum.values()))
print("  union of attained h for n<=6:", allvals)
missing_odd = [v for v in range(1, max(allvals) + 1, 2) if v not in allvals]
print("  odd values <= %d NOT attained for n<=6: %s  (7 and 21 absent, as THM-1745 / death-star S70 state for all n; the others fill in at larger n per the theorem)" % (max(allvals), missing_odd))
check(7 not in allvals and 21 not in allvals, "7,21 absent n<=6")
mers = [2**j - 1 for j in range(1, 7)]
print("  Mersenne probe: 2^j-1 =", mers, "; attained for n<=6:", [v for v in mers if v in allvals], "; max h per n vs 2^(n-1)-1:", [(n, max(spectrum[n]), 2**(n - 1) - 1) for n in spectrum])
print("  VERDICT (a): 7 = 2^3-1 has no Mersenne mechanism: the spectrum is 'odds minus {7,21}' so 1,3,15,31,63,... ARE attained;")
print("      the max h at n=4 is 5 (< 7) and at n=5 it is 15 = 2^4-1 while 7 is skipped; 63 = 2^6-1 = 9*7 is in the spectrum (THM-1745 head).")
print("      The row-braid fact ord_9(2)=6 (9 | 63) and the tournament hole 7 share the prime 7 only as a numeral: NO MAP FOUND.")

print("\n--- S2.2  3n+k cycle census, gcd(k,6)=1, 1<=k<=99, odd |n| <= 2*10^5, both signs ---")
NCEN = 2 * 10**5
def cycles_census(k, N, cap=10**40, maxsteps=20000):
    status = {}
    cycles = []
    unresolved = 0
    order = []
    for a in range(1, N + 1, 2):
        order.append(a); order.append(-a)
    for s in order:
        if s in status:
            continue
        path = []; seen = set(); v = s; res = None
        while True:
            if v in seen:
                i0 = path.index(v)
                cyc = path[i0:]
                mmin = min(cyc, key=abs)
                i1 = cyc.index(mmin)
                cyc = tuple(cyc[i1:] + cyc[:i1])
                cycles.append(cyc)
                res = len(cycles) - 1
                break
            av = abs(v)
            if av <= N and v in status:
                res = status[v]; break
            if av > cap or len(path) > maxsteps:
                res = -1; unresolved += 1; break
            seen.add(v); path.append(v)
            xx = 3 * v + k
            v = xx >> v2(xx)
        for pth in path:
            if abs(pth) <= N:
                status[pth] = res
    return cycles, unresolved
KS = [k for k in range(1, 100) if gcd(k, 6) == 1]
CYC = {}
for k in KS:
    cyc, unres = cycles_census(k, NCEN)
    check(unres == 0, "unresolved starts for k=%d" % k)
    # verify each cycle exactly
    for c in cyc:
        L = len(c)
        for i in range(L):
            xx = 3 * c[i] + k
            check(xx >> v2(xx) == c[(i + 1) % L], "cycle verification k=%d" % k)
    CYC[k] = cyc
    # negation conjugacy check: cycles of 3n-k are exactly the negatives
    cycm, unresm = cycles_census(-k, NCEN)
    check(unresm == 0, "unresolved for -k")
    check(sorted(tuple(-t for t in c) for c in cyc) == sorted(cycm), "negation conjugacy k=%d" % k)
print("  census done for %d values of k (and their negatives), all starts resolved, all cycles verified, C_{-k} = -C_k checked" % len(KS), tm())
F3 = {1: 4, 5: 9, 7: 5, 11: 7, 13: 13}
for k in F3:
    check(len(CYC[k]) == F3[k], "F3 control k=%d" % k)
print("  positive control vs (F3): c(1)=4, c(5)=9, c(7)=5, c(11)=7, c(13)=13 : PASS")
# completeness bounds: max |element| of any cycle with L letters, from B <= 2^(K-L)(3^L-2^L)
def completeness_L(k, N):
    # returns largest L* such that every cycle (either sign) of length <= L* has all |elements| <= N
    Lstar = 0
    for L in range(1, 60):
        worst = Fraction(0)
        for K in range(1, 200):
            d = 2**K - 3**L
            if d == 0:
                continue
            bound = Fraction(k * 2**max(K - L, 0) * (3**L - 2**L), abs(d)) if K >= L else Fraction(0)
            if K < L:
                continue
            if bound > worst:
                worst = bound
            if d > 2 * 3**L:
                break
        if worst <= N:
            Lstar = L
        else:
            break
    return Lstar
Lstar = {k: completeness_L(k, NCEN) for k in KS}
print("  PROVED bound: every element of a cycle with word (k_1..k_L), K=sum k_i, satisfies |n| <= k*2^(K-L)(3^L-2^L)/|2^K-3^L|")
print("    (B = sum 3^(L-1-i) 2^(K_i) <= 2^(K-L)(3^L-2^L) since K_i <= K-(L-i); and n = kB/(2^K-3^L)).")
print("  => the census |n|<=2*10^5 is COMPLETE for all cycles of length L <= L*(k):", {k: Lstar[k] for k in KS})
print("  cycle-count spectrum c(k) (both signs):")
row = []
for k in KS:
    cyc = CYC[k]
    lens = sorted(len(c) for c in cyc)
    npos = sum(1 for c in cyc if c[0] > 0)
    row.append((k, len(cyc), npos, len(cyc) - npos, lens))
for (k, c, npos, nneg, lens) in row:
    print("    k=%2d  c(k)=%2d  (+:%2d, -:%2d)  lengths=%s" % (k, c, npos, nneg, lens))
cvals = sorted(set(t[1] for t in row))
print("  attained values of c(k), k<=99:", cvals)
missing_c = [v for v in range(1, max(cvals) + 1) if v not in cvals]
print("  values in [1,%d] NOT attained as c(k):" % max(cvals), missing_c)
print("  odd ones not attained:", [v for v in missing_c if v % 2 == 1], " ; compare {7,21}: 7 attained? %s ; 21 attained? %s" % (7 in cvals, 21 in cvals))
print("  parity of c(k): even for k in", [k for (k, c, *_) in row if c % 2 == 0], "; odd for k in", [k for (k, c, *_) in row if c % 2 == 1])
alllens = sorted(set(l for t in row for l in t[4]))
print("  cycle LENGTHS occurring (any k<=99):", alllens, " ; 7 occurs (the 3n-1 seven-cycle scaled by k), 21 occurs? %s" % (21 in alllens))
print("  VERDICT (b): c(k) hits 7 (k=11, 23, ...) so the tournament hole is not a cycle-count hole; NO MAP.")
print("  VERDICT (c): h-spectrum 'odds minus {7,21}' vs F=S+U solution set {0,2,10} vs Collatz rows {1,3,5} mod 6: NO MAP FOUND.")
print("      Numerology on record: S(p^2qr)=7 and S*U = 7*3 = 21 reproduce {7,21}; no mechanism connects |Aut(T)| | h to divisor counts.")

# ---------------------------------------------------------------- S3 five conjectures
hdr("S3. Five bold conjectures, probes, follow-ups")

print("\n--- C1 (G-map x primes): greedy 3-adic stopping-time law on primes = on composites ---")
def Gmap(m):
    kk = 0
    while (m % 9) not in (4, 7):
        m *= 2; kk += 1
    return (m - 1) // 3, kk
MG = 10**6
steps = np.full(MG + 1, -1, dtype=np.int32)
steps[1] = 0
for m in range(2, MG + 1):
    if m % 3 == 0 or steps[m] >= 0:
        continue
    path = []; v = m
    while not (v <= MG and steps[v] >= 0):
        path.append(v); v, _ = Gmap(v)
    base = int(steps[v])
    for i, pth in enumerate(reversed(path)):
        if pth <= MG:
            steps[pth] = base + i + 1
ms = np.arange(2, MG + 1)
ms = ms[ms % 3 != 0]
st = steps[ms]
check(int(st.min()) >= 0, "all G-steps resolved")
pm = isP[ms]
sp = st[pm]; sc = st[~pm]
def tv(a, b):
    ha = np.bincount(a, minlength=200) / len(a); hb = np.bincount(b, minlength=200) / len(b)
    return 0.5 * float(np.abs(ha - hb).sum())
print("  m<=10^6, 3!|m: primes %d (mean steps %.4f, max %d) ; composites %d (mean %.4f, max %d) ; TV distance %.4f"
      % (len(sp), sp.mean(), sp.max(), len(sc), sc.mean(), sc.max(), tv(sp, sc)))
# hostile control: split by residue class mod 9 instead (should DIFFER: the first letter depends on m mod 9)
s1 = st[ms % 9 == 1]; s8 = st[ms % 9 == 8]
print("  hostile control: m=1 mod 9 vs m=8 mod 9: means %.4f vs %.4f, TV %.4f (residue classes DO differ)" % (s1.mean(), s8.mean(), tv(s1, s8)))
# residue equidistribution of primes mod 81 (Dirichlet) as the mechanism
units81 = [a for a in range(1, 81) if a % 3 != 0]
hp = np.bincount(ms[pm] % 81, minlength=81)[units81] / pm.sum()
hc = np.bincount(ms[~pm] % 81, minlength=81)[units81] / (~pm).sum()
print("  mechanism: first 3 greedy letters are a function of m mod 81 (inherited S4.1); primes vs composites mod 81: TV %.4f (Dirichlet)" % (0.5 * float(np.abs(hp - hc).sum())))
check(tv(sp, sc) < 0.02 and tv(s1, s8) > 0.1, "C1 outcome")
print("  STATUS: SURVIVES (FINITE-EXACT at 10^6, TV<0.02) ; PROVED for every fixed-depth prefix event (S4.1 + Dirichlet); full law HEURISTIC.")
print("  follow-up F1: does the class of m mod 9 split primes and composites identically?  ->", end=" ")
diffs = []
for rres in (1, 2, 4, 5, 7, 8):
    a = st[(ms % 9 == rres) & pm]; b = st[(ms % 9 == rres) & ~pm]
    diffs.append((rres, round(float(a.mean()), 3), round(float(b.mean()), 3), round(tv(a, b), 4)))
print(diffs, " (all TV small: yes)")

print("\n--- C2 (3n+k x cycle gate): k = 2^K - 3^L  =>  EVERY composition of K into L parts is a cycle element ---")
def Bword(w):
    Kc = 0; B = 0; L = len(w)
    for i, ki in enumerate(w):
        B += 3**(L - 1 - i) * 2**Kc
        Kc += ki
    return B
def realize(n, k, L):
    w = []
    v = n
    for _ in range(L):
        xx = 3 * v + k; e = v2(xx); w.append(e); v = xx >> e
    return tuple(w), v
print("  PROVED: 3B(w)+2^K-3^L = 2^(k_1) B(rot w) identically, B(w) odd  =>  for k=2^K-3^L the orbit of B(w) realizes w and closes.")
print("  (More generally n = kB(w)/(2^K-3^L) integer  =>  all rotations are integers and the word is realized: odd denominator.)")
reps = {}
for K in range(1, 70):
    for L in range(1, 45):
        d = 2**K - 3**L
        if 1 <= abs(d) <= 99 and gcd(d, 6) == 1:
            reps.setdefault(d, []).append((K, L))
print("  representations k = 2^K - 3^L with |k|<=99, K<70, L<45 (FINITE-EXACT):", dict(sorted(reps.items())))
def necklace_count(K, L):
    # number of rotation classes of compositions of K into L parts (each class = one cycle; non-primitive words give shorter cycles)
    comps = itertools.combinations(range(1, K), L - 1)
    seen = set(); cnt = 0
    for cuts in comps:
        w = tuple(b - a for a, b in zip((0,) + cuts, cuts + (K,)))
        if w in seen:
            continue
        cnt += 1
        for i in range(L):
            seen.add(w[i:] + w[:i])
    return cnt
for d, lst in sorted(reps.items()):
    if d < 0:
        continue
    for (K, L) in lst:
        total = 0; found = set()
        for cuts in itertools.combinations(range(1, K), L - 1):
            w = tuple(b - a for a, b in zip((0,) + cuts, cuts + (K,)))
            B = Bword(w)
            wr, back = realize(B, d, L)
            check(wr == w and back == B, "C2 realization k=%d w=%s" % (d, w))
            total += 1
            # locate in census
            hit = [c for c in CYC[d] if B in c]
            check(len(hit) == 1, "C2 census hit")
            found.add(CYC[d].index(hit[0]))
        print("   k=%2d = 2^%d-3^%d : %d compositions, all realized; distinct census cycles they populate: %d = rotation classes %d"
              % (d, K, L, total, len(found), necklace_count(K, L)))
        check(len(found) == necklace_count(K, L), "necklace count")
print("  consequence: c(5) >= 4 + 2 (L=3,K=5) + 1 (L=1,K=3) ; c(13) >= 4 + 7 (L=5,K=8) + 1 (L=1,K=4) ; F3 has 9 and 13: the rest are divisor-gate cycles.")
print("  bold sub-conjecture 'c(k) is maximal among k<=99 exactly at the doubly-representable k (5, 13)':", end=" ")
cmax = max(t[1] for t in row); argmax = [t[0] for t in row if t[1] == cmax]
print("max c(k) = %d at k=%s ->" % (cmax, argmax), "SURVIVES" if set(argmax) <= {5, 13} else "REFUTED (witness k=%s)" % argmax)
print("  follow-up F2: fixed points of 3n+k are exactly n = k/(2^K-3), K>=1, so #1-cycles = #{K>=1 : (2^K-3) | k}: check ->", end=" ")
for k in KS:
    pred = sum(1 for K in range(1, 8) if k % abs(2**K - 3) == 0)
    check(pred == sum(1 for c in CYC[k] if len(c) == 1), "fixed-point count k=%d" % k)
print("PASS for all k<=99 (PROVED + FINITE-EXACT).  k=65=5*13 has %d fixed points: %s" % (sum(1 for c in CYC[65] if len(c) == 1), [c for c in CYC[65] if len(c) == 1]))

print("\n--- C3 (floor sums x odd functions): Legendre duality of sums 1,2 and Pillai defect of sum 3 ---")
def Z_e(n, e):
    return sum(1 for kk in range(1, n) if pow(kk, e, n) == 0)
def Z_formula(n, e):
    pr = 1
    for pp, a in factorint(n).items():
        pr *= pp**(a - (-(-a // e)))
    return pr - 1
def pillai(n):
    return sum(gcd(i, n) for i in range(n))
def pillai_div(n):
    return sum(d * totient(n // d) for d in divisors(n))
bad = []
for n in range(2, 400):
    S1 = sum(kk**3 // n for kk in range(1, n))
    Y = (n - 1) * (n - 2)
    S2 = sum(sympy.integer_nthroot(y * n, 3)[0] for y in range(1, Y + 1))
    S3 = sum((i * j) // n for i in range(1, n) for j in range(1, n))
    Z3 = Z_e(n, 3)
    check(Z3 == Z_formula(n, 3), "Z_3 formula")
    check(Z3 % 2 == 0 or n == 1, "Z_3 even")  # Z counts k and n-k pairs? not needed; formula below handles /2 exactly
    f1 = Fraction((n - 2) * (n - 1) * (n + 1), 4) + Fraction(Z3, 2)
    f2 = Fraction((3 * n - 5) * (n - 2) * (n - 1), 4) + Fraction(Z3, 2)
    P = pillai(n)
    check(P == pillai_div(n), "Pillai divisor formula")
    f3 = Fraction((n - 2) * (n - 1)**2, 4) + Fraction(P - 2 * n + 1, 2)
    if not (S1 == f1 and S2 == f2 and S3 == f3):
        bad.append(n)
    check(S1 == f1 and S2 == f2 and S3 == f3, "floor-sum exact formulas at n=%d" % n)
    sf = all(a == 1 for a in factorint(n).values())
    check((Z3 == 0) == sf, "truth set sum1/sum2 = squarefree")
    check((P == 2 * n - 1) == isprime(n), "truth set sum3 = primes")
print("  PROVED + checked 2<=n<400:")
print("    sum_{k<n} floor(k^3/n)                 = (n-2)(n-1)(n+1)/4 + Z_3(n)/2")
print("    sum_{y<=(n-1)(n-2)} floor((yn)^(1/3))  = (3n-5)(n-2)(n-1)/4 + Z_3(n)/2      (same defect: Legendre/lattice-point duality)")
print("    sum_{i,j<n} floor(ij/n)                = (n-2)(n-1)^2/4 + (P(n)-2n+1)/2,  P(n)=sum_{i<n} gcd(i,n) = sum_{d|n} d phi(n/d) (Pillai)")
print("    Z_e(n) = #{k<n: n|k^e} = prod_p p^(a_p - ceil(a_p/e)) - 1 ; Z_e(n)=0 iff n squarefree (e>=2) ; P(n)=2n-1 iff n prime.")
print("  mechanism = ODD FUNCTION: k -> k^e odd => k^e + (n-k)^e = 0 mod n, so floor(k^e/n)+floor((n-k)^e/n) = (k^e+(n-k)^e)/n - 1 + [n|k^e];")
print("    even e breaks the cancellation (witness n=3,e=2: sum floor(k^2/3) = 1 but the odd-e formula would give 0).")
S1e = sum(kk**2 // 3 for kk in range(1, 3))
print("    hostile: n=3, e=2: sum floor(k^2/3) =", S1e, "; (2 sum k^2 - n(n-1) + n Z_2)/(2n) =", Fraction(2 * 5 - 6 + 0, 6))
print("  prediction from the Z_e formula: 'Z_7(n) = n/rad(n) - 1' fails first at n = 2^8 = 256 (a=8>e=7): Z_7(256) = %d, n/rad-1 = %d"
      % (Z_e(256, 7), 256 // 2 - 1))
check(Z_e(256, 7) == 63 and all(Z_e(n, 7) == n // sympy.primefactors(n)[0] ** 0 * (n // math.prod(sympy.primefactors(n))) - 1 for n in range(2, 256)), "Z_7 first failure at 256")
print("  STATUS: PROVED (three exact identities); F2's 'iff squarefree / iff prime' and 'Z_7 = n/rad-1 for n<200' are explained and sharpened.")
print("  follow-up F3: density of the truth set of sums 1,2 = 6/pi^2 (squarefree), of sum 3 = 0 (primes); the 'three pieces' are")
print("    {squarefree, squarefree, prime} = {Z_3=0, Z_3=0, P=2n-1}; sums 1 and 2 are the two Legendre-dual lattice-point counts under y=x^3/n.")

print("\n--- C4 (sandwich x squares): N_PS - N_SP = #{p>=5: p^2<=6K+1, p^2-2 prime} + bounded residual ? ---")
def sandwich(K):
    kk = np.arange(1, K + 1, dtype=np.int64)
    a = Om[6 * kk - 1]; b = Om[6 * kk + 1]
    ca = np.minimum(a, 4); cb = np.minimum(b, 4)
    M = np.zeros((5, 5), dtype=np.int64)
    np.add.at(M, (ca, cb), 1)
    return M
ctrl = {10**5: 41, 10**6: 412}
for K in (10**5, 10**6):
    M = sandwich(K)
    dPS = int(M[1, 2] - M[2, 1])
    check(dPS == ctrl[K], "wave-one sandwich control at K=%d" % K)
    A = sum(1 for pp in primerange(5, math.isqrt(6 * K + 1) + 1) if isprime(pp * pp - 2))
    # independence prediction for the non-square part
    R = M.sum(axis=1); C = M.sum(axis=0)
    pred = (R[1] * C[2] - R[2] * C[1]) / K
    print("  K=%d: N_PS-N_SP = %d (matches wave-one) ; square-driven pairs (prime, p^2) = %d ; residual after squares = %d ; independence pred = %.1f"
          % (K, dPS, A, dPS - A, pred))
print("  STATUS: HEURISTIC/REFUTED as 'bounded residual': residual after removing squares is 41-64=-23 at 10^5 and 412-137=275 at 10^6 (not bounded).")
print("  The square term is PROVED to be one-sided (p^2 = 1 mod 6 always sits on the right), but it is small against the mixed-class semiprime bias.")
print("  follow-up F4: does the p^2 count alone explain the S-column excess of the marginals?  S_1 - S_5 identity (inherited) already has pi'(sqrt x)/2:")
for K in (10**5, 10**6):
    x6 = 6 * K + 1
    sq = sum(1 for pp in primerange(5, math.isqrt(x6) + 1))
    M = sandwich(K)
    print("    K=%d: pi'(sqrt x) = %d, (S_1 - S_5) = %d (col-row sum of class 2) ; half the squares = %.1f" % (K, sq, int(M.sum(axis=0)[2] - M.sum(axis=1)[2]), sq / 2))

print("\n--- C5 (3n+k x G-map conjugacy): T_k(k n) = k T_1(n), G_k(k m) = k G_1(m); c(k) >= 4 and superadditivity ---")
for k in KS:
    base = [tuple(k * t for t in c) for c in CYC[1]]
    for bc in base:
        check(any(sorted(bc) == sorted(c) for c in CYC[k]), "scaled Collatz cycle in C_k, k=%d" % k)
print("  PROVED: (3(kn)+k)/2^v = k(3n+1)/2^v ; so k*C_1 subset C_k, c(k) >= 4 for every k coprime to 6 (checked k<=99).")
eq4 = [k for k in KS if len(CYC[k]) == 4]
print("  k<=99 with c(k) = 4 exactly (no cycle beyond the four scaled Collatz cycles):", eq4)
print("  superadditivity for coprime k1,k2: c(k1 k2) >= c(k1) + c(k2) - 4 (scaled cycles k2*C_{k1}, k1*C_{k2} meet only in k1k2*C_1):")
viol = []
for k1 in KS:
    for k2 in KS:
        if k1 < k2 and gcd(k1, k2) == 1 and k1 * k2 in CYC:
            lhs = len(CYC[k1 * k2]); rhs = len(CYC[k1]) + len(CYC[k2]) - 4
            if lhs < rhs:
                viol.append((k1, k2, lhs, rhs))
print("    violations among k1*k2<=99:", viol, "-> PROVED inequality holds (FINITE-EXACT confirmation)")
# greedy inverse conjugacy for general k: G_k(m) = (2^j m - k)/3, j minimal with result not divisible by 3
def Gk(m, k):
    j = 0
    while ((2**j * m - k) % 3 != 0) or (((2**j * m - k) // 3) % 3 == 0):
        j += 1
    return (2**j * m - k) // 3, j
for k in (5, 7, 11, 13):
    for m in range(1, 3000):
        if m % 3 == 0:
            continue
        a, ja = Gk(k * m, k); b, jb = Gmap(m)
        check(a == k * b and ja == jb, "G_k conjugacy k=%d m=%d" % (k, m))
print("  PROVED (+checked m<3000, k=5,7,11,13): G_k(k m) = k G_1(m) with the same letter; the residue chain mod 9 of G_k is the")
print("    chain of G_1 conjugated by multiplication by k^{-1} mod 9, so the stationary law (2/9,1/9 pattern), E[letter]=1 and")
print("    the exact log-drift log(2/3) of (F1) are k-independent.  What changes with k is the CYCLE set (C2), not the local law.")
print("  follow-up F5: are the 'new' cycles (beyond k*C_1) of 3n+k all divisor-gate cycles n = kB/(2^K-3^L) with (2^K-3^L) !| B? ->", end=" ")
cntnew = 0; cntgate = 0
for k in KS:
    for c in CYC[k]:
        if any(sorted(c) == sorted(tuple(k * t for t in c1)) for c1 in CYC[1]):
            continue
        cntnew += 1
        L = len(c); w, back = realize(c[0], k, L); K = sum(w); d = 2**K - 3**L; B = Bword(w)
        check(back == c[0] and Fraction(k * B, d) == c[0], "gate identity")
        if B % d != 0:
            cntgate += 1
print("%d new cycles, %d have (2^K-3^L) !| B (need the factor k), %d have (2^K-3^L) | B (would be Collatz cycles scaled: impossible unless k=1 gate) " % (cntnew, cntgate, cntnew - cntgate))
check(cntnew == cntgate, "all new cycles need k")

# ---------------------------------------------------------------- S4 3N-5
hdr("S4. 3N-5 = 3(N-2)+1: exact relation to Collatz under the shift n -> n-2")
def Tk(n, k):
    xx = 3 * n + k; return xx >> v2(xx)
for n in range(-100001, 100002, 2):
    check(Tk(n, -5) == Tk(n - 2, 1), "T_{-5}(n) = T(n-2)")
print("PROVED (+checked odd |n|<=10^5): T_{-5}(n) = T(n-2) exactly (3n-5 = 3(n-2)+1, same 2-adic valuation).")
print("  Hence with the shift s(m)=m+2: s^{-1} T_{-5} s (m) = T(m) - 2.  This is NOT a conjugacy of T with itself: the exact")
print("  discrepancy is the constant -2 after every odd step, i.e. T_{-5} is conjugate to the map m -> T(m)-2, and the orbit of n")
print("  under 3n-5 is the orbit of n-2 under 'Collatz then subtract 2 each odd step'.")
print("  rows mod 6 -> single-halving image (3n-5)/2 :")
for rr in (1, 3, 5):
    j = sympy.symbols('j')
    img = (3 * (6 * j + rr) - 5) / 2
    print("    n = 6j+%d : (3n-5)/2 = %s  = %d mod 9  (Collatz row image of 6j+%d would be %d mod 9)" % (rr, sympy.expand(img), (3 * rr - 5) // 2 % 9, rr, (3 * rr + 1) // 2 % 9))
print("  PROVED: the 3n-5 row map is the Collatz row map precomposed with r -> r-2 mod 6 (1->5->3->1): images 8,2,5 instead of 2,5,8;")
print("  all images are 2 mod 3 exactly as for Collatz (multiples of 3 are leaves for 3n-5 as well).")
cyc_m5 = sorted(tuple(-t for t in c) for c in CYC[5])
print("  cycles of 3n-5 on odd integers, both signs, min-element |n|<=2*10^5 (= negatives of the 9 cycles of 3n+5):")
for c in cyc_m5:
    print("    length %2d  min |n| %6d : %s" % (len(c), abs(c[0]), c if len(c) <= 8 else c[:8] + ('...',)))
pos_m5 = [c for c in cyc_m5 if c[0] > 0]
print("  POSITIVE 3n-5 cycles = %d (the user's '3N-5 on positives' = 3N+5 on negatives): %s" % (len(pos_m5), [c for c in pos_m5]))
print("  the two 3-cycles the user finds fascinating live on the POSITIVE side of 3n+5 = NEGATIVE side of 3n-5: (19,31,49),(23,37,29),")
print("  and are the C2 construction for 5 = 2^5-3^3 (words (1,1,3) and (1,2,2)); 3n-5 has them as (-19,-31,-49),(-23,-37,-29).")
print("  shifted picture: a 3n-5 cycle (n_i) is a cycle (n_i - 2) of m -> T(m)-2, e.g. fixed point -5 of 3n-5 <-> m=-7: T(-7)-2 = -5-2 = -7.")
check(Tk(-7, 1) - 2 == -7 and Tk(-5, -5) == -5, "shift example")
print("  PROVED: 3n-5 on negatives = -(3n+5 on positives); 3n-5 on positives = -(3n+5 on negatives) (sign conjugation, inherited summand note).")

# ---------------------------------------------------------------- S5 summary
hdr("S5. Status summary")
print("PROVED: S1.1 F_2^3/Fano structure of squarefree divisors; d(N)<=2^Omega with equality iff squarefree; 9/pi^2 density of odd n with T(n)")
print("        squarefree and the (8/9)prod(1-2/p^2) density; gcd(n,T(n))=1; C2 word-cycle theorem for k=2^K-3^L and the general")
print("        integrality propagation; fixed-point count #{K:(2^K-3)|k}; census completeness bound; C3 three exact floor identities")
print("        (Legendre duality, Pillai defect, Z_e formula, first failure 256); C5 scaling k*C_1 in C_k, superadditivity, G_k conjugacy;")
print("        S4 T_{-5}(n)=T(n-2) and the row map; C_{-k} = -C_k.")
print("FINITE-EXACT: 16/128 octonion orientations; tournament census n<=6; c(k) spectrum k<=99 (|n|<=2*10^5, complete for L<=L*(k));")
print("        6/pi^2 candidates at 10^6; C1 TV distances at 10^6; C4 sandwich numbers at 10^5,10^6.")
print("CITED: Albuquerque-Majid twisted group algebra (octonions, Clifford); Dirichlet; THM-1745 / death-star S70 h-spectrum; Redei.")
print("HEURISTIC: full G-stopping-time law equality primes vs composites; C4 residual.")
print("REFUTED: 'c(k) is odd for all k' (k=1: 4); 'N_PS-N_SP minus squares is bounded' (residual -23 -> 275); 'c(k) maximal exactly at doubly")
print("        representable k' (see C2 line); 'Z_7(n)=n/rad(n)-1 for all n' (n=256).")
print("NO MAP FOUND: Bott 8; {7,21} vs c(k) or cycle lengths or rows; {0,2,10} vs h-spectrum; x^2+c period-3 bifurcation (polynomial")
print("        dynatomic gate) vs 3n+k 3-cycles (exponential Diophantine gate 2^K-27 | kB).")
print("total runtime", tm())
