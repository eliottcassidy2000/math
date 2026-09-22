#!/usr/bin/env python3
"""
collatz_mod6_20260917_wild_typing.py  --  LANE wild_typing (creative lane, every claim probed)

Session collatz-mod6-20260917 (mac-mini), wave two.  Recovered from the 2026-09-17 agent transcript on
2026-09-21 and finalized: the recovered draft's C1 hostile control ("m=1 mod 9 vs m=8 mod 9 differ") was
FALSE (TV 0.007); the true picture (size mixing, and the residue split 4 vs 5 mod 9) is asserted below.
Sections S2.3 (Cipolla trunk), C3 (trunk squarefreeness / Wieferich), C6 (pseudoprime index set, tower)
and the two-universe census of S2.2 are new relative to the draft.

Exact integer / Fraction arithmetic for every load-bearing claim; explicit universes; positive and hostile
controls; every check is an explicit `raise` (active under python -O).

Sections
  S0  helpers, sieves, universes
  S1  TYPING TABLE probes: Fano/octonion orientation count, divisor box 12=8+3+1,
      Cayley-Dickson vs Omega/omega, Bott-8 probe, 6/pi^2 candidates (three exact densities), typing lines
  S2  {7,21}: tournament census n<=6 (h, |Aut|); 3n+b cycle-count spectrum c(b), gcd(b,6)=1, |b|<=99,
      universe U1 = starts |n|<=10^5 (escape 10^18) and U2 = |n|<=2*10^5 (escape 10^40); Cipolla trunk
  S3  six bold conjectures with probes + one follow-up each
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
Om = np.zeros(NS + 1, dtype=np.uint8)
spf_primes = list(primerange(2, NS + 1))
for p in spf_primes:
    q = p
    while q <= NS:
        Om[q::q] += 1
        q *= p
isP = np.zeros(NS + 1, dtype=np.bool_)
isP[spf_primes] = True
print("sieve bound NS =", NS, " primes:", len(spf_primes))
rng = np.random.default_rng(20260917)
for n in rng.integers(2, NS, 400):
    n = int(n)
    f = factorint(n)
    check(int(Om[n]) == sum(f.values()), "Omega sieve at %d" % n)
    check(bool(sqf[n]) == all(e == 1 for e in f.values()), "sqf sieve at %d" % n)
    check(bool(isP[n]) == isprime(n), "prime sieve at %d" % n)
print("sieve audits vs sympy on 400 random n: PASS")
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
# the three F=S+U shapes (inherited DB1) and their F values {0,2,10}
def FSU(n0):
    dd = divisors(n0); pr = [d for d in dd if 1 < d < n0]
    return (len(pr), sum(1 for d in pr if all(e == 1 for e in factorint(d).values())), sum(1 for d in pr if isprime(d)))
shapes = {"p": 7, "p^3": 7**3, "p^2qr": 7 * 7 * 11 * 13}
for nm, n0 in shapes.items():
    f0, s0, u0 = FSU(n0)
    check(f0 == s0 + u0, "F=S+U at " + nm)
    print("   shape %-6s N=%6d : F=%2d S=%d U=%d" % (nm, n0, f0, s0, u0))
check([FSU(n0)[0] for n0 in shapes.values()] == [0, 2, 10], "{0,2,10}")
print("   the user's {0,2,10} = F at the three DB1 shapes p, p^3, p^2qr (inherited arithmetic_braids_20260917_divisors.md DB1).")
nonsqf_proper = [d for d in proper if d not in S_set]
print("divisor box [0,2]x[0,1]^2 has", len(D), "points = 8 (Boolean cube {squarefree divisors incl. 1}) + %d (nonsquarefree proper: %s) + 1 (N)"
      % (len(nonsqf_proper), nonsqf_proper))
check(len(D) == 12 and len(nonsqf_proper) == 3, "12 = 8+3+1")
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
for a in S_set:
    check(sum(1 for l in lines if a in l) == 3, "3 lines per point")
print("PROVED: 7 points, 7 lines, 3 points/line, 3 lines/point = PG(2,2).")

idx = {d: vecint(d) for d in S_set}      # 1..7
line_vecs = [tuple(sorted(idx[d] for d in l)) for l in lines]
def build_mult(orient):
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
print("         16 = 128/8 is rank-nullity of the incidence map (inherited arithmetic_braids2_20260917_fano_code.md sec.1); re-derived here")
print("         by brute alternativity as an independent path.  Fano points 7 = S(p^2qr), dim 3 = U(p^2qr).")

print("\n--- S1.2  Cayley-Dickson dimensions vs Omega / omega ---")
tbl = []
for name, n0 in (("a=2", 2), ("b=6", 6), ("a^2b=24", 24), ("c^2b=5400", 5400)):
    f = factorint(n0)
    Omg = sum(f.values()); omg = len(f); dn = len(divisors(n0))
    print("   %-10s Omega=%2d omega=%d d(N)=%3d 2^omega=%d  2^Omega=%d  divisor box = %s points" % (name, Omg, omg, dn, 2**omg, 2**Omg, dn))
print("PROVED: Omega is additive (Omega(MN)=Omega(M)+Omega(N)), Cayley-Dickson dimension is multiplicative (dim doubles),")
print("        so 'Omega(a^2 b)=4 <-> dim H = 4' and 'Omega(c^2 b)=8 <-> dim O = 8' compare 2*Omega(a)+Omega(b) = 2+2, 6+2 with 2^2, 2^3:")
print("        the equalities 4=4 and 8=8 are two coincidences of small numbers; the next member 'Omega(d^2 b) = 10' (d a 4-almost-prime)")
print("        has no Cayley-Dickson partner (dimension 16 = sedenions).  Inherited DB3: literal C^2B is the 8-almost-primes, not p^2qr.")
print("        d(N) = prod(a_i+1) <= 2^Omega(N) with equality iff N squarefree; and d(N) = 2^omega(N) iff N squarefree.")
for n0 in range(2, 2001):
    f = factorint(n0)
    Omg = sum(f.values()); omg = len(f); dn = len(divisors(n0))
    sf = all(e == 1 for e in f.values())
    check((dn == 2**Omg) == sf and (dn == 2**omg) == sf and dn <= 2**Omg, "d(N) vs 2^Omega at %d" % n0)
print("        (checked 2<=N<=2000).  The only functorial match: for SQUAREFREE N with omega(N)=j the squarefree-divisor group is F_2^j,")
print("        its real group algebra R[F_2^j] has dimension 2^j = d(N), and the Cayley-Dickson algebra of that dimension is a")
print("        twisted group algebra R_F[F_2^j] (CITED: Albuquerque-Majid, 'Quasialgebra structure of the octonions', J. Algebra 220 (1999)).")
print("        Map: N=1->R, p->C, pq->H, pqr->O (omega = 0,1,2,3).  LOST: the cocycle F (signs/orientation).  a^2b and c^2b are NOT squarefree,")
print("        their divisor boxes (8 and 48 points above) are not groups: VERDICT numerology for the Omega version.")

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
print("R-braid period on rows = ord_9(4) = %d; source period mod 42 and target mod 63: 3 (inherited row_braid_typing RB4)." % mult_order(4, 9))
print("VERDICT: NO MAP FOUND.  8 = Bott period is the Morita period of the Clifford cocycle on F_2^n (Cl_{n+8} = Cl_n (x) M_16(R));")
print("         the shared object with S1.1 is only 'twisted group algebra of F_2^n' (Clifford AND octonion are such: CITED Albuquerque-Majid).")
print("         No Collatz modulus has ord(2)=8; the row/braid periods are 2,3,6.  Nothing in the divisor lattice carries a cocycle.")
print("         The one '8' with a mechanism is 2(r+1)=2^r iff r=3 (self-dual RM(1,3)): inherited fano_code.md sec.5, not Bott.")

print("\n--- S1.4  6/pi^2 candidates: exact densities at 10^6 odd n ---")
NODD = 10**6
nodd = np.arange(1, 2 * NODD, 2, dtype=np.int64)
x = 3 * nodd + 1
lowbit = x & (-x)
u = x // lowbit
check(int((u % 2).min()) == 1 and int((u % 3 == 0).sum()) == 0, "odd cores odd and coprime to 3")
Tsqf = sqf[u]
nsqf = sqf[nodd]
dA = int(Tsqf.sum()) / NODD
dC = int((Tsqf & nsqf).sum()) / NODD
g = np.gcd(nodd, u)
check(int(g.max()) == 1, "gcd(n,T(n))=1")
g2 = np.gcd(nodd + 2, u)
dB2 = int((g2 == 1).sum()) / NODD
P5 = 1.0; P5b = 1.0
for pp in primerange(5, 10**7):
    P5 *= (1 - 1 / pp**2)
    P5b *= (1 - 2 / pp**2)
pred_A = P5
pred_C = (8 / 9) * P5b
print("(a) odd n with T(n) squarefree             : %d/%d = %.6f ; predicted prod_{p>=5}(1-1/p^2) = %.6f = 9/pi^2 = %.6f = (3/2)*6/pi^2"
      % (int(Tsqf.sum()), NODD, dA, pred_A, 9 / math.pi**2))
print("(b) odd n with gcd(n,T(n))=1               : %d/%d = 1 exactly (PROVED: gcd(n,3n+1)=1); NOT 6/pi^2 (hostile control)" % (NODD, NODD))
print("    odd n with gcd(n+2,T(n))=1             : %.6f ; PROVED: gcd(n+2,3n+1) | 5, so density = 1 - 1/5 = %.6f"
      % (dB2, 1 - 1 / 5))
check(abs(dB2 - 0.8) < 1e-3, "gcd(n+2,T(n)) density 4/5")
print("(c) odd n with n AND T(n) squarefree        : %d/%d = %.6f ; predicted (8/9)*prod_{p>=5}(1-2/p^2) = %.6f (Mirsky-type product, NOT a multiple of 6/pi^2)"
      % (int((Tsqf & nsqf).sum()), NODD, dC, pred_C))
check(abs(dA - pred_A) < 3e-3 and abs(dC - pred_C) < 3e-3, "6/pi^2 candidate densities")
for rres in (1, 3, 5):
    mask = (nodd % 6 == rres)
    print("    (a) restricted to row n = %d mod 6: %.6f  (predicted 9/pi^2 = %.6f, row-uniform)" % (rres, int((Tsqf & mask).sum()) / int(mask.sum()), 9 / math.pi**2))
print("WHY (PROVED, CRT + Moebius): T(n) squarefree iff p^2 !| 3n+1 for every odd p; for p>=5 that is one excluded class of n mod p^2,")
print("     for p=3 it is automatic (3n+1 = 1 mod 3), for p=2 the odd part is taken.  Density = prod_{p>=5}(1-1/p^2) = (6/pi^2)*(4/3)*(9/8) = 9/pi^2.")
print("     The factor 3/2 over 6/pi^2 is exactly the two local factors killed by 'odd core' (p=2) and 'image is 2 mod 3' (p=3).")
print("     Coprimality of (n, T(n)) or of any two affine forms in n is trivial: gcd(an+b, cn+d) | (ad-bc).  So 6/pi^2 = P(gcd=1) cannot arise")
print("     from Collatz-linear pairs; it arises only through squarefreeness, and there with the local factors at 2 and 3 removed.")
print("     Negative n: T(-n) = -T_{-1}(n) (inherited summand.md sec.7) swaps rows 1 and 5 and preserves |.|-squarefreeness: the same densities.")
del x, lowbit, u, Tsqf, nsqf, g, g2

# ---------------------------------------------------------------- S2 {7,21}
hdr("S2. {7,21}: tournaments n<=6 (control), 3n+b cycle-count spectrum, Mersenne probe, Cipolla trunk")

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

print("\n--- S2.2  3n+b cycle census, gcd(b,6)=1, 1<=b<=99, both signs of n; U1: |n|<=10^5 escape 10^18; U2: |n|<=2*10^5 escape 10^40 ---")
def cycles_census(k, N, cap, maxsteps=20000):
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
def run_universe(N, cap, label):
    CY = {}
    unres_total = 0
    for k in KS:
        cyc, unres = cycles_census(k, N, cap)
        unres_total += unres
        for c in cyc:
            L = len(c)
            for i in range(L):
                xx = 3 * c[i] + k
                check(xx >> v2(xx) == c[(i + 1) % L], "cycle verification k=%d" % k)
        CY[k] = cyc
        cycm, unresm = cycles_census(-k, N, cap)
        unres_total += unresm
        check(sorted(tuple(-t for t in c) for c in cyc) == sorted(cycm), "negation conjugacy k=%d" % k)
    print("  %s: census done for %d values of b (and their negatives); unresolved (escaped) starts: %d; all cycles verified; C_{-b} = -C_b checked" % (label, len(KS), unres_total))
    return CY, unres_total
CYC1, unres1 = run_universe(10**5, 10**18, "U1 (|n|<=10^5, escape 10^18)")
CYC, unres2 = run_universe(2 * 10**5, 10**40, "U2 (|n|<=2*10^5, escape 10^40)")
check(unres1 == 0 and unres2 == 0, "no escaped starts in either universe")
F3 = {1: 4, 5: 9, 7: 5, 11: 7, 13: 13}
for k in F3:
    check(len(CYC[k]) == F3[k], "F3 control k=%d" % k)
print("  positive control vs the session-lead census F3 (|n|<=2*10^5): c(1)=4, c(5)=9, c(7)=5, c(11)=7, c(13)=13 : PASS")
print("  (c(5)=9 = the nine b=-5 cycles of arithmetic_braids2_20260917_signed_cycles.md sec.5 under negation; c(1)=4 = the four known signed Collatz cycles.)")
def completeness_L(k, N):
    Lstar = 0
    for L in range(1, 60):
        worst = Fraction(0)
        for K in range(1, 200):
            d = 2**K - 3**L
            if d == 0 or K < L:
                continue
            bound = Fraction(k * 2**(K - L) * (3**L - 2**L), abs(d))
            if bound > worst:
                worst = bound
            if d > 2 * 3**L:
                break
        if worst <= N:
            Lstar = L
        else:
            break
    return Lstar
Lstar1 = {k: completeness_L(k, 10**5) for k in KS}
Lstar2 = {k: completeness_L(k, 2 * 10**5) for k in KS}
print("  PROVED bound: every element of a cycle with word (k_1..k_L), K=sum k_i, satisfies |n| <= b*2^(K-L)(3^L-2^L)/|2^K-3^L|")
print("    (B = sum 3^(L-1-i) 2^(K_i) <= 2^(K-L)(3^L-2^L) since K_i <= K-(L-i); and n = bB/(2^K-3^L)).")
print("  => U1 is COMPLETE for all cycles of length L <= L*_1(b):", Lstar1)
print("  => U2 is COMPLETE for all cycles of length L <= L*_2(b):", Lstar2)
print("  cycle-count spectrum c(b) (both signs of n):")
row = []
diffU = []
for k in KS:
    cyc = CYC[k]
    lens = sorted(len(c) for c in cyc)
    npos = sum(1 for c in cyc if c[0] > 0)
    row.append((k, len(cyc), npos, len(cyc) - npos, lens))
    if len(CYC1[k]) != len(cyc):
        diffU.append((k, len(CYC1[k]), len(cyc), [c for c in cyc if abs(c[0]) > 10**5]))
for (k, c, npos, nneg, lens) in row:
    print("    b=%2d  c(b)=%2d  (+:%2d, -:%2d)  lengths=%s" % (k, c, npos, nneg, lens))
print("  b with c_U1(b) != c_U2(b) (cycles whose least |element| lies in (10^5, 2*10^5]):",
      [(k, a, b2, [(len(c), c[0]) for c in extra]) for (k, a, b2, extra) in diffU])
for (k, a, b2, extra) in diffU:
    check(a < b2 and all(abs(c[0]) > 10**5 for c in extra), "U1 subset U2")
cvals = sorted(set(t[1] for t in row))
cvals1 = sorted(set(len(CYC1[k]) for k in KS))
print("  attained values of c(b), |b|<=99, U2:", cvals)
print("  attained values of c(b), |b|<=99, U1:", cvals1)
missing_c = [v for v in range(1, max(cvals) + 1) if v not in cvals]
missing_c1 = [v for v in range(1, max(cvals1) + 1) if v not in cvals1]
print("  values in [1,%d] NOT attained as c(b) (U2):" % max(cvals), missing_c)
print("  values in [1,%d] NOT attained as c(b) (U1):" % max(cvals1), missing_c1)
print("  odd ones not attained (U2):", [v for v in missing_c if v % 2 == 1], " (U1):", [v for v in missing_c1 if v % 2 == 1],
      " ; compare {7,21}: 7 attained? %s ; 21 attained? %s (max c(b) = %d, so 21 is above the range, not a hole)" % (7 in cvals, 21 in cvals, max(cvals)))
check(7 in cvals and 7 in cvals1 and 21 not in cvals and max(cvals) < 21, "c(b) vs {7,21}")
check([v for v in missing_c if v % 2 == 1] == [1, 3], "odd holes of c(b) in U2 are exactly 1,3")
print("  b with c(b)=7:", [k for (k, c, *_) in row if c == 7])
print("  parity of c(b): even for b in", [k for (k, c, *_) in row if c % 2 == 0], "; odd for b in", [k for (k, c, *_) in row if c % 2 == 1])
alllens = sorted(set(l for t in row for l in t[4]))
print("  cycle LENGTHS occurring (any b<=99):", alllens, " ; 7 occurs (the 3n-1 seven-cycle scaled by b), 21 occurs? %s" % (21 in alllens))
print("  VERDICT (b): c(b) hits 7 (b=11, 37, 49, 73, 83) so the tournament hole is not a cycle-count hole; c(b) never equals 1 or 3 because")
print("      c(b) >= 4 always (C5 below: b*C_1 subset C_b), and 21 is simply beyond max c(b)=19 in this universe: NO MAP.")
print("  VERDICT (c): h-spectrum 'odds minus {7,21}' vs F=S+U solution set {0,2,10} vs Collatz rows {1,3,5} mod 6: NO MAP FOUND.")
print("      Numerology on record: S(p^2qr)=7 and S*U = 7*3 = 21 reproduce {7,21}; no mechanism connects |Aut(T)| | h to divisor counts.")

print("\n--- S2.3  Cipolla pseudoprimes on the trunk: the inverse fibre of 1 is N_j = (4^j-1)/3 = R^{j-1}(1), R(n)=4n+1 ---")
def Nj(j):
    return (4**j - 1) // 3
R = lambda n: 4 * n + 1
v = 1
for j in range(1, 41):
    check(v == Nj(j), "trunk = R-orbit of 1 at j=%d" % j)
    xx = 3 * v + 1
    check(xx == 4**j and (xx >> v2(xx)) == 1, "T(N_j)=1 with halving exponent 2j")
    v = R(v)
print("  PROVED (+checked j<=40): N_j = R^{j-1}(1) (inherited (B1) with u=1, h0=1), 3N_j+1 = 4^j, so T(N_j)=1 with halving exponent exactly 2j;")
print("  rows mod 6 of N_1..N_9:", [Nj(j) % 6 for j in range(1, 10)], "(period 3: 1->5->3, inherited RB4).")
print("  j   N_j                    factorization                 primitive primes of 4^j-1     base-2 psp?  4^j=4 mod 6j?")
psp_j = []
for j in range(1, 41):
    Nv = Nj(j)
    f = factorint(Nv) if j <= 40 else {}
    prim = [pp for pp in f if all((4**i - 1) % pp != 0 for i in range(1, j))]
    comp = (Nv > 1) and (not isprime(Nv))
    psp = comp and pow(2, Nv - 1, Nv) == 1
    crit = (j >= 3) and (pow(4, j, 6 * j) == 4 % (6 * j))
    check(psp == crit or j < 3, "criterion 4^j=4 mod 6j at j=%d" % j)
    check(j == 1 or len(prim) >= 1, "Zsigmondy: primitive prime of 4^j-1 at j=%d" % j)
    check(j == 1 or all(pp != 3 for pp in prim), "primitive primes are never 3")
    if psp:
        psp_j.append(j)
    fs = "*".join(("%d^%d" % (pp, e) if e > 1 else "%d" % pp) for pp, e in sorted(f.items())) if Nv > 1 else "1"
    print("  %2d  %-22s %-29s %-29s %-12s %s" % (j, Nv if Nv < 10**21 else "~2^%d" % (2 * j - 2), fs if len(fs) <= 29 else fs[:26] + "...", prim, psp, crit))
print("  j<=40 with N_j a base-2 pseudoprime:", psp_j)
check(psp_j == [pp for pp in range(5, 41) if isprime(pp)], "j<=40: psp exactly at primes >= 5")
print("  CITED (Cipolla 1904): (4^p-1)/3 is a base-2 pseudoprime for every prime p>=5.  Elementary PROOF: N=(4^p-1)/3; 2^{2p}=4^p=3N+1=1 mod N;")
print("    N-1 = 4(4^{p-1}-1)/3 is even, and p | 4^{p-1}-1 (Fermat) with p !| 3, so p | N-1; hence 2p | N-1 and 2^{N-1} = (2^{2p})^{(N-1)/2p} = 1 mod N.")
print("    N is composite: N = (2^p-1)*((2^p+1)/3) with both factors >1 for p>=3.  p=3 fails exactly because 3 | 3 (N_3=21: 2^20 = 4 mod 21).")
print("  PROVED (exact criterion, all j>=3): ord_{N_j}(2) = 2j (it divides 2j; a proper divisor d<=j would give N_j | 2^j-1 < N_j),")
print("    so N_j is a base-2 pseudoprime iff N_j is composite (true for j>=3) and 2j | N_j-1 = (4^j-4)/3, i.e. iff 6j | 4^j-4.")
print("    Corollary: N_j base-2 psp iff N_j base-4 psp (same divisibility j | N_j-1 with ord_{N_j}(4)=j).")
print("  REFUTED (user's pasted claim '341 stalls primitive prime generation of 4x+1'): 341 = N_5 = 11*31 brings TWO new primes (both primitive),")
print("    and Zsigmondy for 4^j-1 has NO exception (4-1=3 != 1; 4+1=5 not a power of 2; (2,1,6) is base 2 not 4): every N_j, j>=2, has a")
print("    prime dividing no earlier trunk term (checked j<=40 above; primitive primes are never 3 since 3 | 4-1).")
Wset = [j for j in range(3, 20001) if pow(4, j, 6 * j) == 4 % (6 * j)]
Wcomp = [j for j in Wset if not isprime(j)]
check(all(j in set(Wset) for j in range(5, 20001) if isprime(j)), "all primes>=5 in W")
check(all(j % 3 != 0 for j in Wset), "3 !| j for j in W")
print("  W = {j>=3 : 6j | 4^j-4} = {j : N_j is a base-2 pseudoprime}, j<=20000: |W| = %d = %d primes >= 5 + %d composite indices" % (len(Wset), len(Wset) - len(Wcomp), len(Wcomp)))
print("    composite members j<=20000:", Wcomp)
print("    even members:", [j for j in Wset if j % 2 == 0])
for j in (85, 91, 341, 946):
    Nv = Nj(j)
    check(not isprime(Nv) and pow(2, Nv - 1, Nv) == 1 and pow(4, Nv - 1, Nv) == 1, "composite-index psp j=%d" % j)
print("    direct verification: N_85, N_91, N_341, N_946 are base-2 (and base-4) pseudoprimes (85=5*17, 91=7*13, 341=11*31, 946=2*11*43).")
print("  PROVED (tower closure): j in W  =>  N_j in W.  Proof: ord_{N_j}(4)=j and 2j | N_j-1 give 4^{N_j-1}=1 mod N_j, i.e. N_j | 4^{N_j}-4;")
print("    N_j is odd (1+4+...+4^{j-1}) and 3 !| N_j (N_j = j mod 3, 3 !| j), and 4^m = 4 mod 6 for all m>=1; CRT gives 6N_j | 4^{N_j}-4.")
print("    Hence every prime p>=5 starts an infinite chain p -> N_p -> N_{N_p} -> ... of base-2 pseudoprimes (Cipolla's construction iterated).")
for pp in (5, 7, 11, 13):
    Nv = Nj(pp)
    check(pow(4, Nv, 6 * Nv) == 4, "N_p in W for p=%d" % pp)
    NN = Nj(Nv) if Nv <= 5461 else None
    if NN is not None:
        check(pow(2, NN - 1, NN) == 1 and not isprime(NN), "second level N_{N_p} psp for p=%d" % pp)
print("    checked: N_5=341, N_7=5461, N_11, N_13 are in W; second level N_341 (682 bits) and N_5461 (10922 bits) verified base-2 pseudoprimes.")

print("\n--- S2.4  TYPING LINES (machine-readable summary of the maps that exist; see the note's table) ---")
typing = [
 ("p^2qr squarefree divisors -> Fano/PG(2,2)", "MAP", "source divisor lattice; target F_2^3; map exponent vector mod 2; preserved S=7 points, 7 square-product lines; lost exponent heights + 112 of 128 orientations; sidecar XOR label; test: lattice selects none of the 16 octonion orientations (S1.1)"),
 ("Omega(a^2b)=4, Omega(c^2b)=8 -> dim H, dim O", "NO MAP (numerology)", "additive count vs multiplicative doubling; next term 10 vs 16 (S1.2)"),
 ("Bott period 8 -> Collatz rows/braid", "NO MAP FOUND", "no Collatz modulus has ord(2)=8; periods 2,3,6 (S1.3)"),
 ("6/pi^2 -> Collatz", "MAP via squarefreeness only", "T(n) squarefree has density 9/pi^2 row-uniformly; gcd-coprimality of affine pairs is trivial (S1.4); floor identities detect squarefreeness (inherited)"),
 ("{7,21} h-holes -> c(b) cycle counts", "NO MAP", "7 = c(11); 1,3 are the only odd holes of c(b) (S2.2)"),
 ("{7,21} -> 63=2^6-1 / Mersenne", "NO MAP", "15,31,63 attained; 7 falls in inter-band gaps (THM-1745) (S2.1)"),
 ("{7,21} -> {0,2,10} / rows {1,3,5}", "NO MAP FOUND", "numerology S=7, S*U=21 only (S2.2)"),
 ("Cipolla trunk (4^j-1)/3 -> inverse braid R(n)=4n+1", "MAP", "N_j = R^{j-1}(1); halving exponent 2j = ord_{N_j}(2); psp iff 6j | 4^j-4 (S2.3)"),
 ("'341 stalls primitive primes' -> Zsigmondy", "REFUTED", "11, 31 primitive at j=5; no exception for base 4 (S2.3)"),
 ("3N-5 -> Collatz", "MAP (shift, not conjugacy)", "T_{-5}(n) = T(n-2); rows precomposed by r -> r-2 (S4); cycles = -(3n+5 cycles) (inherited)"),
 ("odd cycles ~ odd functions ~ odd powers", "MAP (odd functions <-> odd powers), NO MAP to Redei odd cycles", "negation conjugacy and the odd-function floor law (inherited floor lanes); proof template only for tournaments"),
]
for (pair, verdict, detail) in typing:
    print("  [%s] %s :: %s" % (verdict, pair, detail))

# ---------------------------------------------------------------- S3 conjectures
hdr("S3. Six bold conjectures, probes, follow-ups")

print("\n--- C1 (G-map x primes): the greedy 3-adic stopping-time law on primes equals the law on composites ---")
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
    path = []; vv = m
    while not (vv <= MG and steps[vv] >= 0):
        path.append(vv); vv, _ = Gmap(vv)
    base = int(steps[vv])
    for i, pth in enumerate(reversed(path)):
        if pth <= MG:
            steps[pth] = base + i + 1
ms = np.arange(2, MG + 1)
ms = ms[ms % 3 != 0]
st = steps[ms]
check(int(st.min()) >= 0, "all G-steps resolved")
pm = isP[ms]
def tv(a, b):
    ha = np.bincount(a, minlength=200) / len(a); hb = np.bincount(b, minlength=200) / len(b)
    return 0.5 * float(np.abs(ha - hb).sum())
sp = st[pm]; sc = st[~pm]
print("  pooled m<=10^6, 3!|m: primes %d (mean steps %.4f, max %d) ; composites %d (mean %.4f, max %d) ; TV distance %.4f"
      % (len(sp), sp.mean(), sp.max(), len(sc), sc.mean(), sc.max(), tv(sp, sc)))
print("  CORRECTED NEAR MISS (recovered draft): its hostile control 'm=1 mod 9 vs m=8 mod 9 differ' is FALSE:")
mk = (ms >= 2**19) & (ms < 2**20)
print("    TV(1 mod 9, 8 mod 9) = %.4f pooled, %.4f in [2^19,2^20); the two classes have equal means (%.3f vs %.3f in the block)."
      % (tv(st[ms % 9 == 1], st[ms % 9 == 8]), tv(st[mk & (ms % 9 == 1)], st[mk & (ms % 9 == 8)]), st[mk & (ms % 9 == 1)].mean(), st[mk & (ms % 9 == 8)].mean()))
print("  TRUE PICTURE: the step count grows like log m, so any pooled comparison mixes sizes; primes are relatively denser at small m,")
print("    which produces the pooled mean gap %.3f.  Size-controlled comparison (dyadic blocks):" % (sc.mean() - sp.mean()))
blocks = []
for (lo, hi) in ((2**17, 2**18), (2**18, 2**19), (2**19, 2**20)):
    mkb = (ms >= lo) & (ms < hi)
    a = st[mkb & pm]; b = st[mkb & ~pm]
    blocks.append(tv(a, b))
    print("    [%d,%d): primes %6d mean %.4f ; composites %6d mean %.4f ; TV %.4f ; |mean gap| %.4f" % (lo, hi, len(a), a.mean(), len(b), b.mean(), tv(a, b), abs(a.mean() - b.mean())))
check(max(blocks) < 0.03, "C1: size-controlled TV < 0.03")
tv_size = tv(st[(ms >= 2**17) & (ms < 2**18)], st[(ms >= 2**19) & (ms < 2**20)])
tv_res = tv(st[mk & (ms % 9 == 4)], st[mk & (ms % 9 == 5)])
print("  hostile controls (both DIFFER): block [2^17,2^18) vs [2^19,2^20): TV %.4f (size); classes 4 vs 5 mod 9 in [2^19,2^20): TV %.4f"
      % (tv_size, tv_res))
print("    (class means in [2^19,2^20) by m mod 9:", [(rr, round(float(st[mk & (ms % 9 == rr)].mean()), 3)) for rr in (1, 2, 4, 5, 7, 8)], ")")
check(tv_size > 0.1 and tv_res > 0.1, "C1 hostile controls differ")
print("  mechanism: the first J greedy letters are a function of m mod 3^(J+1) (inherited extended_collatz_scc Theorem 4.1), so the")
print("    residue class of m matters (4 mod 9: first step x1/3; 5 mod 9: x8/3) and primes are equidistributed in the unit classes (Dirichlet):")
units81 = [a for a in range(1, 81) if a % 3 != 0]
hp = np.bincount(ms[pm] % 81, minlength=81)[units81] / pm.sum()
hc = np.bincount(ms[~pm] % 81, minlength=81)[units81] / (~pm).sum()
print("    primes vs composites mod 81: TV %.4f." % (0.5 * float(np.abs(hp - hc).sum())))
print("  STATUS: SURVIVES in the size-controlled form (FINITE-EXACT at 10^6: TV < 0.03 in each block); PROVED for every fixed-depth prefix")
print("    event (Theorem 4.1 + Dirichlet); the full-law equality is HEURISTIC.  The pooled statement is REFUTED as stated (size mixing).")
print("  follow-up F1: within a fixed block AND a fixed class mod 9, do primes and composites still agree? ->", end=" ")
diffs = []
for rr in (1, 2, 4, 5, 7, 8):
    a = st[mk & (ms % 9 == rr) & pm]; b = st[mk & (ms % 9 == rr) & ~pm]
    diffs.append((rr, round(float(a.mean()), 3), round(float(b.mean()), 3), round(tv(a, b), 4)))
print(diffs, " (all TV < 0.04 with ~5800 primes per cell: yes)")
check(all(d[3] < 0.04 for d in diffs), "F1")

print("\n--- C2 (signed cycles x cycle gate): b = 2^K - 3^L  =>  EVERY composition of K into L parts is realized as a cycle word ---")
def Bword(w):
    Kc = 0; B = 0; L = len(w)
    for i, ki in enumerate(w):
        B += 3**(L - 1 - i) * 2**Kc
        Kc += ki
    return B
def realize(n, k, L):
    w = []
    vv = n
    for _ in range(L):
        xx = 3 * vv + k; e = v2(xx); w.append(e); vv = xx >> e
    return tuple(w), vv
print("  INHERITED (PROVED in 04-computation/experiments/collatz_mod6_20260917_three_n_plus_k_catalan.py, necklace theorem; and the")
print("  q | b iff of arithmetic_braids2_20260917_signed_cycles.md sec.3): 3B(w)+2^K-3^L = 2^(k_1) B(rot w), so for b=2^K-3^L the orbit of B(w)")
print("  realizes w and closes; distinct necklaces give distinct cycles.  Re-verified here as an independent path (bold form: ALL compositions).")
reps = {}
for K in range(1, 70):
    for L in range(1, 45):
        d = 2**K - 3**L
        if 1 <= abs(d) <= 99 and gcd(d, 6) == 1:
            reps.setdefault(d, []).append((K, L))
print("  representations b = 2^K - 3^L with |b|<=99, K<70, L<45 (FINITE-EXACT):", dict(sorted(reps.items())))
def necklace_count(K, L):
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
            check(wr == w and back == B, "C2 realization b=%d w=%s" % (d, w))
            total += 1
            hit = [c for c in CYC[d] if B in c]
            check(len(hit) == 1, "C2 census hit")
            found.add(CYC[d].index(hit[0]))
        print("   b=%2d = 2^%d-3^%d : %d compositions, all realized; distinct census cycles they populate: %d = rotation classes %d"
              % (d, K, L, total, len(found), necklace_count(K, L)))
        check(len(found) == necklace_count(K, L), "necklace count")
print("  bold sub-conjecture 'c(b) is maximal among b<=99 exactly at the doubly-representable b (5, 13)':", end=" ")
cmax = max(t[1] for t in row); argmax = [t[0] for t in row if t[1] == cmax]
print("max c(b) = %d at b=%s ->" % (cmax, argmax), "SURVIVES" if set(argmax) <= {5, 13} else "REFUTED (witness b=%s, c=%d: 65 = 5*13 inherits both scaled families)" % (argmax, cmax))
check(argmax == [65] and cmax == 19, "C2 sub-conjecture witness")
print("  follow-up F2: fixed points of 3n+b are exactly n = b/(2^K-3), K>=1, so #1-cycles = #{K>=1 : (2^K-3) | b}: check ->", end=" ")
for k in KS:
    pred = sum(1 for K in range(1, 8) if k % abs(2**K - 3) == 0)
    check(pred == sum(1 for c in CYC[k] if len(c) == 1), "fixed-point count b=%d" % k)
print("PASS for all b<=99 (PROVED + FINITE-EXACT).  b=65=5*13 has %d fixed points: %s" % (sum(1 for c in CYC[65] if len(c) == 1), [c for c in CYC[65] if len(c) == 1]))

print("\n--- C3 (floor sums x rows/braid): the squarefree detector Z_3 on the trunk: 'N_j squarefree iff gcd(j, N_j) = gcd(j, 3)' ---")
def Z_e(n, e):
    return sum(1 for kk in range(1, n) if pow(kk, e, n) == 0)
def Z_formula(n, e):
    pr = 1
    for pp, a in factorint(n).items():
        pr *= pp**(a - (-(-a // e)))
    return pr - 1
print("  INHERITED (PROVED, arithmetic_braids2_20260917_floor_reciprocity.md; collatz_mod6_20260917_floor_sums_odd_functions.md S1-S2):")
print("    sum_{k<n} floor(k^3/n) = (n-2)(n-1)(n+1)/4 + Z_3(n)/2 with Z_3(n) = n/n_3 - 1, so the cube floor identity holds iff n is squarefree.")
for j in range(2, 11):
    n = Nj(j)
    S1 = sum(kk**3 // n for kk in range(1, n))
    Z3 = Z_formula(n, 3)
    check(S1 == Fraction((n - 2) * (n - 1) * (n + 1), 4) + Fraction(Z3, 2), "cube floor identity at N_%d" % j)
    if j <= 6:
        check(Z3 == Z_e(n, 3), "Z_3 brute at N_%d" % j)
print("    re-verified on the trunk terms N_2..N_10 (Z_3 by brute force for j<=6).  Detector reading: identity holds at N_j iff N_j squarefree.")
print("  PROVED exact law (LTE): for an odd prime p, v_p(4^j-1) = v_p(4^d-1) + v_p(j/d) when d = ord_p(4) | j (else 0); and v_p(4^d-1) = v_p(2^{p-1}-1).")
print("    Hence p^2 | N_j (p>=5) iff ord_p(4) | j and (p | j or p is a Wieferich prime); 9 | N_j iff 9 | j (v_3(N_j) = v_3(j)).")
print("    So the conjecture 'N_j squarefree iff gcd(j,N_j) = gcd(j,3)' is exactly the statement that no Wieferich prime interferes.")
sqf_status = []
for j in range(1, 41):
    n = Nj(j)
    f = factorint(n)
    is_sf = all(e == 1 for e in f.values())
    conj = (gcd(j, n) == gcd(j, 3))
    check(is_sf == conj, "C3 at j=%d" % j)
    if not is_sf:
        sqf_status.append((j, {pp: e for pp, e in f.items() if e > 1}))
print("  FINITE-EXACT j<=40: conjecture holds; non-squarefree trunk terms and their square factors:", sqf_status)
w1093 = 1093; w3511 = 3511
check(pow(2, w1093 - 1, w1093**2) == 1 and pow(2, w3511 - 1, w3511**2) == 1, "Wieferich primes")
d1093 = mult_order(4, w1093); d3511 = mult_order(4, w3511)
print("  REFUTED in general: Wieferich primes 1093, 3511 (CITED: the only ones below 6.7e15, Dorais-Klyve 2011 / PrimeGrid) with ord_1093(4) = %d, ord_3511(4) = %d:" % (d1093, d3511))
N182 = Nj(182)
check(N182 % (1093**2) == 0 and 182 % 1093 != 0, "1093^2 | N_182")
check(gcd(182, N182) == 1 and gcd(182, 3) == 1, "gcd(182, N_182) = 1 = gcd(182, 3)")
print("    witness j=182 = 2*7*13: 1093^2 | N_182 while 1093 !| 182; gcd(182, N_182) = %d = gcd(182,3) (ord_7(4)=3 and ord_13(4)=6 do not divide 182,"
      % gcd(182, N182))
print("      so the p | j mechanism is silent), i.e. the conjecture predicts 'squarefree' and N_182 is not: REFUTED, Wieferich witness.")
print("    minimal witness: j=182 is the least j with a Wieferich contribution, conditional on no unknown Wieferich prime p with ord_p(4) < 182;")
print("      unconditionally the conjecture is FINITE-EXACT-true for j<=40 and false at j=182.")
print("  STATUS: PROVED law; conjecture FINITE-EXACT for j<=40 and REFUTED at j=182 (Wieferich).  The floor-sum squarefree detector on")
print("    the Collatz trunk sees the Wieferich primes: cross-thread map (rows/braid trunk) x (floor sums) with a concrete sidecar (LTE).")
print("  follow-up F3: is N_j EVER squarefree for j = 0 mod 9? No (9 | N_j).  Fraction of j<=40 with N_j squarefree: %d/40." % (40 - len(sqf_status)))

print("\n--- C4 (sandwich x squares): N_PS - N_SP = #{p>=5: p^2<=6K+1, p^2-2 prime} + bounded residual ? ---")
def sandwich(K):
    kk = np.arange(1, K + 1, dtype=np.int64)
    a = Om[6 * kk - 1]; b = Om[6 * kk + 1]
    ca = np.minimum(a, 4); cb = np.minimum(b, 4)
    M = np.zeros((5, 5), dtype=np.int64)
    np.add.at(M, (ca, cb), 1)
    return M
ctrl = {10**5: 41, 10**6: 412}
resid = {}
for K in (10**5, 10**6):
    M = sandwich(K)
    dPS = int(M[1, 2] - M[2, 1])
    check(dPS == ctrl[K], "sandwich_bias control at K=%d" % K)
    A = sum(1 for pp in primerange(5, math.isqrt(6 * K + 1) + 1) if isprime(pp * pp - 2))
    # independent count of the same pairs straight from the sieve: 6k+1 = p^2 (p>=5 prime) and 6k-1 prime
    kk = np.arange(1, K + 1, dtype=np.int64)
    rt = np.sqrt((6 * kk + 1).astype(np.float64)).round().astype(np.int64)
    sqmask = (rt * rt == 6 * kk + 1) & isP[np.minimum(rt, NS)] & (rt >= 5)
    A2 = int((sqmask & (Om[6 * kk - 1] == 1)).sum())
    check(A == A2, "square-driven pair count at K=%d" % K)
    Rr = M.sum(axis=1); Cc = M.sum(axis=0)
    pred = (Rr[1] * Cc[2] - Rr[2] * Cc[1]) / K
    resid[K] = dPS - A
    print("  K=%d: N_PS-N_SP = %d (matches collatz_mod6_20260917_sandwich_bias.md sec.1) ; square-driven pairs (p^2-2 prime, p^2) = %d ; residual after squares = %d ; independence pred = %.1f"
          % (K, dPS, A, dPS - A, pred))
print("  CORRECTED NEAR MISS (recovered draft): its status line quoted square counts 64 and 137; the true counts of pairs (6k-1, 6k+1) = (p^2-2, p^2)")
print("    with both entries as stated are %d and %d (two independent counts agree)." % (46, 99))
check(resid[10**5] == -5 and resid[10**6] == 313, "C4 residuals")
print("  STATUS: REFUTED as 'bounded residual': residual after removing squares is %d at 10^5 and %d at 10^6 (grows with the drift, not bounded)." % (resid[10**5], resid[10**6]))
print("  The square term is PROVED to be one-sided (p^2 = 1 mod 6 always sits on the right), but it is small against the mixed-class semiprime bias.")
print("  follow-up F4: does the p^2 count alone explain the S-column excess of the marginals?  S_1 - S_5 identity (inherited) already has pi'(sqrt x)/2:")
for K in (10**5, 10**6):
    x6 = 6 * K + 1
    sq = sum(1 for pp in primerange(5, math.isqrt(x6) + 1))
    M = sandwich(K)
    print("    K=%d: pi'(sqrt x) = %d, (S_1 - S_5) = %d (col-row sum of class 2) ; half the squares = %.1f" % (K, sq, int(M.sum(axis=0)[2] - M.sum(axis=1)[2]), sq / 2))

print("\n--- C5 (3n+b x G-map conjugacy): T_b(b n) = b T_1(n), G_b(b m) = b G_1(m); c(b) >= 4 and superadditivity ---")
for k in KS:
    base = [tuple(k * t for t in c) for c in CYC[1]]
    for bc in base:
        check(any(sorted(bc) == sorted(c) for c in CYC[k]), "scaled Collatz cycle in C_b, b=%d" % k)
print("  PROVED: (3(bn)+b)/2^v = b(3n+1)/2^v ; so b*C_1 subset C_b, c(b) >= 4 for every b coprime to 6 (checked b<=99).")
print("  (This is the dilation half of the content theorem, arithmetic_braids2_20260917_signed_cycles.md sec.2: cycles of T_b = disjoint union over g | b of g*(primitive cycles of T_{b/g}).)")
eq4 = [k for k in KS if len(CYC[k]) == 4]
print("  b<=99 with c(b) = 4 exactly (no cycle beyond the four scaled Collatz cycles):", eq4)
print("  superadditivity for coprime b1,b2: c(b1 b2) >= c(b1) + c(b2) - 4 (scaled cycles b2*C_{b1}, b1*C_{b2} meet only in b1b2*C_1):")
viol = []
for k1 in KS:
    for k2 in KS:
        if k1 < k2 and gcd(k1, k2) == 1 and k1 * k2 in CYC:
            lhs = len(CYC[k1 * k2]); rhs = len(CYC[k1]) + len(CYC[k2]) - 4
            if lhs < rhs:
                viol.append((k1, k2, lhs, rhs))
print("    violations among b1*b2<=99:", viol, "-> PROVED inequality holds (FINITE-EXACT confirmation)")
check(viol == [], "superadditivity")
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
        check(a == k * b and ja == jb, "G_b conjugacy b=%d m=%d" % (k, m))
print("  PROVED (+checked m<3000, b=5,7,11,13): G_b(b m) = b G_1(m) with the same letter; the residue chain mod 9 of G_b is the")
print("    chain of G_1 conjugated by multiplication by b^{-1} mod 9, so the stationary law, E[letter]=1 and the drift log(2/3)")
print("    (inherited extended_collatz_scc S4.2) are b-independent.  What changes with b is the CYCLE set (C2), not the local law.")
print("  follow-up F5: are the 'new' cycles (beyond b*C_1) of 3n+b all divisor-gate cycles n = bB/(2^K-3^L) with (2^K-3^L) !| B? ->", end=" ")
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
print("%d new cycles, %d have (2^K-3^L) !| B (need the factor b), %d have (2^K-3^L) | B" % (cntnew, cntgate, cntnew - cntgate))
check(cntnew == cntgate, "all new cycles need b")

print("\n--- C6 (Cipolla trunk x primes): 'N_j is a base-2 pseudoprime iff j is a prime >= 5' ---")
print("  FINITE-EXACT j<=40: holds (S2.3 table).  REFUTED at j=85 (= N_4 = 5*17): N_85 is a base-2 pseudoprime with composite index;")
print("    mechanism: for odd j with 3 !| j, 6j | 4^j-4 iff 4^{j-1} = 1 mod j, i.e. j is a base-4 Fermat (pseudo)prime; 85, 91, 341, ... are base-4 pseudoprimes.")
check(pow(4, 84, 85) == 1 and pow(4, 90, 91) == 1 and pow(4, 14, 15) == 1 and 15 not in set(Wset), "base-4 psp mechanism; 15 excluded by 3 | 15")
print("    hostile: 15 is a base-4 pseudoprime but 3 | 15, so N_15 is NOT a pseudoprime (9 !| 4^15-4 since 3 !| 14): the factor 3 of the trunk matters.")
trunkvals = set(Nj(j) for j in range(1, 12))
check(85 in trunkvals and 91 not in trunkvals and 341 in trunkvals and 5461 in trunkvals, "trunk membership of composite indices")
print("  follow-up F6 (recursion, PROVED in S2.3): W is closed under j -> N_j, so the index set contains p, N_p, N_{N_p}, ... for every prime p>=5;")
print("    the first composite index NOT of the form N_j is 91 = 7*13 (85 = N_4 is of that form; 341 = N_5 and 5461 = N_7 are the tower's second level).")

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
    jj = sympy.symbols('j')
    img = (3 * (6 * jj + rr) - 5) / 2
    print("    n = 6j+%d : (3n-5)/2 = %s  = %d mod 9  (Collatz row image of 6j+%d would be %d mod 9)" % (rr, sympy.expand(img), (3 * rr - 5) // 2 % 9, rr, (3 * rr + 1) // 2 % 9))
print("  PROVED: the 3n-5 row map is the Collatz row map precomposed with r -> r-2 mod 6 (1->5->3->1): images 8,2,5 instead of 2,5,8;")
print("  all images are 2 mod 3 exactly as for Collatz (multiples of 3 are leaves for 3n-5 as well).")
cyc_m5 = sorted(tuple(-t for t in c) for c in CYC[5])
print("  cycles of 3n-5 on odd integers, both signs, min-element |n|<=2*10^5 (= negatives of the 9 cycles of 3n+5):")
for c in cyc_m5:
    print("    length %2d  min |n| %6d : %s" % (len(c), abs(c[0]), c if len(c) <= 8 else c[:8] + ('...',)))
pos_m5 = [c for c in cyc_m5 if c[0] > 0]
print("  POSITIVE 3n-5 cycles = %d (the user's '3N-5 on positives' = 3N+5 on negatives): %s" % (len(pos_m5), [c for c in pos_m5]))
check(len(cyc_m5) == 9 and len(pos_m5) == 3, "nine 3n-5 cycles, three positive")
print("  the two 3-cycles live on the POSITIVE side of 3n+5 = NEGATIVE side of 3n-5: (19,31,49),(23,37,29),")
print("  and are the C2 construction for 5 = 2^5-3^3 (words (1,1,3) and (1,2,2)); 3n-5 has them as (-19,-31,-49),(-23,-37,-29).")
print("  shifted picture: a 3n-5 cycle (n_i) is a cycle (n_i - 2) of m -> T(m)-2, e.g. fixed point -5 of 3n-5 <-> m=-7: T(-7)-2 = -5-2 = -7.")
check(Tk(-7, 1) - 2 == -7 and Tk(-5, -5) == -5, "shift example")
print("  PROVED: 3n-5 on negatives = -(3n+5 on positives); 3n-5 on positives = -(3n+5 on negatives) (sign conjugation, inherited summand note sec.7;")
print("  the nine cycles, their contents d in {1,5} and minimal parameters q in {1,5}: arithmetic_braids2_20260917_signed_cycles.md sec.5).")

# ---------------------------------------------------------------- S5 summary
hdr("S5. Status summary")
print("PROVED: S1.1 F_2^3/Fano structure of squarefree divisors and {0,2,10}; d(N)<=2^Omega with equality iff squarefree; 9/pi^2 density of")
print("        odd n with T(n) squarefree and the (8/9)prod(1-2/p^2) density; gcd(n,T(n))=1, gcd(n+2,T(n)) | 5; census completeness bound;")
print("        S2.3 Cipolla (elementary proof), criterion 6j | 4^j-4, tower closure of W, Zsigmondy exception-freeness for base 4;")
print("        C2 re-verification of the necklace theorem (inherited) and the fixed-point count #{K:(2^K-3)|b}; C3 LTE law for squares in N_j;")
print("        C5 scaling b*C_1 in C_b, superadditivity, G_b conjugacy; S4 T_{-5}(n)=T(n-2) and the row map; C_{-b} = -C_b.")
print("FINITE-EXACT: 16/128 octonion orientations; tournament census n<=6; c(b) spectrum |b|<=99 in U1 and U2 (complete for L<=L*(b));")
print("        6/pi^2 candidates at 10^6; C1 block TV distances at 10^6; C4 sandwich numbers at 10^5,10^6; W to 20000; C3 j<=40.")
print("CITED: Cipolla 1904; Albuquerque-Majid twisted group algebra; Dirichlet; THM-1745 / death-star S70 h-spectrum; Redei; Wieferich list.")
print("HEURISTIC: full G-stopping-time law equality primes vs composites (size-controlled); C4 residual.")
print("REFUTED: '341 stalls primitive primes' (11,31 at j=5); 'N_j psp iff j prime' (j=85); 'N_j squarefree iff gcd(j,N_j)=gcd(j,3)' (j=182, 1093^2);")
print("        draft hostile 'classes 1 vs 8 mod 9 differ' (TV 0.007); pooled C1 (size mixing); 'N_PS-N_SP minus squares bounded' (-5 -> 313);")
print("        'c(b) maximal exactly at doubly representable b' (b=65, c=19).")
print("NO MAP FOUND: Bott 8; {7,21} vs c(b) or cycle lengths or rows or {0,2,10}; Omega(a^2b)=4 / Omega(c^2b)=8 vs H, O (numerology);")
print("        x^2+c period-3 bifurcation vs 3n+b 3-cycles (exponential Diophantine gate 2^K-27 | bB).")
print("total runtime", tm())
