#!/usr/bin/env python3
"""
collatz_mod6_20260917_wild_typing_audit.py -- independent adversarial recompute for lane wild_typing
(session collatz-mod6-20260917, mac-mini; audit written 2026-09-21).

Every key number of collatz_mod6_20260917_wild_typing.{md,out} is recomputed here with code written
from scratch (different cycle-census algorithm, different tournament DP, different Fano test, direct
pow() pseudoprime tests), plus the boundary/quantifier probes the audit brief asks for:
  A  trunk N_j=(4^j-1)/3: R-orbit, ord_{N_j}(2)=2j, criterion 6j | 4^j-4, W to 20000, base-4-only indices
     (the 2-adic gap v_2(N_j-1)=2), tower closure, Cipolla for primes p<=200, Zsigmondy j<=60, j=15 hostile
  B  LTE law for v_p(N_j) vs factorint (p<=200, j<=60); C3 conjecture j<=60; Wieferich witness 182; 3511
  C  densities at 10^6 odd n: T(n) squarefree, joint, gcd(n+2,T(n)), rows
  D  tournament census n<=6 (per-tournament bit DP; Aut via numpy over permutations)
  E  3n+b cycle census |b|<=99, |n|<=2*10^5, escape 10^40 (own algorithm); completeness L*(b); spectrum;
     F2 fixed points; C5 scaling, superadditivity (true form with c(1)); F5; necklace theorem; doubly
     representable witness b=-1
  F  greedy 3-adic G-map stopping-time law: primes vs composites, blocks, residue classes (m<=10^6)
  G  sandwich N_PS-N_SP at 10^5, 10^6; square-driven pairs; S_1-S_5; pi'(sqrt x)
  H  3n-5 = T(n-2); row images
  I  Fano: 16/128 orientations by norm multiplicativity (independent of alternativity)
  J  divisor-box facts; ord_m(2)=8 moduli
All checks are explicit raises (survive python -O).  RAM < 1 GB, runtime ~1-2 min.
"""
import sys, time, math, itertools
from fractions import Fraction
from math import gcd
import numpy as np
from sympy import factorint, isprime, primerange, n_order

T0 = time.time()
def tm():
    return "[t=%.1fs]" % (time.time() - T0)
def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)
def hdr(s):
    print(); print("=" * 78); print(s); print("=" * 78)
def v2(x):
    return (x & -x).bit_length() - 1
def Nj(j):
    return (4**j - 1) // 3

# ------------------------------------------------------------------ A trunk / Cipolla
hdr("A. Trunk N_j = (4^j-1)/3: orders, criterion, W, base-4-only gap, tower, Cipolla, Zsigmondy")
v = 1
for j in range(1, 61):
    check(v == Nj(j), "R-orbit at j=%d" % j)
    check(3 * v + 1 == 4**j, "3N_j+1=4^j at j=%d" % j)
    v = 4 * v + 1
print("N_j = R^(j-1)(1) and 3N_j+1 = 4^j: checked j<=60")
for j in range(2, 41):
    N = Nj(j)
    check(n_order(2, N) == 2 * j, "ord_{N_j}(2) = 2j at j=%d" % j)
    check(n_order(4, N) == j, "ord_{N_j}(4) = j at j=%d" % j)
print("ord_{N_j}(2) = 2j and ord_{N_j}(4) = j (sympy n_order): checked 2<=j<=40")
# composite for j>=3, prime at j=2
check(isprime(Nj(2)) and all(not isprime(Nj(j)) for j in range(3, 41)), "N_j composite iff j>=3 (j<=40)")
# criterion: base-2 psp iff 6j | 4^j-4 (j>=3), direct pow test j<=400
psp_direct = []
for j in range(3, 401):
    N = Nj(j)
    d2 = pow(2, N - 1, N) == 1
    d4 = pow(4, N - 1, N) == 1
    crit = (4**j - 4) % (6 * j) == 0
    check(d2 == crit, "criterion 6j | 4^j-4 vs direct base-2 test at j=%d" % j)
    check(d4 == (((N - 1) % j) == 0), "base-4 psp iff j | N_j-1 at j=%d" % j)
    if d2:
        psp_direct.append(j)
print("direct base-2 test == criterion 6j | 4^j-4 for all 3<=j<=400; psp indices j<=40:", [j for j in psp_direct if j <= 40])
check([j for j in psp_direct if j <= 40] == [j for j in range(5, 41) if isprime(j)], "j<=40 psp exactly primes >= 5")
check(psp_direct[:12] == [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43], "first 12")
# v_2(N_j - 1) = 2 for all j>=2 : base-4 psp needs j | N_j-1, base-2 needs 2j | N_j-1
for j in range(2, 200):
    check(v2(Nj(j) - 1) == 2, "v_2(N_j-1)=2 at j=%d" % j)
print("PROVED+checked j<200: v_2(N_j-1) = 2 exactly (N_j-1 = 4(4^(j-1)-1)/3 with 4^(j-1)-1 odd).")
W = [j for j in range(3, 20001) if (4**j - 4) % (6 * j) == 0]
# faster equivalent via pow for the base-4-only search
W4only = []
for j in range(3, 20001):
    if pow(4, j, 3 * j) == 4 % (3 * j):      # j | (4^j-4)/3  <=> 3j | 4^j-4
        if (4**j - 4) % (6 * j) != 0:
            W4only.append(j)
print("W = {3<=j<=20000 : 6j | 4^j-4}: |W| = %d, primes>=5 in W: %d, composites: %d" % (len(W), sum(1 for j in W if isprime(j)), sum(1 for j in W if not isprime(j))))
Wc = [j for j in W if not isprime(j)]
print("  composite members (first 13):", Wc[:13], " even members:", [j for j in W if j % 2 == 0])
check(len(W) == 2321 and len(Wc) == 61, "|W| = 2321, 61 composite")
check(Wc[:13] == [85, 91, 341, 451, 703, 946, 1105, 1247, 1271, 1387, 1729, 1891, 2047], "composite members")
check([j for j in W if j % 2 == 0] == [946, 2926, 8614], "even members")
check(all(j % 3 != 0 for j in W), "3 !| j in W")
check(all(j in set(W) for j in primerange(5, 20001)), "all primes >= 5 in W")
print("  base-4-only indices (j | N_j-1 but 2j !| N_j-1, i.e. N_j base-4 psp but NOT base-2 psp), j<=20000:", W4only)
print("  => these are exactly the j = 4 mod 8 with j | N_j-1; the note's 'base-2 iff base-4' corollary needs 4 !| j.")
for j in W4only[:3]:
    N = Nj(j)
    check(pow(4, N - 1, N) == 1 and pow(2, N - 1, N) != 1 and not isprime(N), "base-4-only witness j=%d" % j)
    print("     witness j=%d: N_j is a base-4 pseudoprime (4^(N-1)=1 mod N) and 2^(N-1) != 1 mod N; %d bits" % (j, N.bit_length()))
# tower closure
for p in (5, 7, 11, 13):
    N = Nj(p)
    check((4**N - 4) % (6 * N) == 0 if N < 6000 else pow(4, N, 6 * N) == 4, "N_p in W, p=%d" % p)
for p in (5, 7):
    NN = Nj(Nj(p))
    check(pow(2, NN - 1, NN) == 1 and NN % 3 != 0, "second level N_{N_p} base-2 psp p=%d" % p)
print("tower: N_5, N_7, N_11, N_13 in W; N_341 (%d bits) and N_5461 (%d bits) are base-2 pseudoprimes" % (Nj(341).bit_length(), Nj(5461).bit_length()))
# base-4-only tower question: does j in W4only give N_j in W?  (only the base-4 part is inherited)
for j in W4only[:2]:
    N = Nj(j)
    check(pow(4, N, 6 * N) == 4, "N_j in W for base-4-only j=%d" % j)
print("  tower closure also holds from base-4-only indices: j | N_j-1 already gives N_j | 4^{N_j}-4, so N_j in W for j=%s" % W4only[:2])
# Cipolla direct for primes p<=200
for p in primerange(5, 201):
    N = Nj(p)
    check(pow(2, N - 1, N) == 1 and N % (2**p - 1) == 0 and N // (2**p - 1) == (2**p + 1) // 3, "Cipolla p=%d" % p)
check(pow(2, 20, 21) == 4, "N_3=21 not psp")
print("Cipolla: (4^p-1)/3 base-2 psp for all primes 5<=p<=200 (direct); N_3=21: 2^20 = 4 mod 21")
# Zsigmondy primitive primes for 4^j-1, j<=60
prev = set()
for j in range(1, 61):
    f = factorint(4**j - 1)
    prim = [p for p in f if p not in prev]
    check(len(prim) >= 1, "primitive prime at j=%d" % j)
    check(j == 1 or 3 not in prim, "3 primitive only at j=1")
    if j == 5:
        check(sorted(prim) == [11, 31], "j=5 primitive primes 11, 31")
    prev |= set(f)
print("Zsigmondy 4^j-1: a primitive prime at every j<=60; j=5 brings {11,31}: '341 stalls' REFUTED")
check(pow(4, 14, 15) == 1 and 15 not in set(W) and pow(2, Nj(15) - 1, Nj(15)) != 1, "j=15 hostile")
print("hostile j=15: base-4 psp index, 3 | 15, N_15 not a pseudoprime")

# ------------------------------------------------------------------ B LTE / C3
hdr("B. LTE law v_p(N_j) and the C3 conjecture; Wieferich witnesses")
def vp(n, p):
    e = 0
    while n % p == 0:
        n //= p; e += 1
    return e
for p in primerange(5, 200):
    d = n_order(4, p)
    base = vp(4**d - 1, p)
    check(base == vp(2**(p - 1) - 1, p), "v_p(4^d-1) = v_p(2^(p-1)-1) p=%d" % p)
    for j in range(1, 61):
        pred = (base + vp(j // d, p)) if j % d == 0 else 0
        check(vp(Nj(j), p) == pred, "LTE at p=%d j=%d" % (p, j))
for j in range(1, 61):
    check(vp(Nj(j), 3) == vp(j, 3), "v_3(N_j)=v_3(j) at j=%d" % j)
print("LTE law verified against direct valuations for 5<=p<200, j<=60; v_3(N_j) = v_3(j) for j<=60")
nonsf = []
for j in range(1, 61):
    N = Nj(j); f = factorint(N)
    sf = all(e == 1 for e in f.values())
    conj = (gcd(j, N) == gcd(j, 3))
    check(sf == conj, "C3 conjecture at j=%d" % j)
    if not sf:
        nonsf.append((j, sorted((p, e) for p, e in f.items() if e > 1)))
print("C3 conjecture 'N_j squarefree iff gcd(j,N_j)=gcd(j,3)' holds for j<=60; non-squarefree j<=40:", [t for t in nonsf if t[0] <= 40])
check([t[0] for t in nonsf if t[0] <= 40] == [9, 10, 18, 20, 21, 27, 30, 36, 40], "non-squarefree indices j<=40")
check(sum(1 for j in range(1, 41) if all(e == 1 for e in factorint(Nj(j)).values())) == 31, "31/40 squarefree")
for p in (1093, 3511):
    check(pow(2, p - 1, p * p) == 1, "Wieferich %d" % p)
d1093 = n_order(4, 1093); d3511 = n_order(4, 3511)
print("ord_1093(4) = %d, ord_3511(4) = %d" % (d1093, d3511))
check(d1093 == 182 and d3511 == 1755, "orders")
N182 = Nj(182)
check(N182 % 1093**2 == 0 and gcd(182, N182) == 1 and gcd(182, 3) == 1, "witness j=182")
print("j=182: 1093^2 | N_182, gcd(182,N_182) = 1 = gcd(182,3): conjecture REFUTED at 182")
# least j failing, under the LTE law, using only the two known Wieferich primes:
fails = [j for j in range(1, 2000) if (j % 182 == 0 and j % 1093 != 0 and j % 9 != 0
         and all(not (j % n_order(4, p) == 0 and j % p == 0) for p in primerange(5, 200)))]
print("  under the LTE law the failing j (only Wieferich primes 1093, 3511 known) are multiples of 182 with 1093 !| j, 9 !| j and no p|j with ord_p(4)|j; first ones:", fails[:6])
check(fails[0] == 182, "least witness 182 (conditional)")
print("  3511 never produces a witness: 9 | ord_3511(4) = 1755 = 3^3*5*13, so any multiple has 9 | j and the conjecture's two sides are both false.")
check(1755 % 9 == 0, "9 | 1755")
# cube-floor identity on trunk with brute R3
for j in range(2, 8):
    n = Nj(j)
    S = sum(k**3 // n for k in range(1, n))
    R3 = sum(1 for k in range(n) if pow(k, 3, n) == 0)
    check(Fraction(S) == Fraction((n - 2) * (n - 1) * (n + 1), 4) + Fraction(R3 - 1, 2), "cube floor at N_%d" % j)
    check((R3 == 1) == all(e == 1 for e in factorint(n).values()), "R3=1 iff squarefree at N_%d" % j)
print("cube floor identity sum floor(k^3/n) = (n-2)(n-1)(n+1)/4 + (R3(n)-1)/2 with brute R3: N_2..N_7")

# ------------------------------------------------------------------ C densities
hdr("C. Densities at 10^6 odd n")
NS = 6 * 10**6 + 1
sqf = np.ones(NS + 1, dtype=np.bool_); sqf[0] = False
for p in primerange(2, math.isqrt(NS) + 1):
    sqf[p * p::p * p] = False
nodd = np.arange(1, 2 * 10**6, 2, dtype=np.int64)
x = 3 * nodd + 1
u = x // (x & -x)
cA = int(sqf[u].sum()); cC = int((sqf[u] & sqf[nodd]).sum())
g2 = np.gcd(nodd + 2, u)
check(int(np.gcd(nodd, u).max()) == 1, "gcd(n,T(n))=1")
check(set(np.unique(g2).tolist()) <= {1, 5}, "gcd(n+2,T(n)) in {1,5}")
dB = int((g2 == 1).sum())
print("T(n) squarefree: %d/10^6 = %.6f (9/pi^2 = %.6f); joint: %d = %.6f; gcd(n+2,T(n))=1: %d = %.6f" % (cA, cA / 1e6, 9 / math.pi**2, cC, cC / 1e6, dB, dB / 1e6))
check(cA == 911902 and cC == 737515 and dB == 800000, "density counts")
P5 = 1.0; P5b = 1.0
for p in primerange(5, 10**7):
    P5 *= 1 - 1 / p**2; P5b *= 1 - 2 / p**2
print("  products: prod_{p>=5}(1-1/p^2) = %.6f ; (8/9) prod_{p>=5}(1-2/p^2) = %.6f" % (P5, 8 / 9 * P5b))
check(abs(P5 - 9 / math.pi**2) < 1e-6 and abs(8 / 9 * P5b - 0.737449) < 2e-6, "predicted values")
for r in (1, 3, 5):
    m = nodd % 6 == r
    print("  row %d: %.6f" % (r, int((sqf[u] & m).sum()) / int(m.sum())))
rows = [int((sqf[u] & (nodd % 6 == r)).sum()) / int((nodd % 6 == r).sum()) for r in (1, 3, 5)]
check(all(abs(t - 9 / math.pi**2) < 5e-4 for t in rows), "row uniform")
del x, u, g2

# ------------------------------------------------------------------ D tournaments
hdr("D. Tournament census n=3..6 (independent bit DP + numpy Aut)")
def census(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    hs = np.zeros(2**m, dtype=np.int64)
    # adjacency bitmasks per tournament: out[v] = bitmask of w with v->w
    codes = np.arange(2**m, dtype=np.int64)
    out = [np.zeros(2**m, dtype=np.int64) for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        bit = (codes >> b) & 1
        out[i] |= bit << j
        out[j] |= (1 - bit) << i
    full = (1 << n) - 1
    # DP over subsets vectorised: f[S][v] = number of paths covering S ending at v
    f = {}
    for v in range(n):
        f[(1 << v, v)] = np.ones(2**m, dtype=np.int64)
    for S in range(1, 1 << n):
        if bin(S).count("1") < 2:
            continue
        for v in range(n):
            if not (S >> v) & 1:
                continue
            Sp = S & ~(1 << v)
            acc = np.zeros(2**m, dtype=np.int64)
            for w in range(n):
                if (Sp >> w) & 1:
                    acc += f[(Sp, w)] * ((out[w] >> v) & 1)
            f[(S, v)] = acc
    for v in range(n):
        hs += f[(full, v)]
    aut = np.zeros(2**m, dtype=np.int64)
    for perm in itertools.permutations(range(n)):
        ok = np.ones(2**m, dtype=np.bool_)
        for (i, j) in pairs:
            a, b = perm[i], perm[j]
            ok &= (((out[a] >> b) & 1) == ((out[i] >> j) & 1))
        aut += ok
    return hs, aut
union = set()
for n in range(3, 7):
    hs, aut = census(n)
    vals = sorted(set(hs.tolist())); av = sorted(set(aut.tolist()))
    check(all(h % 2 == 1 for h in vals) and int((hs % aut).max()) == 0, "Redei / Aut | h at n=%d" % n)
    check(int(hs.sum()) == math.factorial(n) * 2**((n - 1) * (n - 2) // 2), "sum h = n! 2^T(n-2) at n=%d" % n)
    print("  n=%d: h values %s ; |Aut| values %s" % (n, vals, av))
    union |= set(vals)
    if n == 6:
        check(vals == [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45] and av == [1, 3, 5, 9], "n=6 spectrum")
missing = [v for v in range(1, 46, 2) if v not in union]
print("  odd values <= 45 not attained for n<=6:", missing)
check(missing == [7, 21, 35, 39], "missing 7,21,35,39")
print("  (63 is attained at n=8 and every odd value in [1,609] except 7, 21 by n=8: THM-1370, FINITE-EXACT there; 7, 21 omitted for ALL n: THM-1370 / THM-200.)")

# ------------------------------------------------------------------ E cycle census
hdr("E. 3n+b cycle census, gcd(b,6)=1, 1<=b<=99, |n|<=2*10^5 both signs, escape 10^40 (own algorithm)")
NLIM = 2 * 10**5; CAP = 10**40
def Tb(n, b):
    y = 3 * n + b
    return y >> v2(y)
def census_b(b):
    """Returns list of cycles (as tuples starting at min |element|) with some element |n|<=NLIM."""
    seen_cycle_elems = {}     # element -> cycle id, for all cycle elements (any size)
    done = {}                 # start -> True for |start|<=NLIM already classified
    cycles = []
    escaped = 0
    for a in range(1, NLIM + 1, 2):
        for s in (a, -a):
            if s in done:
                continue
            path = []; pos = {}; v = s; found = None
            while True:
                if v in pos:
                    cyc = path[pos[v]:]
                    m = min(cyc, key=abs); i = cyc.index(m)
                    cyc = tuple(cyc[i:] + cyc[:i])
                    cycles.append(cyc)
                    cid = len(cycles) - 1
                    for e in cyc:
                        seen_cycle_elems[e] = cid
                    found = cid
                    break
                if v in seen_cycle_elems:
                    found = seen_cycle_elems[v]; break
                if abs(v) <= NLIM and v in done:
                    found = done[v]; break
                if abs(v) > CAP:
                    escaped += 1; found = -1; break
                pos[v] = len(path); path.append(v)
                v = Tb(v, b)
            for e in path:
                if abs(e) <= NLIM:
                    done[e] = found
    return cycles, escaped
KS = [b for b in range(1, 100) if gcd(b, 6) == 1]
CY = {}
for b in KS:
    cyc, esc = census_b(b)
    check(esc == 0, "escape b=%d" % b)
    for c in cyc:
        L = len(c)
        for i in range(L):
            check(Tb(c[i], b) == c[(i + 1) % L], "cycle verify b=%d" % b)
    CY[b] = cyc
cm5, esc = census_b(-5)
check(esc == 0 and sorted(tuple(-t for t in c) for c in CY[5]) == sorted(cm5), "C_{-5} = -C_5")
print("census done %s; no escapes; all cycles verified; C_{-5} = -C_5" % tm())
c = {b: len(CY[b]) for b in KS}
print("  c(b):", c)
expect = {1: 4, 5: 9, 7: 5, 11: 7, 13: 13, 17: 8, 19: 6, 23: 15, 25: 12, 29: 9, 31: 6, 35: 12, 37: 7, 41: 5, 43: 5, 47: 11, 49: 7, 53: 5,
          55: 15, 59: 11, 61: 6, 65: 19, 67: 5, 71: 11, 73: 7, 77: 12, 79: 8, 83: 7, 85: 15, 89: 5, 91: 17, 95: 15, 97: 8}
check(c == expect, "c(b) table matches the lane .out")
lens = {b: sorted(len(t) for t in CY[b]) for b in KS}
check(lens[5] == [1, 1, 1, 2, 3, 3, 7, 17, 17] and lens[65] == [1, 1, 1, 1, 2, 3, 3, 5, 5, 5, 5, 5, 5, 5, 7, 12, 15, 17, 17], "length lists")
vals = sorted(set(c.values()))
print("  attained c(b):", vals, "; holes in [1,19]:", [v for v in range(1, 20) if v not in vals], "; b with c=7:", [b for b in KS if c[b] == 7])
check(vals == [4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 17, 19] and [b for b in KS if c[b] == 7] == [11, 37, 49, 73, 83], "spectrum")
alll = sorted(set(l for b in KS for l in lens[b]))
print("  cycle lengths:", alll)
check(alll == [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18, 20, 22, 24, 26, 32, 36, 41, 56] and all(7 in lens[b] for b in KS), "lengths")
# completeness bound recomputed from scratch
def Lstar(b, N):
    L = 0
    while True:
        L += 1
        worst = Fraction(0)
        K = L
        while True:
            d = 2**K - 3**L
            if d != 0:
                worst = max(worst, Fraction(b * 2**(K - L) * (3**L - 2**L), abs(d)))
            if d > 3 * 3**L:
                break
            K += 1
        if worst > N:
            return L - 1
Ls = {b: Lstar(b, NLIM) for b in KS}
print("  L*(b) at N=2*10^5:", Ls)
check(Ls[1] == 23 and Ls[5] == 21 and Ls[7] == 21 and all(Ls[b] == 16 for b in (11, 13, 17, 19)) and all(Ls[b] == 11 for b in KS if b >= 23), "L*")
# bound sanity on a brute universe: every found cycle element satisfies the bound
def word(n, b, L):
    w = []; vv = n
    for _ in range(L):
        y = 3 * vv + b; e = v2(y); w.append(e); vv = y >> e
    return tuple(w), vv
def Bw(w):
    L = len(w); K = 0; B = 0
    for i, k in enumerate(w):
        B += 3**(L - 1 - i) * 2**K; K += k
    return B, K
for b in KS:
    for cyc in CY[b]:
        L = len(cyc)
        for n0 in cyc:
            w, back = word(n0, b, L); B, K = Bw(w)
            check(back == n0 and Fraction(b * B, 2**K - 3**L) == n0, "gate identity b=%d" % b)
            check(abs(n0) <= Fraction(b * 2**(K - L) * (3**L - 2**L), abs(2**K - 3**L)), "bound on element")
print("  gate identity n = bB/(2^K-3^L) and the element bound hold on every element of every found cycle")
# F2 fixed points
for b in KS:
    pred = sum(1 for K in range(1, 8) if b % abs(2**K - 3) == 0)
    check(pred == sum(1 for t in CY[b] if len(t) == 1), "fixed points b=%d" % b)
check(sorted(t[0] for t in CY[65] if len(t) == 1) == [-65, 5, 13, 65], "b=65 fixed points")
print("  F2: #fixed points = #{K>=1 : (2^K-3) | b} for all b<=99; b=65: [-65, 5, 13, 65]")
# C5 scaling and superadditivity (true form)
for b in KS:
    for c1 in CY[1]:
        sc = sorted(b * t for t in c1)
        check(any(sorted(t) == sc for t in CY[b]), "b*C_1 in C_b, b=%d" % b)
check([b for b in KS if c[b] == 4] == [1], "only b=1 has c=4")
viol = [(b1, b2) for b1 in KS for b2 in KS if b1 < b2 and gcd(b1, b2) == 1 and b1 * b2 <= 99 and c[b1 * b2] < c[b1] + c[b2] - c[1]]
check(viol == [], "superadditivity with c(1)=4 in-universe")
print("  b*C_1 subset C_b (all b); only b=1 has c(b)=4; superadditivity c(b1b2) >= c(b1)+c(b2)-c(1): no violation, b1b2<=99")
# F5
newc = 0; need_b = 0
for b in KS:
    for cyc in CY[b]:
        if any(sorted(cyc) == sorted(b * t for t in c1) for c1 in CY[1]):
            continue
        newc += 1
        w, _ = word(cyc[0], b, len(cyc)); B, K = Bw(w)
        if B % (2**K - 3**len(cyc)) != 0:
            need_b += 1
print("  F5: %d cycles beyond b*C_1, %d have (2^K-3^L) !| B" % (newc, need_b))
check(newc == 175 and need_b == 175, "F5")
# max and doubly representable
reps = {}
for K in range(1, 70):
    for L in range(1, 45):
        d = 2**K - 3**L
        if 1 <= abs(d) <= 99 and gcd(d, 6) == 1:
            reps.setdefault(d, []).append((K, L))
mx = max(c.values()); arg = [b for b in KS if c[b] == mx]
print("  max c(b) = %d at b=%s ; representations of 65: %s, of -65: %s ; doubly representable |b|: %s" % (mx, arg, reps.get(65), reps.get(-65), sorted(set(abs(d) for d, l in reps.items() if len(l) >= 2))))
check(mx == 19 and arg == [65] and reps.get(65) is None and reps[-65] == [(4, 4)], "b=65")
check(len(reps[-1]) == 2 and c[1] == 4, "b=-1 doubly representable with the MINIMUM count 4 (sharper witness against 'doubly representable => maximal')")
# necklace theorem re-verification
def necklaces(K, L):
    seen = set(); cnt = 0
    for cuts in itertools.combinations(range(1, K), L - 1):
        w = tuple(b2 - a2 for a2, b2 in zip((0,) + cuts, cuts + (K,)))
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
        hit = set(); tot = 0
        for cuts in itertools.combinations(range(1, K), L - 1):
            w = tuple(b2 - a2 for a2, b2 in zip((0,) + cuts, cuts + (K,)))
            B, KK = Bw(w)
            ww, back = word(B, d, L)
            check(ww == w and back == B, "realization b=%d w=%s" % (d, w))
            ids = [i for i, cyc in enumerate(CY[d]) if B in cyc]
            check(len(ids) == 1, "census hit")
            hit.add(ids[0]); tot += 1
        check(len(hit) == necklaces(K, L), "necklace count b=%d" % d)
        if d in (13, 47, 5):
            print("  b=%d = 2^%d-3^%d: %d compositions realized on %d distinct cycles = necklace count" % (d, K, L, tot, len(hit)))
print("  necklace theorem re-verified for the 12 positive representations b<=99")

# ------------------------------------------------------------------ F G-map
hdr("F. Greedy 3-adic G-map stopping time: primes vs composites (m<=10^6)")
MG = 10**6
def G(m):
    k = 0
    while m % 9 not in (4, 7):
        m *= 2; k += 1
    return (m - 1) // 3
steps = np.full(MG + 1, -1, dtype=np.int32); steps[1] = 0
for m in range(2, MG + 1):
    if m % 3 == 0 or steps[m] >= 0:
        continue
    path = []; vv = m
    while not (vv <= MG and steps[vv] >= 0):
        path.append(vv); vv = G(vv)
        check(len(path) < 10**5, "runaway")
    b0 = int(steps[vv])
    for i, q in enumerate(reversed(path)):
        if q <= MG:
            steps[q] = b0 + i + 1
isP = np.zeros(MG + 1, dtype=np.bool_); isP[list(primerange(2, MG + 1))] = True
ms = np.arange(2, MG + 1); ms = ms[ms % 3 != 0]; st = steps[ms]; pm = isP[ms]
check(int(st.min()) >= 0, "resolved")
def tv(a, b):
    ha = np.bincount(a, minlength=200) / len(a); hb = np.bincount(b, minlength=200) / len(b)
    return 0.5 * float(np.abs(ha - hb).sum())
print("  primes %d composites %d ; pooled means %.4f %.4f ; pooled TV %.4f" % (pm.sum(), (~pm).sum(), st[pm].mean(), st[~pm].mean(), tv(st[pm], st[~pm])))
check(int(pm.sum()) == 78497 and int((~pm).sum()) == 588169, "counts")
check(abs(tv(st[pm], st[~pm]) - 0.0174) < 6e-4 and abs(st[pm].mean() - 28.7373) < 1e-3 and abs(st[~pm].mean() - 29.0404) < 1e-3, "pooled")
blk = []
for lo, hi in ((2**17, 2**18), (2**18, 2**19), (2**19, 2**20)):
    mk = (ms >= lo) & (ms < hi)
    t = tv(st[mk & pm], st[mk & ~pm]); blk.append(t)
    print("  block [2^%d,2^%d): TV %.4f ; means %.4f %.4f" % (round(math.log2(lo)), round(math.log2(hi)), t, st[mk & pm].mean(), st[mk & ~pm].mean()))
check(all(abs(a - b) < 6e-4 for a, b in zip(blk, (0.0187, 0.0143, 0.0122))), "block TVs")
mk = (ms >= 2**19) & (ms < 2**20)
t18 = tv(st[ms % 9 == 1], st[ms % 9 == 8]); t18b = tv(st[mk & (ms % 9 == 1)], st[mk & (ms % 9 == 8)])
t45 = tv(st[mk & (ms % 9 == 4)], st[mk & (ms % 9 == 5)])
tsz = tv(st[(ms >= 2**17) & (ms < 2**18)], st[mk])
print("  classes 1 vs 8 mod 9: TV %.4f pooled, %.4f in top block ; 4 vs 5 mod 9 top block: TV %.4f ; size control TV %.4f" % (t18, t18b, t45, tsz))
check(abs(t18 - 0.0072) < 6e-4 and abs(t18b - 0.0088) < 6e-4 and abs(t45 - 0.2953) < 6e-4 and abs(tsz - 0.2085) < 6e-4, "residue / size controls")
f1 = [tv(st[mk & (ms % 9 == r) & pm], st[mk & (ms % 9 == r) & ~pm]) for r in (1, 2, 4, 5, 7, 8)]
print("  F1 per-class TV in top block:", [round(t, 4) for t in f1])
check(max(f1) < 0.04, "F1")
units = [a for a in range(1, 81) if a % 3]
hp = np.bincount(ms[pm] % 81, minlength=81)[units] / pm.sum(); hc = np.bincount(ms[~pm] % 81, minlength=81)[units] / (~pm).sum()
print("  primes vs composites mod 81: TV %.4f" % (0.5 * float(np.abs(hp - hc).sum())))

# ------------------------------------------------------------------ G sandwich
hdr("G. Sandwich N_PS - N_SP, square-driven pairs, S_1 - S_5")
Om = np.zeros(NS + 1, dtype=np.uint8)
for p in primerange(2, NS + 1):
    q = p
    while q <= NS:
        Om[q::q] += 1; q *= p
for K in (10**5, 10**6):
    k = np.arange(1, K + 1, dtype=np.int64)
    a = Om[6 * k - 1]; b = Om[6 * k + 1]
    dPS = int(((a == 1) & (b == 2)).sum() - ((a == 2) & (b == 1)).sum())
    sq = [p for p in primerange(5, math.isqrt(6 * K + 1) + 1)]
    A = sum(1 for p in sq if isprime(p * p - 2))
    x6 = 6 * K + 1
    S1 = int((Om[np.arange(1, x6 + 1, 6)] == 2).sum()); S5 = int((Om[np.arange(5, x6 + 1, 6)] == 2).sum())
    print("  K=%d: N_PS-N_SP = %d ; pairs (p^2-2 prime, p^2) = %d ; residual %d ; S_1-S_5 = %d ; pi'(sqrt x) = %d" % (K, dPS, A, dPS - A, S1 - S5, len(sq)))
    check((K, dPS, A, S1 - S5, len(sq)) in ((10**5, 41, 46, 66, 135), (10**6, 412, 99, 244, 361)), "sandwich numbers at K=%d" % K)
print("  (41, 412 match collatz_mod6_20260917_sandwich_bias.md sec.1; 66, 244 and pi'/2 = 67.5, 180.5 match its Theorem 2.2 table)")

# ------------------------------------------------------------------ H 3n-5
hdr("H. 3n-5 = T(n-2)")
for n in range(-100001, 100002, 2):
    check(Tb(n, -5) == Tb(n - 2, 1), "T_{-5}(n)=T(n-2)")
for r, img, col in ((1, 8, 2), (3, 2, 5), (5, 5, 8)):
    for jj in range(-50, 51):
        n = 6 * jj + r
        check(((3 * n - 5) // 2) % 9 == img and ((3 * n + 1) // 2) % 9 == col, "row images")
print("T_{-5}(n)=T(n-2) for odd |n|<=10^5; row images (3n-5)/2 = 8,2,5 mod 9 vs Collatz 2,5,8")
pos = sorted(c for c in cm5 if c[0] > 0)
check(len(cm5) == 9 and len(pos) == 3 and [len(t) for t in pos] == [1, 2, 7], "nine 3n-5 cycles, three positive")
print("  3n-5 cycles: 9 (3 positive: lengths 1,2,7), = negatives of the b=5 cycles")

# ------------------------------------------------------------------ I Fano / octonions by norm multiplicativity
hdr("I. Fano orientations giving a composition algebra (norm test, independent of the alternativity test)")
p, q, r = 2, 3, 5
S7 = [d for d in (2, 3, 5, 6, 10, 15, 30)]
def vec(d):
    return (d % 2 == 0) * 4 + (d % 3 == 0) * 2 + (d % 5 == 0)
lines = sorted(set(tuple(sorted((vec(a), vec(b), vec(a) ^ vec(b)))) for a in S7 for b in S7 if a < b))
check(len(lines) == 7 and all(all(1 <= t <= 7 for t in l) for l in lines), "7 lines")
sqlines = sorted(set(tuple(sorted((a, b, [cc for cc in S7 if vec(cc) == vec(a) ^ vec(b)][0]))) for a in S7 for b in S7 if a < b))
check(all(math.isqrt(a * b * cc)**2 == a * b * cc for a, b, cc in sqlines), "line products are squares")
print("  square-product lines:", sqlines)
rng = np.random.default_rng(3)
pairs = [(rng.integers(-5, 6, 8), rng.integers(-5, 6, 8)) for _ in range(12)]
good = 0
for orient in itertools.product((0, 1), repeat=7):
    mult = {}
    for i in range(8):
        mult[(0, i)] = (1, i); mult[(i, 0)] = (1, i)
    for i in range(1, 8):
        mult[(i, i)] = (-1, 0)
    for (a, b, cc), o in zip(lines, orient):
        cyc = [(a, b, cc), (b, cc, a), (cc, a, b)] if o == 0 else [(a, cc, b), (cc, b, a), (b, a, cc)]
        for (xx, yy, zz) in cyc:
            mult[(xx, yy)] = (1, zz); mult[(yy, xx)] = (-1, zz)
    ok = True
    for X, Y in pairs:
        Z = np.zeros(8, dtype=np.int64)
        for i in range(8):
            for j in range(8):
                s, kk = mult[(i, j)]
                Z[kk] += s * X[i] * Y[j]
        if int((Z * Z).sum()) != int((X * X).sum()) * int((Y * Y).sum()):
            ok = False; break
    good += ok
print("  orientations with multiplicative norm on 12 random integer pairs: %d of 128" % good)
check(good == 16, "16 composition orientations")

# ------------------------------------------------------------------ J divisor facts
hdr("J. Divisor-box facts and ord_m(2)=8")
from sympy import divisors
def FSU(N):
    pr = [d for d in divisors(N) if 1 < d < N]
    return len(pr), sum(1 for d in pr if all(e == 1 for e in factorint(d).values())), sum(1 for d in pr if isprime(d))
check(FSU(60) == (10, 7, 3) and FSU(7) == (0, 0, 0) and FSU(343) == (2, 1, 1) and FSU(7007) == (10, 7, 3), "F,S,U profiles")
for N in range(2, 5001):
    f = factorint(N); dn = len(divisors(N)); sf = all(e == 1 for e in f.values())
    check(dn <= 2**sum(f.values()) and (dn == 2**sum(f.values())) == sf and (dn == 2**len(f)) == sf, "d(N) vs 2^Omega at %d" % N)
print("  F=S+U profiles (0,0,0),(2,1,1),(10,7,3) at p, p^3, p^2qr; d(N)<=2^Omega with equality iff squarefree, N<=5000")
check(len(divisors(24)) == 8 and len(divisors(5400)) == 48 and sum(factorint(5400).values()) == 8, "24, 5400")
o8 = [m for m in range(3, 2000, 2) if n_order(2, m) == 8]
check(o8 == [17, 51, 85, 255] and [n_order(2, m) for m in (3, 9, 7, 63, 21, 27)] == [2, 6, 3, 6, 6, 18], "orders")
print("  ord_m(2)=8, odd m<2000:", o8, "; orders at 3,9,7,63,21,27:", [n_order(2, m) for m in (3, 9, 7, 63, 21, 27)])
print()
print("AUDIT COMPLETE: all checks passed", tm())
