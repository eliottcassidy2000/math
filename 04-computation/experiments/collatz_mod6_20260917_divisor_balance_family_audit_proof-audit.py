#!/usr/bin/env python3
"""Independent adversarial audit (proof-audit lens) of lane divisor_balance_family.

Does NOT import the explorer's script.  Recomputes every load-bearing number by a
different construction:
  * p^2qr and p^3 indicators to 10^7 are built by DIRECT GENERATION of the triples
    (p, q<r) and of cubes, not by an omega/Omega sieve;
  * primes by an independent sieve;
  * the 16 linear cells by a brute-force profile window (r<=8, exponents<=20) AND by an
    independent multiplicative-partition enumeration of T(r) on the proved support bounds;
  * the k-free cells by brute force for k=1..8;
  * partitions per Omega, sandwich census, side law.
Also runs the hostile probes the audit found: Omega-injectivity of the INFINITE cell (2,0),
and the k=1 shape count.
Every check is an explicit RuntimeError (active under -O).  RAM < 400 MB, ~30 s.
"""
from __future__ import annotations

import math
import sys
import time
from itertools import combinations_with_replacement, product
from math import isqrt, prod

import numpy as np

T0 = time.time()


def require(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def hdr(s):
    print("\n" + "=" * 72 + "\n" + s + "\n" + "=" * 72)


# ---------------------------------------------------------------- profile formulas (DB1)
def F_of(prof):
    return prod(a + 1 for a in prof) - 2


def S_of(prof):
    return (1 << len(prof)) - 1 - (1 if all(a == 1 for a in prof) else 0)


def U_of(prof):
    return len(prof) - (1 if tuple(prof) == (1,) else 0)


def Sk_of(prof, k):
    return prod(min(a + 1, k) for a in prof) - 1 - (1 if all(a < k for a in prof) else 0)


def D_of(prof):
    return F_of(prof) - S_of(prof) - U_of(prof)


# independent direct-divisor audit of the formulas on a smaller window (different code path)
def factor_trial(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


hdr("0. Direct divisor audit of F, S_k, U formulas on 2<=N<=30000 (independent code)")
for n in range(2, 30001):
    fac = factor_trial(n)
    prof = tuple(sorted(fac.values(), reverse=True))
    divs = [d for d in range(2, n) if n % d == 0]
    F = len(divs)
    U = sum(1 for d in divs if d in fac)
    require(F == F_of(prof) and U == U_of(prof), f"F/U at {n}")
    for k in range(1, 8):
        Sk = sum(1 for d in divs if all(d % (p ** k) for p in fac))
        require(Sk == Sk_of(prof, k), f"S_{k} at {n}")
print("F, U, S_1..S_7 formulas confirmed by direct enumeration for 2<=N<=30000")
print(f"[t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- 1. density by direct generation
hdr("1. Shape indicators to 10^7 by DIRECT GENERATION (not by an omega/Omega sieve)")
X = 10_000_000
sieve = np.ones(X + 1, dtype=bool)
sieve[:2] = False
for i in range(2, isqrt(X) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
primes = np.flatnonzero(sieve)
plist = primes.tolist()
print(f"pi(10^7) = {len(plist)}")
require(len(plist) == 664579, "pi(10^7) must be 664579")

# prime cubes
is_cube = np.zeros(X + 1, dtype=bool)
for p in plist:
    if p ** 3 > X:
        break
    is_cube[p ** 3] = True

# p^2 q r, q<r, q,r != p: generate every triple exactly once
is_p2qr = np.zeros(X + 1, dtype=bool)
ntriples = 0
for p in plist:
    pp = p * p
    if pp * 2 * 3 > X:
        break
    for qi, q in enumerate(plist):
        rmax = X // (pp * q)
        if rmax <= q:
            break
        if q == p:
            continue
        hi = int(np.searchsorted(primes, rmax, side="right"))
        lo = qi + 1
        if hi <= lo:
            continue
        rs = primes[lo:hi]
        vals = pp * q * rs
        is_p2qr[vals] = True
        cnt = hi - lo
        if q < p <= rmax:
            is_p2qr[pp * q * p] = False   # r = p is not allowed; p^3 q is never a p^2qr number
            cnt -= 1
        ntriples += cnt
require(int(is_p2qr.sum()) == ntriples, "each p^2qr generated exactly once")
print(f"generated {ntriples} numbers p^2qr <= 10^7 (each exactly once)")

cum_p = np.cumsum(sieve, dtype=np.int64)
cum_c = np.cumsum(is_cube, dtype=np.int64)
cum_4 = np.cumsum(is_p2qr, dtype=np.int64)

expected = {10 ** 5: (9592, 14, 9346), 10 ** 6: (78498, 25, 87338), 10 ** 7: (664579, 47, 804249)}
P2 = 0.4522474200410654985
print("x        primes   p^3   p^2qr   ratio   asym    asym/exact")
for x, (e1, e3, e4) in expected.items():
    n1, n3, n4 = int(cum_p[x]), int(cum_c[x]), int(cum_4[x])
    require((n1, n3, n4) == (e1, e3, e4), f"shape counts at {x}: {(n1,n3,n4)} vs claimed {(e1,e3,e4)}")
    asym = P2 * x * math.log(math.log(x)) / math.log(x)
    print(f"{x:<8} {n1:>7} {n3:>5} {n4:>7}  {n4/n1:.4f}  {asym:>8.0f}  {asym/n4:.4f}")
print("CONFIRMED: C1a shape counts at 10^5, 10^6, 10^7")

# crossover and sign changes of the strict lead
diff = cum_4 - cum_p
lead = diff > 0
first = int(np.argmax(lead))
require(lead[first], "no crossover")
require(first == 145119, f"first crossover {first}")
require(3 * 13 * 61 ** 2 == 145119 and bool(is_p2qr[145119]), "145119 = 3*13*61^2 is p^2qr")
require((int(cum_4[first]), int(cum_p[first])) == (13433, 13432), "counts at crossover")
flips = np.flatnonzero(lead[1:] != lead[:-1]) + 1
print(f"first n with #p2qr > #primes: {first}; counts {int(cum_4[first])} vs {int(cum_p[first])}")
print(f"strict-lead flip points: {flips.tolist()}  (count {len(flips)})")
require(flips.tolist() == [145119, 145121, 145132, 145133, 145138, 145139, 145148, 147709, 147725, 147727, 147908, 147919, 147925], "flip list")
require(bool(lead[147925:].all()), "p^2qr leads on [147925, 10^7]")
require(int(diff[X]) == 139670, f"lead at 10^7 = {int(diff[X])}")
# what happens at the flip points (hostile: each flip must be a prime or a p^2qr)
kinds = ["P" if sieve[n] else ("Q" if is_p2qr[n] else "?") for n in flips.tolist()]
print(f"flip-point kinds (P=prime, Q=p^2qr): {kinds}")
require("?" not in kinds, "every flip point is a prime or a p^2qr number")
naive = math.exp(math.exp(1 / P2))
print(f"naive P(2) loglog x = 1 scale: x ~ {naive:.0f}, ratio {first/naive:.2f}")
print("CONFIRMED: C1c crossover, 13 flips, lead 139670 at 10^7")

# the prime-cube count is pi(floor(x^(1/3)))
for x in (10 ** 5, 10 ** 6, 10 ** 7):
    c = int(round(x ** (1 / 3)))
    while c ** 3 > x:
        c -= 1
    while (c + 1) ** 3 <= x:
        c += 1
    require(int(cum_c[x]) == int(cum_p[c]), "cube count = pi(x^(1/3))")

# numerical sanity for Theorem 1's dominator: S_p(x)/M(x) vs 1/p^2 at x=10^7 for small p
M = P2 * X * math.log(math.log(X)) / math.log(X) / P2   # x loglog x / log x
print("head terms S_p(10^7)/M(10^7) vs 1/p^2 (proof of Theorem 1 says ratio -> 1/p^2, dominated by C'/p^2):")
for p in (2, 3, 5, 7, 11):
    pp = p * p
    # count directly from the generated indicator: n = p^2 q r with v_p(n) = 2 exactly
    idx = np.flatnonzero(is_p2qr[::pp])  # multiples of p^2
    n_mult = idx * pp
    Sp = int(np.count_nonzero((n_mult // pp) % p != 0))
    print(f"  p={p}: S_p/M = {Sp/M:.4f}   1/p^2 = {1/pp:.4f}   ratio {Sp/M*pp:.3f}")
print(f"[t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- 2. linear cells
hdr("2. Linear cells F = alpha S + beta U: brute-force window + independent T(r) factorization")


def mult_partitions(T, r, lo=2):
    """All nondecreasing r-tuples of integers >= lo with product T."""
    if r == 1:
        return [(T,)] if T >= lo else []
    out = []
    f = lo
    while f ** r <= T:
        if T % f == 0:
            out.extend((f,) + rest for rest in mult_partitions(T // f, r - 1, f))
        f += 1
    return out


claimed = {
    (0, 0): {(1,)},
    (0, 1): {(1,), (2,), (1, 1)},
    (0, 2): {(1,), (3,), (2, 1), (1, 1, 1)},
    (0, 3): {(1,), (4,), (3, 1)},
    (1, 1): {(1,), (3,), (2, 1, 1)},
    (1, 2): {(1,), (4,), (2, 2)},
    (1, 3): {(1,), (5,), (2, 2, 1), (2, 1, 1, 1, 1)},
    (2, 1): {(1,), (4,), (4, 1), (2, 2, 1, 1)},
    (2, 2): {(1,), (5,), (3, 2), (5, 1), (4, 1, 1, 1)},
    (2, 3): {(1,), (6,), (6, 1)},
    (3, 0): {(1,), (4,)},
    (3, 1): {(1,), (5,)},
    (3, 2): {(1,), (6,), (4, 2)},
    (3, 3): {(1,), (7,), (3, 3, 1), (7, 1, 1)},
}
RMAX, EMAX = 8, 20
window = {(a, b): set() for a in range(4) for b in range(4)}
nprof = 0
for r in range(1, RMAX + 1):
    for prof in combinations_with_replacement(range(EMAX, 0, -1), r):
        nprof += 1
        F, S, U = F_of(prof), S_of(prof), U_of(prof)
        for a in range(4):
            for b in range(4):
                if F == a * S + b * U:
                    window[(a, b)].add(prof)
print(f"window: {nprof} profiles, r<=8, exponents<=20")
for cell, sols in claimed.items():
    require(window[cell] == sols, f"cell {cell}: window {sorted(window[cell])} vs claimed {sorted(sols)}")
# infinite cells: (1,0) = squarefree or (2); (2,0) = (1) or (3,1^(r-1))
for prof in window[(1, 0)]:
    require(all(a == 1 for a in prof) or prof == (2,), f"(1,0) hostile {prof}")
require(all((tuple([1] * r) in window[(1, 0)]) for r in range(1, RMAX + 1)) and (2,) in window[(1, 0)], "(1,0) content")
require(window[(2, 0)] == {(1,)} | {tuple([3] + [1] * (r - 1)) for r in range(1, RMAX + 1)}, "(2,0) content")
print("CONFIRMED: all 16 cells agree with the claimed profile lists on the window")

# independent enumeration on the PROVED support bounds via multiplicative partitions of T(r)
def T_of(a, b, r):
    return 2 + a * ((1 << r) - 1) + b * r


rbound = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 2, (1, 0): 1, (1, 1): 3, (1, 2): 4, (1, 3): 5}
for a in (2, 3):
    for b in range(4):
        if (a, b) != (2, 0):
            rbound[(a, b)] = 7
# re-derive the alpha in {0,1} bounds from 3*2^(r-1) <= T(r)
for (a, b), rb in rbound.items():
    if a in (0, 1):
        rr = max([r for r in range(1, 40) if 3 * (1 << (r - 1)) <= T_of(a, b, r)] + [0])
        require(rr == rb, f"support bound {(a,b)}: {rr} vs {rb}")
# re-derive the alpha in {2,3} bound: r - k <= v2(T), k <= 3, and check r=8..12 empty
def v2(n):
    n = abs(n)
    return (n & -n).bit_length() - 1 if n else 10 ** 9


for a in (2, 3):
    for b in range(4):
        if (a, b) == (2, 0):
            continue
        for r in range(4, 13):
            T = T_of(a, b, r)
            c = b * r + 2 - a
            require(c != 0 and abs(c) < (1 << r), f"c_r nonzero and < 2^r at {(a,b,r)}")
            require(v2(T) == v2(c), f"v2(T)=v2(c_r) at {(a,b,r)}")
            require(T / (1 << r) <= 4.5, f"T/2^r <= 4.5 at {(a,b,r)}")
        # r <= 3 + log2(3r+1) fails for r >= 8
        require(all(r > 3 + math.log2(3 * r + 1) for r in range(8, 40)), "r<=7 bound")
for cell, sols in claimed.items():
    if cell == (0, 0):
        continue
    a, b = cell
    found = {(1,)}
    for r in range(1, rbound[cell] + 1):
        for fac in mult_partitions(T_of(a, b, r), r):
            if max(fac) >= 3:
                found.add(tuple(sorted((f - 1 for f in fac), reverse=True)))
    for r in range(2, 40):   # squarefree composites
        if (1 - a) * ((1 << r) - 2) == b * r:
            found.add(tuple([1] * r))
    require(found == sols, f"T(r)-enumeration for {cell}: {sorted(found)} vs {sorted(sols)}")
print("CONFIRMED: multiplicative-partition enumeration on the proved r-bounds reproduces every claimed cell")

# prime-power law
for cell, sols in claimed.items():
    a, b = cell
    pp = sorted(p for p in sols if len(p) == 1)
    require(pp == ([(1,)] if a + b == 0 else [(1,), (a + b + 1,)]), f"prime-power law {cell}")
print("CONFIRMED: prime-power law p^(alpha+beta+1)")

# HOSTILE: Omega-injectivity.  The explorer labels INFINITE cells as 'collides' by fiat.
print("\nHOSTILE: Omega values of the (2,0) family p^3 q_1...q_(r-1), r=1..12:",
      [3 + (r - 1) for r in range(1, 13)], "plus 1 for the prime -> pairwise DISTINCT")
omegas_20 = [1] + [sum(tuple([3] + [1] * (r - 1))) for r in range(1, 13)]
require(len(set(omegas_20)) == len(omegas_20), "(2,0) is Omega-injective")
print("  => the claim '(1,1) is the unique prime-cube cell with pairwise distinct Omega' is FALSE:")
print("     (2,0) also has pairwise distinct Omega (Omega(p^3 m) = omega(m)+3 = r+2, injective in r).")
print("     Only (0,2) among {(0,2),(1,1),(2,0)} has an Omega collision (p^3, p^2q, pqr all have Omega=3).")
inj = {}
for cell, sols in claimed.items():
    om = [sum(p) for p in sols]
    inj[cell] = len(set(om)) == len(om)
require(sorted(c for c in inj if inj[c]) == [(0, 0), (1, 1), (2, 1), (2, 3), (3, 0), (3, 1)], "finite injective cells")
require(sorted(c for c in inj if inj[c] and len(claimed[c]) == 3) == [(1, 1), (2, 3)], "3-profile injective")
print("CONFIRMED (finite cells only): Omega-injective finite cells = (0,0),(1,1),(2,1),(2,3),(3,0),(3,1); 3-profile ones (1,1),(2,3)")

# tight profiles (2,1^(r-1)) across the grid
tight = [((a, b), r) for a in range(4) for b in range(4) for r in range(1, 30) if 3 * (1 << (r - 1)) == T_of(a, b, r)]
require(tight == [((0, 1), 1), ((0, 2), 2), ((1, 0), 1), ((1, 1), 3), ((1, 3), 5)], f"tight list {tight}")
print(f"tight profiles in the grid: {tight}")

# Theorem 2 (general finiteness) numeric probe: cells alpha,beta in 4..8, brute force r<=8, exp<=40
print("\nTheorem 2 probe: cells with 4<=alpha,beta<=8, window r<=8, exponents<=40")
big = {}
for r in range(1, 9):
    for prof in combinations_with_replacement(range(40, 0, -1), r):
        F, S, U = F_of(prof), S_of(prof), U_of(prof)
        for a in range(4, 9):
            for b in range(0, 9):
                if F == a * S + b * U:
                    big.setdefault((a, b), set()).add(prof)
for a in range(4, 9):
    for b in range(0, 9):
        sols = big.get((a, b), set())
        require((a + b + 1,) in sols, f"prime-power law fails in {(a,b)}")
        # Theorem 2 predicts r <= k + v2(alpha-2) for beta=0, k <= log_{3/2}(alpha+beta+2)
        kmax = int(math.log(a + b + 2) / math.log(1.5))
        if b == 0:
            rb = kmax + v2(a - 2)
        else:
            rb = max(r for r in range(1, 200) if r <= kmax + math.log2(b * r + 2))
        require(all(len(p) <= rb for p in sols), f"support bound violated in {(a,b)}: {sols}")
print("  all 45 probed cells: prime-power law holds and every solution respects the Theorem 2 support bound")
print(f"[t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- 3. lattice and k-free
hdr("3. Lattice identity and k-free refinement")
# F - S - U = |Q| - (r+1) for every nonsquarefree profile in the window (identity check)
for r in range(1, 7):
    for prof in combinations_with_replacement(range(8, 0, -1), r):
        if all(a == 1 for a in prof):
            continue
        box = list(product(*[range(a + 1) for a in prof]))
        Q = [e for e in box if max(e) >= 2]
        require(D_of(prof) == len(Q) - (r + 1), f"lattice identity at {prof}")
print("CONFIRMED: F-S-U = |Q|-(r+1) for all nonsquarefree profiles, r<=6, exp<=8")
for r in range(1, 12):
    prof = tuple([2] + [1] * (r - 1))
    q = 1 << (r - 1)
    require(D_of(prof) == q - r - 1, "tight profile defect")
    require((q == r + 1) == (r == 3), "|Q|=r+1 iff r=3")
    if r >= 3:
        require(D_of(prof) == sum(math.comb(r - 1, j) for j in range(2, r - 1)), "unreached count")
print("CONFIRMED: tight-profile |Q| = 2^(r-1), = r+1 iff r=3, defect = #{d|m: 2<=omega(d)<=r-2}, r<=11")
# N=60 explicit
Q60 = [d for d in range(1, 61) if 60 % d == 0 and any(d % (p * p) == 0 for p in (2, 3, 5))]
require(Q60 == [4, 12, 20, 60], "Q(60)")
require(sorted(2 * math.lcm(2, l) for l in (2, 3, 5)) == [4, 12, 20], "atom map at 60")
require(sorted(d // 2 for d in Q60) == [2, 6, 10, 30], "d/p hostile")
print("CONFIRMED: N=60 example (Q, atom map, apex, d/p hostile)")

# k-free cells, brute force k=1..8
def sk_claimed(k):
    if k == 1:
        return {(1,), (2,), (1, 1)}
    s = {(1,), (k + 1,), (k, 1, 1)}
    if k >= 3:
        s.add((k, 2))
    return s


for k in range(1, 9):
    hits = set()
    for r in range(1, 8):
        for prof in combinations_with_replacement(range(16, 0, -1), r):
            if F_of(prof) == Sk_of(prof, k) + U_of(prof):
                hits.add(prof)
    require(hits == sk_claimed(k), f"k={k}: {sorted(hits)} vs {sorted(sk_claimed(k))}")
print("CONFIRMED: F = S_k + U shapes for k=1..8 (window r<=7, exp<=16)")
print("  NOTE: k=1 ALSO has exactly three shapes (p, p^2, pq); 'k=2 is the exceptional 3-shape case' is only true among k>=2")
# Theorem 3 boundary values
for prof, dk in (((3,), -1), ((2, 1, 1), -3)):
    require(F_of(prof) - Sk_of(prof, 3) - U_of(prof) == dk, "k=3 boundary")
# direct integer counts of solutions N<=200000 for k=1..6 (claimed 63214, 36372, 25045, 21096, 19442, 18684)
spf = list(range(200001))
for i in range(2, isqrt(200000) + 1):
    if spf[i] == i:
        for j in range(i * i, 200001, i):
            if spf[j] == j:
                spf[j] = i
counts = [0] * 7
cellcounts = {c: 0 for c in claimed}
c10 = c20 = 0
for n in range(2, 200001):
    m, f = n, {}
    while m > 1:
        p = spf[m]
        while m % p == 0:
            f[p] = f.get(p, 0) + 1
            m //= p
    prof = tuple(sorted(f.values(), reverse=True))
    for k in range(1, 7):
        counts[k] += prof in sk_claimed(k)
    for c in claimed:
        cellcounts[c] += prof in claimed[c]
    c10 += (all(a == 1 for a in prof) or prof == (2,))
    c20 += (prof == (1,) or prof == tuple([3] + [1] * (len(prof) - 1)))
require(counts[1:] == [63214, 36372, 25045, 21096, 19442, 18684], f"k-counts {counts[1:]}")
claimed_counts = {(0, 0): 17984, (0, 1): 63214, (0, 2): 68867, (0, 3): 22140, (1, 1): 36372, (1, 2): 18122,
                  (1, 3): 21629, (2, 1): 21840, (2, 2): 20334, (2, 3): 18493, (3, 0): 17992, (3, 1): 17989,
                  (3, 2): 18042, (3, 3): 18565}
for c, v in claimed_counts.items():
    require(cellcounts[c] == v, f"count {c}: {cellcounts[c]} vs {v}")
require(c10 == 121666 and c20 == 32827, f"(1,0),(2,0) counts {c10},{c20}")
print("CONFIRMED: #N<=2e5 columns for all 16 cells and for k=1..6")
print(f"[t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- 4. Omega partitions + census
hdr("4. Partitions per Omega, D witness, sandwich census")


def partitions(n, m=None):
    m = n if m is None else m
    if n == 0:
        yield ()
        return
    for first in range(min(n, m), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


pk = [len(list(partitions(k))) for k in range(1, 9)]
require(pk == [1, 2, 3, 5, 7, 11, 15, 22], "p(k)")
nsol = [sum(1 for p in partitions(k) if D_of(p) == 0) for k in range(1, 9)]
require(nsol == [1, 0, 1, 1, 0, 0, 0, 0], "solutions per Omega")
require(sorted(set(D_of(p) for p in partitions(4))) == [-4, 0, 1, 2], "D over Omega=4")
require(D_of((2, 2)) == 2 and D_of((3, 1)) == 1, "witness")
print("CONFIRMED: p(k) table, 1,0,1,1,0,0,0,0 solutions, witness p^2q^2 (D=2) vs p^3q (D=1)")

K = 1_000_000
ks = np.arange(1, K + 1, dtype=np.int64)
L, R = 6 * ks - 1, 6 * ks + 1
res = {}
for name, arr in (("prime", sieve), ("p^3", is_cube), ("p^2qr", is_p2qr)):
    res[name] = (int(arr[L].sum()), int(arr[R].sum()))
    print(f"  {name:<6} left {res[name][0]:>7} right {res[name][1]:>7}")
require(res == {"prime": (206502, 206345), "p^3": (21, 19), "p^2qr": (35061, 34876)}, f"census {res}")
sol = sieve | is_cube | is_p2qr
both = int((sol[L] & sol[R]).sum())
twin = int((sieve[L] & sieve[R]).sum())
require((both, twin) == (55613, 37915), f"both/twin {(both, twin)}")
# side law, including the boundary primes 2 and 3 (neither cube is an endpoint)
for p in plist:
    c = p ** 3
    if c > 6 * K + 1:
        break
    require((c % 6 == 5) == (p % 6 == 5) and (c % 6 == 1) == (p % 6 == 1), f"side law at p={p}")
require(8 % 6 == 2 and 27 % 6 == 3, "p=2,3 cubes are not endpoints")
cl = [p ** 3 for p in plist if p ** 3 <= 6 * K - 1 and p % 6 == 5]
cr = [p ** 3 for p in plist if p ** 3 <= 6 * K + 1 and p % 6 == 1]
require(len(cl) == 21 and len(cr) == 19 and cl[:4] == [125, 1331, 4913, 12167] and cr[:4] == [343, 2197, 6859, 29791], "cube lists")
print(f"CONFIRMED: census, both={both}, twin={twin}, mixed={both-twin}; side law for all p with p^3 <= 6*10^6+1 (p=2,3 excluded automatically)")
print(f"[t={time.time()-T0:.1f}s]")
print("\nALL AUDIT CHECKS PASSED")
