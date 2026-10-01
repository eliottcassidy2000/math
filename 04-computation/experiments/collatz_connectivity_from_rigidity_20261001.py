#!/usr/bin/env python3
"""Trying to prove Collatz connectivity from rigidity (opus S15, eighteenth note, 2026-10-01).

Non-shortcut maps T_{a,b}(n) = n/2 (even), a n + b (odd); Collatz = (3, 1); the minus sheet = (3, -1).
Checks:
  A. barrier: the rigidity conclusion (depth-D backward trees identify vertices) holds for 3n-1 and 5n+1,
     which have three components each (so rigidity cannot imply connectivity)
  B. two-place density: explicit n in every class mod 2^a 3^b whose orbit reaches 1 (Terras reduction +
     the family 2^f (4^j - 1)/3), verified by direct iteration
  C. the quantifier exchange on the minus sheet: delta_k(n) = least element of the basin of {1} with
     n's parity and n's class mod 3^k grows without bound for n = 5, 17 (genuine second components)
  D. the backward sieve on a minimal counterexample: fraction of odd 3-adic classes whose backward tree
     has no predecessor with multiplier 2^e/3^j < 1 within depth D
  E. the two-sided trap is realised on the minus sheet: 5 and 17 are the minima of their components
Reproduce: python 04-computation/experiments/collatz_connectivity_from_rigidity_20261001.py
"""
import sys
import math
from functools import lru_cache
from collections import defaultdict

sys.setrecursionlimit(1000000)
FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def vp(n, p):
    if n == 0:
        return 10 ** 9
    n = abs(n)
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


_canon = {}


def cid(children):
    key = tuple(sorted(children))
    v = _canon.get(key)
    if v is None:
        v = len(_canon)
        _canon[key] = v
    return v


LEAF = cid(())

# ---------------- A. the barrier: rigidity also holds for 3n-1 and 5n+1 ----------------
print("== A. backward trees identify vertices for 3n-1 and 5n+1 (both have three cycles) ==")


def make_int_tree(a, b):
    @lru_cache(maxsize=None)
    def tree(n, d):
        if d == 0:
            return LEAF
        ch = [tree(2 * n, d - 1)]
        if n % 2 == 0 and (n - b) % a == 0 and (n - b) // a >= 1 and ((n - b) // a) % 2 == 1:
            ch.append(tree((n - b) // a, d - 1))
        return cid(ch)
    return tree


def cycles_of(a, b, N=3000, cap=10 ** 30):
    found = set()
    for n0 in range(1, N + 1):
        seen = {}
        n = n0
        while n < cap and n not in seen:
            seen[n] = 1
            n = n // 2 if n % 2 == 0 else a * n + b
        if n < cap:
            c = [n]
            m = n // 2 if n % 2 == 0 else a * n + b
            while m != n:
                c.append(m)
                m = m // 2 if m % 2 == 0 else a * m + b
            found.add(min(c))
    return sorted(found)


for (a, b, NMAX, D) in ((3, -1, 2000, 25), (5, 1, 624, 33)):
    tree = make_int_tree(a, b)
    ids = {}
    dup = []
    for n in range(1, NMAX + 1):
        if n % a == 0:
            continue
        c = tree(n, D)
        if c in ids:
            dup.append((ids[c], n))
        ids[c] = n
    cyc = cycles_of(a, b)
    print(f"  {a}n{b:+d}: depth-{D} backward trees of the {len(ids)} vertices n <= {NMAX} prime to {a}: pairwise distinct = {not dup}; cycles (least elements) = {cyc}")
    check(not dup and len(cyc) >= 3, f"A: {a}n{b:+d} is backward-rigid on n <= {NMAX} yet has {len(cyc)} cycles (>= 3 components)")

# ---------------- B. two-place density of the trivial component ----------------
print("\n== B. every class mod 2^a 3^b contains a number reaching 1 (explicit construction) ==")


def reaches_one(n, cap=10 ** 6):
    k = 0
    while n != 1 and k < cap:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        k += 1
    return n == 1


def shortcut(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


_logtab = {}


def log4_table(T):
    """j in [1, T] indexed by 4^j mod 3T (4 has order T = 3^c modulo 3T = 3^(c+1))"""
    if T not in _logtab:
        tab = {}
        y = 1
        for j in range(1, T + 1):
            y = y * 4 % (3 * T)
            tab[y] = j
        _logtab[T] = tab
    return _logtab[T]


def construct(r, a, b):
    """n = r (mod 2^a 3^b) with T1^a(n) = (4^j - 1)/3 (which reaches 1), built by the Terras reduction."""
    M2, M3 = 2 ** a, 3 ** b
    r0 = r % M2
    if r0 == 0:
        r0 = M2
    # the first a shortcut steps are fixed on the class r0 mod 2^a: T1^a(r0 + 2^a t) = A0 + 3^p t
    x, p = r0, 0
    for _ in range(a):
        if x % 2:
            p += 1
        x = shortcut(x)
    A0 = x
    t0 = ((r - r0) * pow(M2, -1, M3)) % M3 if b > 0 else 0   # n = r0 + 2^a t = r (mod 3^b) iff t = t0
    T = 3 ** (p + b)
    lower = A0 + 3 ** p * t0
    target = lower % T                      # need (4^j - 1)/3 = target (mod T), i.e. 4^j = 3 target + 1 (mod 3T)
    j = log4_table(T)[(3 * target + 1) % (3 * T)]
    m = (4 ** j - 1) // 3
    while m < lower:
        j += T
        m = (4 ** j - 1) // 3
    t = (m - A0) // 3 ** p
    return r0 + M2 * t, m, j


ok = True
count = 0
for a in range(0, 6):
    for b in range(0, 4):
        M = 2 ** a * 3 ** b
        for r in range(M):
            n, m, j = construct(r, a, b)
            x = n
            for _ in range(a):
                x = shortcut(x)
            # m = (4^j - 1)/3 is odd and 3m + 1 = 4^j, so m reaches 1 in 2j + 1 steps
            if n <= 0 or n % M != r or x != m or m % 2 != 1 or 3 * m + 1 != 4 ** j:
                ok = False
                print("  failed", r, a, b)
                break
            count += 1
check(ok, f"B: for every residue class mod 2^a 3^b (a <= 5, b <= 3; {count} classes) the constructed n lies in the class and T1^a(n) = (4^j-1)/3, which reaches 1")

# ---------------- C. the quantifier exchange on the minus sheet ----------------
print("\n== C. delta_k(n) on the minus sheet 3n-1: least element of the basin of {1} with n's parity and class mod 3^k ==")
N = 3_000_000
basin = bytearray(N + 1)       # 1, 2, 3 = basin of the cycle through 1, 5, 17
cyc_min = {1: 1, 2: 1, 5: 2, 14: 2, 7: 2, 20: 2, 10: 2}
c17 = [17]
y = 17
while True:
    y = y // 2 if y % 2 == 0 else 3 * y - 1
    if y == 17:
        break
    c17.append(y)
for v in c17:
    cyc_min[v] = 3
for v, lab in cyc_min.items():
    if v <= N:
        basin[v] = lab
for n0 in range(1, N + 1):
    if basin[n0]:
        continue
    path = []
    n = n0
    while True:
        if n <= N and basin[n]:
            lab = basin[n]
            break
        if n in cyc_min:
            lab = cyc_min[n]
            break
        path.append(n)
        n = n // 2 if n % 2 == 0 else 3 * n - 1
    for v in path:
        if v <= N:
            basin[v] = lab
dens = [basin.count(L) / N for L in (1, 2, 3)]
print(f"  basin counting densities on [1, {N}]: {{1}}: {dens[0]:.5f}, {{5,..}}: {dens[1]:.5f}, {{17,..}}: {dens[2]:.5f} (HYP-9165: 0.327, 0.325, 0.348)")


def delta(n, k, lab=1):
    step = 2 * 3 ** k
    m = n % step
    if m == 0:
        m = step
    while m <= N:
        if basin[m] == lab:
            return m
        m += step
    return None


rows = []
for n in (5, 17, 7, 25, 41):
    row = [delta(n, k) for k in range(1, 12)]
    rows.append((n, row))
    print(f"  n = {n:3d} (basin {basin[n]}): delta_k for k = 1..11:", row)
# None means no basin-{1} element of the class below N, i.e. delta_k(n) > N
grow = all((r[1][-1] if r[1][-1] is not None else N + 1) > 3 ** 9 for r in rows)
check(grow, f"C: for vertices of the two non-trivial components, delta_11(n) exceeds 3^9 (None = beyond N = {N})")
ratio = [round(r[1][k - 1] / 3 ** k, 2) for r in rows[:2] for k in (8, 9, 10, 11) if r[1][k - 1]]
print("  delta_k(n)/3^k for n = 5, 17 at k = 8..11:", ratio, "(bounded ratio: delta_k grows like 3^k)")

# ---------------- D. the backward sieve on a minimal counterexample ----------------
print("\n== D. backward sieve: odd n with no predecessor of multiplier 2^e/3^j < 1 within depth D ==")
L2, L3 = math.log(2), math.log(3)
P3 = [3 ** i for i in range(80)]
INF = float('inf')


@lru_cache(maxsize=None)
def Me(x, d):
    """min log-multiplier over non-root vertices of the depth-d truncation of e(x); x mod 3^ceil(d/2)"""
    if d == 0:
        return INF
    best = L2 + min(0.0, Me((2 * x) % P3[d // 2], d - 1))
    if x % 3 == 1:
        best = min(best, -L3 + min(0.0, Mo(((x - 1) // 3) % P3[(d - 1) // 2], d - 1)))
    return best


@lru_cache(maxsize=None)
def Mo(x, d):
    if d == 0:
        return INF
    return L2 + min(0.0, Me((2 * x) % P3[d // 2], d - 1))


print("   D   classes (odd, 3 not dividing, mod 3^floor(D/2))   surviving fraction")
surv = {}
for D in range(4, 23, 2):
    m = D // 2
    tot = 0
    good = 0
    for x in range(P3[m]):
        if x % 3 == 0:
            continue
        tot += 1
        if Mo(x, D) >= 0:
            good += 1
    surv[D] = good / tot
    print(f"  {D:2d}   {tot:8d}   {good / tot:.5f}")
check(Mo(2, 6) < 0 and all(Mo(x, 10) < 0 for x in range(4, P3[5], 9)),
      "D: odd n = 2 mod 3 and odd n = 4 mod 9 always have a smaller predecessor ((2n-1)/3, (8n-5)/9)")
print("  multiples of 3 have the bare backward tree, so the backward sieve never excludes them")

# ---------------- E. the two-sided trap on the minus sheet ----------------
print("\n== E. the minima of the minus-sheet components ==")
mins = {}
for v in range(1, N + 1):
    lab = basin[v]
    if lab not in mins:
        mins[lab] = v
print("  least element of each basin in [1, N]:", mins)
check(mins == {1: 1, 2: 5, 3: 17}, "E: 5 and 17 are the minima of their components: the two-sided no-descent trap is realised by integers on the minus sheet")

print("\n" + ("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED: {FAILS}"))
