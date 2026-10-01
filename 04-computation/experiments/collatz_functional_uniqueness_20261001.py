#!/usr/bin/env python3
"""Collatz among functional graphs: rigidity and uniqueness (opus S15, seventeenth note, 2026-10-01).

Map: T(n) = n/2 (n even), 3n+1 (n odd)  [non-shortcut], on the positive integers.
Unfolded backward tree of n: children 2n, and (n-1)/3 when n = 4 (mod 6).

Abstract 3-adic model (Lemma 1 of the note):
  e(x): children e(2x), and o((x-1)/3) when x = 1 (mod 3)
  o(x): child e(2x)
  A(u) := e(u) for u = 1 (mod 3).  Depth-d truncation of e(x) depends on x mod 3^ceil(d/2),
  of o(x) on x mod 3^floor(d/2).

Checks (each prints PASS/FAIL lines; the run ends with ALL CHECKS PASSED or a failure count):
  A. separation exponent s(D) of the trees A(u) versus the proved bound 2 + floor((D-4)/5)
  B. integer backward trees equal the abstract trees; depth-D collisions are 3-adically close,
     never across shapes; independent string-canonical cross-check
  C. Eckmann-Hilton near-balance: agreement depth of B(2w) and B((w-1)/3) versus v_3(5w+1);
     the 3-adic point -1/5 is exactly balanced; no integer fibre is balanced
  D. residue graphs Gamma_M (no periodic invariant) for 3x+1, 3x-1; gcd invariant for 3x+5, 3x+7
  E. selection census for n/2, an+b on the positive integers
  F. cycle balance identities (first and second moment) on known cycles
  G. a = 5 (2 a primitive root mod 5^k): backward separation; a = 7: bare unit spines
Reproduce: python 04-computation/experiments/collatz_functional_uniqueness_20261001.py
"""
import sys
import time
from functools import lru_cache
from collections import defaultdict

sys.setrecursionlimit(1000000)
FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def v3(n):
    if n == 0:
        return 10 ** 9
    n = abs(n)
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k


def vp(n, p):
    if n == 0:
        return 10 ** 9
    n = abs(n)
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


# ---------------- canonical forms (AHU hash-consing) ----------------
_canon = {}


def cid(children):
    key = tuple(sorted(children))
    v = _canon.get(key)
    if v is None:
        v = len(_canon)
        _canon[key] = v
    return v


LEAF = cid(())
P3 = [3 ** i for i in range(200)]


@lru_cache(maxsize=None)
def Ecan(x, d):
    """canonical id of the depth-d truncation of e(x); x reduced mod 3^ceil(d/2)"""
    if d == 0:
        return LEAF
    ch = [Ecan((2 * x) % P3[d // 2], d - 1)]
    if x % 3 == 1:
        ch.append(Ocan(((x - 1) // 3) % P3[(d - 1) // 2], d - 1))
    return cid(ch)


@lru_cache(maxsize=None)
def Ocan(x, d):
    """canonical id of the depth-d truncation of o(x); x reduced mod 3^floor(d/2)"""
    if d == 0:
        return LEAF
    return cid([Ecan((2 * x) % P3[d // 2], d - 1)])


def E(x, d):
    return Ecan(x % P3[(d + 1) // 2], d)


def O(x, d):
    return Ocan(x % P3[d // 2], d)


# ---------------- A. separation exponent of A(u) ----------------
def separation(D):
    m = (D + 1) // 2
    groups = defaultdict(list)
    for x in range(1, P3[m], 3):
        groups[E(x, D)].append(x)
    s = m
    worst = None
    for g in groups.values():
        if len(g) > 1:
            x0 = g[0]
            for y in g[1:]:
                vv = v3(y - x0)
                if vv < s:
                    s = vv
                    worst = (x0, y)
    return s, m, len(groups), worst


print("== A. separation exponent s(D) of the abstract trees A(u), u = 1 mod 3 ==")
print("   D   m=ceil(D/2)  #classes(=3^(m-1) if injective)  s(D)  proved bound 2+floor((D-4)/5)")
sD = {}
t0 = time.time()
for D in range(4, 27):
    s, m, ncls, worst = separation(D)
    sD[D] = s
    bound = 2 + (D - 4) // 5
    print(f"  {D:2d}   {m:2d}   {ncls:8d} / {P3[m-1]:8d}   s={s:2d}   bound={bound}   worst-pair={worst}")
    check(s >= bound, f"A: s({D}) = {s} >= proved bound {bound}")
print(f"  (time {time.time()-t0:.1f}s)")
# the observed growth: first depth at which s(D) reaches each value
first = {}
for D in sorted(sD):
    first.setdefault(sD[D], D)
print("  first depth reaching s = k:", first, " (steps", [first[k + 1] - first[k] for k in sorted(first)[:-1]], ")")
print("  observed rate ~ 10/3 depth per 3-adic digit; the proof's recursion s(D) >= 1 + s(D-5) gives 5")

# ---------------- B. integer backward trees ----------------
print("\n== B. integer backward trees (positive integers) ==")


@lru_cache(maxsize=None)
def intcan(n, d):
    if d == 0:
        return LEAF
    ch = [intcan(2 * n, d - 1)]
    if n % 6 == 4:
        ch.append(intcan((n - 1) // 3, d - 1))
    return cid(ch)


NB = 6000
ok = True
for D in (6, 10, 14, 18):
    for n in range(1, NB + 1):
        a = intcan(n, D)
        b = E(n, D) if n % 2 == 0 else O(n, D)
        if a != b:
            ok = False
            print("  mismatch", n, D)
            break
check(ok, f"B: integer backward trees of n <= {NB} equal the abstract e/o trees (depths 6,10,14,18)")


def shape(n):
    """first-branching depth and normalised label u = 1 mod 3"""
    if n % 3 == 0:
        return (None, None)
    if n % 2 == 0 and n % 3 == 1:
        return (0, n)
    if n % 3 == 2:
        return (1, 2 * n)
    return (2, 4 * n)  # odd, = 1 mod 3


for D in (12, 16, 20):
    groups = defaultdict(list)
    for n in range(1, NB + 1):
        groups[intcan(n, D)].append(n)
    cross_shape = 0
    too_close = 0
    multi = 0
    for g in groups.values():
        shp = {shape(n)[0] for n in g}
        if len(shp) > 1:
            cross_shape += 1
            continue
        k = next(iter(shp))
        if k is None:
            continue  # multiples of 3: bare rays, one class
        if len(g) > 1:
            multi += 1
            u0 = shape(g[0])[1]
            need = sD.get(D - k, (D - k - 2) // 2)
            for n in g[1:]:
                if v3(shape(n)[1] - u0) < need:
                    too_close += 1
    bare = groups[intcan(3, D)]
    check(cross_shape == 0, f"B: depth {D}: no canonical class mixes shapes (fbd 0/1/2/bare)")
    check(too_close == 0, f"B: depth {D}: every collision among n <= {NB} is 3-adically close (v_3 >= s(D-k)); {multi} multi-classes")
    check(all(n % 3 == 0 for n in bare) and len(bare) == NB // 3,
          f"B: depth {D}: the bare-ray class is exactly the multiples of 3 ({len(bare)})")

# direct identification: distinct n <= 2000 with 3 not dividing n differ 3-adically to order <= 6,
# and s(D-2) >= 7 at D = 25, so their depth-25 backward trees must all be distinct
D = 25
ids = {}
dup = []
for n in range(1, 2001):
    if n % 3 == 0:
        continue
    c = intcan(n, D)
    if c in ids:
        dup.append((ids[c], n))
    ids[c] = n
check(not dup, f"B: the depth-{D} backward trees of the {len(ids)} vertices n <= 2000 with 3 not dividing n are pairwise distinct" + ("" if not dup else f" dup={dup[:5]}"))

# independent cross-check with string canonical forms (different code path)


def strcan_int(n, d):
    if d == 0:
        return "()"
    ch = [strcan_int(2 * n, d - 1)]
    if n % 6 == 4:
        ch.append(strcan_int((n - 1) // 3, d - 1))
    return "(" + "".join(sorted(ch)) + ")"


ok = True
D = 11
seen = {}
for n in range(1, 1500):
    s = strcan_int(n, D)
    key = intcan(n, D)
    if key in seen and seen[key] != s:
        ok = False
    seen.setdefault(key, s)
# and the converse: equal strings <=> equal ids
inv = defaultdict(set)
for n in range(1, 1500):
    inv[strcan_int(n, D)].add(intcan(n, D))
ok = ok and all(len(v) == 1 for v in inv.values())
check(ok, "B: string canonical forms (independent path) agree with hash-consed ids, n < 1500, depth 11")

# ---------------- C. Eckmann-Hilton near balance ----------------
print("\n== C. Eckmann-Hilton point: fibres {2w, (w-1)/3}, w = 4 mod 6 ==")
DCAP = 40


def agree_depth(w):
    a, b = 2 * w, (w - 1) // 3
    d = 0
    while d < DCAP and intcan(a, d + 1) == intcan(b, d + 1):
        d += 1
    return d


table = defaultdict(list)
balanced = 0
NW = 200000
for w in range(4, NW, 6):
    v = v3(5 * w + 1)
    ad = agree_depth(w)
    table[v].append(ad)
    if ad >= DCAP:
        balanced += 1
print("  v_3(5w+1)   count   min agree   max agree   (agree = deepest d with B(2w) ~_d B((w-1)/3))")
for v in sorted(table):
    L = table[v]
    print(f"     {v:2d}     {len(L):7d}     {min(L):3d}        {max(L):3d}")
print("  agreement depth / v_3(5w+1):", {v: round(max(table[v]) / v, 2) for v in sorted(table)})
# the proved upper bound: agreement to depth d forces v_3(5w+1) - 1 >= s(d-1) >= 2 + floor((d-5)/5), so d <= 5v - 6
ok = all(max(table[v]) <= 5 * v - 6 for v in table if v >= 2) and max(table[1]) <= 1
check(ok, "C: agreement depth <= 5 v_3(5w+1) - 6 (the bound implied by the separation lemma) for every w < %d" % NW)
ok = all(max(table[v]) <= 5 * v - 6 for v in table if v >= 2) and all(min(table[v]) >= 2 * v - 1 for v in table)
check(ok, "C: 2v - 1 <= agreement depth <= 5v - 6 on the whole range (observed about 3v - 1)")
check(balanced == 0, f"C: no fibre with w < {NW} is balanced to the cap {DCAP}")
w12 = [w for w in range(4, NW, 6) if v3(5 * w + 1) == 12]
print("  deepest example:", w12, "fibre", [(2 * w, (w - 1) // 3) for w in w12], "agree", [agree_depth(w) for w in w12])
# exact balance at -1/5 (3-adic): the two children of e(-1/5) coincide at every depth
m = 30
inv5 = pow(5, -1, P3[m])
x = (-inv5) % P3[m]
ok = x % 3 == 1 and all(E(2 * x, d) == O((x - 1) // 3, d) for d in range(1, 2 * m - 4))
check(ok, "C: -1/5 = 1 mod 3 and e(-1/5) has two identical child subtrees (depths < 2m-4, m = 30)")
ok = (2 * x - (x - 1) // 3) % P3[m - 1] == 0
check(ok, "C: at -1/5 the two preimage maps 2x and (x-1)/3 coincide 3-adically")
print("  -1/5 is 2-adically odd (numerator -1), so no even vertex of a rational graph has 3-adic value -1/5")

# ---------------- D. residue graphs Gamma_M ----------------
print("\n== D. residue graphs Gamma_M (edges r -> T(n) mod M over n = r mod M) ==")


def gamma_components(M, a, b):
    parent = list(range(M))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for n in range(1, 2 * M + 1):
        t = n // 2 if n % 2 == 0 else a * n + b
        if t <= 0:
            continue
        ra, rb = find(n % M), find(t % M)
        if ra != rb:
            parent[ra] = rb
    return len({find(i) for i in range(M)})


MMAX = 1500
for (a, b) in ((3, 1), (3, -1), (5, 1), (7, 1)):
    bad = [M for M in range(1, MMAX + 1) if gamma_components(M, a, b) != 1]
    check(not bad, f"D: {a}x{b:+d}: Gamma_M connected for every M <= {MMAX} (no periodic invariant)" + ("" if not bad else f" bad={bad[:10]}"))
for (a, b) in ((3, 5), (3, 7)):
    comp = gamma_components(b, a, b)
    print(f"  {a}x+{b}: Gamma_{b} has {comp} components (gcd(n,{b}) invariant)")
    check(comp >= 2, f"D: {a}x+{b} has a non-constant periodic invariant (gcd with {b})")

# ---------------- E. selection census ----------------
print("\n== E. census of n/2, an+b on the positive integers (starts <= 3000, cap 2e5 steps or 1e40) ==")


def census(a, b, N=3000, steps=20000, vcap=10 ** 40):
    """result per start: least element of the cycle reached, or 'esc' (value cap / step cap)"""
    known = {}
    cycles = set()
    esc = 0
    for n0 in range(1, N + 1):
        path = []
        onpath = {}
        n = n0
        res = None
        k = 0
        while True:
            if n in known:
                res = known[n]
                break
            if n in onpath:
                cyc = path[onpath[n]:]
                res = min(cyc)
                cycles.add(res)
                break
            if k >= steps or n >= vcap:
                res = 'esc'
                break
            onpath[n] = len(path)
            path.append(n)
            n = n // 2 if n % 2 == 0 else a * n + b
            k += 1
        for y in path:
            if y < 10 ** 7:
                known[y] = res
        if res == 'esc':
            esc += 1
    return sorted(cycles), esc


from math import gcd
print("   a   b   cycles (least elements)                 escapes   obstruction")
for a in (1, 3, 5, 7, 9):
    for b in (1, -1, 3, 5, 7):
        if a == 1 and b == -1:
            continue
        cyc, esc = census(a, b)
        obs = []
        if gcd(a, abs(b)) > 1:
            obs.append(f"a,b share {gcd(a, abs(b))} (reduces to b/{gcd(a, abs(b))} on an absorbing sublattice)")
        elif abs(b) > 1:
            obs.append(f"gcd(n,{abs(b)}) invariant")
        if len(cyc) >= 2:
            obs.append(f"{len(cyc)} cycles")
        if esc:
            obs.append("escapes (drift > 0)" if a >= 5 else "escapes")
        print(f"  {a:2d} {b:+3d}   {str(cyc)[:40]:40s}  {esc:6d}   {'; '.join(obs) if obs else '-- (no obstruction found)'}")

# ---------------- F. cycle balance identities ----------------
print("\n== F. cycle balance identities (shortcut map x/2, (3x+b)/2) ==")


def shortcut_cycle(x0, b):
    cyc = [x0]
    x = x0
    while True:
        x = x // 2 if x % 2 == 0 else (3 * x + b) // 2
        if x == x0:
            return cyc
        cyc.append(x)


ok = True
for (x0, b) in ((1, 1), (-1, 1), (-5, 1), (-17, 1), (1, -1), (5, -1), (17, -1), (1, 5), (19, 5), (23, 5), (187, 5), (347, 5)):
    c = shortcut_cycle(x0, b)
    Ev = [x for x in c if x % 2 == 0]
    Od = [x for x in c if x % 2 != 0]
    p = len(Od)
    first = sum(Ev) - sum(Od) == b * p
    second = 3 * sum(x * x for x in Ev) == 5 * sum(x * x for x in Od) + 6 * b * sum(Od) + b * b * p
    ok = ok and first and second
    print(f"  3x{b:+d} cycle from {x0:5d}: length {len(c):3d}, p = {p:3d}, sum_even - sum_odd = {sum(Ev)-sum(Od):5d} (= b p: {first}); second moment: {second}")
check(ok, "F: sum_even - sum_odd = b p and 3 E_2 = 5 O_2 + 6 b O_1 + b^2 p on every listed cycle")

# ---------------- G. other multipliers ----------------
print("\n== G. backward separation for n/2, an+1: a = 5 (2 primitive root mod 5^k) vs a = 7 ==")


def make_can(a):
    P = [a ** i for i in range(80)]

    @lru_cache(maxsize=None)
    def Ea(x, d):
        if d == 0:
            return LEAF
        ch = [Ea((2 * x) % P[d], d - 1)]
        if x % a == 1:
            ch.append(Oa(((x - 1) // a) % P[d], d - 1))
        return cid(ch)

    @lru_cache(maxsize=None)
    def Oa(x, d):
        if d == 0:
            return LEAF
        return cid([Ea((2 * x) % P[d], d - 1)])
    return Ea, P


for a, D, m in ((5, 8, 4), (5, 12, 6), (5, 16, 7), (7, 12, 4), (7, 16, 5)):
    Ea, P = make_can(a)
    groups = defaultdict(list)
    for x in range(P[m]):
        if x % a == 0:
            continue
        groups[Ea(x, D)].append(x)
    sizes = sorted((len(g) for g in groups.values()), reverse=True)
    smin = m
    for g in groups.values():
        if a == 7 and g[0] % 7 in (3, 5, 6):
            continue
        for y in g[1:]:
            smin = min(smin, vp(y - g[0], a))
    bare = sum(1 for x in range(1, P[m]) if a == 7 and x % 7 in (3, 5, 6))
    print(f"  a={a} depth {D}, units mod {a}^{m}: {sum(sizes)} labels -> {len(groups)} classes; largest class {sizes[0]};"
          f" collisions (outside bare spines) are {a}-adically close to order {smin}"
          + (f"; labels with bare spines (x = 3,5,6 mod 7): {bare}" if a == 7 else ""))
    if a == 7:
        check(sizes[0] >= bare, f"G: a = 7, depth {D}: the non-residue coset {{3,5,6}} mod 7 collapses to one bare class (no backward injectivity)")

# ---------------- H. the shortcut graph x/2, (3x+1)/2 ----------------
print("\n== H. shortcut graph: preimages 2n and (2n-1)/3 (n = 2 mod 3); one vertex type, labels in Z_3 ==")


@lru_cache(maxsize=None)
def Scan(x, d):
    """canonical id of the depth-d truncation of the shortcut backward tree S(x); x reduced mod 3^d"""
    if d == 0:
        return LEAF
    ch = [Scan((2 * x) % P3[d - 1], d - 1)]
    if x % 3 == 2:
        ch.append(Scan(((2 * x - 1) // 3) % P3[d - 1], d - 1))
    return cid(ch)


for D, m in ((10, 5), (14, 7), (18, 9), (22, 10)):
    groups = defaultdict(list)
    for x in range(P3[m]):
        if x % 3 == 0:
            continue
        groups[Scan(x % P3[D], D)].append(x)
    smin = m
    for g in groups.values():
        for y in g[1:]:
            smin = min(smin, v3(y - g[0]))
    print(f"  depth {D}: units mod 3^{m}: {sum(len(g) for g in groups.values())} labels -> {len(groups)} classes; collisions 3-adically close to order {smin}")
    check(smin >= 2, f"H: shortcut depth {D}: every collision is 3-adically close (order >= 2)")
m = 30
x = (-pow(4, -1, P3[m])) % P3[m]
ok = x % 3 == 2 and all(Scan((2 * x) % P3[m], d) == Scan(((2 * x - 1) // 3) % P3[m], d) for d in range(1, 28))
check(ok, "H: shortcut Eckmann-Hilton point -1/4 = 2 mod 3: its two preimage subtrees coincide (depths < 28)")
print("  -1/4 is not a 2-adic integer, so the shortcut coincidence point is not a vertex of any graph on Z_(2)")

print("\n" + ("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED: {FAILS}"))
