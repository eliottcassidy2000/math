"""Independent exhaustive census for THM-5040 at n <= 7 (all labelled graphs; numpy for d, audit code for parity),
plus checks of the corollaries and of the Remark on graphs without a universal vertex.

usage: python indep_census.py [NMAX]   (default 7)
"""
import sys, os, json, math, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_lib import *

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
rng = random.Random(7)
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'egcp_20261109')


def census(n):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    m = len(pairs)
    G = np.arange(1 << m, dtype=np.uint32)
    closed = []
    for v in range(n):
        c = np.full(1 << m, 1 << v, dtype=np.uint8)
        for e, (a, b) in enumerate(pairs):
            if a == v:
                c |= (((G >> e) & 1) << b).astype(np.uint8)
            elif b == v:
                c |= (((G >> e) & 1) << a).astype(np.uint8)
        closed.append(c)
    full = (1 << n) - 1
    d = np.zeros(1 << m, dtype=np.uint16)
    for S in range(1, 1 << n):
        dom = np.zeros(1 << m, dtype=np.uint8)
        for v in range(n):
            if (S >> v) & 1:
                dom |= closed[v]
        d += (dom == full)
    univ = np.zeros(1 << m, dtype=bool)
    for v in range(n):
        univ |= (closed[v] == full)
    return pairs, d, univ


def decode(n, pairs, code):
    return graph(n, [pairs[e] for e in range(len(pairs)) if (int(code) >> e) & 1])


def M0(n):
    return 2 ** n - n - 1 if n % 2 == 0 else 2 ** n - n - 4


def even_values(n, pairs, d, mask=None):
    """set of odd k such that some even graph (restricted to mask) on exactly n vertices has d = k;
    for every k NOT realised, every graph with d = k (in mask) was checked and found odd"""
    res = set()
    checked_odd = {}
    for k in range(1, 2 ** n, 2):
        sel = (d == k) if mask is None else ((d == k) & mask)
        idx = np.flatnonzero(sel)
        if len(idx) == 0:
            continue
        order = list(idx)
        rng.shuffle(order)
        found = False
        for code in order:
            if is_even(decode(n, pairs, code)):
                found = True
                break
        if found:
            res.add(k)
        else:
            checked_odd[k] = len(idx)
    return res, checked_odd


print("Exhaustive census of all labelled graphs, n <= %d" % NMAX)
for n in range(1, NMAX + 1):
    pairs, d, univ = census(n)
    # sanity: every d odd (Brouwer-Csorba-Schrijver)
    assert np.all(d % 2 == 1)
    E, oddonly = even_values(n, pairs, d)
    E0, oddonly0 = even_values(n, pairs, d, ~univ)
    print("n=%d: %d labelled graphs, all d odd; E(n) == O(n): %s" % (n, len(d), E == O(n)))
    print("     values realised only by odd graphs (all realisers checked): %s" % (
        {k: c for k, c in sorted(oddonly.items())}))
    if n >= 5:
        claim = set(range(1, M0(n) + 1, 2)) - {3, 7, 13, M0(n) - 2}
        print("     no-universal-vertex values E0(n) == Remark's set: %s  (M0 = %d)" % (E0 == claim, M0(n)))
        if E0 != claim:
            print("       E0 - claim:", sorted(E0 - claim), " claim - E0:", sorted(claim - E0))
    sys.stdout.flush()

# Remark (no universal vertex) against the packet's own census data for n <= 9
js = json.load(open(os.path.join(EXP, 'base_sets.json')))
E0b = {int(k): set(v) for k, v in js['E0'].items()}
for n in range(5, 10):
    claim = set(range(1, M0(n) + 1, 2)) - {3, 7, 13, M0(n) - 2}
    print("base_sets.json E0(%d) == Remark's set: %s" % (n, E0b[n] == claim))

# ---------------- Corollaries 2-3 arithmetic ----------------
print("Corollary 2/3 arithmetic")


def V(k):
    n = 1
    while k > M(n):
        n += 1
    return n + 1 if k == X(n) else n


bad = 0
for k in range(1, 2 ** 18, 2):
    if k in (3, 7, 13):
        continue
    v = V(k)
    N = math.ceil(math.log2(k + 1))
    # V(k) == least n with k in O(n)
    least = next(n for n in range(1, 40) if k in O(n)) if k < 2 ** 12 else v
    bad += not (N <= v <= N + 1) or least != v
print("  ceil(log2(k+1)) <= V(k) <= ceil(log2(k+1)) + 1 and V(k) = least n with k in O(n), odd k < 2^18: failures", bad)
for m in range(2, 21):
    extra = [k for k in range(1, 2 ** (m - 1), 2) if k not in (3, 7, 13) and V(k) == math.ceil(math.log2(k + 1)) + 1
             and math.ceil(math.log2(k + 1)) + 1 == m]
    formula = 1 + (m - 2) // 2
    flag = "" if len(extra) == formula else "   <-- differs"
    print("  m=%2d: values needing ceil(log2(k+1)) + 1 = m vertices: %3d, Corollary 3 formula 1 + floor((m-2)/2) = %2d%s" % (
        m, len(extra), formula, flag))
for m in range(1, 9):
    k = 2 ** m - 1
    print("  2^%d - 1 = %d: attainable %s, V = %s (claimed m + 1 = %d)" % (
        m, k, k not in (3, 7, 13), V(k) if k not in (3, 7, 13) else None, m + 1))
