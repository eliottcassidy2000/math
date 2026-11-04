#!/usr/bin/env python3
"""
primegame_deletion_tower_analysis_20261104.py  (opus-2026-11-04-S72)

Analysis for the note 05-knowledge/results/primegame_deletion_tower_path_weighted_decks_20261104.md.

Input: the node files tower_n.txt (n <= 8) written by
    deletion_tower_census DATA_DIR-run:  deletion_tower_census 8 DATA_DIR
Usage:
    python3 primegame_deletion_tower_analysis_20261104.py DATA_DIR

Sections:
  A  the order-4 and order-5 face rows (THM-796 / session129 table) and Lemma W worked on S_4
  B  line classes: orbits of <kappa, sigma> on fixed-path tilings, Burnside 2^(m-2) + 2^(h-2), checked n <= 7
  C  Fibonacci-valued rows and the source-extension sum rule (Lemma S), checked on every node n <= 8
  D  PRIMEGAME's ten primes as the coordinates of order 6: integer codes of the 34 order-6 nodes, and a calibrated
     test of whether PRIMEGAME's fractions act as moves between order-6 node codes
  E  axis-profile nodes
"""
import sys, random, itertools
from fractions import Fraction
from math import gcd
from functools import reduce

D = sys.argv[1] if len(sys.argv) > 1 else "."
NMAX = 8


def load(n):
    rows = []
    for line in open(f"{D}/tower_{n}.txt"):
        if line.startswith('#'):
            continue
        a, W, R, code = line.split('|')
        f = a.split()
        rows.append(dict(id=int(f[0]), sc=int(f[1]), aut=int(f[2]), H=int(f[3]), til=int(f[4]),
                         score=tuple(int(x) for x in f[5].split(',')),
                         W={int(k): int(v) for k, v in (t.split(':') for t in W.split())},
                         R={int(k): int(v) for k, v in (t.split(':') for t in R.split())},
                         code=int(code)))
    return rows


T = {n: load(n) for n in range(2, NMAX + 1)}


def decode(n, code):
    """adjacency out[i] bitmask from the census code (bits ordered by q = 1..n-1, p = 0..q-1, MSB first)"""
    tot = n * (n - 1) // 2
    b = tot - 1
    out = [0] * n
    for q in range(1, n):
        for p in range(q):
            bit = (code >> b) & 1
            b -= 1
            if bit:
                out[p] |= 1 << q
            else:
                out[q] |= 1 << p
    return out


def ham_counts(n, out):
    full = (1 << n) - 1
    def ends(g):
        dp = [[0] * n for _ in range(full + 1)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1, full + 1):
            for v in range(n):
                c = dp[S][v]
                if not c:
                    continue
                for w in range(n):
                    if (g[v] >> w) & 1 and not (S >> w) & 1:
                        dp[S | 1 << w][w] += c
        return [dp[full][v] for v in range(n)]
    op = [0] * n
    for i in range(n):
        for j in range(n):
            if (out[i] >> j) & 1:
                op[j] |= 1 << i
    return ends(op), ends(out)   # (h_s, h_e)


print("=" * 100)
print("A. Face rows at orders 4 and 5 (raw R over the lower nodes, in census order) and Lemma W on S_4")
for n in (4, 5):
    low = T[n - 1]
    print(f"  order {n}: lower nodes " + ", ".join(f"[{r['id']}] H={r['H']}" for r in low))
    for r in T[n]:
        print(f"    node {r['id']}: H={r['H']:>3} |Aut|={r['aut']} {'self-converse' if r['sc'] else 'converse pair'} "
              f"tilings={r['til']:>2} score={r['score']} R={[r['R'].get(k, 0) for k in range(len(low))]} W={[r['W'].get(k, 0) for k in range(len(low))]}")
S4 = [r for r in T[4] if r['H'] == 5][0]
out = decode(4, S4['code'])
hs, he = ham_counts(4, out)
print(f"  S_4 (code {S4['code']}): h_s = {hs}, h_e = {he}, h_s + h_e = {[a + b for a, b in zip(hs, he)]};"
      f" W = {S4['W']} = 2|Aut| R = 2 x {S4['R']}")

print("=" * 100)
print("B. Line classes: orbits of <kappa, sigma> on fixed-path tilings")


def line_classes(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    m = len(pairs)
    pos = {p: k for k, p in enumerate(pairs)}
    refl = [pos[(n - 1 - j, n - 1 - i)] for (i, j) in pairs]
    full = (1 << m) - 1
    seen = set()
    orbits = 0
    for t in range(1 << m):
        if t in seen:
            continue
        s = 0
        for k in range(m):
            if (t >> k) & 1:
                s |= 1 << refl[k]
        orb = {t, t ^ full, s, s ^ full}
        seen |= orb
        orbits += 1
    h = (n - 1) ** 2 // 4
    return orbits, 2 ** (m - 2) + 2 ** (h - 2) if m >= 2 else None, m, h


for n in range(3, 8):
    o, f, m, h = line_classes(n)
    formula = Fraction(2 ** m + 2 ** h, 4)
    print(f"  n={n}: m=binom(n-1,2)={m}, h=floor((n-1)^2/4)={h}: line classes = {o}; Burnside (2^m + 2^h)/4 = {formula}")
for n in range(8, 11):
    m = (n - 1) * (n - 2) // 2
    h = (n - 1) ** 2 // 4
    print(f"  n={n}: Burnside (2^m + 2^h)/4 = {(2 ** m + 2 ** h) // 4}")

print("=" * 100)
print("C. Fibonacci-valued rows, and the source-extension sum rule")
fib = [1, 2]
while fib[-1] < 10 ** 7:
    fib.append(fib[-1] + fib[-2])
fibidx = {v: i for i, v in enumerate(fib)}
for n in range(4, NMAX + 1):
    hits = []
    for r in T[n]:
        vals = sorted(r['R'].values())
        if len(vals) >= 2 and len(set(vals)) == len(vals) and all(v in fibidx for v in vals):
            idx = [fibidx[v] for v in vals]
            if all(idx[i + 1] - idx[i] == 1 for i in range(len(idx) - 1)):
                hits.append(r)
    nmulti = sum(1 for r in T[n] if len(r['R']) >= 2)
    print(f"  order {n}: {len(hits)} of {nmulti} multi-entry rows are runs of consecutive Fibonacci numbers")
    for r in hits:
        print(f"     node {r['id']} H={r['H']} |Aut|={r['aut']} {'sc' if r['sc'] else 'pair'} tilings={r['til']} R={dict(sorted(r['R'].items()))}")

# Lemma S check: for T' of order n-1 and T = s -> T' (source added), W(T) = H(T') e_[T'] + sum_x h_e^{T'}(x) e_[s -> T'-x]
# verified on every merged node of order n-1 <= 7 by recomputing W(T) from scratch in this script.
canon_cache = {}


def canon(n, out):
    key = (n, tuple(out))
    if key in canon_cache:
        return canon_cache[key]
    best = None
    for perm in itertools.permutations(range(n)):
        code = 0
        for q in range(1, n):
            for p in range(q):
                code = (code << 1) | ((out[perm[p]] >> perm[q]) & 1)
        if best is None or code < best:
            best = code
    canon_cache[key] = best
    return best


def merged_id_table(n):
    tab = {}
    for r in T[n]:
        out = decode(n, r['code'])
        c = canon(n, out)
        op = [0] * n
        for i in range(n):
            for j in range(n):
                if (out[i] >> j) & 1:
                    op[j] |= 1 << i
        tab[c] = r['id']
        tab[canon(n, op)] = r['id']
    return tab


def delete(n, out, x):
    res = []
    for i in range(n):
        if i == x:
            continue
        r = 0
        kk = 0
        for j in range(n):
            if j == x:
                continue
            if (out[i] >> j) & 1:
                r |= 1 << kk
            kk += 1
        res.append(r)
    return res


def conv_out(n, out):
    op = [0] * n
    for i in range(n):
        for j in range(n):
            if (out[i] >> j) & 1:
                op[j] |= 1 << i
    return op


def join(a, A, b, B):
    """ordinal sum A => B: vertices 0..a-1 carry A, a..a+b-1 carry B, every A-vertex beats every B-vertex"""
    n = a + b
    out = [0] * n
    for i in range(a):
        out[i] = A[i] | (((1 << b) - 1) << a)
    for j in range(b):
        out[a + j] = B[j] << a
    return out


def reps_with_converses(k):
    if k == 1:
        return [[0]]
    res = []
    for r in T[k]:
        o = decode(k, r['code'])
        res.append(o)
        if not r['sc']:
            res.append(conv_out(k, o))
    return res


# Lemma J (join rule): for T = A => B every Hamiltonian path runs through A and then B, so
#   W(T) = H(B) sum_{x in A} h_s^A(x) e_[(A - x) => B] + H(A) sum_{y in B} h_e^B(y) e_[A => (B - y)];
# Lemma S (source extension) is the case |A| = 1.  Checked on every join of ordinary classes with |A| + |B| <= 6.
checked = 0
fails = 0
tables = {k: merged_id_table(k) for k in range(2, 6)}
for n in range(3, 7):
    tab_low = tables[n - 1]
    for a in range(1, n):
        b = n - a
        for A in reps_with_converses(a):
            for B in reps_with_converses(b):
                t = join(a, A, b, B)
                hs, he = ham_counts(n, t)
                Wd = {}
                for x in range(n):
                    v = tab_low[canon(n - 1, delete(n, t, x))]
                    Wd[v] = Wd.get(v, 0) + hs[x] + he[x]
                Wd = {k: v for k, v in Wd.items() if v}
                hsA, heA = ham_counts(a, A) if a > 1 else ([1], [1])
                hsB, heB = ham_counts(b, B) if b > 1 else ([1], [1])
                HA, HB = sum(heA), sum(heB)
                Wp = {}
                for x in range(a):
                    sub = join(a - 1, delete(a, A, x), b, B) if a > 1 else B
                    v = tab_low[canon(n - 1, sub)]
                    Wp[v] = Wp.get(v, 0) + HB * hsA[x]
                for y in range(b):
                    sub = join(a, A, b - 1, delete(b, B, y)) if b > 1 else A
                    v = tab_low[canon(n - 1, sub)]
                    Wp[v] = Wp.get(v, 0) + HA * heB[y]
                Wp = {k: v for k, v in Wp.items() if v}
                checked += 1
                if Wp != Wd:
                    fails += 1
print(f"  Lemma J (join rule; Lemma S = the case |A| = 1) checked on all {checked} joins A => B of ordinary classes with |A| + |B| <= 6: {fails} failures")

# where the order-5 Fibonacci rows go under a source extension (rows at order 6, over the order-5 nodes)
tab6 = merged_id_table(6)
byid6 = {r['id']: r for r in T[6]}
for r5 in T[5]:
    if sorted(r5['R'].values()) in ([2, 3], [5, 8]):
        o = decode(5, r5['code'])
        for label, t in (("s -> T'", join(1, [0], 5, o)), ("T' -> s", join(5, o, 1, [0]))):
            u = tab6[canon(6, t)]
            print(f"  order-5 node {r5['id']} (H={r5['H']}, row {dict(sorted(r5['R'].items()))}): {label} is order-6 node {u}, row {dict(sorted(byid6[u]['R'].items()))}")

print("=" * 100)
print("D. PRIMEGAME's ten primes as the coordinates of order 6")
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
fracs = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11), (15, 2), (1, 7), (55, 1)]


def val(x):
    v = [0] * 10
    for i, p in enumerate(primes):
        while x % p == 0:
            x //= p
            v[i] += 1
    assert x == 1
    return v


fvec = [tuple(a - b for a, b in zip(val(p), val(q))) for p, q in fracs]
low5 = T[5]
order5 = sorted(low5, key=lambda r: (r['H'], r['til'], sorted(r['R'].items())))
assign = {r['id']: primes[i] for i, r in enumerate(order5)}
print("  canonical assignment (order-5 nodes sorted by H, tilings, row):")
for r in order5:
    print(f"     prime {assign[r['id']]:>2} <- order-5 node {r['id']} (H={r['H']}, tilings={r['til']}, R={dict(sorted(r['R'].items()))})")
codes6 = []
for r in T[6]:
    N = 1
    for k, v in r['R'].items():
        N *= assign[k] ** v
    g = reduce(gcd, r['R'].values())
    Np = 1
    for k, v in r['R'].items():
        Np *= assign[k] ** (v // g)
    codes6.append((r, N, Np))
print(f"  the 34 order-6 nodes as integers over PRIMEGAME's primes (raw row exponents; primitive in brackets):")
for r, N, Np in codes6:
    print(f"     node {r['id']:>2} H={r['H']:>2} {'sc' if r['sc'] else 'pr'}: {N}" + ("" if N == Np else f"  [{Np}]"))
print(f"  distinct raw codes: {len(set(N for _, N, _ in codes6))} / 34; distinct primitive codes: {len(set(Np for _, _, Np in codes6))} / 34")


def move_count(assignment):
    vecs = []
    for r in T[6]:
        v = [0] * 10
        for k, e in r['R'].items():
            v[primes.index(assignment[k])] += e
        vecs.append(tuple(v))
    fset = {}
    for i, f in enumerate(fvec):
        fset.setdefault(f, []).append(i)
    hits = 0
    for a in range(len(vecs)):
        for b in range(len(vecs)):
            if a == b:
                continue
            d = tuple(x - y for x, y in zip(vecs[b], vecs[a]))
            if d in fset:
                hits += len(fset[d])
    return hits


h0 = move_count(assign)
random.seed(20261104)
ids = [r['id'] for r in low5]
samples = []
for _ in range(20000):
    perm = primes[:]
    random.shuffle(perm)
    samples.append(move_count({k: perm[i] for i, k in enumerate(ids)}))
samples.sort()
ge = sum(1 for s in samples if s >= h0)
print(f"  PRIMEGAME fractions as moves between order-6 node codes (N(u') = f N(u)): {h0} (pair, fraction) incidences under the canonical assignment")
print(f"  random assignments (20,000): min {samples[0]}, median {samples[len(samples) // 2]}, mean {sum(samples) / len(samples):.3f}, max {samples[-1]}; P(random >= canonical) = {ge / len(samples):.4f}")

print("=" * 100)
print("E. Axis-profile nodes (face row supported on one lower node)")
for n in range(3, NMAX + 1):
    ax = [r for r in T[n] if len(r['R']) == 1]
    print(f"  order {n}: {len(ax)} axis nodes; H values {sorted(r['H'] for r in ax)}")
