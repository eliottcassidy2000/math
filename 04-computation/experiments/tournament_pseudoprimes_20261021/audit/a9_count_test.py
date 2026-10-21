# AUDIT a9 (independent auditor, 2026-10-21): THM-4746 statement 6 (tournament-count test) and its table.
from sympy import isprime, n_order, factorint
direct = {}
for line in open("a2_direct_4_2000.out"):
    n, ok, sv, cv = line.split(); direct[int(n)] = int(ok)
def count_test(n): return pow(2, n * (n - 1) // 2, n) == pow(2, (n - 1) // 2, n)
def psp2(n): return pow(2, n - 1, n) == 1
bad = 0; passes = []
for n in range(9, 20001, 2):
    if isprime(n): continue
    a = count_test(n); b = ((n - 1) ** 2 // 2) % n_order(2, n) == 0
    if a != b: bad += 1; print("MISMATCH", n)
    if a: passes.append(n)
    if psp2(n) and not a: bad += 1; print("psp not in count test", n)
print("odd composite n <= 20000: count test <=> ord_n(2) | (n-1)^2/2 mismatches:", bad)
print("odd composite count-test passes <= 1200:", [n for n in passes if n <= 1200])
print("   of these NOT base-2 pseudoprimes:", [n for n in passes if n <= 1200 and not psp2(n)])
for n in [85, 133, 481, 511, 949, 341, 561, 645, 1105, 21, 1023, 1729]:
    print(f"n={n} {factorint(n)}: count test {count_test(n)}, base-2 psp {psp2(n)}, ord_n(2)={n_order(2, n)}, "
          f"strong (Camion) test {bool(direct.get(n))}")
# universal residues of the count test: odd r, p > r, p | 2^((r-1)^2/2) - 1
for r in [1, 3, 5, 7]:
    v = 2 ** ((r - 1) ** 2 // 2) - 1
    print(f"r={r}: 2^((r-1)^2/2)-1 = {v} = {factorint(v) if v else 0}")
# Burnside reading at prime m (brute force m = 3, 5, 7): rotation classes of all / strong labelled tournaments
import itertools
def strong(n, adj):
    full = (1 << n) - 1
    def reach(a):
        seen = 1; fr = 1
        while fr:
            nx = 0
            for v in range(n):
                if fr >> v & 1: nx |= a[v]
            fr = nx & ~seen; seen |= nx
        return seen == full
    radj = [0] * n
    for v in range(n):
        for w in range(n):
            if adj[v] >> w & 1: radj[w] |= 1 << v
    return reach(adj) and reach(radj)
for m in [3, 5, 7]:
    pairs = list(itertools.combinations(range(m), 2)); idx = {pq: i for i, pq in enumerate(pairs)}
    rot = []
    for (i, j) in pairs:  # bit for pair (i,j) set means i->j ; rotation x->x+1
        a, b = (i + 1) % m, (j + 1) % m
        rot.append((idx[(min(a, b), max(a, b))], a < b))
    seen = set(); orb_all = 0; orb_strong = 0
    N = 1 << len(pairs)
    for bits in range(N):
        if bits in seen: continue
        orbit = []; x = bits
        while x not in orbit:
            orbit.append(x)
            y = 0
            for k, (t, same) in enumerate(rot):
                bit = x >> k & 1
                if not same: bit ^= 1
                y |= bit << t
            x = y
        for o in orbit: seen.add(o)
        orb_all += 1
        adj = [0] * m
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
        orb_strong += strong(m, adj)
    T = 2 ** (m * (m - 1) // 2); c = 2 ** ((m - 1) // 2)
    s = {3: 2, 5: 544, 7: 1677488}[m]
    print(f"m={m}: rotation classes all {orb_all} vs (T+(m-1)c)/m = {(T + (m - 1) * c) // m}; strong {orb_strong} vs (s+(m-1)c)/m = {(s + (m - 1) * c) // m}")
