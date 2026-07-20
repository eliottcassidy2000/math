#!/usr/bin/env python3
"""
h_spectrum_n7_exhaustive_boxeph_S152.py  (HYP-8220)

OWNER S152: exhaust the n=7 trivial-Aut stratum for h = 7 and 21.

METHOD (covers ALL strata, not just trivial-Aut, so the result is total):
every 7-tournament minus a vertex is a 6-tournament, so augmenting canonical
representatives of ALL 56 six-vertex classes by one new vertex in all 2^6
orientation patterns reaches EVERY 7-vertex isomorphism class.  3584 candidate
tournaments -> the COMPLETE n=7 h-spectrum (h is an isomorphism invariant).
This settles h=7 and h=21 at n=7 definitively (S151 already excluded every
|Aut|>1 stratum; this removes the sampling caveat for |Aut|=1 as well).

Also computed: the complete STRONG n=7 h-spectrum (strongness is iso-invariant)
and its minimum -- the WHY behind any gap: gaps at n survive iff the strong
spectra of sizes <= n cannot tile the value multiplicatively (condensation
monoid, S146).  Canonical n=7 representatives (456 classes, A000568 check) are
saved to 05-knowledge/results/n7_class_reps_boxeph_S152.txt for the n=8 stage.

boxeph-2026-07-20-S152.
"""

from itertools import permutations

N6, N7 = 6, 7

def hpaths(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c or not (S >> v) & 1: continue
            m = adj[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[(1 << n) - 1][v] for v in range(n))

def is_strong(adj, n):
    full = (1 << n) - 1
    for s in range(n):
        seen = 1 << s
        stack = [s]
        while stack:
            u = stack.pop()
            m = adj[u] & ~seen
            while m:
                w = (m & -m).bit_length() - 1
                seen |= 1 << w
                stack.append(w)
                m &= m - 1
        if seen != full: return False
    return True

# ---- pair-bit machinery -----------------------------------------------------
def pair_index(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    return pairs, idx

def perm_tables(n):
    """for each permutation: list of (src_bit, dst_bit, flip)"""
    pairs, idx = pair_index(n)
    tabs = []
    for p in permutations(range(n)):
        t = []
        for k, (i, j) in enumerate(pairs):
            a, b = p[i], p[j]
            if a < b: t.append((k, idx[(a, b)], 0))
            else:     t.append((k, idx[(b, a)], 1))
        tabs.append(t)
    return pairs, tabs

def canon_int(bits, tabs, npairs):
    best = None
    for t in tabs:
        v = 0
        for (s, d, f) in t:
            if ((bits >> s) & 1) ^ f: v |= 1 << d
        if best is None or v < best: best = v
    return best

def bits_to_adj(bits, n, pairs):
    adj = [0] * n
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: adj[i] |= 1 << j
        else: adj[j] |= 1 << i
    return adj

# ---- stage 1: canonical n=6 reps -------------------------------------------
print("stage 1: n=6 classes...")
pairs6, tabs6 = perm_tables(N6)
reps6 = set()
for bits in range(1 << 15):
    reps6.add(canon_int(bits, tabs6, 15))
reps6 = sorted(reps6)
print("  n=6 classes: %d (expect 56)" % len(reps6))
assert len(reps6) == 56

# ---- stage 2: augment to n=7, full h-spectrum ------------------------------
print("stage 2: augment 56 x 64 -> full n=7 spectrum...")
pairs7, idx7 = pair_index(N7)
h_all, h_strong = {}, {}
cands = []
for r6 in reps6:
    adj6 = bits_to_adj(r6, N6, pairs6)
    for pat in range(1 << N6):
        adj = [adj6[v] | ((1 << N6) if not (pat >> v) & 1 else 0) for v in range(N6)]
        adj.append(pat)          # new vertex 6: beats v iff pat bit set
        h = hpaths(adj, N7)
        h_all[h] = h_all.get(h, 0) + 1
        if is_strong(adj, N7):
            h_strong[h] = h_strong.get(h, 0) + 1
        cands.append(adj)
print("  candidates: %d" % len(cands))
print("\nCOMPLETE n=7 h-spectrum (%d distinct values):" % len(h_all))
print(" ", sorted(h_all))
print("\nh = 7 attained at n=7: %s" % (7 in h_all))
print("h = 21 attained at n=7: %s" % (21 in h_all))
print("h = 35 attained at n=7: %s" % (35 in h_all))
print("\nSTRONG n=7 spectrum (%d values): %s" % (len(h_strong), sorted(h_strong)))
print("MIN strong h at n=7: %d" % min(h_strong))
missing_odds = [v for v in range(1, max(h_all) + 1, 2) if v not in h_all]
print("missing odd values up to max: %s" % missing_odds)

# ---- stage 3: canonicalize to 456 reps for the n=8 stage -------------------
print("\nstage 3: canonicalize 3584 -> n=7 class reps (A000568 check)...")
pairs7l, tabs7 = perm_tables(N7)
reps7 = set()
for adj in cands:
    bits = 0
    for k, (i, j) in enumerate(pairs7l):
        if (adj[i] >> j) & 1: bits |= 1 << k
    reps7.add(canon_int(bits, tabs7, 21))
print("  n=7 classes reached: %d (expect 456)" % len(reps7))
assert len(reps7) == 456
with open("05-knowledge/results/n7_class_reps_boxeph_S152.txt", "w") as f:
    for r in sorted(reps7): f.write("%d\n" % r)
print("  reps saved.  DONE.")
