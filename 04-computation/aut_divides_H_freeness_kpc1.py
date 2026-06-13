#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aut_divides_H_freeness_kpc1.py
kind-pasteur-2026-06-10-S1, Thread B (HYP-2369 / LEM-003 candidate)

CLAIM (universal freeness): for ANY digraph D, Aut(D) acts FREELY on the set
of directed Hamiltonian paths (a directed Ham path's arc set determines its
vertex sequence: unique in-degree-0 start vertex, then induct along the path;
an automorphism fixing the arc set therefore fixes every vertex). Hence every
orbit has size exactly |Aut(D)| and |Aut(D)| | H(D).

CYCLE CAVEAT: the statement FAILS for directed Hamiltonian CYCLES (no
distinguished start; rotations can fix a cycle's arc set). Demonstrated below.

This script verifies, with EXACT integer arithmetic:
  PART 0: self-test of the bit-table isomorphism transform vs naive (n=5).
  PART 1: ALL 2^10 = 1024 labeled n=5 tournaments: H by Held-Karp bitmask DP,
          |Aut| by brute force over S_5; assert |Aut| | H for every one.
  PART 2: ALL 2^15 = 32768 labeled n=6 tournaments: same with S_6.
          Burnside cross-check: sum(|Aut|) = n! * #iso-classes (A000568).
  PART 3: EXPLICIT ORBIT CONSTRUCTION (freeness, not just divisibility):
          for EVERY labeled n=5 and n=6 tournament with |Aut| > 1, and for
          100 random n=6 tournaments, build the Aut-action on Ham-path
          sequences and assert every orbit has size exactly |Aut|
          (equivalently: no non-identity automorphism fixes any Ham path),
          and #orbits * |Aut| == H.
  PART 4: CYCLE CAVEAT, explicit: cyclic triangle C3 (1 Ham cycle, |Aut|=3,
          3 does not divide 1, rotation fixes the cycle) and the regular
          circulant RQ5 (i->i+1,i+2 mod 5): rotation-fixed cycles exist,
          orbit sizes are NOT all |Aut|, divisibility fails.

Pure integers throughout. No numpy/sympy. Runtime target: a few minutes.
"""

import itertools
import random
import time

T0 = time.time()


def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def adj_from_mask(mask, n, pairs):
    """adj[v] = bitmask of out-neighbours. bit k of mask: 1 means i->j for pairs[k]=(i,j)."""
    adj = [0] * n
    for k, (i, j) in enumerate(pairs):
        if (mask >> k) & 1:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return adj


def mask_from_adj(adj, n, pairs):
    m = 0
    for k, (i, j) in enumerate(pairs):
        if (adj[i] >> j) & 1:
            m |= 1 << k
    return m


def naive_transform(mask, perm, n, pairs, pidx):
    """Image of tournament `mask` under vertex permutation perm (naive)."""
    out = 0
    for k, (i, j) in enumerate(pairs):
        bit = (mask >> k) & 1
        a, b = perm[i], perm[j]
        if a < b:
            t, flip = pidx[(a, b)], 0
        else:
            t, flip = pidx[(b, a)], 1
        if bit ^ flip:
            out |= 1 << t
    return out


def build_tables(n):
    """For each permutation sigma of [n], the map mask -> sigma(mask) equals
    TLO[s][mask & 255] ^ THI[s][mask >> 8] ^ F[s] (disjoint bit-permutation
    tables XOR a constant flip mask). Returns (perms, TLO, THI, F, pairs, pidx)."""
    pairs = pairs_of(n)
    m = len(pairs)
    pidx = {p: k for k, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    TLO, THI, F = [], [], []
    lo_bits = min(8, m)
    for perm in perms:
        tgt = [0] * m
        f = 0
        for k, (i, j) in enumerate(pairs):
            a, b = perm[i], perm[j]
            if a < b:
                tgt[k] = pidx[(a, b)]
            else:
                tgt[k] = pidx[(b, a)]
                f |= 1 << tgt[k]
        bitmap = [1 << t for t in tgt]
        tl = [0] * (1 << lo_bits)
        for x in range(1, 1 << lo_bits):
            lsb = x & -x
            tl[x] = tl[x ^ lsb] ^ bitmap[lsb.bit_length() - 1]
        hi_n = m - lo_bits
        th = [0] * (1 << max(hi_n, 0))
        for x in range(1, 1 << hi_n):
            lsb = x & -x
            th[x] = th[x ^ lsb] ^ bitmap[lo_bits + lsb.bit_length() - 1]
        TLO.append(tl)
        THI.append(th)
        F.append(f)
    return perms, TLO, THI, F, pairs, pidx


def ham_count(adj, n):
    """Held-Karp bitmask DP: number of directed Hamiltonian paths (exact int)."""
    full = (1 << n) - 1
    dp = [0] * ((full + 1) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for S in range(1, full + 1):
        b = S * n
        for v in range(n):
            if not (S >> v) & 1:
                continue
            c = dp[b + v]
            if not c:
                continue
            t = adj[v] & ~S & full
            while t:
                w = t & -t
                t ^= w
                dp[(S | w) * n + (w.bit_length() - 1)] += c
    return sum(dp[full * n + v] for v in range(n))


def ham_paths(adj, n):
    """All directed Hamiltonian paths as vertex sequences (tuples)."""
    full = (1 << n) - 1
    out = []

    def rec(v, vis, seq):
        if vis == full:
            out.append(tuple(seq))
            return
        t = adj[v] & ~vis
        while t:
            w = t & -t
            t ^= w
            wi = w.bit_length() - 1
            seq.append(wi)
            rec(wi, vis | w, seq)
            seq.pop()

    for s in range(n):
        rec(s, 1 << s, [s])
    return out


def ham_cycles_arcsets(adj, n):
    """All directed Hamiltonian cycles as ARC SETS (frozensets). Each cycle
    counted once (canonical traversal starts at vertex 0)."""
    full = (1 << n) - 1
    res = []

    def rec(v, vis, seq):
        if vis == full:
            if (adj[v] >> seq[0]) & 1:
                res.append(frozenset(zip(seq, seq[1:] + [seq[0]])))
            return
        t = adj[v] & ~vis
        while t:
            w = t & -t
            t ^= w
            wi = w.bit_length() - 1
            rec(wi, vis | w, seq + [wi])

    rec(0, 1, [0])
    return res


def aut_list(adj, n, perms):
    """Brute-force automorphism list (as permutation tuples)."""
    out = []
    for perm in perms:
        ok = True
        for u in range(n):
            au = adj[u]
            pu = perm[u]
            for v in range(n):
                if u != v and ((au >> v) & 1) != ((adj[pu] >> perm[v]) & 1):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(perm)
    return out


# ===================================================================
print("=" * 70)
print("PART 0: SELF-TEST — bit-table transform vs naive (n=5, 300 random)")
print("=" * 70)
n = 5
perms5, TLO5, THI5, F5, pairs5, pidx5 = build_tables(5)
rng = random.Random(20260610)
for _ in range(300):
    mask = rng.randrange(1 << 10)
    s = rng.randrange(120)
    fast = TLO5[s][mask & 255] ^ THI5[s][mask >> 8] ^ F5[s]
    slow = naive_transform(mask, perms5[s], 5, pairs5, pidx5)
    assert fast == slow, (mask, s)
print("  300/300 random (mask, sigma) pairs agree. Self-test PASSED.")

# ===================================================================
print()
print("=" * 70)
print("PART 1: ALL 2^10 LABELED n=5 TOURNAMENTS — assert |Aut| | H")
print("=" * 70)
t1 = time.time()
N5 = 1 << 10
autcnt5 = [1] * N5  # identity counted
for s in range(1, 120):
    tl, th, f = TLO5[s], THI5[s], F5[s]
    for m in range(N5):
        if tl[m & 255] ^ th[m >> 8] ^ f == m:
            autcnt5[m] += 1
H5 = [0] * N5
fail5 = 0
h0_5 = 0
all_odd5 = True
for m in range(N5):
    adj = adj_from_mask(m, 5, pairs5)
    h = ham_count(adj, 5)
    H5[m] = h
    if h == 0:
        h0_5 += 1
    if h % 2 == 0:
        all_odd5 = False
    if h % autcnt5[m] != 0:
        fail5 += 1
        print("  FAIL: mask=%d |Aut|=%d H=%d" % (m, autcnt5[m], h))
burn5 = sum(autcnt5)
print("  tournaments checked : %d" % N5)
print("  divisibility failures: %d" % fail5)
print("  H = 0 cases          : %d  (Redei: impossible for tournaments)" % h0_5)
print("  all H odd (Redei)    : %s" % all_odd5)
print("  |Aut| values seen    : %s" % sorted(set(autcnt5)))
print("  Burnside check: sum |Aut| = %d = 120 * %d  (A000568(5) = 12: %s)"
      % (burn5, burn5 // 120, burn5 == 120 * 12))
assert fail5 == 0 and burn5 == 120 * 12
print("  PART 1 PASSED in %.1fs" % (time.time() - t1))

# ===================================================================
print()
print("=" * 70)
print("PART 2: ALL 2^15 LABELED n=6 TOURNAMENTS — assert |Aut| | H")
print("=" * 70)
t2 = time.time()
perms6, TLO6, THI6, F6, pairs6, pidx6 = build_tables(6)
N6 = 1 << 15
autcnt6 = [1] * N6
R6 = range(N6)
for s in range(1, 720):
    tl, th, f = TLO6[s], THI6[s], F6[s]
    fixed = [m for m in R6 if tl[m & 255] ^ th[m >> 8] ^ f == m]
    for m in fixed:
        autcnt6[m] += 1
print("  Aut counts done in %.1fs" % (time.time() - t2))
t2b = time.time()
H6 = [0] * N6
fail6 = 0
h0_6 = 0
all_odd6 = True
hmax = 0
for m in R6:
    adj = adj_from_mask(m, 6, pairs6)
    h = ham_count(adj, 6)
    H6[m] = h
    if h == 0:
        h0_6 += 1
    if h % 2 == 0:
        all_odd6 = False
    if h > hmax:
        hmax = h
    if h % autcnt6[m] != 0:
        fail6 += 1
        print("  FAIL: mask=%d |Aut|=%d H=%d" % (m, autcnt6[m], h))
burn6 = sum(autcnt6)
nontriv6 = sum(1 for a in autcnt6 if a > 1)
print("  Held-Karp H for all masks done in %.1fs" % (time.time() - t2b))
print("  tournaments checked : %d" % N6)
print("  divisibility failures: %d" % fail6)
print("  H = 0 cases          : %d  (Redei: impossible for tournaments)" % h0_6)
print("  all H odd (Redei)    : %s" % all_odd6)
print("  max H at n=6         : %d" % hmax)
print("  |Aut| values seen    : %s" % sorted(set(autcnt6)))
print("  masks with |Aut| > 1 : %d" % nontriv6)
print("  Burnside check: sum |Aut| = %d = 720 * %d  (A000568(6) = 56: %s)"
      % (burn6, burn6 // 720, burn6 == 720 * 56))
assert fail6 == 0 and burn6 == 720 * 56
print("  PART 2 PASSED in %.1fs total" % (time.time() - t2))

# ===================================================================
print()
print("=" * 70)
print("PART 3: EXPLICIT ORBITS — freeness (orbit size == |Aut| exactly)")
print("=" * 70)
t3 = time.time()


def check_orbits(mask, n, pairs, perms, expected_aut, expected_H):
    """Build the Aut-action on Ham-path sequences explicitly.
    Returns (ok, |Aut|, H, #orbits). ok = every orbit has size exactly |Aut|."""
    adj = adj_from_mask(mask, n, pairs)
    A = aut_list(adj, n, perms)
    if expected_aut is not None:
        assert len(A) == expected_aut, (mask, len(A), expected_aut)
    paths = ham_paths(adj, n)
    if expected_H is not None:
        assert len(paths) == expected_H, (mask, len(paths), expected_H)
    pathset = set(paths)
    orbits = set()
    ok = True
    for p in paths:
        imgs = set()
        for sg in A:
            q = tuple(sg[v] for v in p)
            if q not in pathset:   # sanity: automorphisms map Ham paths to Ham paths
                return (False, len(A), len(paths), -1)
            imgs.add(q)
        if len(imgs) != len(A):    # some non-identity sigma fixed p
            ok = False
        orbits.add(frozenset(imgs))
    if len(orbits) * len(A) != len(paths):
        ok = False
    return (ok, len(A), len(paths), len(orbits))


# 3a: EVERY labeled n=5 tournament with |Aut| > 1
bad = 0
cnt = 0
for m in range(N5):
    if autcnt5[m] > 1:
        ok, a, h, no = check_orbits(m, 5, pairs5, perms5, autcnt5[m], H5[m])
        cnt += 1
        if not ok:
            bad += 1
            print("  n=5 FREENESS FAIL: mask=%d |Aut|=%d H=%d" % (m, a, h))
print("  3a: all %d labeled n=5 tournaments with |Aut|>1: freeness failures = %d" % (cnt, bad))
assert bad == 0

# 3b: EVERY labeled n=6 tournament with |Aut| > 1 (exhaustive where non-vacuous)
bad = 0
cnt = 0
aut_orbit_profile = {}
for m in R6:
    if autcnt6[m] > 1:
        ok, a, h, no = check_orbits(m, 6, pairs6, perms6, autcnt6[m], H6[m])
        cnt += 1
        if not ok:
            bad += 1
            print("  n=6 FREENESS FAIL: mask=%d |Aut|=%d H=%d" % (m, a, h))
        aut_orbit_profile.setdefault((a, h, no), 0)
        aut_orbit_profile[(a, h, no)] = aut_orbit_profile[(a, h, no)] + 1
print("  3b: all %d labeled n=6 tournaments with |Aut|>1: freeness failures = %d" % (cnt, bad))
print("      (|Aut|, H, #orbits) profile (orbits = H/|Aut| exactly):")
for k in sorted(aut_orbit_profile):
    print("        |Aut|=%d H=%2d orbits=%2d : %d tournaments" % (k[0], k[1], k[2], aut_orbit_profile[k]))
assert bad == 0

# 3c: 100 random n=6 tournaments (any |Aut|), explicit action
rng = random.Random(20260610)
sample = [rng.randrange(N6) for _ in range(100)]
bad = 0
ntriv = 0
for m in sample:
    ok, a, h, no = check_orbits(m, 6, pairs6, perms6, autcnt6[m], H6[m])
    if a > 1:
        ntriv += 1
    if not ok:
        bad += 1
        print("  RANDOM FAIL: mask=%d" % m)
print("  3c: 100 random n=6 masks: freeness failures = %d (%d had |Aut|>1)" % (bad, ntriv))
assert bad == 0
print("  PART 3 PASSED in %.1fs" % (time.time() - t3))

# ===================================================================
print()
print("=" * 70)
print("PART 4: THE CYCLE CAVEAT — freeness FAILS for Hamiltonian CYCLES")
print("=" * 70)


def cycle_report(name, adj, n, perms):
    A = aut_list(adj, n, perms)
    cycles = ham_cycles_arcsets(adj, n)
    HC = len(cycles)
    cycset = set(cycles)
    fixed_pairs = 0   # (non-identity sigma, cycle) with sigma(cycle) == cycle
    orbits = set()
    for c in cycles:
        imgs = set()
        for sg in A:
            img = frozenset((sg[u], sg[v]) for (u, v) in c)
            assert img in cycset
            imgs.add(img)
            if img == c and sg != tuple(range(n)):
                fixed_pairs += 1
        orbits.add(frozenset(imgs))
    orbit_sizes = sorted(len(o) for o in orbits)
    # NOTE: orbits as sets of arc-sets; a cycle fixed by k automorphisms sits
    # in an orbit of size |Aut|/(k+1) by orbit-stabilizer.
    div = (HC % len(A) == 0)
    print("  %s:" % name)
    print("    |Aut| = %d, #HamCycles (arc sets) = %d" % (len(A), HC))
    print("    |Aut| divides #HC ? %s  (%d %% %d = %d)"
          % (div, HC, len(A), HC % len(A)))
    print("    (non-identity sigma, fixed cycle) incidences: %d  -> action %s"
          % (fixed_pairs, "NOT FREE" if fixed_pairs else "free"))
    print("    orbit sizes: %s  (free would force all = %d)" % (orbit_sizes, len(A)))
    return div, fixed_pairs


# C3: cyclic triangle 0->1->2->0
perms3 = list(itertools.permutations(range(3)))
adjC3 = [0] * 3
for (u, v) in [(0, 1), (1, 2), (2, 0)]:
    adjC3[u] |= 1 << v
div3, fp3 = cycle_report("C3 cyclic triangle (0->1->2->0)", adjC3, 3, perms3)
assert not div3 and fp3 > 0

# RQ5: circulant on Z_5, i -> i+1, i+2
adjR5 = [0] * 5
for i in range(5):
    for d in (1, 2):
        adjR5[i] |= 1 << ((i + d) % 5)
div5, fp5 = cycle_report("RQ5 circulant (i -> i+1, i+2 mod 5)", adjR5, 5, perms5)
assert fp5 > 0
# contrast: paths on the same tournaments ARE free
okP, aP, hP, noP = check_orbits(mask_from_adj(adjR5, 5, pairs_of(5)), 5, pairs5, perms5, None, None)
print("  CONTRAST (paths on RQ5): |Aut|=%d, H=%d, orbits=%d, all of size |Aut|: %s"
      % (aP, hP, noP, okP))
assert okP

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("  LEM (universal freeness): Aut acts FREELY on directed Ham PATHS;")
print("  |Aut| | H verified for ALL 1024 labeled n=5 and ALL 32768 labeled")
print("  n=6 tournaments; orbit sizes == |Aut| verified explicitly for every")
print("  tournament with |Aut|>1 at n=5,6 plus 100 random n=6 masks.")
print("  CYCLE caveat confirmed: C3 has 1 Ham cycle, |Aut|=3, 3 does not")
print("  divide 1 (rotation fixes the cycle). Paths are order-rigid; cycles")
print("  are not. Zero Paley/QR/Eisenstein content anywhere in the proof.")
print("  Total runtime: %.1fs" % (time.time() - T0))
