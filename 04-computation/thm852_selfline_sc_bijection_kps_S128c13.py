#!/usr/bin/env python3
"""
thm852_selfline_sc_bijection_kps_S128c13.py
===========================================
kind-pasteur-2026-07-15-S128 (cont.13, overnight).  THE BLACK SELF-LINE LAW (opus S306/S310/S311):

    #{non-grid-sym tilings t : T(kappa t) iso T(t)}  =  SC(n)     (n >= 5)

kappa = complement tiling (flip all m tile bits).  Verified 8, 12, 88 at n = 5, 6, 7 by opus.
THIS SCRIPT: (1) structure data for the bijection hunt at n = 5, 6 (and 7 verify):
    - the Klein-four action {1, g, kappa, g kappa} on the quasi-fixed set X (g = grid reflection)
    - per-tiling: gridsym?, cls, cls SC?, |Aut|, the iso-witness sigmas, g-orbit structure
    - WEIGHTED counts W = sum |Aut| over X (the Burnside-affine dual)
 (2) THE n=8 CHECK by invariant-filter + exact iso:
    - all 2^21 tilings: compute (sorted scores, c3) invariant of T(t) and T(kappa t); survivors
      (matching invariants) get exact iso test (permutation search with pruning); count
      non-grid-sym quasi-fixed tilings; TARGET = SC(8) = 176 (from merged=3528 = (6880+SC)/2).
 (3) Burnside-affine cross-check at n = 5, 6, 7, 8: W(n) = sum_sigma #{t : sigma T(t) = T(kappa t)}
      via affine GF(2) systems -- fast, exact, independent of (2).
"""
import sys
from math import comb
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def setup(n):
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    m = len(tiles)
    gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
    return tiles, tidx, m, gmap

def beats_from_bits(tv, n, tidx):
    B = [[False] * n for _ in range(n)]
    for k in range(2, n + 1):
        B[k - 1][k - 2] = True
    for (x, y), i in tidx.items():
        if tv[i]:
            B[x - 1][y - 1] = True
        else:
            B[y - 1][x - 1] = True
    return B

def inv_pair(B, n):
    d = [sum(B[u]) for u in range(n)]
    c3 = comb(n, 3) - sum(comb(x, 2) for x in d)
    return (tuple(sorted(d)), c3)

def iso(B1, B2, n):
    """exact iso test with degree-class pruning."""
    d1 = [sum(B1[u]) for u in range(n)]
    d2 = [sum(B2[u]) for u in range(n)]
    if sorted(d1) != sorted(d2):
        return False
    buckets = defaultdict(list)
    for v in range(n):
        buckets[d2[v]].append(v)
    order = sorted(range(n), key=lambda u: (d1[u], u))
    assign = [-1] * n
    used = [False] * n
    def bt(i):
        if i == n:
            return True
        u = order[i]
        for v in buckets[d1[u]]:
            if used[v]:
                continue
            ok = True
            for j in range(i):
                w = order[j]
                if B1[u][w] != B2[assign[u] if False else v][assign[w]]:
                    ok = False
                    break
                if B1[w][u] != B2[assign[w]][v]:
                    ok = False
                    break
            if ok:
                assign[u] = v
                used[v] = True
                if bt(i + 1):
                    return True
                used[v] = False
                assign[u] = -1
        return False
    return bt(0)

# ---------- (1)+(3) small n structure + weighted counts ----------
print("=" * 96)
print("(1) structure of the quasi-fixed set X, n=5,6 (exact, full enumeration)")
print("=" * 96)
for n in [5, 6]:
    tiles, tidx, m, gmap = setup(n)
    X = []
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        kv = [1 - b for b in tv]
        B1 = beats_from_bits(tv, n, tidx)
        B2 = beats_from_bits(kv, n, tidx)
        if inv_pair(B1, n) == inv_pair(B2, n) and iso(B1, B2, n):
            gs = all(tv[i] == tv[gmap[i]] for i in range(m))
            X.append((t, gs))
    ngs = [t for t, gs in X if not gs]
    print("n=%d: |X|=%d  gridsym=%d  NON-gridsym=%d  (SC(n) target: %s)"
          % (n, len(X), len(X) - len(ngs), len(ngs), {5: 8, 6: 12, 7: 88}.get(n, '?')))
    # Klein-four orbit structure on the non-gridsym part
    orb = set()
    reps = []
    for t in ngs:
        if t in orb:
            continue
        tv = [(t >> i) & 1 for i in range(m)]
        gt = sum(tv[gmap[i]] << i for i in range(m))
        kt = t ^ ((1 << m) - 1)
        gkt = gt ^ ((1 << m) - 1)
        o = {t, gt, kt, gkt}
        orb |= o
        reps.append(sorted(o))
    print("   Klein-4 orbits on non-gridsym X: %d orbits, sizes %s"
          % (len(reps), Counter(len(set(o)) for o in reps)))

# ---------- (3) Burnside-affine weighted counts ----------
print()
print("=" * 96)
print("(3) Burnside-affine: W(n) = sum_sigma #{t : sigma T(t) = T(kappa t)} = sum_{t in X} |Aut(T(t))|")
print("=" * 96)
def burnside_affine(n, split_gridsym=True):
    tiles, tidx, m, gmap = setup(n)
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    # arc slot -> (tile index or None base, orientation meaning)
    # T(t): arc between a<b: if b-a==1: b beats a. else tile (b+1?, ...) -- vertices 1..n: pair (a,b), a<b
    # tile (x,y): x>y+1, bit=1 => x beats y.
    def slot(a, b):
        # a<b in 0-index -> vertices a+1<b+1
        if b - a == 1:
            return None, True   # b beats a always (base)
        return tidx[(b + 1, a + 1)], None
    Wtot = 0
    Wgs = 0
    from itertools import permutations as perms
    for sig in perms(range(n)):
        # constraints: for each pair (u,v): orientation of sigma T(t) on (sig u, sig v) == T(kappa t) on (sig u, sig v)
        # sigma T(t) arc sig(u)->sig(v) iff T(t) has u->v.
        # Build linear system over GF(2) with variables t_0..t_{m-1}; equation per unordered pair (p<q):
        # [T(t) orientation on (u,v) with (sig u, sig v)=(p,q)] == [T(kappa t) orientation on (p,q)]
        # encode orientation as bit: 1 iff higher-index vertex beats lower? use: arrow p->q (p<q)? define f(p,q)=1 iff p beats q.
        rows = []
        rhs = []
        ok = True
        for (u, v) in pairs:  # u<v
            p, q = sig[u], sig[v]
            if p > q:
                p, q = q, p
                uu, vv = v, u
            else:
                uu, vv = u, v
            # f1 = [p beats q in sigma T(t)] = [uu beats vv in T(t)]
            s1, base1 = slot(uu, vv) if uu < vv else slot(vv, uu)
            # orientation: for pair (a<b): "a beats b"? base: b beats a => f(a,b)=0. tile (b+1,a+1) bit=1 => b+1 beats a+1 => f(a,b)=0 when bit=1? f(a,b)= [a beats b] = bit==0 => f = 1 + t_e mod 2.
            # careful with uu<vv or vv<uu:
            a, b = (uu, vv) if uu < vv else (vv, uu)
            e, isbase = slot(a, b)
            if isbase:
                f1_const, f1_var = (0, None)   # f(a,b) = 0 (b beats a)
            else:
                f1_const, f1_var = (1, e)       # f(a,b) = 1 + t_e
            if uu > vv:   # we want [uu beats vv] = 1 - f(a,b)
                f1_const ^= 1
            # f2 = [p beats q in T(kappa t)]: pair (p<q): slot e2: kappa flips bits: f(p,q) = 1 + (1+t_e2) = t_e2; base unchanged.
            e2, isbase2 = slot(p, q)
            if isbase2:
                f2_const, f2_var = 0, None
            else:
                f2_const, f2_var = 0, e2        # f = 1 + kbit = 1 + (1 + t) = t
            # equation: f1 == f2
            row = [0] * m
            c = f1_const ^ f2_const
            if f1_var is not None:
                row[f1_var] ^= 1
            if f2_var is not None:
                row[f2_var] ^= 1
            if any(row):
                rows.append((row, c))
            elif c:
                ok = False
                break
        if not ok:
            continue
        # gauss over GF(2)
        mat = [r[:] + [c] for r, c in rows]
        rank = 0
        cols = m
        for col in range(cols):
            piv = None
            for r in range(rank, len(mat)):
                if mat[r][col]:
                    piv = r
                    break
            if piv is None:
                continue
            mat[rank], mat[piv] = mat[piv], mat[rank]
            for r in range(len(mat)):
                if r != rank and mat[r][col]:
                    mat[r] = [a ^ b for a, b in zip(mat[r], mat[rank])]
            rank += 1
        consistent = all(not (all(x == 0 for x in row[:-1]) and row[-1]) for row in mat)
        if consistent:
            Wtot += 1 << (m - rank)
    return Wtot

for n in [5, 6, 7]:
    W = burnside_affine(n)
    print("n=%d: W = sum_{t in X} |Aut| = %d" % (n, W))
print("\n(n=8 Burnside-affine and the direct n=8 filter run in part (2) -- see below)")

# ---------- (2) n=8 exact check ----------
print()
print("=" * 96)
print("(2) n=8 CHECK: count non-gridsym quasi-fixed tilings; TARGET SC(8) = 176")
print("=" * 96)
n = 8
tiles, tidx, m, gmap = setup(n)
print("m=%d, tilings=%d; pass 1 = invariant filter (scores, c3)..." % (m, 1 << m))
surv = []
for t in range(1 << m):
    tv = [(t >> i) & 1 for i in range(m)]
    B1 = beats_from_bits(tv, n, tidx)
    kv = [1 - b for b in tv]
    B2 = beats_from_bits(kv, n, tidx)
    if inv_pair(B1, n) == inv_pair(B2, n):
        surv.append(t)
    if t % (1 << 18) == 0:
        print("   ...%d/%d, survivors so far %d" % (t, 1 << m, len(surv)), flush=True)
print("pass 1 survivors: %d; pass 2 = exact iso..." % len(surv))
cnt_qf = 0
cnt_ngs = 0
for idx, t in enumerate(surv):
    tv = [(t >> i) & 1 for i in range(m)]
    kv = [1 - b for b in tv]
    B1 = beats_from_bits(tv, n, tidx)
    B2 = beats_from_bits(kv, n, tidx)
    if iso(B1, B2, n):
        cnt_qf += 1
        if not all(tv[i] == tv[gmap[i]] for i in range(m)):
            cnt_ngs += 1
    if idx % 5000 == 0:
        print("   ...iso %d/%d, qf=%d ngs=%d" % (idx, len(surv), cnt_qf, cnt_ngs), flush=True)
print("n=8 RESULT: quasi-fixed=%d  gridsym-qf=%d  NON-gridsym-qf=%d  TARGET SC(8)=176: %s"
      % (cnt_qf, cnt_qf - cnt_ngs, cnt_ngs, "MATCH" if cnt_ngs == 176 else "** MISMATCH **"))
print("DONE")
