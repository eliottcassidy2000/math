#!/usr/bin/env python3
"""
opus-2026-07-05-S82 -- HYP-4119: THE SHADOW DICHOTOMY AT LEVEL 169 (exhaustive).

The 169-grid: for a primitive 12-family with residues pinned {1..12}, the strict witness
question at t = a/169 (want min_i dist(v_i a/169) >= 14/169 > 1/13) depends only on the
PATTERN kappa : r in {1..12} -> Z/13, where v_r = r + 13*k_r and kappa_r = k_r mod 13
(kappa_r = 0 <=> unlifted OR invisible lift k = 0 mod 13).

CELLS: a in [1,168], a != 0 mod 13   (156 cells; a = 13a' gives dist exactly 1/13, never strict).
KILL:  position r with digit kappa kills cell a  <=>  a*(r+13*kappa) mod 169 in [0,13] u [156,168]
       (equivalently the mod-13 digit of a*v lands in {0,12}; exactly 2 kappa per (a,r)).
FULL COVER (no strict witness): union of the 12 kill sets = all 156 cells.

KNOWN full covers: the DILATION SHADOWS P_lambda, lambda in units(169): the mod-169 shadow
of the tight family lambda*{1..12} (position r holds lambda*i, i = lambda^{-1} r mod 13).
Dilated standards are non-primitive (gcd = lambda) -- excluded upstream by the primitivity
split -- but any family SHARING the shadow mod 169 is also witness-free at this level.

QUESTION (exhaustive, not sampled): are the full covers with >= 7 visible (nonzero kappa)
EXACTLY the dilation shadows?  A DFS over patterns, cell-driven (branch on the least-covered
uncovered cell), prunes: (i) capacity -- remaining positions * 24 < uncovered; (ii) zeros
budget <= 5 (>= 7 visible).  Plus: bottom-6 vs top-6 kill anatomy, per-position coverage.
"""
import sys, time
from itertools import product

MOD = 169
CELLS = [a for a in range(1, 169) if a % 13 != 0]
CIDX = {a: i for i, a in enumerate(CELLS)}
FULL = (1 << 156) - 1

def bad(x):  # x = v*a mod 169 in the closed bad zone (dist <= 13/169, i.e. NOT >= 14/169)
    return x <= 13 or x >= 156

# kill masks: KM[r][kappa] = bitmask of cells killed
KM = [[0] * 13 for _ in range(13)]  # index r = 1..12
for r in range(1, 13):
    for kap in range(13):
        v = r + 13 * kap
        m = 0
        for a in CELLS:
            if bad((a * v) % MOD):
                m |= 1 << CIDX[a]
        KM[r][kap] = m

sizes = sorted({bin(KM[r][k]).count('1') for r in range(1, 13) for k in range(13)})
print(f"kill-set sizes (should be all 24): {sizes}")

# ---- dilation shadows
def shadow(lam):
    """pattern of the dilated standard lambda*{1..12} mod 169: position r holds lambda*i,
    i = lambda^{-1} r mod 13; kappa_r = ((lambda*i mod 169) - r)/13 mod 13."""
    laminv = pow(lam % 13, -1, 13)
    pat = {}
    for r in range(1, 13):
        i = (laminv * r) % 13
        val = (lam * i) % MOD
        assert val % 13 == r % 13
        pat[r] = ((val - r) // 13) % 13
    return tuple(pat[r] for r in range(1, 13))

shadows = {}
for lam in range(1, MOD):
    if lam % 13 != 0:
        shadows.setdefault(shadow(lam), []).append(lam)
print(f"distinct dilation shadows: {len(shadows)} (from 156 units)")

def cover_mask(pat):
    m = 0
    for r in range(1, 13):
        m |= KM[r][pat[r - 1]]
    return m

sh_full = {p: ls for p, ls in shadows.items() if cover_mask(p) == FULL}
print(f"shadows that are FULL covers: {len(sh_full)}/{len(shadows)} (tightness check: must be all)")
vis_hist = {}
for p in sh_full:
    nv = sum(1 for x in p if x != 0)
    vis_hist[nv] = vis_hist.get(nv, 0) + 1
print(f"shadow visible-count histogram: {dict(sorted(vis_hist.items()))}")

# ---- exhaustive DFS for ALL full covers (column-deficit pruning, cell-driven)
import time as _t
t0 = _t.time()
nodes = 0
found = []

killers = {a: [] for a in CELLS}
for r in range(1, 13):
    for kap in range(13):
        m = KM[r][kap]
        for a in CELLS:
            if m >> CIDX[a] & 1:
                killers[a].append((r, kap))
for a in CELLS:
    assert len(killers[a]) == 24

# per-column masks: column of cell a=(b+13c) is b = a % 13
COLMASK = [0]*13
for a in CELLS:
    COLMASK[a % 13] |= 1 << CIDX[a]

import sys
sys.setrecursionlimit(10000)

def dfs(assigned, covered):
    global nodes
    nodes += 1
    if nodes % 2000000 == 0:
        print(f"  ...{nodes} nodes, {len(found)} covers, {_t.time()-t0:.0f}s", flush=True)
    if covered == FULL:
        found.append(dict(assigned))
        return
    nrem = 12 - len(assigned)
    # column-deficit prune: each unassigned position adds <= 2 cells per column
    for b in range(1, 13):
        deficit = 13 - bin(covered & COLMASK[b]).count("1")
        if deficit > 2 * nrem:
            return
    # branch on the uncovered cell with fewest viable killers
    best_a, best_opts = None, None
    for a in CELLS:
        if covered >> CIDX[a] & 1:
            continue
        opts = [(r, kap) for (r, kap) in killers[a] if r not in assigned]
        if best_opts is None or len(opts) < len(best_opts):
            best_a, best_opts = a, opts
            if len(opts) <= 2:
                break
    if not best_opts:
        return
    for (r, kap) in best_opts:
        assigned[r] = kap
        dfs(assigned, covered | KM[r][kap])
        del assigned[r]

dfs({}, 0)
print(f"DFS complete: {nodes} nodes, {_t.time()-t0:.1f}s; cover-CORES found: {len(found)}", flush=True)

# classify cores: a core (some positions assigned, rest FREE) is a family of patterns.
# Deduplicate into full patterns only where needed: check whether every completion
# pattern that is a full cover is a dilation shadow.
sh_set = set(shadows.keys())
core_patterns = set()
for assigned in found:
    free = [r for r in range(1, 13) if r not in assigned]
    if len(free) > 4:
        print(f"  WIDE core (free={len(free)}): {assigned} -- classifying by free expansion")
    from itertools import product as _prod
    for combo in _prod(range(13), repeat=len(free)):
        pat = [0]*12
        for r in range(1, 13):
            pat[r-1] = assigned.get(r, 0)
        for r, kv in zip(free, combo):
            pat[r-1] = kv
        core_patterns.add(tuple(pat))
print(f"distinct full-cover patterns (expanded): {len(core_patterns)}", flush=True)
really_full = {p for p in core_patterns if cover_mask(p) == FULL}
print(f"  of which actually full: {len(really_full)}")
non_shadow = really_full - sh_set
print(f"  NON-SHADOW full covers: {len(non_shadow)}")
import collections as _c
hist = _c.Counter(sum(1 for x in p if x) for p in really_full)
print(f"  visible histogram of ALL full covers: {dict(sorted(hist.items()))}")
if non_shadow:
    hist2 = _c.Counter(sum(1 for x in p if x) for p in non_shadow)
    print(f"  visible histogram of NON-SHADOW covers: {dict(sorted(hist2.items()))}")
    for p in sorted(non_shadow)[:10]:
        print(f"    non-shadow: kappa={p} visible={sum(1 for x in p if x)}")
else:
    print("  *** SHADOW DICHOTOMY HOLDS: full covers at level 169 = the 156 dilation shadows, exactly")

# ---- bottom-6 anatomy (alignment theme)
print()
print("bottom-6 vs top-6 kill anatomy:")
import collections
overlap = collections.Counter()
for a in CELLS:
    bl = sum(1 for (r, k) in killers[a] if r <= 6)
    overlap[bl] += 1
print(f"  per-cell killer count from bottom-6 positions (of 24): {dict(sorted(overlap.items()))}")
# which positions are FORCED in covers: appear in every found core?
pos_freq = collections.Counter()
for assigned, _ in found:
    for r in assigned:
        pos_freq[r] += 1
print(f"  position frequency in cover cores: {dict(sorted(pos_freq.items()))}")
print(f"TOTAL runtime {time.time()-t0:.1f}s")
