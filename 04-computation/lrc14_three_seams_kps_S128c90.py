#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c90 -- HYP-7955: three creative seams investigated
concurrently (owner: DFS through past threads, propose + test new hypotheses).

SEAM G -- THE QUANTIZED SEVEN-TILE WALL.  c(p) := minimum number of folded
danger-APs A(w) = {fold(jw) : j <= dk}, dk = floor(p/14), needed to cover the
folded line P = [1..(p-1)/2].  Density: |P|/|A(w)| -> 7 exactly as p grows --
the r=7 crown's zero-coefficient wall (THM-1153/1155) in discrete costume.
S128c88 data gave c(p) = 9 at p = 43, 61, 71.  Here: c(p) across more p, the
excess c(p) - ceil(h/dk), and a witness minimum cover's structure.

SEAM H -- MOMENT BLINDNESS / THE PTE BRIDGE.  The cage (HYP-7940) pins
families by even power sums to degree 26; the Vitali wall (claude 06-03,
phase-1 era) says finite moments cannot decide tightness.  Quantify the
blindness PROFILE: build (S2,S4)-fibers (and (S2,S4,S6)) of 13-sets from
PTE-style collisions of 3-subsets (4-subsets), and measure how much M varies
within a fiber.  Target: a triage-law axis "moment degree >= d* or blind."

SEAM J -- THE CRT INTERFERENCE WITNESS.  boxeph-S130's near-miss: no multiple
of 29 => M >= 2/29 = 96.6% of 1/14 (the 1/406 sliver, MISTAKE-093).  Post
THM-1289 the micro-gap is unconditionally empty but INEFFECTIVE; effective
witnesses still matter.  Test: does the two-prime grid q = 29*31 = 899 (and
the single primes) recover >= 1/14 witnesses on the no-mult-29 class?
"""
import sys, time, random
from math import gcd
from fractions import Fraction as F
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build
from lrc14_ladder_realization_crossN_kps_S128c86 import M_exact_arg

random.seed(90)

# ---------------- seam G ----------------
def cover_exists(p, size_cap, budget=30_000_000):
    """Does a cover of P by <= size_cap folded danger-APs exist? MRV search."""
    h, dk, maskA, cand, FULL, fold, inv = build(p)
    nodes = 0
    usedset = set()
    found = [None]
    def rec(U, banned):
        nonlocal nodes
        if found[0]: return
        nodes += 1
        if nodes > budget: raise TimeoutError
        if U == 0:
            found[0] = sorted(usedset); return
        rem = size_cap - len(usedset)
        if rem <= 0: return
        if U.bit_count() > dk*rem: return
        V = U; seen = 0; bestcs = None
        while V and seen < 24:
            i = (V & -V).bit_length() - 1
            x = i + 1
            cs = [w for w in cand[x] if w not in usedset and w not in banned]
            if not cs: return
            if bestcs is None or len(cs) < len(bestcs):
                bestcs = cs
                if len(cs) <= 1: break
            V &= V - 1; seen += 1
        newban = set()
        for w in bestcs:
            if found[0]: return
            usedset.add(w)
            rec(U & ~maskA[w], banned | newban)
            usedset.discard(w)
            newban.add(w)
    # WLOG 1 in W? NO -- a size-capped cover need not contain 1 (scaling maps
    # covers to covers but changes WHICH w's; the class of covers is
    # scaling-invariant so testing 1-in-W suffices for EXISTENCE).
    usedset.add(1)
    try:
        rec(FULL & ~maskA[1], frozenset())
    except TimeoutError:
        return None, nodes  # budget exceeded: unknown
    return found[0], nodes

def seamG():
    print("== SEAM G: the quantized seven-tile wall ==", flush=True)
    print("  c(p) = min # folded danger-APs covering P; density bound ceil(h/dk)", flush=True)
    for p in (29, 43, 61, 71, 101, 113, 127, 151):
        h = (p-1)//2; dk = p//14
        triv = -(-h // dk)
        row = f"  p={p:3d} h={h:2d} dk={dk:2d} density-bound={triv}"
        cp = None
        for size in range(max(triv, 7), 12):
            t0 = time.time()
            W, nodes = cover_exists(p, size)
            if W is None and nodes > 0 and W is None and nodes >= 30_000_000:
                row += f"  size<={size}: BUDGET"; break
            if W:
                cp = size
                row += f"  c(p)={size} witness w-set={W} [{nodes:,}n {time.time()-t0:.1f}s]"
                break
            else:
                row += f"  no<={size}[{nodes:,}n]"
        print(row, flush=True)

# ---------------- seam H ----------------
def seamH():
    print("\n== SEAM H: moment blindness within (S2,S4)- and (S2,S4,S6)-fibers ==", flush=True)
    from itertools import combinations
    # degree-4 fibers: 3-subset collisions on (sum sq, sum 4th)
    UNIV = list(range(1, 46))
    buckets = {}
    for X in combinations(UNIV, 3):
        key = (sum(v*v for v in X), sum(v**4 for v in X))
        buckets.setdefault(key, []).append(X)
    pairs = []
    for key, xs in buckets.items():
        if len(xs) >= 2:
            for i in range(len(xs)):
                for j in range(i+1, len(xs)):
                    if set(xs[i]).isdisjoint(xs[j]): pairs.append((xs[i], xs[j]))
    print(f"  degree-4: {len(pairs)} disjoint colliding 3-subset pairs in [1,45]", flush=True)
    random.shuffle(pairs)
    dmax = (F(0), None); straddle = []
    n = 0
    for X, Y in pairs[:250]:
        used = set(X) | set(Y)
        C = [v for v in range(1, 100) if v not in used][:10]
        V = sorted(C + list(X)); W = sorted(C + list(Y))
        MV, _ = M_exact_arg(V); MW, _ = M_exact_arg(W)
        d = abs(MV - MW); n += 1
        if d > dmax[0]: dmax = (d, (V, MV, W, MW))
        lo, hi = min(MV, MW), max(MV, MW)
        for rung in (F(1,13), F(1,12), F(2,25), F(1,11)):
            if lo < rung <= hi: straddle.append((V, MV, W, MW, rung)); break
    print(f"  tested {n} fiber pairs: max dM = {dmax[0]} ~ {float(dmax[0]):.4f}", flush=True)
    if dmax[1]:
        V, MV, W, MW = dmax[1]
        print(f"    extreme pair: M={MV} vs M={MW}", flush=True)
        print(f"      V={V}", flush=True)
        print(f"      W={W}  (same S2, S4)", flush=True)
    print(f"  rung-straddling pairs (blind to a certificate rung): {len(straddle)}", flush=True)
    for V, MV, W, MW, rung in straddle[:3]:
        print(f"    rung {rung}: {MV} vs {MW}", flush=True)
    # degree-6 fibers: 4-subset collisions on (S2, S4, S6)
    buckets6 = {}
    for X in combinations(range(1, 41), 4):
        key = (sum(v*v for v in X), sum(v**4 for v in X), sum(v**6 for v in X))
        buckets6.setdefault(key, []).append(X)
    pairs6 = []
    for key, xs in buckets6.items():
        if len(xs) >= 2:
            for i in range(len(xs)):
                for j in range(i+1, len(xs)):
                    if set(xs[i]).isdisjoint(xs[j]): pairs6.append((xs[i], xs[j]))
    print(f"  degree-6: {len(pairs6)} disjoint colliding 4-subset pairs in [1,40]", flush=True)
    dmax6 = (F(0), None)
    for X, Y in pairs6[:120]:
        used = set(X) | set(Y)
        C = [v for v in range(1, 100) if v not in used][:9]
        V = sorted(C + list(X)); W = sorted(C + list(Y))
        MV, _ = M_exact_arg(V); MW, _ = M_exact_arg(W)
        d = abs(MV - MW)
        if d > dmax6[0]: dmax6 = (d, (V, MV, W, MW))
    print(f"  max dM within degree-6 fibers: {dmax6[0]} ~ {float(dmax6[0]):.4f}"
          f"{' (no pairs)' if dmax6[1] is None else ''}", flush=True)

# ---------------- seam J ----------------
def bestgrid(V, q):
    best = 0
    for a in range(1, q//2 + 1):
        mn = q
        for v in V:
            r = (v*a) % q
            if r > q - r: r = q - r
            if r < mn: mn = r
        if mn > best: best = mn
    return F(best, q)

def seamJ():
    print("\n== SEAM J: the CRT interference witness (mod 29, 31, 899) ==", flush=True)
    target = F(1, 14)
    for name, V in (("AP13", list(range(1,14))), ("GW", list(range(1,12))+[13,24])):
        b29, b31, b899 = bestgrid(V,29), bestgrid(V,31), bestgrid(V,899)
        print(f"  {name}: best29={b29} best31={b31} best899={b899} "
              f"(1/14={float(target):.5f}; 899-grid {'>=' if b899 >= target else '<'} 1/14)", flush=True)
    wins29 = wins899 = 0; n = 0; worst899 = (F(1), None)
    for _ in range(200):
        V = sorted(random.sample([v for v in range(1, 81) if v % 29 != 0], 13))
        if gcd(*V[:2]) > 1: pass
        n += 1
        b29 = bestgrid(V, 29); b899 = bestgrid(V, 899)
        if b29 >= target: wins29 += 1
        if b899 >= target: wins899 += 1
        if b899 < worst899[0]: worst899 = (b899, V)
    print(f"  {n} random no-mult-29 families: mod-29 grid >= 1/14 on {wins29}, "
          f"mod-899 grid >= 1/14 on {wins899}", flush=True)
    print(f"  worst mod-899 value seen: {worst899[0]} ~ {float(worst899[0]):.5f}", flush=True)

if __name__ == "__main__":
    seamG()
    seamH()
    seamJ()
