#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c88 -- HYP-7930 part 2: MINIMAL covers of the
level-1 sieve and the ANSATZ-COMPLETENESS question.

S-T name the k=13 bottleneck as "a better understanding of speed tuples that
do not have a witness time in an ansatz."  The level-1 improper structures
are exactly the MINIMAL covering w-sets (size <= 13) of the folded line
(see lrc14_I13p1_acceptance_test_kps_S128c88.py for the formulation).  Here
we enumerate them exactly at completed primes and classify each against the
repo's known-family ansatz:

  AP    c*{1..13}                     (the tight AP)
  GW    c*{1..11,13,24}              (Goddyn-Wong second tight family)
  K2    c*{1..12,26}                 (2/27 slack-1 D=2, S128c86)
  F3    c*{1..11,13,36}              (3/41 mediant)
  K3    c*{1..12,39}, F4 c*{1..11,13,48}, K-ladder c*{1..12, 13m} m<=6,
  F-ladder c*{1..11,13,12m} m<=6, prefix APs c*{1..j} j<=13 (sub-AP cores)

A minimal cover whose speed-set matches NO ansatz entry is an "unclassified
improper structure" -- the objects S-T say need understanding.  Their count
and shape as p grows is the deliverable.

Method: MRV-dual branching (every minimal cover appears as a complete path:
at each pivot, a minimal cover supplies a candidate, and no proper prefix of
a minimal cover can already cover).  Collect path-sets at U==0, dedupe,
filter to minimal (drop-one test), classify.  Sizes < 13 allowed (more
improper freedom).
"""
import sys, time
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build

def minimal_covers(p, max_size=13, budget=None):
    h, dk, maskA, cand, FULL, fold, inv = build(p)
    found = set()
    nodes = 0
    usedset = set()
    def rec(U, banned):
        nonlocal nodes
        nodes += 1
        if budget and nodes > budget: raise TimeoutError
        if U == 0:
            found.add(frozenset(usedset))
            return
        if len(usedset) >= max_size: return
        rem = max_size - len(usedset)
        if U.bit_count() > dk*rem: return
        # pivot: least-candidates uncovered point (sample lowest 24)
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
            usedset.add(w)
            rec(U & ~maskA[w], banned | newban)
            usedset.discard(w)
            newban.add(w)
    # NOTE: unlike the acceptance test we do NOT fix 1 in W here -- minimal
    # covers are classified up to scaling at the CLASSIFY stage instead, so
    # we enumerate all classes with 1 included WLOG (scaling maps any cover
    # to one containing 1) -- same normalization, stated for clarity.
    usedset.add(1)
    rec(FULL & ~maskA[1], frozenset())
    usedset.discard(1)
    # minimality filter
    mins = []
    for W in found:
        m = 0
        ok = True
        for w in W:
            rest = 0
            for u in W:
                if u != w: rest |= maskA[u]
            if rest == FULL: ok = False; break
        if ok: mins.append(sorted(W))
    return mins, nodes, (h, dk, maskA, cand, FULL, fold, inv)

def ansatz_library(p, ctx):
    """Precompute {frozenset(folded speed-set) : family tag} for all dilates."""
    h, dk, maskA, cand, FULL, fold, inv = ctx
    base_families = [("AP", list(range(1, 14))), ("GW", [1,2,3,4,5,6,7,8,9,10,11,13,24])]
    for m in range(2, 7):
        base_families.append((f"K{m}", list(range(1, 13)) + [13*m]))
        base_families.append((f"F{m}", list(range(1, 12)) + [13, 12*m]))
    for j in range(2, 14):
        base_families.append((f"pAP{j}", list(range(1, j+1))))
    lib = {}
    for c in range(1, h+1):
        for name, fam in base_families:
            FV = frozenset(fold(c*v % p) for v in fam if c*v % p != 0)
            lib.setdefault(FV, name)
    return lib

def classify(p, W, ctx, lib):
    h, dk, maskA, cand, FULL, fold, inv = ctx
    V = frozenset(fold(inv[w]) for w in W)       # the speed-set (folded)
    return lib.get(V)

def main():
    for p in (43, 61, 71):
        t0 = time.time()
        try:
            mins, nodes, ctx = minimal_covers(p, budget=60_000_000)
        except TimeoutError:
            print(f"p={p}: BUDGET EXCEEDED, no conclusion", flush=True); continue
        lib = ansatz_library(p, ctx)
        sizes = {}
        for W in mins: sizes[len(W)] = sizes.get(len(W), 0) + 1
        uncls = []
        tags = {}
        for W in mins:
            tag = classify(p, W, ctx, lib)
            if tag is None: uncls.append(W)
            else: tags[tag] = tags.get(tag, 0) + 1
        h, dk = ctx[0], ctx[1]
        print(f"p={p} (h={h}, dk={dk}): minimal covers = {len(mins)}, "
              f"sizes {sizes}, nodes={nodes:,}, {time.time()-t0:.1f}s", flush=True)
        print(f"   classified: {tags}", flush=True)
        print(f"   UNCLASSIFIED: {len(uncls)}", flush=True)
        for W in uncls[:12]:
            V = sorted(ctx[5](ctx[6][w]) for w in W)
            print(f"     w-set {W}  ->  speeds(folded) {V}", flush=True)
        if len(uncls) > 12: print(f"     ... and {len(uncls)-12} more", flush=True)

if __name__ == "__main__":
    main()
