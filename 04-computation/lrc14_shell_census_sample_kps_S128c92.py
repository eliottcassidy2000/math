#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c92 -- HYP-7995 Lane A (v2, SAMPLING census).

The direct enumeration of minimal covers is infeasible at p ~ 300 (40M DFS
nodes reach zero complete covers -- the c91 finding again: exhaustive order
cannot find).  v2 samples the minimal-cover space instead:
  sample = randomized-greedy cover (top-t choice) -> DROP-ONE MINIMALIZATION
           (remove redundant sets in random order until minimal)
  abundance = capture-recapture on distinct samples (collisions c among k
           draws => N-hat ~ k^2/(2c); zero collisions => N >> k^2/2 lower bd)
  structure = exact-ansatz hits + NEAR-SHELL Hamming distance (min
           substitutions of the speed-set to any ansatz dilate) + size dist.
Collapse verdict: few distinct covers, repeated collisions, small near-shell
distances = SHELL COLLAPSE; endless distinct spread covers = generic sea.
Sample min size also upper-bounds c(p).
"""
import sys, time, random
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build
from lrc14_I13p1_minimal_covers_kps_S128c88 import ansatz_library

random.seed(92)

def sample_minimal(p, ctx, topt):
    h, dk, maskA, cand, FULL, fold, inv = ctx
    ws = list(range(1, h+1))
    U = FULL; chosen = []
    for _ in range(13):
        best = []
        for w in ws:
            c = (maskA[w] & U).bit_count()
            if c: best.append((c, w))
        if not best: return None
        best.sort(reverse=True)
        c, w = random.choice(best[:topt])
        chosen.append(w); U &= ~maskA[w]
        if U == 0: break
    if U != 0: return None
    # minimalize: drop redundant sets in random order
    random.shuffle(chosen)
    i = 0
    while i < len(chosen):
        rest = 0
        for j, u in enumerate(chosen):
            if j != i: rest |= maskA[u]
        if rest == FULL: chosen.pop(i)
        else: i += 1
    return frozenset(chosen)

def main():
    HARD = {349, 457, 467, 683, 719, 797, 809, 823, 839, 907, 983, 1103, 1123, 1153, 1163}
    NS = 800
    print(f"== Lane A v2: SAMPLING shell census ({NS} greedy draws/prime) ==", flush=True)
    for p in (307, 331, 349, 397, 457, 467):
        t0 = time.time()
        ctx = build(p)
        h, dk, fold, inv = ctx[0], ctx[1], ctx[5], ctx[6]
        lib = ansatz_library(p, ctx)
        libsets = list(lib.keys())
        seen = {}
        fails = 0
        for _ in range(NS):
            W = sample_minimal(p, ctx, topt=random.choice((2, 3, 4)))
            if W is None: fails += 1; continue
            seen[W] = seen.get(W, 0) + 1
        k = sum(seen.values()); distinct = len(seen)
        coll = k - distinct
        est = (k*k)/(2*coll) if coll else None
        sizes = {}; hits = {}; near = {}
        for W in seen:
            sizes[len(W)] = sizes.get(len(W), 0) + 1
            V = frozenset(fold(inv[w]) for w in W)
            tag = lib.get(V)
            if tag: hits[tag] = hits.get(tag, 0) + 1
            d = min(len(V.symmetric_difference(F))//2 for F in libsets)
            near[d] = near.get(d, 0) + 1
        tagp = "thin-sea-HARD" if p in HARD else "greedy-easy"
        print(f"  p={p} ({tagp}, dk={dk}, h={h}):", flush=True)
        print(f"    draws ok={k} fail={fails}; DISTINCT minimal covers = {distinct}; "
              f"collisions = {coll}; abundance {'N~%.0f' % est if est else 'N >> %d (no collisions)' % (k*k//2)}", flush=True)
        print(f"    sizes {dict(sorted(sizes.items()))}  (min size = upper bound on c(p))", flush=True)
        print(f"    exact-ansatz hits: {hits or 'NONE'}", flush=True)
        print(f"    near-shell Hamming histogram: {dict(sorted(near.items()))}  [{time.time()-t0:.0f}s]", flush=True)

if __name__ == "__main__":
    main()
