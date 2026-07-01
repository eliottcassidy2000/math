#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S81 -- the CONFIGURATION RULE (explicit fixing) for kappa(n), and an n=7 upper-bound search.

Part 1: for n=4 and n=5, print an explicit optimal config: the free arcs (shape) + the fixed arcs and THEIR
        orientations (the 'configuration rule'), and the resulting 2^kappa tournaments -> all iso classes.
Part 2: n=7 -- randomized/greedy search for a covering config at k=11,12 (upper bound on kappa(7)).
"""
import itertools, random
from math import comb
from collections import Counter

A000568 = {3:2, 4:4, 5:12, 6:56, 7:456}

def setup(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    def relabel(mask, perm):
        out = 0
        for (i, j), k in idx.items():
            bit = (mask >> k) & 1
            pi, pj = perm[i], perm[j]
            a, b = (pi, pj) if pi < pj else (pj, pi)
            out |= (bit if pi < pj else 1 - bit) << idx[(a, b)]
        return out
    return pairs, idx, len(pairs), perms, relabel

def canon(mask, perms, relabel):
    return min(relabel(mask, p) for p in perms)

# ---------- Part 1: explicit configuration rule for n=4,5 ----------
def find_config(n):
    pairs, idx, P, perms, relabel = setup(n)
    canon_of = {}
    for mask in range(1 << P):
        canon_of[mask] = canon(mask, perms, relabel)
    nclasses = len(set(canon_of.values())); assert nclasses == A000568[n]
    k = 1 + comb(n-2, 2)
    for F in itertools.combinations(range(P), k):
        Fbits = list(F); rest = [b for b in range(P) if b not in F]
        subs = []
        for s in range(1 << k):
            mm = 0
            for t in range(k):
                if (s >> t) & 1: mm |= (1 << Fbits[t])
            subs.append(mm)
        for fixs in range(1 << len(rest)):
            fixmask = 0
            for t, b in enumerate(rest):
                if (fixs >> t) & 1: fixmask |= (1 << b)
            seen = set(canon_of[mm | fixmask] for mm in subs)
            if len(seen) == nclasses:
                return pairs, F, rest, fixmask, k
    return None

for n in [4, 5]:
    pairs, F, rest, fixmask, k = find_config(n)
    free_arcs = [pairs[b] for b in F]
    fixed = [(pairs[b], (fixmask >> b) & 1) for b in rest]   # (arc, 1 => i->j)
    print(f"=== n={n}: kappa={k} config ===")
    print(f"  FREE arcs (the subcube, shape): {free_arcs}")
    print(f"  FIXED arcs + orientation (1 = i->j, 0 = j->i)  [the CONFIGURATION RULE]:")
    print(f"     {fixed}")
    # show the 2^k tournaments hit distinct classes
    print(f"  => the {1<<k} subcube tournaments realize all {A000568[n]} iso classes.\n")

# ---------- Part 2: n=7 randomized upper-bound search ----------
print("=== n=7: randomized search for a covering config (upper bound on kappa(7); formula predicts 11) ===")
pairs, idx, P, perms, relabel = setup(7)
def coverage(Fbits, fixmask):
    seen = set()
    for s in range(1 << len(Fbits)):
        mm = fixmask
        for t in range(len(Fbits)):
            if (s >> t) & 1: mm |= (1 << Fbits[t])
        seen.add(canon(mm, perms, relabel))
        if len(seen) == 456: return 456
    return len(seen)

random.seed(7)
best = {}
for k in [11, 12, 13]:
    bestcov = 0; besttup = None; tries = 0
    # try: random free-arc set + a few random fixings, keep best
    import time
    t0 = None
    for attempt in range(18):
        F = random.sample(range(P), k); Fbits = list(F); rest = [b for b in range(P) if b not in F]
        for _ in range(3):
            fixmask = 0
            for b in rest:
                if random.random() < 0.5: fixmask |= (1 << b)
            cov = coverage(Fbits, fixmask); tries += 1
            if cov > bestcov: bestcov = cov; besttup = ([pairs[b] for b in F], )
            if cov == 456: break
        if bestcov == 456: break
    best[k] = bestcov
    print(f"  k={k}: best coverage found = {bestcov}/456  ({tries} configs tried)  {'COVERS ALL' if bestcov==456 else ''}")
print(f"\n  info-floor kappa(7) >= ceil(log2 456) = 9;  formula predicts kappa(7)=11.")
print(f"  (search gives an UPPER-bound witness only if it reaches 456; else inconclusive.)")
