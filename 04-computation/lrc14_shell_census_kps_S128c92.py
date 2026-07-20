#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c92 -- HYP-7995 Lane A: the SHELL-COLLAPSE CENSUS
at p ~ 300-500.

Question (the c91 reframe): does I(13,p,1) collapse to the tight-family shell
as p grows -- minimal covers few, near the ansatz -- reversing the S128c88
generic sea (1,280 / 25,711 / 260,568 minimal covers at p = 43/61/71, 99%+
unclassified)?  Census primes mix greedy-easy (307, 331, 397) and thin-sea-
hard (349, 457, 467) to test the thin-sea ranking against actual cover counts.

Per prime: minimal covers (soft node budget; truncation REPORTED), size
distribution, exact-ansatz hits, NEAR-SHELL Hamming distance histogram
(min substitutions to any ansatz dilate), and p-1 factorization.

Lane B / S6 rider: across all 130 triage primes (c91), correlate greedy
difficulty with the largest prime factor of p-1 (safe primes p = 2q+1 have
maximally non-smooth multiplicative structure -- do they cluster greedy-hard?).
"""
import sys, time
sys.path.insert(0, '04-computation')
from lrc14_I13p1_acceptance_test_kps_S128c88 import build
from lrc14_I13p1_minimal_covers_kps_S128c88 import ansatz_library

def minimal_covers_soft(p, max_size=13, budget=40_000_000):
    h, dk, maskA, cand, FULL, fold, inv = build(p)
    found = set(); nodes = 0; truncated = False
    usedset = set()
    def rec(U, banned):
        nonlocal nodes, truncated
        nodes += 1
        if nodes > budget: truncated = True; return
        if U == 0:
            found.add(frozenset(usedset)); return
        if len(usedset) >= max_size: return
        rem = max_size - len(usedset)
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
            if truncated: return
            usedset.add(w)
            rec(U & ~maskA[w], banned | newban)
            usedset.discard(w)
            newban.add(w)
    usedset.add(1)
    rec(FULL & ~maskA[1], frozenset())
    usedset.discard(1)
    mins = []
    for Wc in found:
        ok = True
        for w in Wc:
            rest = 0
            for u in Wc:
                if u != w: rest |= maskA[u]
            if rest == FULL: ok = False; break
        if ok: mins.append(sorted(Wc))
    return mins, nodes, truncated, (h, dk, maskA, cand, FULL, fold, inv)

def largest_pf(n):
    f = 1; d = 2
    while d*d <= n:
        while n % d == 0: f = max(f, d); n //= d
        d += 1
    return max(f, n) if n > 1 else f

def main():
    print("== Lane A: shell-collapse census ==", flush=True)
    HARD = {349, 457, 467, 683, 719, 797, 809, 823, 839, 907, 983, 1103, 1123, 1153, 1163}
    for p in (307, 331, 349, 397, 457, 467):
        t0 = time.time()
        mins, nodes, trunc, ctx = minimal_covers_soft(p)
        h, dk, fold, inv = ctx[0], ctx[1], ctx[5], ctx[6]
        lib = ansatz_library(p, ctx)
        libsets = list(lib.keys())
        sizes = {}; hits = {}; near = {}
        for Wc in mins: sizes[len(Wc)] = sizes.get(len(Wc), 0) + 1
        for Wc in mins:
            V = frozenset(fold(inv[w]) for w in Wc)
            tag = lib.get(V)
            if tag: hits[tag] = hits.get(tag, 0) + 1
            d = min(len(V.symmetric_difference(F))//2 for F in libsets)
            near[d] = near.get(d, 0) + 1
        tagp = "thin-sea-HARD" if p in HARD else "greedy-easy"
        print(f"  p={p} ({tagp}, dk={dk}, h={h}, largest-pf(p-1)={largest_pf(p-1)}):", flush=True)
        print(f"    minimal covers = {len(mins)}{' (TRUNCATED at budget -- lower bound)' if trunc else ' (complete)'}"
              f"  nodes={nodes:,}  [{time.time()-t0:.0f}s]", flush=True)
        print(f"    sizes {dict(sorted(sizes.items()))}  exact-ansatz hits {hits or 'NONE'}", flush=True)
        print(f"    near-shell Hamming histogram {dict(sorted(near.items()))}", flush=True)

    print("\n== Lane B rider (S6): greedy difficulty vs multiplicative structure of p-1 ==", flush=True)
    import re, io
    tri = io.open('05-knowledge/results/lrc14_extinction_triage_kps_S128c91.out', encoding='utf-8').read()
    easy = [int(m) for m in re.findall(r"p=(\d+) ALIVE \(greedy cover", tri)]
    hard = sorted(HARD)
    def stats(ps):
        lps = sorted(largest_pf(p-1) / (p-1) for p in ps)
        safe = sum(1 for p in ps if largest_pf(p-1) == (p-1)//2)
        return lps[len(lps)//2], safe, len(ps)
    me, se, ne = stats([p for p in easy if p not in HARD])
    mh, sh, nh = stats(hard)
    print(f"  greedy-easy (n={ne}): median largest-pf(p-1)/(p-1) = {me:.3f}, safe primes = {se} ({100*se/ne:.0f}%)", flush=True)
    print(f"  thin-sea-hard (n={nh}): median = {mh:.3f}, safe primes = {sh} ({100*sh/nh:.0f}%)", flush=True)
    print(f"  S6 verdict: {'CORRELATED (hard primes are multiplicatively unsmooth)' if mh > me and sh/nh > se/ne else 'NOT clearly correlated'}", flush=True)

if __name__ == "__main__":
    main()
