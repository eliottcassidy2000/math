#!/usr/bin/env python3
"""
lrc14_S_per_height_boxeph_S139.py  (HYP-8050)

|S| PER TOWER HEIGHT — the experiment three sessions of gate laws point at
(first-rung collapse conjecture, HYP-8025 §4), run at the largest affordable p.

OBJECT.  At threshold 1/14, "improper mod 7^a 2^b p" IS level-1 improperness of
the composite modulus Q (band floor(Q/14)) — the tower is one family I(13,Q,1/14)
over 14-smooth multiples of p.  The countable proxy per height is kind-pasteur's
c88 object generalized to composite Q, enumerated PRIMALLY (their w = v^{-1} dual
needs units; at composite Q the forced 2-/7-multiples are non-units):

  N_mc(Q) = # MINIMAL sets R of folded residues r in [1, Q/2] whose danger sets
            D(r) = {scales b <= Q/2 : dist_Q(r b) <= floor(Q/14)} cover all scales,
            |R| <= 13, minimality checked exactly at each leaf.

Improper 13-tuples = minimal cover + free slots, so N_mc and the min-size
distribution are the shape of |S| per height.  The first-rung collapse predicts:
the size-vs-13 slack TIGHTENS at the first 7-rung (Q: p -> 7p) and RELAXES in
the tower interior (7p -> 49p), increasingly so as p grows.

VALIDATION: N_mc(43)=1,280, N_mc(61)=25,711, N_mc(71)=260,568 (kind-pasteur c88,
counted in the dual; A=B=C gates make the primal count equal).  Budget-capped
enumerations report lower bounds, never fabricated exacts.

boxeph-2026-07-19-S139.  Pure Python, exact integers.
"""

import sys
from math import gcd

def dist(x, m):
    r = x % m
    return r if r <= m - r else m - r

def min_covers(Q, node_budget, count_cap=300_000, max_size=13):
    dk = Q // 14
    H = Q // 2
    danger = {}
    for r in range(1, H + 1):
        danger[r] = frozenset(b for b in range(1, H + 1) if dist(r * b, Q) <= dk)
    cand = {b: [r for r in range(1, H + 1) if b in danger[r]] for b in range(1, H + 1)}
    found = set()
    sizes = {}
    nodes = 0
    complete = True
    chosen = []

    def minimal(Rset):
        for r in Rset:
            rest = set()
            for x in Rset:
                if x != r:
                    rest |= danger[x]
            if len(rest) == H:
                return False
        return True

    def rec(uncovered, banned):
        nonlocal nodes, complete
        if nodes > node_budget or len(found) >= count_cap:
            complete = False
            return
        nodes += 1
        if not uncovered:
            R = frozenset(chosen)
            if R not in found and minimal(R):
                found.add(R)
                sizes[len(R)] = sizes.get(len(R), 0) + 1
            return
        if len(chosen) >= max_size:
            return
        it = iter(uncovered)
        us = [next(it) for _ in range(min(20, len(uncovered)))]
        b = min(us, key=lambda x: sum(1 for r in cand[x] if r not in banned))
        opts = [r for r in cand[b] if r not in banned]
        # prune: enough residues must remain to cover
        for r in opts:
            chosen.append(r)
            rec(uncovered - danger[r], banned)
            chosen.pop()
            banned = banned | {r}
    rec(frozenset(range(1, H + 1)), frozenset())
    return len(found), sizes, nodes, complete

def run(pairs):
    print("%-8s %-6s %-12s %-9s %-30s %s" %
          ("Q", "p·h", "N_mc", "complete", "size histogram", "nodes"))
    for (p, h, Q, budget) in pairs:
        n, sizes, nodes, comp = min_covers(Q, budget)
        hist = " ".join("%d:%d" % (k, sizes[k]) for k in sorted(sizes))
        print("%-8d %-6s %-12s %-9s %-30s %d" %
              (Q, "%dx%s" % (p, h), ("%d" % n) if comp else (">=%d" % n),
               comp, hist, nodes))
        results.append((p, h, Q, n, comp, sizes))

results = []
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "a"
    if mode == "a":
        # p=43 full tower + validation anchors 61, 71 at height 1
        run([(43, "1", 43, 2_000_000),
             (61, "1", 61, 4_000_000),
             (43, "2", 86, 4_000_000),
             (43, "7", 301, 6_000_000)])
    elif mode == "b":
        run([(71, "1", 71, 14_000_000),
             (43, "14", 602, 8_000_000)])
    elif mode == "c":
        run([(61, "2", 122, 6_000_000),
             (61, "7", 427, 8_000_000)])
    print("\nREADING KEY: 'p·h' = prime x height-multiplier. Watch (i) N_mc growth vs Q,")
    print("(ii) the size histogram's shift toward 13 (slack = 13 - minsize) at the first")
    print("7-rung vs the 2-rung, (iii) whether slack relaxes at greater heights. The")
    print("first-rung collapse conjecture (HYP-8025) predicts the slack squeeze at p->7p")
    print("sharpens as p grows.")
