#!/usr/bin/env python3
"""
lrc_nerve_factor_s521.py   claudebox-2026-06-01-S521

Two more complex-vertex encodings of LRC (reflection:
07-reflections/lrc-vertex-taxonomy-capstone-s521.md).

(1) NERVE / WINDING (topological combinatorics): the sum(v_i) danger sub-arcs
    (runner i, cluster k: center k/v_i, radius 1/(n v_i)) cover S^1  <=>  NOT
    strictly lonely <=> strict-LRC fails. Covering of a circle by arcs <=> their
    overlap structure winds once around. Vertices = sub-arcs; edges = overlaps.
(2) FACTOR COMPLEXITY (symbolic dynamics): the crossing word (which runner crosses
    a sector boundary, in time order) is the cutting sequence of the torus line;
    its factor complexity p(k) measures orbit complexity (bounded/periodic for
    rational lines).
"""
from fractions import Fraction as F

def subarcs(sp, n):
    A = []
    for v in sp:
        for k in range(v):
            c = F(k, v); r = F(1, n*v)
            A.append(((c-r) % 1, (c+r) % 1))
    return A
def covers_circle(sp, n):
    A = subarcs(sp, n)
    pts = sorted(set(a for arc in A for a in arc))
    aug = pts + [pts[0] + 1]
    for x, y in zip(aug, aug[1:]):
        mid = ((x + y) / 2) % 1
        if not any((lo <= mid < hi) if lo < hi else (mid >= lo or mid < hi) for lo, hi in A):
            return False
    return True
def crossing_word(sp, n):
    ev = []
    for idx, v in enumerate(sp):
        for k in range(n*v): ev.append((F(k, n*v), idx))
    ev.sort()
    return [idx for _, idx in ev]
def factor_complexity(word, k):
    L = len(word); fac = set()
    for i in range(L): fac.add(tuple(word[(i+j) % L] for j in range(k)))
    return len(fac)

def main():
    print("(1) Nerve/winding: do danger sub-arcs cover S^1? (cover <=> strict-LRC fails)")
    for sp in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11)]:
        n = len(sp)+1
        cov = covers_circle(list(sp), n)
        print(f"   {sp}: {sum(sp)} sub-arcs, cover circle? {cov}  -> "
              f"{'strictly non-lonely (tight/boundary)' if cov else 'LONELY (uncovered gap)'}")
    print("\n(2) Factor complexity p(k) of the crossing word (symbolic dynamics):")
    for sp in [(1,2,3),(1,2,4,7)]:
        n = len(sp)+1; W = crossing_word(list(sp), n)
        ps = [factor_complexity(W, k) for k in range(1, 7)]
        print(f"   {sp} (word len {len(W)}): p(1..6)={ps}  (plateaus: rational line => periodic)")

if __name__ == "__main__":
    main()
