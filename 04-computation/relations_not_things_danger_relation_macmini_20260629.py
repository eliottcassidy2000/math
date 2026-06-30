#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relations, not things: the LRC danger relation and the metagraph dominance relation as
the PRIMARY objects; the only invariant is the relation's self-composition (the 2nd moment);
a coboundary = a SEPARABLE (rank-1) relation. (mac-mini-2026-06-29-S24)

Extends the coboundary lens (HYP-3562) + klein THM-588 (no thing-invariant, only the quadratic
3-cycle RELATION) + klein-S5 (the reference-collapse relation, S_n clean vs Z_14 dirty).

THESIS (relations not things):
 - a TOURNAMENT is the dominance RELATION (a->b); its only S_n-invariant is the relational content
   (cyclicity / 3-cycle, THM-588); the transitive (a linear ORDER, a 'thing') is the trivial point.
 - the LRC danger is the RELATION D(v,t) = [||v t|| < 1/14], essentially BILINEAR (the product v*t),
   NOT separable.
 - the SECOND MOMENT = the relation composed with its transpose (D D^T = the sheet/pair correlation,
   THM-589/579) -- the ONLY invariant (THM-588).
 - ESSENTIAL = the relation has rank > 1 (not separable); COBOUNDARY = rank-1 = f(v) g(t) (trivializable
   by a 'thing'/potential). Disproof = the relation is a coboundary; proof = essential.
"""
from __future__ import annotations
import functools, math, itertools
print = functools.partial(print, flush=True)
try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False


def frac(x): return x - math.floor(x)
def dist(x): f = frac(x); return min(f, 1 - f)


def danger_relation(speeds, n, T):
    """D[i][a] = 1 if speed v_i is 'safe' (||v_i * a/T|| >= 1/n) at time a/T, else 0 (DANGER).
    We use the SAFE relation (complement of danger) so a lonely time = a column of all-1s."""
    return [[1 if dist(v * a / T) >= 1/n - 1e-12 else 0 for a in range(T)] for v in speeds]


def main():
    print("=" * 80)
    print("Relations, not things: the danger relation as the primary object (mac-mini-S24)")
    print("=" * 80)

    # ---- (1) the LRC SAFE relation, its rank (essentiality) and Gram (the 2nd moment) ----
    print("\n[1] LRC SAFE relation S[i][a] = [speed v_i safe at time a/T] (a lonely time = an all-safe")
    print("    column). The relation is the OBJECT; a lonely point is DERIVED (a full column).")
    for speeds, n, name in [([1, 2, 3], 4, "{1,2,3} LRC4 extremal"),
                            ([2, 3, 4], 4, "{2,3,4} LRC4 covering"),
                            ([1, 2, 3, 4, 5], 6, "{1..5} LRC6 extremal")]:
        T = 12 * n   # fine enough lattice of times a/T
        R = danger_relation(speeds, n, T)
        lonely_cols = sum(1 for a in range(T) if all(R[i][a] for i in range(len(speeds))))
        if HAVE_NP:
            M = np.array(R, dtype=float)
            rank = int(np.linalg.matrix_rank(M, tol=1e-9))
        else:
            rank = "(numpy NA)"
        print(f"    {name:24s}: relation {len(speeds)}x{T}, rank={rank} (>1 => ESSENTIAL, not separable/"
              f"coboundary); lonely (all-safe) columns: {lonely_cols}")
    print("    => the SAFE relation has rank > 1 (essential): it is NOT a coboundary f(v)g(t).")
    print("    A disproof would need the relation to TRIVIALIZE (separable => covers); the bilinear")
    print("    product v*t in ||v t|| forbids it (the multiplicative structure = the essentiality).")

    # ---- (2) the 2nd moment = the relation composed with its transpose (the ONLY invariant) ----
    print("\n[2] The SECOND MOMENT = the relation's SELF-COMPOSITION (Gram), the pair-correlation:")
    speeds, n = [2, 3, 4], 4; T = 12 * n
    R = danger_relation(speeds, n, T)
    if HAVE_NP:
        M = np.array(R, dtype=float)
        gram_time = M.T @ M       # T x T: G[a,b] = #speeds safe at BOTH a and b (sheet-pair overlap)
        gram_speed = M @ M.T      # k x k: G[i,j] = #times both v_i,v_j safe (speed-pair overlap)
        print(f"    {speeds}: D D^T (speed-pair overlap, {len(speeds)}x{len(speeds)}):")
        for row in gram_speed.astype(int): print(f"      {list(row)}")
        print(f"    diag = #safe-times per speed (the 1st moment); off-diag = the PAIR RELATION (2nd")
        print(f"    moment). klein THM-588: the off-diag (pair) is the ONLY invariant -- no 1st-order.")
    else:
        print("    (numpy unavailable; the Gram D D^T is the pair-overlap = the 2nd moment.)")

    # ---- (3) the metagraph dominance relation: only the 3-cycle (relation) is invariant ----
    print("\n[3] Metagraph: a tournament IS the dominance RELATION; only the 3-cycle (a RELATION among")
    print("    3 vertices = failure of transitivity) is S_n-invariant (THM-588, mult(1)=0,mult(2)=1).")
    print("    The transitive tournament = a linear ORDER = the trivial 'thing' (0 cyclicity).")
    for n in (4, 5, 6):
        # count 3-cycles per iso class is the unique level-2 invariant; show range over iso classes
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        perms = list(itertools.permutations(range(n)))
        c3vals = set()
        seen = set()
        for bits in range(1 << len(pairs)):
            arc = [[False]*n for _ in range(n)]
            for b,(i,j) in enumerate(pairs):
                if (bits>>b)&1: arc[i][j]=True
                else: arc[j][i]=True
            canon = min(tuple(1 if arc[s[i]][s[j]] else 0 for i in range(n) for j in range(n) if i!=j) for s in perms)
            if canon in seen: continue
            seen.add(canon)
            c3 = sum(1 for a in range(n) for b2 in range(n) for c in range(n)
                     if a<b2<c and ((arc[a][b2] and arc[b2][c] and arc[c][a]) or (arc[b2][a] and arc[c][b2] and arc[a][c])))
            c3vals.add(c3)
        print(f"    n={n}: 3-cycle counts over iso classes = {sorted(c3vals)} (0 = transitive/thing;")
        print(f"          >0 = relational content). This single relation-invariant carries the metagraph.")

    print("\n" + "=" * 80)
    print("RELATIONS NOT THINGS: the object is the RELATION (dominance / danger), the only invariant")
    print("is its SELF-COMPOSITION (the 2nd-moment pair-correlation, THM-588's quadratic), a lonely")
    print("point / SC tournament is DERIVED, and a disproof = the relation being a COBOUNDARY")
    print("(separable / rank-1 / trivializable by a thing). The bilinear v*t (LRC) and the cyclic")
    print("3-cycle (metagraph) are the essentialities that forbid trivialization.")
    print("=" * 80)


if __name__ == "__main__":
    main()
