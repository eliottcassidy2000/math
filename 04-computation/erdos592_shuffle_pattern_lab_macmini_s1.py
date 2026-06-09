#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erdős Problem 592 lab — the shuffle-pattern calculus for ordinal partition relations.
mac-mini-2026-06-09-S1  (T768, HYP-2344..2346, THM-453)

PROBLEM (Erdős, $1000, OPEN): characterise countable ordinals alpha with
    alpha -> (alpha, 3)^2
i.e. every red/blue colouring of [alpha]^2 has a red set of full type alpha or a blue triangle.

REFRAME (THM-453 part 1): alpha -> (alpha,3)^2  <=>  every TRIANGLE-FREE graph on alpha
has an independent set of order type alpha.  Negative witness = triangle-free graph
killing all full-type independent sets.

MODEL: omega^n = W(n) = strictly increasing n-tuples from omega, ordered lexicographically
(Larson's model; order type omega^n).

PATTERN of a pair {x,y} of increasing n-tuples: the order/equality pattern of the 2n
values (rank-collapsed).  Every known witness colouring is PATTERN-MEASURABLE (HYP-2345).

A pattern set S (blue patterns) defines graph G_S on W(n).
  * Triangle-freeness of G_S = a FINITE check (no triple of n-tuples with all three
    pairwise patterns in S; triples need <= 3n values, rank-collapse exhausts all cases).
  * "Unavoidable" = every subset of type omega^n contains an S-pair.
    FACT (proved below in comments, used as the miniature): X subseteq W(n) has type
    omega^n iff X contains a full grid = an omega-branching tree of height n with values
    strictly increasing along branches (leaves = the tuples).  So:
       S avoidable <=> exists an S-free full (infinite) grid.
    FINITE SHADOW: max branching t such that an S-free t-grid exists in value-universe [V].
    If max t stays bounded as V grows -> S behaves unavoidable (negative-witness shaped);
    if t grows -> avoidable.  (Finite side exact; infinite extrapolation labelled as such.)

KNOWN ANSWERS THE LAB MUST REPRODUCE (literature; Specker):
  n=2: omega^2 -> (omega^2, 3): EVERY triangle-free S must be avoidable.
  n=3: omega^3 -/-> (omega^3, 3): SOME triangle-free S is unavoidable.

Part 5 records the tile-pair dictionary at n=2 (HYP-2344): patterns <-> the repo's
staircase tile-pair classes.

All sections print exhaustive counts; nothing is sampled unless labelled SEARCH.
"""

import itertools, sys
from functools import lru_cache

# ----------------------------------------------------------------------------------
# Pattern machinery
# ----------------------------------------------------------------------------------

def pattern(x, y):
    """Canonical pattern of unordered pair {x,y} of increasing tuples, x != y.
    Returns (px, py) with px the lex-smaller tuple's rank-tuple."""
    a, b = (x, y) if x < y else (y, x)
    vals = sorted(set(a) | set(b))
    rk = {v: i for i, v in enumerate(vals)}
    return (tuple(rk[v] for v in a), tuple(rk[v] for v in b))

def all_patterns(n):
    """All patterns of pairs of increasing n-tuples (values may be shared)."""
    pats = set()
    rng = range(2 * n)
    for x in itertools.combinations(rng, n):
        for y in itertools.combinations(rng, n):
            if x != y:
                pats.add(pattern(x, y))
    return sorted(pats)

def strict_patterns(n):
    """Patterns with all 2n values distinct (the 'spread' ones)."""
    return [p for p in all_patterns(n) if len(set(p[0]) | set(p[1])) == 2 * n]

def triangle_table(n):
    """All triples of pairwise patterns realizable by x <lex y <lex z increasing n-tuples.
    Exhaustive: any triple uses <= 3n distinct values; enumerate over [0, 3n)."""
    triples = set()
    rng = range(3 * n)
    tuples = list(itertools.combinations(rng, n))
    for x in tuples:
        for y in tuples:
            if y <= x:
                continue
            pxy = pattern(x, y)
            for z in tuples:
                if z <= y:
                    continue
                triples.add((pxy, pattern(y, z), pattern(x, z)))
    return triples

def is_triangle_free(S, triples):
    Sset = set(S)
    for (p, q, r) in triples:
        if p in Sset and q in Sset and r in Sset:
            return False
    return True

def maximal_triangle_free_sets(pats, triples):
    """All maximal triangle-free subsets of pats (exact, via recursive extension)."""
    pats = list(pats)
    results = set()

    def extend(S, candidates):
        extended = False
        for i, p in enumerate(candidates):
            S2 = S | {p}
            if is_triangle_free(S2, triples):
                extended = True
                extend(S2, candidates[i + 1:])
        if not extended:
            # maximal w.r.t. candidates ordering — verify true maximality
            if all(not is_triangle_free(S | {q}, triples) for q in pats if q not in S):
                results.add(frozenset(S))

    extend(set(), pats)
    return sorted(results, key=lambda s: (-len(s), sorted(s)))

# ----------------------------------------------------------------------------------
# Grids.  A t-grid of height n in universe [V]: tree, branching t at every level,
# values strictly increasing along each branch, children of a node have t distinct
# values (> parent value).  Leaves = increasing n-tuples.  X(grid) = set of leaves.
#
# WHY GRIDS: X ⊆ W(n) has order type omega^n  <=>  X ⊇ leaves of a full
# (omega-branching) grid.  Proof sketch (induction on n): otp(X) = Σ_{a} otp(X_a)
# over first coordinates a; each otp(X_a) ≤ omega^{n-1}; the ordered sum is omega^n
# iff infinitely many sections have full type omega^{n-1}.  [Used for the miniature.]
# ----------------------------------------------------------------------------------

def exists_S_free_grid(n, t, V, S, all_pairs_blue=False, return_witness=False):
    """Exact backtracking: does an S-free t-grid of height n exist with values < V?
    If all_pairs_blue: require every leaf-pair pattern IN S instead (for sanity checks).
    Leaves are built in DFS order; each completed leaf is checked against all previous."""
    Sset = set(S)
    # nodes of the tree in DFS leaf order: we realize branches one at a time.
    # state: stack of (level, value) along current path + list of completed leaves
    leaves = []
    # children values chosen per internal node, to enforce branching t with distinct vals:
    # We generate the tree recursively.

    sys.setrecursionlimit(100000)

    def ok(leaf):
        for prev in leaves:
            if prev == leaf:
                return False
            p = pattern(prev, leaf)
            inS = p in Sset
            if all_pairs_blue and not inS:
                return False
            if (not all_pairs_blue) and inS:
                return False
        return True

    def build(prefix, lo):
        """Complete a subtree under 'prefix' (tuple of values so far), choosing t children
        with values in (lo, V).  Returns True if subtree completed consistently."""
        level = len(prefix)
        if level == n:
            if ok(prefix):
                leaves.append(prefix)
                return True
            return False

        def choose(children, start):
            if len(children) == t:
                return True
            for v in range(start, V):
                children.append(v)
                if build(prefix + (v,), v):
                    if choose(children, v + 1):
                        return True
                    # undo this child's subtree leaves
                    sub = t ** (n - level - 1)
                    del leaves[-sub:]
                children.pop()
            return False

        return choose([], lo + 1)

    found = build((), -1)
    if return_witness:
        return (found, list(leaves))
    return found

def max_grid_branching(n, V, S, tmax=8):
    """Largest t <= tmax with an S-free t-grid in [V] (exact)."""
    best = 0
    for t in range(1, tmax + 1):
        if t ** n > V * 4 and t > 2:  # hopeless quick cut (need t distinct values/level anyway)
            pass
        if exists_S_free_grid(n, t, V, S):
            best = t
        else:
            break
    return best

# ----------------------------------------------------------------------------------
# Pretty names for n=2 patterns (the tile dictionary, HYP-2344)
# ----------------------------------------------------------------------------------

NAME2 = {
    ((0, 1), (2, 3)): "DISJOINT   (Larson form 0)",
    ((0, 2), (1, 3)): "CROSS      (Larson form 1)",
    ((0, 3), (1, 2)): "NEST       (Larson form 2)",
    ((0, 1), (0, 2)): "SAME-LEFT  (Larson form 3; tiles sharing upper vertex)",
    ((0, 2), (1, 2)): "SAME-RIGHT (tiles sharing lower vertex)",
    ((0, 1), (1, 2)): "ADJACENT   (chain b=c; consecutive tiles)",
}

def fmt(p):
    return NAME2.get(p, str(p))

# ----------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------

def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

def main():
    # ---------------- n = 2 ----------------
    section("PART 1 — n=2 pattern algebra (omega^2; Specker POSITIVE level)")
    pats2 = all_patterns(2)
    print(f"all patterns of pairs of increasing 2-tuples: {len(pats2)}")
    for p in pats2:
        print("   ", p, "=", fmt(p))
    tri2 = triangle_table(2)
    print(f"realizable (P_xy,P_yz,P_xz) triples: {len(tri2)}")

    # self-composition / single-pattern triangles
    print("\nsingle patterns P with a monochromatic-P triangle (P o P ∋ P realizable):")
    for p in pats2:
        mono = any(a == b == c == p for (a, b, c) in tri2)
        print(f"    {fmt(p):55s} mono-triangle: {mono}")

    mtf2 = maximal_triangle_free_sets(pats2, tri2)
    print(f"\nmaximal triangle-free pattern sets at n=2: {len(mtf2)}")
    for S in mtf2:
        print("    {" + ", ".join(fmt(p) for p in sorted(S)) + "}")

    section("PART 2 — n=2 grid avoidability (Specker positive shadow): for EVERY "
            "triangle-free S there should be S-free grids of growing branching")
    for S in mtf2:
        row = []
        for V in (6, 9, 12, 16, 20):
            row.append((V, max_grid_branching(2, V, S, tmax=6)))
        print("   S={" + ",".join(fmt(p).split()[0] for p in sorted(S)) + "}:",
              " ".join(f"V={V}:t<= {t}" for (V, t) in row))

    # ---------------- n = 3 ----------------
    section("PART 3 — n=3 pattern algebra (omega^3; Specker NEGATIVE level)")
    pats3 = all_patterns(3)
    spats3 = strict_patterns(3)
    print(f"all patterns: {len(pats3)}   strict (12 distinct values... actually 6): {len(spats3)}")
    print("building triangle table over [0,9) (exhaustive; C(9,3)=84 tuples)...")
    tri3 = triangle_table(3)
    print(f"realizable pattern-triples: {len(tri3)}")

    print("\nmono-triangle status of the", len(spats3), "strict patterns:")
    for p in spats3:
        mono = any(a == b == c == p for (a, b, c) in tri3)
        print(f"    {p}   mono-triangle: {mono}")

    print("\nenumerating maximal triangle-free subsets of the STRICT patterns "
          "(witness colourings can be taken strict-measurable after thinning)...")
    mtf3 = maximal_triangle_free_sets(spats3, tri3)
    print(f"maximal triangle-free strict sets: {len(mtf3)}")

    section("PART 4 — n=3 grid hunt: which triangle-free S kill grids "
            "(Specker negative shadow: expect SOME S with bounded t)")
    print("(exact backtracking; t-grid has t^3 leaves; V up to 14)")
    interesting = []
    for S in mtf3:
        caps = []
        for V in (8, 11, 14):
            caps.append(max_grid_branching(3, V, S, tmax=4))
        tag = "  <-- GRID-KILLER candidate" if caps[-1] <= 1 else ""
        if caps[-1] <= 2:
            interesting.append((S, caps))
        print(f"   |S|={len(S)} S={sorted(S)} maxt(V=8,11,14)={caps}{tag}")

    section("PART 5 — the tile dictionary (HYP-2344): n=2 patterns = staircase tile-pair classes")
    print("""
A tile in the repo = an arc (a,b), a>=b+2, i.e. an ordered pair — the staircase cell.
Two tiles t1=(x,y), t2=(x',y') stand in exactly one of the relations:
  share upper vertex / share lower vertex / chain (y=x' or y'=x) / cross / nest / disjoint.
The n=2 pair patterns above are EXACTLY these classes (with Larson's forms 0,1,2,3 the
ones surviving her thinning).  This is the content of HYP-2344 at level n=2: the pattern
calculus of omega^2-partitions IS the staircase tile-pair geometry.
""")
    # mechanical check: classification function on tile pairs matches pattern classes
    def tile_rel(t1, t2):
        (a, b), (c, d) = t1, t2  # tiles as (low, high) increasing pairs here
        s1, s2 = set(t1), set(t2)
        if t1 == t2: return "equal"
        if a == c: return "SAME-LEFT"
        if b == d: return "SAME-RIGHT"
        if b == c or d == a: return "ADJACENT"
        if b < c or d < a: return "DISJOINT"
        if (a < c < b < d) or (c < a < d < b): return "CROSS"
        return "NEST"
    mism = 0
    for x in itertools.combinations(range(8), 2):
        for y in itertools.combinations(range(8), 2):
            if x >= y: continue
            pn = fmt(pattern(x, y)).split()[0]
            tn = tile_rel(x, y)
            if pn != tn:
                mism += 1
                if mism < 5:
                    print("   MISMATCH", x, y, pn, tn)
    print(f"tile-relation vs pattern-class mechanical check over [0,8)^2 pairs: "
          f"{'0 mismatches — DICTIONARY VERIFIED' if mism == 0 else f'{mism} mismatches'}")

    return mtf3, tri3

if __name__ == "__main__":
    main()
