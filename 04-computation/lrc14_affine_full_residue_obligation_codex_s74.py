#!/usr/bin/env python3
"""
codex-2026-06-21-S74

Affine/full-residue proof-obligation scout for LRC(14).

This is a synthesis script, not a new atlas.  The user asked to use the CJJ
complete LP hierarchy, the full-residue reduction, and tournament extremality
as proof inspiration.  Incoming work suggests the right CJJ fix is an affine
quotient rather than a linear-code quotient.  This script checks the first
possible pitfall: Aff(1,F_7) is 2-transitive, so residue-only affine moments may
forget too much on the HYP-2749 full-residue stratum.

Tournament Analysis is included with proof lenses as vertices.  The pairwise
observable is a score vector:

    (pins_consec, avoids_circularity, preserves_gap_word, formal_ready, reusable)

The switch is lexicographic comparison of that vector; ties are broken by name.
The Hamiltonian path count is then a fingerprint of how rigid the proof-order
tournament is.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
import importlib.util
import os


P = 7


def load_measS7():
    here = os.path.dirname(__file__)
    path = os.path.join(here, "lrc_q108_relation_code_mds_kps.py")
    spec = importlib.util.spec_from_file_location("relation_code_mds", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.measS7


measS7 = load_measS7()


def aff_maps():
    for a in range(1, P):
        for b in range(P):
            yield a, b


def canonical_affine_multiset(residues):
    """Canonical Aff(1,F_7)-orbit representative of a residue multiset."""
    reps = []
    for a, b in aff_maps():
        reps.append(tuple(sorted((a * r + b) % P for r in residues)))
    return min(reps)


def affine_subset_orbits(size):
    """Orbit partition of ordinary subsets of F_7 under Aff(1,F_7)."""
    seen = set()
    orbits = []
    for S in combinations(range(P), size):
        if S in seen:
            continue
        orb = {tuple(sorted((a * x + b) % P for x in S)) for a, b in aff_maps()}
        seen |= orb
        orbits.append(sorted(orb))
    return sorted(orbits, key=lambda o: (len(o), o[0]))


def residue_counts(E):
    return tuple(sorted(Counter(e % P for e in E).values(), reverse=True))


def affine_pair_signature(E):
    res = [e % P for e in E]
    collision = 0
    distinct = 0
    for i, j in combinations(range(len(res)), 2):
        if res[i] == res[j]:
            collision += 1
        else:
            distinct += 1
    return distinct, collision


def affine_r_profile(E, r):
    """Profile of r-index subsets after quotienting residues by Aff(1,F_7)."""
    res = [e % P for e in E]
    prof = Counter()
    for idxs in combinations(range(len(res)), r):
        prof[canonical_affine_multiset([res[i] for i in idxs])] += 1
    return tuple(sorted(prof.items()))


def full_residue_shapes_k8(max_e=15):
    for rest in combinations(range(1, max_e + 1), 7):
        E = (0,) + rest
        if set(e % P for e in E) == set(range(P)):
            yield E


def gap_word(E):
    E = tuple(sorted(E))
    return tuple(E[i + 1] - E[i] for i in range(len(E) - 1))


def is_ap_shape(E):
    g = gap_word(E)
    return len(set(g)) == 1


def summarize_full_residue_k8():
    shapes = list(full_residue_shapes_k8())
    rows = []
    for E in shapes:
        rows.append(
            {
                "E": E,
                "measure": measS7(E),
                "res_count": residue_counts(E),
                "pair": affine_pair_signature(E),
                "triple": affine_r_profile(E, 3),
                "quad": affine_r_profile(E, 4),
                "gap": gap_word(E),
                "is_ap": is_ap_shape(E),
            }
        )
    max_measure = max(r["measure"] for r in rows)
    argmax = [r for r in rows if r["measure"] == max_measure]
    consec = tuple(range(8))
    consec_row = next(r for r in rows if r["E"] == consec)

    def class_count(key):
        return len({r[key] for r in rows})

    def sharing_ap(key):
        return sum(1 for r in rows if r[key] == consec_row[key])

    return {
        "rows": rows,
        "total": len(rows),
        "consec": consec_row,
        "max_measure": max_measure,
        "argmax": argmax,
        "class_counts": {
            "residue_count": class_count("res_count"),
            "affine_pair": class_count("pair"),
            "affine_triple_profile": class_count("triple"),
            "affine_quad_profile": class_count("quad"),
            "integer_gap_word": class_count("gap"),
        },
        "sharing_ap": {
            "residue_count": sharing_ap("res_count"),
            "affine_pair": sharing_ap("pair"),
            "affine_triple_profile": sharing_ap("triple"),
            "affine_quad_profile": sharing_ap("quad"),
            "integer_gap_word": sharing_ap("gap"),
        },
    }


LENSES = {
    "full_residue_localization": (3, 3, 1, 3, 3),
    "relation_pair_marginals": (3, 3, 2, 2, 3),
    "affine_gap_word": (3, 3, 3, 1, 3),
    "signed_delsarte_dual": (2, 2, 1, 3, 3),
    "theta_prime_ceiling": (1, 3, 1, 3, 2),
    "affine_residue_pairs": (0, 3, 0, 2, 2),
    "boolean_mobius_integrality": (1, 1, 1, 2, 2),
    "signed_tanner_dessin": (1, 3, 0, 2, 1),
    "raw_tanner_expansion": (0, 2, 0, 1, 1),
}


def tournament_from_scores(scores):
    names = list(scores)
    adj = {u: set() for u in names}
    for i, u in enumerate(names):
        for v in names[i + 1 :]:
            if (scores[u], u) >= (scores[v], v):
                adj[u].add(v)
            else:
                adj[v].add(u)
    return names, adj


def count_three_cycles(names, adj):
    c = 0
    for a, b, c0 in combinations(names, 3):
        edges = 0
        for u, v in ((a, b), (b, c0), (c0, a)):
            if v in adj[u]:
                edges += 1
        if edges in (0, 3):
            c += 1
    return c


def sccs(names, adj):
    rev = {u: set() for u in names}
    for u in names:
        for v in adj[u]:
            rev[v].add(u)

    seen = set()
    order = []

    def dfs1(u):
        seen.add(u)
        for v in adj[u]:
            if v not in seen:
                dfs1(v)
        order.append(u)

    for u in names:
        if u not in seen:
            dfs1(u)

    comps = []
    seen.clear()

    def dfs2(u, comp):
        seen.add(u)
        comp.append(u)
        for v in rev[u]:
            if v not in seen:
                dfs2(v, comp)

    for u in reversed(order):
        if u not in seen:
            comp = []
            dfs2(u, comp)
            comps.append(tuple(sorted(comp)))
    return comps


def hamiltonian_paths(names, adj):
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    outmask = [0] * n
    for u in names:
        i = idx[u]
        for v in adj[u]:
            outmask[i] |= 1 << idx[v]

    @lru_cache(None)
    def dp(mask, last):
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if prev_mask & (1 << prev) and (outmask[prev] & (1 << last)):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def print_tournament():
    names, adj = tournament_from_scores(LENSES)
    scores = {u: len(adj[u]) for u in names}
    hist = Counter(scores.values())
    comps = sccs(names, adj)
    hp = hamiltonian_paths(names, adj)
    c3 = count_three_cycles(names, adj)
    ranked = sorted(names, key=lambda u: (LENSES[u], u), reverse=True)
    print("TOURNAMENT ANALYSIS: proof lenses as vertices")
    print("  pair observable: (pins, noncircular, preserves_gap, formal_ready, reusable)")
    print("  switch: lexicographic larger vector wins; ties by name")
    print("  ranking:")
    for u in ranked:
        print(f"    {u:28s} score_vector={LENSES[u]} outdegree={scores[u]}")
    print(f"  score histogram: {dict(sorted(hist.items()))}")
    print(f"  directed 3-cycles: {c3}")
    print(f"  SCC sizes: {[len(c) for c in comps]}")
    print(f"  Hamiltonian path count: {hp}")


def main():
    print("=" * 78)
    print("LRC14 affine/full-residue proof-obligation scout (codex S74)")
    print("=" * 78)
    print("Source readings:")
    print("  arXiv:2211.01248 = Exact Completeness of LP Hierarchies for Linear Codes.")
    print("  arXiv:2112.09221 = A Complete Linear Programming Hierarchy for Linear Codes.")
    print("  repo signal: HYP-2749 reduces consec-max to the full-Z/7-residue stratum.")
    print("  repo signal: HYP-2750 proposes an affine CJJ quotient for the AP extremizer.")
    print()

    print("AFFINE ORBITS ON ORDINARY SUBSETS OF F_7")
    for r in range(0, 8):
        orbs = affine_subset_orbits(r)
        sizes = sorted(len(o) for o in orbs)
        reps = [o[0] for o in orbs]
        print(f"  r={r}: orbit_count={len(orbs)} orbit_sizes={sizes} reps={reps}")
    print("  Consequence: level <=2 residue-only affine data is forced to collapse.")
    print()

    summary = summarize_full_residue_k8()
    consec = summary["consec"]
    print("K=8 FULL-RESIDUE STRATUM IN [0,15]")
    print(f"  total shapes: {summary['total']}")
    print(f"  consec: E={consec['E']} measS7={consec['measure']} ({float(consec['measure']):.6f})")
    print(
        f"  max measS7={summary['max_measure']} ({float(summary['max_measure']):.6f}); "
        f"argmax_count={len(summary['argmax'])}"
    )
    print(f"  argmax shapes: {[r['E'] for r in summary['argmax']]}")
    print("  quotient class counts over this stratum:")
    for key, val in summary["class_counts"].items():
        print(f"    {key:24s}: {val:4d} classes; AP-class size={summary['sharing_ap'][key]}")
    print()
    print("INTERPRETATION")
    print("  The residue-count and affine-pair quotients see every full-residue k=8 shape")
    print("  as the same object: one duplicated residue and all seven residues present.")
    print("  Even affine triple/quad residue profiles remain far coarser than the integer")
    print("  gap word.  Thus the affine CJJ route must retain generated gap-word or")
    print("  relation-code marginal data after quotienting; residue-only affine moments")
    print("  cannot be the final proof certificate.")
    print()
    print_tournament()
    print()
    print("PROOF OBLIGATION PROPOSED")
    print("  P0: use HYP-2749 to discard non-full residue shapes.")
    print("  P1: normalize the remaining shapes by translation+dilation, but retain the")
    print("      integer gap/generated-word data rather than only residues.")
    print("  P2: prove AP maximizes the level-2 relation-code pair-marginal bound on this")
    print("      affine-normalized full-residue quotient.")
    print("  P3: splice that extremality into the signed THM-534 Delsarte dual; use the")
    print("      theta-prime/tournament bridge as a ceiling comparison, not as a proof by")
    print("      itself.")


if __name__ == "__main__":
    main()
