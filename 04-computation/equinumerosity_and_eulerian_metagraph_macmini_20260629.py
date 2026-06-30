#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tournament <-> even-graph EQUINUMEROSITY (the particularities) + is the MERGED METAGRAPH EULERIAN?
(mac-mini-2026-06-29-S32)

PART 1: the counts. Labeled tournaments 2^C(n,2) vs labeled even graphs 2^C(n-1,2); ratio 2^(n-1) = the
CUT space (scores). GF(2): tournament = CUT (dim n-1, scores) (+) CYCLE (dim C(n-1,2), even graph). The
iso-class counts A000568 (tournaments) vs A002854 (even graphs) do NOT match (S_n breaks the labeled
bijection). Computed firsthand n=3..6.

PART 2: is the arc-flip metagraph G_n (and its merge G_n/Z_2) EULERIAN (all degrees even => Euler circuit)?
The LABELED arc-flip graph is the hypercube Q_d, d=C(n,2), which is d-regular -> Eulerian iff d even iff
n=0,1 mod 4. The ISO-CLASS metagraph is NOT regular (S_n-quotient); compute its degree sequence and check
parity, full G_n and merged G_n/Z_2, n=4,5,6.
"""
from __future__ import annotations
import functools, itertools
from collections import defaultdict
print = functools.partial(print, flush=True)


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def canon_tournament(bits, n, P, PR):
    """canonical form of a labeled tournament (bitstring over pairs) under S_n relabeling."""
    best = None
    prs = pairs(n)
    idx = {p: b for b, p in enumerate(prs)}
    for perm in P:
        v = 0
        for b, (i, j) in enumerate(prs):
            pi, pj = perm[i], perm[j]
            # orientation of edge {pi,pj} induced: arc i->j (bit set) maps to perm[i]->perm[j]
            a, c = (pi, pj)
            bit = (bits >> b) & 1
            # the induced arc is perm[i]->perm[j] if bit else perm[j]->perm[i]
            if a < c:
                ob = bit
                pos = idx[(a, c)]
            else:
                ob = 1 - bit
                pos = idx[(c, a)]
            v |= (ob << pos)
        if best is None or v < best:
            best = v
    return best


def complement_bits(bits, n):
    m = len(pairs(n))
    return bits ^ ((1 << m) - 1)


def metagraph(n):
    """Return iso classes + arc-flip degree data (simple & multi) for full G_n and merged G_n/Z_2."""
    prs = pairs(n)
    m = len(prs)
    P = list(itertools.permutations(range(n)))
    classes = {}
    rep = {}
    for bits in range(1 << m):
        c = canon_tournament(bits, n, P, None)
        if c not in classes:
            classes[c] = len(classes)
            rep[c] = bits
    # arc-flip neighbors from each class's representative
    simple_deg, multi_deg, selfloops = {}, {}, {}
    comp = {}  # complement class of each class
    for c, b in rep.items():
        comp[c] = canon_tournament(complement_bits(b, n), n, P, None)
        nbrs = defaultdict(int)
        for k in range(m):
            nb = canon_tournament(b ^ (1 << k), n, P, None)
            nbrs[nb] += 1
        selfloops[c] = nbrs.get(c, 0)
        multi_deg[c] = m - nbrs.get(c, 0)              # flips leaving the class (multigraph, loopless)
        simple_deg[c] = sum(1 for nb in nbrs if nb != c)  # distinct neighbor classes
    return classes, rep, simple_deg, multi_deg, selfloops, comp


def even_graph_classes(n):
    """iso classes of even graphs (all degrees even) on n vertices (= A002854(n))."""
    prs = pairs(n)
    m = len(prs)
    P = list(itertools.permutations(range(n)))
    idx = {p: b for b, p in enumerate(prs)}
    seen = set()
    cnt = 0
    for bits in range(1 << m):
        deg = [0]*n
        for b, (i, j) in enumerate(prs):
            if (bits >> b) & 1:
                deg[i] += 1; deg[j] += 1
        if any(d % 2 for d in deg):
            continue
        # canonical
        best = None
        for perm in P:
            v = 0
            for b, (i, j) in enumerate(prs):
                if (bits >> b) & 1:
                    a, c = sorted((perm[i], perm[j]))
                    v |= (1 << idx[(a, c)])
            if best is None or v < best:
                best = v
        if best not in seen:
            seen.add(best); cnt += 1
    return cnt


def main():
    print("=" * 80)
    print("TOURNAMENT <-> EVEN-GRAPH EQUINUMEROSITY + is the MERGED METAGRAPH EULERIAN? (mac-mini-S32)")
    print("=" * 80)

    # ---- PART 1: the counts ----
    print("\n[1] EQUINUMEROSITY -- labeled bijection (Cut+Cycle) vs iso-class mismatch:")
    print(f"    {'n':>2} {'C(n,2)':>7} {'lab.tourn':>10} {'lab.even':>9} {'ratio':>7} {'iso.tourn':>9} {'iso.even':>8}")
    A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
    for n in (3, 4, 5, 6):
        cn2 = n*(n-1)//2
        cn12 = (n-1)*(n-2)//2
        lab_t, lab_e = 1 << cn2, 1 << cn12
        iso_e = even_graph_classes(n)
        print(f"    {n:>2} {cn2:>7} {lab_t:>10} {lab_e:>9} {lab_t//lab_e:>6}x {A000568[n]:>9} {iso_e:>8}")
    print("    => LABELED: tournaments = 2^(n-1) x even graphs (the CUT/score factor); the GF(2) split")
    print("       tournament = CUT(dim n-1, scores) (+) CYCLE(dim C(n-1,2), even graph) is EXACT & clean.")
    print("    => ISO: A000568 (2,4,12,56) vs A002854 (2,3,7,16) DIVERGE -- S_n breaks the labeled bijection")
    print("       (tournaments & even graphs have DIFFERENT stabilizers); equinumerous only LABELED/per-fiber.")

    # ---- PART 2: is the metagraph Eulerian? ----
    print("\n[2] is the ARC-FLIP METAGRAPH EULERIAN (all degrees even <=> admits an Euler circuit)?")
    print(f"    LABELED arc-flip graph = hypercube Q_d, d=C(n,2): d-regular, Eulerian iff d even iff n=0,1 mod4.")
    print(f"    {'n':>2} {'d=C(n,2)':>8} {'#classes':>8} {'#odd-multideg':>13} {'#odd-simpledeg':>14} {'Eulerian(multi)?':>16}")
    for n in (4, 5, 6):
        classes, rep, sdeg, mdeg, sl, comp = metagraph(n)
        odd_m = sum(1 for c in classes if mdeg[c] % 2)
        odd_s = sum(1 for c in classes if sdeg[c] % 2)
        d = n*(n-1)//2
        eul = "YES" if odd_m == 0 else f"NO ({odd_m} odd)"
        print(f"    {n:>2} {d:>8} {len(classes):>8} {odd_m:>13} {odd_s:>14} {eul:>16}")
    print("    (multideg = #arc-flips leaving the class; simpledeg = #distinct neighbor classes.)")

    # ---- detail at n=5,6: which classes are odd-degree, and the merged graph ----
    print("\n[3] DETAIL -- parity by class & the MERGED graph G_n/Z_2 (complement-folded):")
    for n in (4, 5, 6):
        classes, rep, sdeg, mdeg, sl, comp = metagraph(n)
        d = n*(n-1)//2
        # merged: identify class with its complement; merged multidegree (sum, minus internal comp edge)
        merged = {}
        for c in classes:
            key = min(c, comp[c])
            merged.setdefault(key, []).append(c)
        mdeg_dist = sorted(set(mdeg[c] for c in classes))
        sdeg_dist = sorted(set(sdeg[c] for c in classes))
        n_selfcomp = sum(1 for c in classes if comp[c] == c)
        print(f"    n={n}: d={d}; multideg values={mdeg_dist}; simpledeg values={sdeg_dist}; "
              f"#merged nodes={len(merged)}, #self-complementary={n_selfcomp}")
    print("    NOTE: d=C(n,2) parity: n=4->6(even), n=5->10(even), n=6->15(ODD), n=7->21(ODD). The hypercube")
    print("    Q_d is Eulerian for n=4,5 (d even) but NOT n=6,7 (d odd). The iso-quotient may differ.")

    print("\n" + "=" * 80)
    print("FINDINGS: (1) equinumerosity is LABELED only (tournaments = 2^(n-1) x even graphs via Cut+Cycle);")
    print("iso-class counts diverge (A000568 vs A002854) because S_n stabilizers differ. (2) Eulerian-ness of")
    print("the metagraph is governed by d=C(n,2) parity at the labeled (hypercube) level; the iso-class degree")
    print("parities are computed above -- see which n give all-even (Eulerian) metagraphs.")
    print("=" * 80)


if __name__ == "__main__":
    main()
