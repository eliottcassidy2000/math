#!/usr/bin/env python3
"""Auditor bridges-note, script A6 (independent code): explicit isomorphism search between the Cayley-Dickson sign
tournaments CD_n (conventions of THM-5133 (D), (G)) and the cores T_(2^n - 1) of the THM-448 tower, and their
converses.  Backtracking with vertex colours (H of the out-link and of the in-link)."""
import sys
from collections import Counter
import numpy as np

sys.path.insert(0, ".")
from a5_cayley_dickson_vs_tower import cd_tournament, tower_tournament, H_of  # noqa: E402


def colours(m, arcs):
    out = {v: [w for w in range(m) if (v, w) in arcs] for v in range(m)}
    inn = {v: [w for w in range(m) if (w, v) in arcs] for v in range(m)}
    if m <= 31:
        return {v: (H_of(out[v], arcs), H_of(inn[v], arcs)) for v in range(m)}
    return {v: 0 for v in range(m)}


def isomorphisms(m, A, B, colA, colB, limit=None):
    """all (or up to limit) isomorphisms A -> B respecting the colours"""
    adjA = [[(u, v) in A for v in range(m)] for u in range(m)]
    adjB = [[(u, v) in B for v in range(m)] for u in range(m)]
    # order the vertices of A: rarest colour first
    cnt = Counter(colA.values())
    order = sorted(range(m), key=lambda v: (cnt[colA[v]], v))
    cand = {a: [b for b in range(m) if colB[b] == colA[a]] for a in range(m)}
    found = []
    img = {}
    used = [False] * m
    nodes = [0]

    def rec(k):
        if limit is not None and len(found) >= limit:
            return
        if k == m:
            found.append(dict(img))
            return
        a = order[k]
        for b in cand[a]:
            if used[b]:
                continue
            ok = True
            for t in range(k):
                a2 = order[t]
                if adjA[a2][a] != adjB[img[a2]][b]:
                    ok = False
                    break
            if ok:
                nodes[0] += 1
                img[a] = b
                used[b] = True
                rec(k + 1)
                used[b] = False
                del img[a]

    rec(0)
    return found, nodes[0]


def main():
    for n in (3, 4, 5):
        m, cd = cd_tournament(n)
        m2, tw = tower_tournament(n)
        assert m == m2
        twop = {(v, u) for (u, v) in tw}
        ccd, ctw, ctwop = colours(m, cd), colours(m, tw), colours(m, twop)
        aut_cd, nd0 = isomorphisms(m, cd, cd, ccd, ccd)
        aut_tw, nd00 = isomorphisms(m, tw, tw, ctw, ctw)
        iso1, nd1 = isomorphisms(m, cd, tw, ccd, ctw)
        iso2, nd2 = isomorphisms(m, cd, twop, ccd, ctwop)
        anti, nd3 = isomorphisms(m, tw, twop, ctw, ctwop)
        print("n = %d (%d vertices): |Aut(CD_n)| = %d, |Aut(T)| = %d; isomorphisms CD_n -> T: %d (search nodes %d); "
              "CD_n -> converse of T: %d (search nodes %d); T -> converse of T (anti-automorphisms): %d"
              % (n, m, len(aut_cd), len(aut_tw), len(iso1), nd1, len(iso2), nd2, len(anti)))
        if iso1:
            f = iso1[0]
            assert all(((f[u], f[v]) in tw) for (u, v) in cd)
            print("    explicit isomorphism CD_%d -> T_%d (vertex i -> f(i), 0-based): %s" % (n, m, [f[i] for i in range(m)]))
        if iso2:
            f = iso2[0]
            assert all(((f[u], f[v]) in twop) for (u, v) in cd)
            print("    explicit isomorphism CD_%d -> converse of T_%d: %s" % (n, m, [f[i] for i in range(m)]))


if __name__ == "__main__":
    main()
