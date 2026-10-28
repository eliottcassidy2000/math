"""Referee census of oriented excess -1 carriers, n = 3..NMAX (default 7).

For every oriented carrier D (oriented graph, no 2-cycles, n-1 arcs) up to isomorphism:
  A = Aut(D); decide LINEAR trijectivity by direct isomorphism test of signed A-sets
  (multiset of (A x Z2)-orbit types) against every A-set X of size n-1 built from coset spaces;
  decide WEAK trijectivity by comparing marks of the A-sets {T >= D} and Tour(X) at every subgroup;
  decide criterion (C) chi(G^K)=1 at every subgroup (and cyclic subgroups, and via the
  intersection-closure of element fixed sets).  Independent of the author's code.
Usage: python census.py [NMAX]
"""
import itertools
import json
import sys
import time
from collections import Counter

import numpy as np

from core import (PermGroup, signed_pair_types, tournament_count, tournament_marks,
                  disjoint_union, chi_induced, fixed_set, automorphisms, cycles)


def graph_reps(n, m):
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    idxarr = np.full((n, n), -1, dtype=np.int64)
    for (i, j), k in pid.items():
        idxarr[i, j] = k
        idxarr[j, i] = k
    pm = np.empty((len(perms), len(pairs)), dtype=np.int64)
    for k, (i, j) in enumerate(pairs):
        pm[:, k] = idxarr[perms[:, i], perms[:, j]]
    one = np.int64(1)

    def bits_all(el):
        if not el:
            return np.zeros(len(perms), dtype=np.int64)
        return (one << pm[:, list(el)]).sum(axis=1)

    level = {0: ()}
    for _ in range(m):
        nxt = {}
        for _, el in level.items():
            for k in range(len(pairs)):
                if k in el:
                    continue
                el2 = tuple(sorted(el + (k,)))
                c = int(bits_all(el2).min())
                if c not in nxt:
                    nxt[c] = el2
        level = nxt
    graphs = []
    for c, el in sorted(level.items()):
        b = bits_all(el)
        own = sum(1 << k for k in el)
        aut = [tuple(int(x) for x in perms[t]) for t in np.nonzero(b == own)[0]]
        graphs.append(([pairs[k] for k in el], aut))
    return graphs


def oriented_reps(n, edges, autG):
    m = len(edges)
    reps = {}
    for mask in range(2 ** m):
        arcs = [(a, b) if (mask >> t) & 1 == 0 else (b, a) for t, (a, b) in enumerate(edges)]
        best = None
        for p in autG:
            code = sum(1 << (p[a] * n + p[b]) for a, b in arcs)
            if best is None or code < best:
                best = code
        if best not in reps:
            reps[best] = arcs
    return list(reps.values())


def kind_flags(n, arcs):
    outd = Counter(a for a, b in arcs)
    ind = Counter(b for a, b in arcs)
    func = sum(1 for v in range(n) if outd[v] == 0) == 1 and all(outd[v] <= 1 for v in range(n))
    cofunc = sum(1 for v in range(n) if ind[v] == 0) == 1 and all(ind[v] <= 1 for v in range(n))
    # tree: connected with n-1 edges
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in arcs:
        par[f(a)] = f(b)
    tree = len({f(v) for v in range(n)}) == 1
    return tree, func, cofunc


def analyse(n, arcs, full=True):
    E = set(tuple(sorted(a)) for a in arcs)
    autA = automorphisms(n, arcs)
    G = PermGroup(autA)
    V_act = G.elems  # action on V is the permutation itself
    rec = {"arcs": arcs, "order": G.N}
    rec["tree"], rec["functional"], rec["cofunctional"] = kind_flags(n, arcs)
    # criterion via intersection-closure of element fixed sets (no subgroup enumeration)
    elem_fix = {fixed_set(G, [g]) for g in range(G.N)}
    closure = set(elem_fix)
    frontier = set(elem_fix)
    while frontier:
        new = set()
        for F1 in frontier:
            for F2 in elem_fix:
                F3 = F1 & F2
                if F3 not in closure:
                    new.add(F3)
        closure |= new
        frontier = new
    rec["crit_intersection"] = all(chi_induced(F, E)[0] == 1 for F in closure)
    rec["crit_cyclic"] = all(chi_induced(F, E)[0] == 1 for F in elem_fix)
    VA = fixed_set(G, range(G.N))
    rec["A_fixes_vertex"] = len(VA) > 0
    # involution / odd defects (elementwise)
    rec["involution_defect"] = any(G.order_of[g] == 2 and chi_induced(fixed_set(G, [g]), E)[0] != 1
                                   for g in range(G.N))
    rec["odd_elem_defect"] = any(G.order_of[g] % 2 == 1 and g != G.e and
                                 chi_induced(fixed_set(G, [g]), E)[0] != 1 for g in range(G.N))
    Pi_types = signed_pair_types(G, V_act, n, E)
    rec["count"] = tournament_count(G, V_act, n, E)
    # cheap linear test: X = V \ {r}
    lin_cheap = False
    for r in sorted(VA):
        keep = [v for v in range(n) if v != r]
        pos = {v: i for i, v in enumerate(keep)}
        actX = [tuple(pos[g[v]] for v in keep) for g in G.elems]
        if signed_pair_types(G, actX, n - 1) == Pi_types:
            lin_cheap = True
            break
    rec["linear_cheap"] = lin_cheap
    if not full:
        rec["linear"] = True if lin_cheap else None
        return rec
    subs = G.all_subgroups()
    classes = G.subgroup_classes(subs)
    reps = [H for H, _ in classes]
    rec["n_subgroups"] = len(subs)
    rec["crit_all"] = all(chi_induced(fixed_set(G, H), E)[0] == 1 for H in subs)
    rec["odd_sub_defect"] = any(len(H) % 2 == 1 and len(H) > 1 and
                                chi_induced(fixed_set(G, H), E)[0] != 1 for H in subs)
    # transitive A-sets of size <= n-1
    trans = []
    for H in reps:
        idx = G.N // len(H)
        if idx <= n - 1:
            trans.append((idx, H))
    trans.sort(key=lambda t: (t[0], sorted(t[1])))
    # all multisets of transitive sets with total size n-1
    cands = []

    def rec_build(start, remaining, chosen):
        if remaining == 0:
            cands.append(list(chosen))
            return
        for t in range(start, len(trans)):
            sz = trans[t][0]
            if sz <= remaining:
                chosen.append(t)
                rec_build(t, remaining - sz, chosen)
                chosen.pop()
    rec_build(0, n - 1, [])
    cached = {}
    for t, (sz, H) in enumerate(trans):
        cached[t] = G.coset_action(H)
    Pi_marks = tournament_marks(G, V_act, n, E, reps, subs)
    linear_X = []
    weak_X = []
    count_X = set()
    for cand in cands:
        actX, m = disjoint_union([cached[t] for t in cand])
        if m != n - 1:
            raise RuntimeError("size")
        desc = [trans[t][0] for t in cand]
        if signed_pair_types(G, actX, m) == Pi_types:
            linear_X.append(desc)
        if tournament_marks(G, actX, m, frozenset(), reps, subs) == Pi_marks:
            weak_X.append(desc)
        count_X.add(tournament_count(G, actX, m, frozenset()))
    rec["linear"] = len(linear_X) > 0
    rec["weak"] = len(weak_X) > 0
    rec["n_candidates"] = len(cands)
    rec["count_matches_some_X"] = rec["count"] in count_X
    rec["linear_X_orbit_sizes"] = linear_X[:3]
    if rec["linear"] != lin_cheap:
        rec["note_cheap_vs_full"] = True
    return rec


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    nmin = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    out = {}
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        graphs = graph_reps(n, n - 1)
        carriers = []
        for edges, autG in graphs:
            for arcs in oriented_reps(n, edges, autG):
                carriers.append(arcs)
        recs = []
        for arcs in carriers:
            big = len(automorphisms(n, arcs)) > 200
            recs.append(analyse(n, arcs, full=not big))
        out[n] = recs
        lin = sum(1 for r in recs if r["linear"])
        weak_only = sum(1 for r in recs if r.get("weak") and not r["linear"])
        neither = sum(1 for r in recs if r.get("weak") is False and not r["linear"])
        big = sum(1 for r in recs if "crit_all" not in r)
        print(f"n={n}: underlying graphs {len(graphs)}, oriented carriers {len(recs)}, "
              f"linear {lin}, weak-only {weak_only}, neither {neither}, big-group(skipped full) {big}, "
              f"{time.time() - t0:.1f}s", flush=True)
    with open(f"census_n{nmax}.json" if nmin == 3 else f"census_n{nmin}_{nmax}.json", "w") as fh:
        json.dump({str(k): v for k, v in out.items()}, fh)


if __name__ == "__main__":
    main()
