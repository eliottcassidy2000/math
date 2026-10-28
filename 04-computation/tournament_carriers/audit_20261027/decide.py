"""Referee: decide linear / weak trijectivity of an oriented carrier for a GIVEN group of
automorphisms (not necessarily the full Aut(D)).  Independent of the author's code."""
from core import (PermGroup, signed_pair_types, tournament_count, tournament_marks,
                  disjoint_union, chi_induced, fixed_set)


def decide(n, arcs, elems, jdrop=1, want_counts=True):
    """jdrop = j: compare with A-sets X of size n-j (j=1 trijection; j>=2 Maxwell ladder)."""
    E = set(tuple(sorted(a)) for a in arcs)
    G = PermGroup(elems)
    for g in G.elems:  # sanity: every element is an automorphism of D
        assert set((g[a], g[b]) for a, b in arcs) == set(map(tuple, arcs))
    V_act = G.elems
    subs = G.all_subgroups()
    classes = G.subgroup_classes(subs)
    reps = [H for H, _ in classes]
    res = {"order": G.N, "n_subgroups": len(subs)}
    res["crit"] = all(chi_induced(fixed_set(G, H), E)[0] == 1 for H in subs)
    res["defects"] = sorted({(len(H), chi_induced(fixed_set(G, H), E)[0] - 1) for H in subs
                             if chi_induced(fixed_set(G, H), E)[0] != 1})
    # Maxwell count  |E^K| = j|V^K| - C(j+1,2)
    j = jdrop
    res["maxwell_count"] = all(
        chi_induced(fixed_set(G, H), E)[1] == j * len(fixed_set(G, H)) - j * (j + 1) // 2
        for H in subs)
    res["A_fixed"] = len(fixed_set(G, range(G.N)))
    Pi_types = signed_pair_types(G, V_act, n, E)
    Pi_marks = tournament_marks(G, V_act, n, E, reps, subs)
    m = n - jdrop
    trans = sorted(((G.N // len(H), H) for H in reps if G.N // len(H) <= m),
                   key=lambda t: (t[0], sorted(t[1])))
    cached = [G.coset_action(H) for _, H in trans]
    cands = []

    def build(start, rem, ch):
        if rem == 0:
            cands.append(list(ch))
            return
        for t in range(start, len(trans)):
            if trans[t][0] <= rem:
                ch.append(t)
                build(t, rem - trans[t][0], ch)
                ch.pop()
    build(0, m, [])
    lin, weak, counts = [], [], set()
    for c in cands:
        actX, mm = disjoint_union([cached[t] for t in c]) if c else ([tuple()] * G.N, 0)
        desc = [(trans[t][0], len(trans[t][1])) for t in c]  # (orbit size, stabiliser order)
        if signed_pair_types(G, actX, mm) == Pi_types:
            lin.append(desc)
        if tournament_marks(G, actX, mm, frozenset(), reps, subs) == Pi_marks:
            weak.append(desc)
        if want_counts:
            counts.add(tournament_count(G, actX, mm, frozenset()))
    res["linear"] = bool(lin)
    res["linear_X"] = lin[:4]
    res["weak"] = bool(weak)
    res["weak_X"] = weak[:4]
    if want_counts:
        res["count"] = tournament_count(G, V_act, n, E)
        res["count_matches"] = res["count"] in counts
    return res
