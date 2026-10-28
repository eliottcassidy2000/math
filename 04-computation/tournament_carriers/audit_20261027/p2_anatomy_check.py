"""Referee: subgroup-level anatomy of the two j=3, n=5 P2 disagreements (author's listed arcs)."""
from math import comb
from core import automorphisms, PermGroup, chi_induced, fixed_set
from decide import decide

for arcs in ([(0, 1), (0, 3), (0, 4), (2, 1), (1, 3), (4, 1), (2, 3), (2, 4), (3, 4)],):
    n, j = 5, 3
    E = set(tuple(sorted(a)) for a in arcs)
    G = PermGroup(automorphisms(n, arcs))
    subs = G.all_subgroups()
    print("|A| =", G.N, " element orders:", sorted(G.order_of))
    for H in sorted(subs, key=len):
        F = fixed_set(G, H)
        e = chi_induced(F, E)[1]
        ok = e == j * len(F) - comb(j + 1, 2)
        print(f"  |K|={len(H)} |V^K|={len(F)} |E^K|={e} count_ok={ok}")
    r = decide(n, arcs, G.elems, jdrop=j, want_counts=False)
    print("  linear 3-jective:", r["linear"], "X:", r["linear_X"], " count criterion:", r["maxwell_count"])
