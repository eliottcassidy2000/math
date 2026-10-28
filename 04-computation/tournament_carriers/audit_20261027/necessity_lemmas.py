"""Referee: numerical check of the lemmas in the referee's proof of FULL necessity
(linear => chi(G^K) = 1 for every subgroup K of A = Aut(D)), on every oriented carrier n <= NMAX.

Lemma A  (twins): if two distinct components are isomorphic as oriented graphs and have chi != 0,
         the transposition swapping them is an involution with delta = -2 chi != 0.
Lemma B  (sign):  if every component that has an isomorphic twin has chi = 0, then delta(K) >= 0
         for every subgroup K.
Lemma E  (cyclic => all): if delta(<g>) = 0 for every element g, then delta(K) = 0 for all K.
Lemma F  (fixed vertex): linear => A fixes a vertex.
Usage: python necessity_lemmas.py NMAX"""
import sys
import time

from census import graph_reps, oriented_reps, analyse
from core import automorphisms, PermGroup, chi_induced, fixed_set


def components(n, arcs):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in arcs:
        par[f(a)] = f(b)
    comp = {}
    for v in range(n):
        comp.setdefault(f(v), []).append(v)
    return list(comp.values())


def iso_components(arcs, c1, c2):
    if len(c1) != len(c2):
        return False
    s1, s2 = set(c1), set(c2)
    a1 = [(a, b) for a, b in arcs if a in s1]
    a2 = [(a, b) for a, b in arcs if a in s2]
    if len(a1) != len(a2):
        return False
    idx1 = {v: i for i, v in enumerate(c1)}
    idx2 = {v: i for i, v in enumerate(c2)}
    m = len(c1)
    A1 = [(idx1[a], idx1[b]) for a, b in a1]
    A2 = set((idx2[a], idx2[b]) for a, b in a2)
    # brute force small components
    import itertools
    for p in itertools.permutations(range(m)):
        if set((p[a], p[b]) for a, b in A1) == A2:
            return True
    return False


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    nmin = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        stats = dict(carriers=0, twin_nonzero=0, twin_nonzero_and_linear=0, B_hyp=0, B_viol=0,
                     E_viol=0, F_viol=0)
        for edges, autG in graph_reps(n, n - 1):
            for arcs in oriented_reps(n, edges, autG):
                stats["carriers"] += 1
                E = set(tuple(sorted(a)) for a in arcs)
                comps = components(n, arcs)
                chis = [chi_induced(set(c), E)[0] for c in comps]
                twin_nonzero = any(chis[i] != 0 and iso_components(arcs, comps[i], comps[j])
                                   for i in range(len(comps)) for j in range(i + 1, len(comps)))
                elems = automorphisms(n, arcs)
                if len(elems) > 200:
                    continue
                G = PermGroup(elems)
                subs = G.all_subgroups()
                deltas = [chi_induced(fixed_set(G, H), E)[0] - 1 for H in subs]
                cyc_ok = all(chi_induced(fixed_set(G, [g]), E)[0] == 1 for g in range(G.N))
                all_ok = all(d == 0 for d in deltas)
                lin = analyse(n, arcs)["linear"]
                if twin_nonzero:
                    stats["twin_nonzero"] += 1
                    stats["twin_nonzero_and_linear"] += lin
                else:
                    stats["B_hyp"] += 1
                    if any(d < 0 for d in deltas):
                        stats["B_viol"] += 1
                if cyc_ok != all_ok:
                    stats["E_viol"] += 1
                if lin and len(fixed_set(G, range(G.N))) == 0:
                    stats["F_viol"] += 1
        print(f"n={n}: {stats} [{time.time() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    main()
