"""Referee: independent recomputation of the pre-registered P2 test (Maxwell ladder).
For j=2 (n=4,5,6; 2n-3 arcs) and j=3 (n=5,6; 3n-6 arcs): every oriented graph up to iso;
linear j-jectivity = signed A-set iso Pi ~ pairs(X), |X| = n-j (search over all X);
count criterion |E^K| = j|V^K| - C(j+1,2) for all K; weak = A-set iso {T>=D} ~ Tour(X).
Also the anatomy of every disagreement.  Run: python maxwell.py"""
import json
import time
from collections import Counter
from math import comb

from census import graph_reps, oriented_reps
from core import automorphisms, PermGroup, chi_induced, fixed_set
from decide import decide


def anatomy(n, arcs, j):
    E = set(tuple(sorted(a)) for a in arcs)
    G = PermGroup(automorphisms(n, arcs))
    subs = G.all_subgroups()
    bad = []
    for H in subs:
        F = fixed_set(G, H)
        e = chi_induced(F, E)[1]
        if e != j * len(F) - comb(j + 1, 2):
            # alternative fixed counts x with C(x,2) = C(v,2) - e
            xs = [x for x in range(0, n) if comb(x, 2) == comb(len(F), 2) - e]
            bad.append((len(H), len(F), e, tuple(xs)))
    return sorted(set(bad)), len(fixed_set(G, range(G.N)))


def main():
    out = {}
    for j, ns in ((2, (4, 5, 6)), (3, (5, 6))):
        for n in ns:
            t0 = time.time()
            m = j * n - comb(j + 1, 2)
            carriers = []
            for edges, autG in graph_reps(n, m):
                carriers += oriented_reps(n, edges, autG)
            tally = Counter()
            dis = []
            for arcs in carriers:
                r = decide(n, arcs, automorphisms(n, arcs), jdrop=j, want_counts=False)
                key = (r["linear"], r["maxwell_count"], r["weak"])
                tally[key] += 1
                if r["linear"] != r["maxwell_count"]:
                    bad, fixA = anatomy(n, arcs, j)
                    dis.append({"arcs": arcs, "order": r["order"], "linear": r["linear"],
                                "count": r["maxwell_count"], "A_fixed_vertices": fixA,
                                "count_failures(|K|,|V^K|,|E^K|,x-options)": bad})
            out[f"j={j},n={n}"] = {"tally": {str(k): v for k, v in tally.items()}, "disagreements": dis}
            print(f"j={j} n={n} ({m} arcs): carriers {len(carriers)}; (linear, count, weak) {dict(tally)}; "
                  f"disagreements {len(dis)} [{time.time() - t0:.1f}s]", flush=True)
            for d in dis:
                print(f"    lin={d['linear']} count={d['count']} |A|={d['order']} A-fixed={d['A_fixed_vertices']} "
                      f"failures={d['count_failures(|K|,|V^K|,|E^K|,x-options)']}")
    json.dump(out, open("maxwell_referee.json", "w"), indent=1)


if __name__ == "__main__":
    main()
