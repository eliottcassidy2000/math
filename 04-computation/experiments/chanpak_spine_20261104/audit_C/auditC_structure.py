#!/usr/bin/env python3
"""auditC_structure.py -- audit C: structural verification of THM-4876 (2) and (3) on the (T, v) pairs with increment 2 or 4
listed by auditC_census del (files del<n>.out, lines "P <tour> <v> <H> <H(T-v)>").

For increment 2 (k = 1): H(T-v) in {1,3}; H(T) in {3,5}; exactly one odd cycle through v; it is a triangle (v,a,b);
T - {v,a,b} transitive; the only arc from A = N+(v) to B = N-(v) is a -> b; NOT both (x in A' with x -> a) and
(y in B' with b -> y); T - v has at most one odd cycle; the realisations by C3 and the strong 4-tournament.
For increment 4 (k = 2): exactly two odd cycles through v, both with transitive complement; H(T-v) in {1,5}.
"""
import sys, os
from collections import Counter
from auditC_lib import *

HERE = os.path.dirname(os.path.abspath(__file__))


def check_order(n):
    fn = os.path.join(HERE, "del%d.out" % n)
    res = {1: Counter(), 2: Counter()}
    viol = Counter()
    pairs = {1: 0, 2: 0}
    hv_sets = {1: set(), 2: set()}
    with open(fn) as f:
        for line in f:
            if not line.startswith("P "):
                continue
            _, tstr, v, Ht, Hv = line.split()
            v, Ht, Hv = int(v), int(Ht), int(Hv)
            adj = parse_gentourng_line(n, tstr)
            k = (Ht - Hv) // 2
            pairs[k] += 1
            hv_sets[k].add(Hv)
            full = (1 << n) - 1
            E = dp_end(n, adj)
            hs = hsub_table(n, E)
            assert hs[full] == Ht and hs[full & ~(1 << v)] == Hv, "re-derived counts disagree with the C census"
            cyc = odd_cycles_through(n, adj, v)
            comps = [hs[full & ~S] for S, L in cyc]
            if k == 1:
                if Hv not in (1, 3):
                    viol["k1: H(T-v) not in {1,3}"] += 1
                if Ht not in (3, 5):
                    viol["k1: H(T) not in {3,5}"] += 1
                if len(cyc) != 1:
                    viol["k1: #odd cycles through v != 1"] += 1
                    continue
                S, L = cyc[0]
                if L != 3:
                    viol["k1: the cycle is not a triangle"] += 1
                if not is_transitive(n, adj, full & ~S):
                    viol["k1: complement of the triangle not transitive"] += 1
                A = adj[v]
                B = full & ~A & ~(1 << v)
                # the triangle (v, a, b): a in A, b in B
                tri = [u for u in range(n) if (S >> u) & 1 and u != v]
                a = [u for u in tri if (A >> u) & 1]
                b = [u for u in tri if (B >> u) & 1]
                if len(a) != 1 or len(b) != 1:
                    viol["k1: triangle not of the form (v, a, b) with a in A, b in B"] += 1
                    continue
                a, b = a[0], b[0]
                arcsAB = [(x, y) for x in range(n) if (A >> x) & 1 for y in range(n) if (B >> y) & 1 and (adj[x] >> y) & 1]
                if arcsAB != [(a, b)]:
                    viol["k1: arcs from A to B other than a -> b"] += 1
                Ap = A & ~(1 << a)
                Bp = B & ~(1 << b)
                in_a = any((adj[x] >> a) & 1 for x in range(n) if (Ap >> x) & 1)
                out_b = any((adj[b] >> y) & 1 for y in range(n) if (Bp >> y) & 1)
                if in_a and out_b:
                    viol["k1: both an in-neighbour of a in A' and an out-neighbour of b in B' exist"] += 1
                n2, adj2 = delete_vertex(n, adj, v)
                nodd = len(all_odd_cycles(n2, adj2))
                if nodd > 1:
                    viol["k1: T - v has more than one odd cycle"] += 1
                res[1][(Hv, Ht, len(cyc), tuple(sorted(comps)), nodd)] += 1
            else:
                if Hv not in (1, 5):
                    viol["k2: H(T-v) not in {1,5}"] += 1
                if len(cyc) != 2:
                    viol["k2: #odd cycles through v != 2"] += 1
                if any(not is_transitive(n, adj, full & ~S) for S, L in cyc):
                    viol["k2: a complement is not transitive"] += 1
                res[2][(Hv, Ht, len(cyc), tuple(sorted(comps)), tuple(sorted(L for S, L in cyc)))] += 1
    print("order %d: increment-2 pairs %d (H(T-v) values %s), increment-4 pairs %d (H(T-v) values %s)" % (n, pairs[1], sorted(hv_sets[1]), pairs[2], sorted(hv_sets[2])))
    print("   increment-2 profiles (H(T-v), H(T), #cycles through v, complement counts, #odd cycles of T-v) -> count: %s" % dict(res[1]))
    print("   increment-4 profiles (H(T-v), H(T), #cycles, complement counts, cycle lengths) -> count: %s" % dict(res[2]))
    print("   violations: %s" % (dict(viol) if viol else "none"))
    return viol


if __name__ == "__main__":
    tot = Counter()
    for n in range(3, 11):
        if os.path.exists(os.path.join(HERE, "del%d.out" % n)):
            tot += check_order(n)
    print("TOTAL violations over orders 3..10:", dict(tot) if tot else "none")
