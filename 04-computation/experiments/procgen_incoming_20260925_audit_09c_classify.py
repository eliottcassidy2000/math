"""Seed pre-screen, larger range and classification of control hits.
Plus: odd starts <= N, first edge vs all later edges, edges touching 1 excluded.
Minus / 5n+1 / 3n+5: classify hits by whether the later (or earlier) edge lies on a cycle."""
import sys, collections
from procgen_incoming_20260925_audit_09_seed_berggren_antichain import orbit_edges, comparable, U

def cycles_of(q, s, bound=10**6, cap=400):
    cyc = set()
    for x in range(1, 20001, 2):
        seen = {}
        y = x; i = 0
        while y not in seen and i < cap and y < bound:
            seen[y] = i; y = U(y, q, s); i += 1
        if y in seen:
            # collect cycle
            z = y
            while True:
                cyc.add(z); z = U(z, q, s)
                if z == y: break
    return cyc

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200001
for (q, s, name) in ((3, 1, "PLUS"), (3, -1, "MINUS"), (5, 1, "5n+1"), (3, 5, "3n+5")):
    cyc = cycles_of(q, s)
    hits = collections.Counter(); ex = []
    lim = N if name == "PLUS" else 20001
    for x in range(3, lim + 1, 2):
        if q == 3 and s == 5 and x % 5 == 0 and False: pass
        E = [e for e in orbit_edges(x, q, s) if e[1] > 1]
        if not E: continue
        P = E[0]
        for d, Qe in enumerate(E[1:], start=1):
            c = comparable(P, Qe)
            if c:
                oncyc = (Qe[0] in cyc and Qe[1] in cyc) or (P[0] in cyc and P[1] in cyc)
                key = (c, "cycle-edge" if oncyc else ("consecutive" if d == 1 else "other"))
                hits[key] += 1
                if key[1] == "other" and len(ex) < 8: ex.append((x, d, P, Qe, c))
    print(name, "starts<=%d" % lim, "cycle elements(min):", sorted(cyc)[:12], "hits:", dict(hits))
    for e in ex: print("    other:", e)
