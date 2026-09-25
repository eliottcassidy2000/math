from procgen_incoming_20260925_audit_09_seed_berggren_antichain import orbit_edges, comparable
import collections
for (q, s, name) in ((3, 1, "PLUS"), (3, -1, "MINUS"), (5, 1, "5n+1")):
    hits = collections.Counter(); ex = []
    smallQ = collections.Counter()
    for x in range(3, 20002, 2):
        E = [e for e in orbit_edges(x, q, s) if e[1] > 1]  # drop edges touching 1
        if not E: continue
        P = E[0]
        for d, Q in enumerate(E[1:], start=1):
            c = comparable(P, Q)
            if c:
                hits[c] += 1
                smallQ[Q if c == 'desc' else P] += 1
                if len(ex) < 6: ex.append((x, d, P, Q, c))
    print(name, "hits excluding edges at 1:", dict(hits), "most common ancestor nodes:", smallQ.most_common(6))
    for e in ex: print("   ", e)
