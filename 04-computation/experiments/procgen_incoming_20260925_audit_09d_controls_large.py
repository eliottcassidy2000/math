import collections
from procgen_incoming_20260925_audit_09_seed_berggren_antichain import orbit_edges, comparable, U
from procgen_incoming_20260925_audit_09c_classify import cycles_of
for (q, s, name, N) in ((3, -1, "MINUS", 200001), (5, 1, "5n+1", 100001), (7, 1, "7n+1", 50001)):
    cyc = cycles_of(q, s)
    hits = collections.Counter(); spor = collections.Counter()
    for x in range(3, N + 1, 2):
        E = [e for e in orbit_edges(x, q, s, cap=300) if e[1] > 1]
        if not E: continue
        P = E[0]
        for d, Qe in enumerate(E[1:], start=1):
            c = comparable(P, Qe)
            if c:
                oncyc = (Qe[0] in cyc and Qe[1] in cyc) or (P[0] in cyc and P[1] in cyc)
                hits[(c, "cycle-edge" if oncyc else "other")] += 1
                if not oncyc: spor[(Qe if c == "desc" else P)] += 1
    print(name, "starts<=%d" % N, "cycles(min elems):", sorted(cyc)[:10], "hits:", dict(hits), "sporadic ancestor nodes:", spor.most_common(8), flush=True)
