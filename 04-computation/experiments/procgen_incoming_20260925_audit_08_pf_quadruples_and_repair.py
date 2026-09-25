#!/usr/bin/env python3
"""(1) Link check: THM-4472 AM-fair quadruples vs codex's four-core Pfaffian switching
invariant.  For T_b (b=+-1), odd s_o with |s_o|>=5, s_e=s_o+b, V={s_o,s_e,T(s_o),T(s_e)},
map arcs x->T(x), other pairs smaller->larger.  Record (b*s_o sign, scores, |Pf|, H).
(2) Independent audit of decoder_pair_repair: minimum number of arc reversals making
some perfect matching of Q[C3,C3,C3,1] a family of pair modules, over all 64 cores.
Claimed distribution {3:24, 6:24, 7:12, 8:4}.
"""
from itertools import permutations, combinations


def T(n, b):
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def pf_and_H(verts, arcs):
    idx = {v: i for i, v in enumerate(verts)}
    S = [[0] * 4 for _ in range(4)]
    for (x, y) in arcs:
        S[idx[x]][idx[y]] = 1
        S[idx[y]][idx[x]] = -1
    pf = S[0][1] * S[2][3] - S[0][2] * S[1][3] + S[0][3] * S[1][2]
    H = sum(1 for p in permutations(range(4)) if all(S[p[i]][p[i + 1]] == 1 for i in range(3)))
    scores = tuple(sorted(sum(1 for j in range(4) if S[i][j] == 1) for i in range(4)))
    return abs(pf), H, scores


def quadruples():
    table = {}
    for b in (1, -1):
        for so in list(range(-301, -4, 2)) + list(range(5, 302, 2)):
            se = so + b
            V = [so, se, T(so, b), T(se, b)]
            if len(set(V)) < 4:
                continue
            arcs = set()
            maparcs = {(so, T(so, b)), (se, T(se, b))}
            for (x, y) in maparcs:
                arcs.add((x, y))
            for x, y in combinations(V, 2):
                if (x, y) in maparcs or (y, x) in maparcs:
                    continue
                arcs.add((min(x, y), max(x, y)))
            apf, H, sc = pf_and_H(V, arcs)
            key = (b, 1 if b * so > 0 else -1, sc, apf, H)
            table[key] = table.get(key, 0) + 1
    for k, v in sorted(table.items()):
        print("  b=%+d sgn(b*s_o)=%+d scores=%s |Pf|=%d H=%d : %d quadruples" % (k + (v,)))
    ok = all(k[3] == k[4] == 2 + k[1] for k in table)
    print("  |Pf| = H = 2+sgn(b s_o) on all quadruples:", ok)


def repair_census():
    # vertices: triangles {0,1,2},{3,4,5},{6,7,8} cyclic, singleton 9; core vertex i -> block i
    block = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3]
    cyc = {(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (6, 7), (7, 8), (8, 6)}
    core_pairs = list(combinations(range(4), 2))
    # all perfect matchings of 10 vertices
    def matchings(rem):
        if not rem:
            yield []
            return
        a = rem[0]
        for i in range(1, len(rem)):
            b = rem[i]
            rest = rem[1:i] + rem[i + 1:]
            for m in matchings(rest):
                yield [(a, b)] + m
    Ms = list(matchings(list(range(10))))
    assert len(Ms) == 945
    dist = {}
    for mask in range(64):
        core = {}
        for i, (a, b) in enumerate(core_pairs):
            core[(a, b)] = 1 if (mask >> i) & 1 else -1
            core[(b, a)] = -core[(a, b)]
        def arc(u, v):  # +1 if u->v
            bu, bv = block[u], block[v]
            if bu == bv:
                return 1 if (u, v) in cyc else -1
            return core[(bu, bv)]
        best = None
        for M in Ms:
            cost = 0
            for (p, r) in combinations(range(5), 2):
                P, R = M[p], M[r]
                k = sum(1 for u in P for v in R if arc(u, v) == 1)
                cost += min(k, 4 - k)
            if best is None or cost < best:
                best = cost
        dist[best] = dist.get(best, 0) + 1
    print("  minimum reversal distribution over 64 cores:", dict(sorted(dist.items())))


if __name__ == "__main__":
    print("(1) THM-4472 quadruples and the Pfaffian switching class:")
    quadruples()
    print("(2) pair-repair census:")
    repair_census()
