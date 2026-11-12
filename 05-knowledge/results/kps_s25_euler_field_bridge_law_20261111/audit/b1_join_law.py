"""b1_join_law.py - audit of Proposition B1 by brute force (auditor's own code, tt.py).

T = X => C => Y (C transitive, i >= 1 vertices), arc between the roots reversed (v -> u).  Claims:

    H(T)     = H(X) H(Y) + 2^i W+_X(u) W-_Y(v)
    W+_T(v)  = H(X) W+_Y(v) + 2^i W+_X(u) tau(Y - v)

Tests
  J1  ALL pairs of labeled rooted tournaments with at most 4 vertices each (285 x 285), i = 1
  J2  random pairs (orders 1..4, random root), i = 2 and i = 3: 3000 cases each
  J3  random pairs with orders up to 5, i = 1..3: 600 cases
  J4  W+_X(u) >= H(X) and W-_X(u) >= H(X) for ALL labeled rooted tournaments with at most 6 vertices; equality cases
  J5  i = 0 (no chain): the formula of B1 is FALSE in general (the parity split is needed) - count of failures
"""
import random

from tt import tables, invariants, ordinal_sum, reverse_arc, transitive, all_tournaments


def rooted_data(n, adj):
    tabs = tables(n, adj)
    return [invariants(n, adj, u, tabs) for u in range(n)]


def check(X, u, invX, Y, v, invY, i):
    nX, aX = X
    nY, aY = Y
    blocks = [(nX, aX)] + ([transitive(i)] if i >= 1 else []) + [(nY, aY)]
    n, adj, offs = ordinal_sum(blocks)
    gu = offs[0] + u
    gv = offs[-1] + v
    reverse_arc(adj, gu, gv)
    tabs = tables(n, adj)
    invT = invariants(n, adj, gv, tabs)
    predH = invX["H"] * invY["H"] + (1 << i) * invX["Wp"] * invY["Wm"]
    predW = invX["H"] * invY["Wp"] + (1 << i) * invX["Wp"] * invY["tau_minus_root"]
    return invT["H"] == predH, invT["Wp"] == predW


def main():
    random.seed(20261111)
    pool = {}
    for n in range(1, 6):
        pool[n] = []
        for adj in all_tournaments(n):
            pool[n].append((n, adj, rooted_data(n, adj)))
    small = [t for n in range(1, 5) for t in pool[n]]
    rooted = [(n, adj, u, inv[u]) for (n, adj, inv) in small for u in range(n)]
    print("rooted labeled tournaments with at most 4 vertices:", len(rooted))

    tot = badH = badW = 0
    for (nX, aX, u, iX) in rooted:
        for (nY, aY, v, iY) in rooted:
            okH, okW = check((nX, aX), u, iX, (nY, aY), v, iY, 1)
            tot += 1
            badH += (not okH)
            badW += (not okW)
    print("J1  all pairs, i = 1: %d cases; H formula failures: %d; W+ formula failures: %d" % (tot, badH, badW))

    for i in (2, 3):
        tot = badH = badW = 0
        for _ in range(3000):
            (nX, aX, u, iX) = random.choice(rooted)
            (nY, aY, v, iY) = random.choice(rooted)
            okH, okW = check((nX, aX), u, iX, (nY, aY), v, iY, i)
            tot += 1
            badH += (not okH)
            badW += (not okW)
        print("J2  random pairs, i = %d: %d cases; H formula failures: %d; W+ formula failures: %d" % (i, tot, badH, badW))

    tot = badH = badW = 0
    for _ in range(600):
        nX = random.randint(3, 5)
        nY = random.randint(3, 5)
        (_, aX, invX) = random.choice(pool[nX])
        (_, aY, invY) = random.choice(pool[nY])
        u = random.randrange(nX)
        v = random.randrange(nY)
        i = random.randint(1, 3)
        okH, okW = check((nX, aX), u, invX[u], (nY, aY), v, invY[v], i)
        tot += 1
        badH += (not okH)
        badW += (not okW)
    print("J3  random pairs, orders 3..5, i = 1..3: %d cases; H formula failures: %d; W+ formula failures: %d" % (tot, badH, badW))

    # J4
    tot = bad = eq = eq_nonsink = 0
    for n in range(1, 7):
        for adj in all_tournaments(n):
            tabs = tables(n, adj)
            for u in range(n):
                inv = invariants(n, adj, u, tabs)
                tot += 1
                if inv["Wp"] < inv["H"] or inv["Wm"] < inv["H"]:
                    bad += 1
                if inv["Wp"] == inv["H"]:
                    eq += 1
                    if bin(adj[u]).count("1") != 0:
                        eq_nonsink += 1
    print("J4  W+_X(u) >= H(X) and W-_X(u) >= H(X) on all labeled rooted tournaments with at most 6 vertices: %d cases, %d violations; W+ = H in %d cases (%d of them with a root that is not a sink, e.g. any vertex of the 3-cycle)" % (tot, bad, eq, eq_nonsink))

    # J5: without a chain the 2 x 2 law fails
    tot = bad = 0
    for (nX, aX, u, iX) in rooted[:120]:
        for (nY, aY, v, iY) in rooted[:120]:
            okH, okW = check((nX, aX), u, iX, (nY, aY), v, iY, 0)
            tot += 1
            bad += (not okH)
    print("J5  i = 0 (direct join) with the same formula: %d cases, %d failures of the H formula (expected: many; the parity split is needed)" % (tot, bad))


if __name__ == "__main__":
    main()
