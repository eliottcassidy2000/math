"""cd_lines_aut.py - (i) the line rule for every generalized Cayley-Dickson algebra of dimension <= 32;
(ii) automorphism group orders of the sign tournaments CD_3, CD_4 and of the type A example; converse test;
(iii) block recursion G_(n+1) = [[G, G^T], [G D, -G^T D]] of the sign table, n <= 4."""
import itertools, sys
from cd_more import table, c3_of
from cd_signs import sign_table


def line_rule(gs):
    n = len(gs)
    N = 1 << n
    g = table(list(gs))
    V = range(1, N)
    sq = [g[i][i] for i in range(N)]
    if not all(g[i][j] == -g[j][i] for i in V for j in V if i != j):
        return False, 0
    arc = lambda i, j: g[i][j] == 1
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    ok = True
    cyc = 0
    for L in lines:
        neg = [x for x in L if sq[x] == -1]
        is_c = c3_of(arc, L) == 1
        cyc += is_c
        if len(neg) == 3:
            ok &= is_c
        elif len(neg) == 1:
            m = neg[0]
            o = [x for x in L if x != m]
            ok &= (not is_c) and ((arc(o[0], m) and arc(m, o[1])) or (arc(o[1], m) and arc(m, o[0])))
        else:
            ok = False
    # the units of square -1 together with 0 form a subgroup
    negset = {0} | {i for i in V if sq[i] == -1}
    ok &= all((a ^ b) in negset for a in negset for b in negset)
    return ok, cyc


def automorphisms(adj):
    n = len(adj)
    out_sig = [tuple(sorted(sum(adj[v][w] and adj[w][u] for w in range(n)) for u in range(n) if adj[v][u])) for v in range(n)]
    count = 0
    perm = [-1] * n
    used = [False] * n

    def rec(k):
        nonlocal count
        if k == n:
            count += 1
            return
        for img in range(n):
            if used[img] or out_sig[img] != out_sig[k]:
                continue
            if all(adj[j][k] == adj[perm[j]][img] for j in range(k)):
                perm[k] = img
                used[img] = True
                rec(k + 1)
                used[img] = False
        perm[k] = -1
    rec(0)
    return count


def canon_key(adj):
    n = len(adj)
    return tuple(sorted((sum(adj[v]), tuple(sorted(sum(adj[v][w] and adj[w][u] for w in range(n)) for u in range(n) if adj[v][u]))) for v in range(n)))


if __name__ == "__main__":
    for n in (2, 3, 4, 5):
        res = {}
        for gs in itertools.product((-1, 1), repeat=n):
            ok, cyc = line_rule(gs)
            res.setdefault((ok, cyc), 0)
            res[(ok, cyc)] += 1
        print("dimension %d: (line rule holds, number of cyclic lines) -> number of parameter choices: %s" % (1 << n, sorted(res.items())))
    for n in (2, 3, 4):
        g = sign_table(n)
        g2 = sign_table(n + 1)
        N = 1 << n
        D = [1] + [-1] * (N - 1)
        ok = True
        for i in range(N):
            for j in range(N):
                ok &= g2[i][j] == g[i][j] and g2[i][j + N] == g[j][i] and g2[i + N][j] == g[i][j] * D[j] and g2[i + N][j + N] == -g[j][i] * D[j]
        print("block recursion for the sign table, n = %d -> %d: %s" % (n, n + 1, ok))
    for n in (3, 4):
        g = sign_table(n)
        N = 1 << n
        adj = [[1 if (i != j and g[i + 1][j + 1] == 1) else 0 for j in range(N - 1)] for i in range(N - 1)]
        conv = [[adj[j][i] for j in range(N - 1)] for i in range(N - 1)]
        print("CD_%d: automorphisms %d; same vertex invariants as its converse: %s" % (n, automorphisms(adj), canon_key(adj) == canon_key(conv)))
