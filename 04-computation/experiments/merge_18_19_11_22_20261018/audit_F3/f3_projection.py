#!/usr/bin/env python3
"""Audit F3, item (b), optional structured heuristic: line-projection colourings of G_q = Cay(F_(q^2), mu_(q+1)).

For q prime the proper additive subgroups of F_(q^2) are the q + 1 lines L through 0.  A colouring constant on the
cosets of L is a colouring of the circulant Cay(F_(q^2)/L, image of mu), possible iff L meets mu only in {} (i.e.
N(L - 0) consists of non-squares).  rhombus_overlap_elliptic_20261011.md G2's 5-colouring of the q = 19 residue
graph is of this type ((x, y) -> c(x + y)).  Here: for each admissible line, the exact chromatic number and
independence number of the quotient circulant (own backtracking), so
  kappa(q) <= min_L chi(circ_L)   and   alpha(G_q) >= q * max_L alpha(circ_L).
Also tests whether the session's colouring files are constant on the cosets of some line.
Usage: python3 f3_projection.py
"""
import os, sys, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)


def qnr(q):
    return next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)


def circle(q, r):
    return [(a, b) for a in range(q) for b in range(q) if (a * a - r * b * b) % q == 1]


def lines(q):
    # direction vectors (1, m) and (0, 1)
    return [(1, m) for m in range(q)] + [(0, 1)]


def functional_for_line(q, v):
    # a linear functional with kernel F_q v: l(a, b) = v1 * a - v0 * b
    return lambda a, b: (v[1] * a - v[0] * b) % q


def chromatic_circulant(q, S, kmax=12):
    """Exact chromatic number of Cay(Z_q, S) by backtracking (vertex 0 gets colour 0)."""
    S = set(s % q for s in S)
    nbr = [[(v + s) % q for s in S] for v in range(q)]
    for k in range(1, kmax + 1):
        col = [-1] * q
        order = list(range(q))

        def bt(i, used):
            if i == q:
                return True
            v = order[i]
            forb = {col[u] for u in nbr[v] if col[u] >= 0}
            for c in range(min(k, used + 1)):
                if c not in forb:
                    col[v] = c
                    if bt(i + 1, max(used, c + 1)):
                        return True
                    col[v] = -1
            return False

        if bt(0, 0):
            return k, col[:]
    return None, None


def alpha_circulant(q, S):
    S = set(s % q for s in S)
    best = 0
    # branch and bound on vertex sets containing 0 (vertex-transitive)
    adj = [set((v + s) % q for s in S) for v in range(q)]

    def rec(cands, size):
        nonlocal best
        if size + len(cands) <= best:
            return
        if not cands:
            best = max(best, size)
            return
        v = min(cands)
        rec([u for u in cands if u != v and u not in adj[v]], size + 1)
        rec([u for u in cands if u != v], size)

    rec([u for u in range(1, q) if u not in adj[0]], 1)
    return best


def main():
    out = []
    pr = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    for q in [11, 13, 17, 19, 23, 29, 31, 37]:
        r = qnr(q)
        mu = circle(q, r)
        res = []
        for v in lines(q):
            l = functional_for_line(q, v)
            S = sorted({l(a, b) for (a, b) in mu})
            if 0 in S:
                continue
            k, col = chromatic_circulant(q, S)
            al = alpha_circulant(q, S)
            res.append((k, al, v, len(S)))
        kmin = min(x[0] for x in res)
        amax = max(x[1] for x in res)
        pr('q = %2d: %d admissible lines (L cap mu empty); quotient circulant |S| = %s; chi over lines = %s; '
           'alpha over lines = %s => kappa(q) <= %d, alpha(G_q) >= %d (q^2/alpha >= %.3f)'
           % (q, len(res), sorted({x[3] for x in res}), sorted({x[0] for x in res}), sorted({x[1] for x in res}),
              kmin, q * amax, q * q / (q * amax)))
    pr()
    for name in ['kappa_19_colouring_k5.txt', 'kappa_23_colouring_k8.txt', 'kappa_29_colouring_k7.txt',
                 'kappa_31_colouring_k8.txt']:
        lines_ = open(os.path.join(PARENT, name)).read().split('\n')
        q = int(lines_[0].split('q=')[1].split()[0])
        col = [int(c) for c in lines_[1].split()]
        inv = []
        for v in lines(q):
            if all(col[a * q + b] == col[((a + v[0]) % q) * q + (b + v[1]) % q] for a in range(q) for b in range(q)):
                inv.append(v)
        sizes = sorted([col.count(c) for c in set(col)], reverse=True)
        pr('%s: invariant under translation by lines %s ; colour-class sizes %s' % (name, inv, sizes))
    open(os.path.join(HERE, 'f3_projection.out'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
