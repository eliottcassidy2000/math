#!/usr/bin/env python3
"""Seed pre-screen: extend the four-vertex AM-fair quadruple (THM-4472) to the
six-vertex 'two-step' tournament V={s_o,s_e,T s_o,T s_e,T^2 s_o,T^2 s_e}
(map arcs x->T(x) for the first two steps; every other pair smaller->larger),
and ask whether its switching invariant |Pf| (6x6 skew sign matrix) and H mod 4
are functions of local data only (sign of b*s_o and the parity pattern of the
images), for b=+1,-1 and both sides of 0.  If yes, the invariant is a repackaged
sign law (no new input).  Also records H (Hamiltonian paths)."""
from itertools import permutations
from collections import defaultdict


def T(n, b):
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def pfaffian(S):
    n = len(S)
    if n == 0:
        return 1
    tot = 0
    for j in range(1, n):
        if S[0][j] == 0:
            continue
        rest = [k for k in range(1, n) if k != j]
        sub = [[S[a][c] for c in rest] for a in rest]
        sign = (-1) ** (j - 1)
        tot += sign * S[0][j] * pfaffian(sub)
    return tot


def ham_paths(S):
    n = len(S)
    return sum(1 for p in permutations(range(n)) if all(S[p[i]][p[i + 1]] == 1 for i in range(n - 1)))


def build(b, so):
    se = so + b
    xs = [so, se, T(so, b), T(se, b), T(T(so, b), b), T(T(se, b), b)]
    if len(set(xs)) < 6:
        return None
    maparcs = {(so, xs[2]), (se, xs[3]), (xs[2], xs[4]), (xs[3], xs[5])}
    idx = {v: i for i, v in enumerate(xs)}
    S = [[0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            x, y = xs[i], xs[j]
            if (x, y) in maparcs:
                a, c = x, y
            elif (y, x) in maparcs:
                a, c = y, x
            else:
                a, c = (x, y) if x < y else (y, x)
            S[idx[a]][idx[c]] = 1
            S[idx[c]][idx[a]] = -1
    return xs, S


def main():
    table = defaultdict(lambda: defaultdict(int))
    for b in (1, -1):
        for so in list(range(-4001, -8, 2)) + list(range(9, 4002, 2)):
            r = build(b, so)
            if r is None:
                continue
            xs, S = r
            pf = abs(pfaffian(S))
            H = ham_paths(S)
            # local data: side sign, parities of first images
            local = (1 if b * so > 0 else -1, xs[2] % 2, xs[3] % 2)
            table[(b,) + local][(pf, H % 4)] += 1
    for k in sorted(table):
        print("b=%+d sgn(b s_o)=%+d parities(T s_o,T s_e)=(%d,%d):" % k, dict(table[k]))
    det = all(len(v) == 1 for v in table.values())
    print("(|Pf|, H mod 4) determined by (b, sgn(b s_o), image parities):", det)


if __name__ == "__main__":
    main()
