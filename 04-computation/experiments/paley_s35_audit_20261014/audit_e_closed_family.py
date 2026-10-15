#!/usr/bin/env python3
"""Audit (e2) of S35 / HYP-9277: closed families HP(T) cap Ext(Q), their balance, and the gold-partition test;
Paley P_19 by a subset DP; and whether strong 8-tournaments reach delta = 4/11 (the session's annealing reports
nothing below 64/161).

Definitions as stated in HYP-9277 / the atlas note (re-implemented, nothing imported):
  family F(D) = { HP(D) cap Ext(Q) : Q a poset }, generated from HP(D) by intersecting with half-spaces {x before y};
  delta(S) = max_pairs min(p, 1-p);
  gold(S): exists a first comparison (a,b) such that for each nonempty answer S1 (t1 = |S1| >= 2) there is a second
           comparison (c,d) with |S| >= t1 + max(|S1 cap {c<d}|, |S1 cap {d<c}|)   (Peczarski's weak t0 >= t1 + t2,
           adaptive second comparison; any pairs allowed).
Isomorphism classes by orbit marking (each class counted once, as in the session).
"""
import itertools, time
from fractions import Fraction
import numpy as np
import numba as nb


def ham_paths(n, arc):
    out = []

    def rec(path, used):
        if len(path) == n:
            out.append(tuple(path))
            return
        u = path[-1]
        for v in range(n):
            if not used >> v & 1 and arc[u][v]:
                path.append(v)
                rec(path, used | 1 << v)
                path.pop()
    for s in range(n):
        rec([s], 1 << s)
    return out


def family_stats(n, arc):
    HP = ham_paths(n, arc)
    if len(HP) < 2:
        return 0, Fraction(1, 2), 0, None
    pos = [[0] * n for _ in HP]
    for k, p in enumerate(HP):
        for i, v in enumerate(p):
            pos[k][v] = i
    h = {}
    for x in range(n):
        for y in range(n):
            if x != y:
                h[(x, y)] = sum(1 << k for k in range(len(HP)) if pos[k][x] < pos[k][y])
    full = (1 << len(HP)) - 1
    fam = {full}
    stack = [full]
    while stack:
        S = stack.pop()
        for m in h.values():
            T = S & m
            if T and T not in fam:
                fam.add(T)
                stack.append(T)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    members = 0
    worst = Fraction(1, 2)
    worst_S = None
    gold_fail = 0
    for S in fam:
        t0 = bin(S).count('1')
        if t0 < 2:
            continue
        members += 1
        d = max(Fraction(min(bin(S & h[(x, y)]).count('1'), bin(S & h[(y, x)]).count('1')), t0) for x, y in pairs)
        if d < worst:
            worst, worst_S = d, S
        ok_any = False
        for a, b in pairs:
            ok = True
            for S1 in (S & h[(a, b)], S & h[(b, a)]):
                t1 = bin(S1).count('1')
                if t1 == 0 or t1 == 1:
                    continue
                if not any(t0 >= t1 + max(bin(S1 & h[(c, d_)]).count('1'), bin(S1 & h[(d_, c)]).count('1')) for c, d_ in pairs):
                    ok = False
                    break
            if ok:
                ok_any = True
                break
        if not ok_any:
            gold_fail += 1
    ws = None if worst_S is None else [HP[k] for k in range(len(HP)) if worst_S >> k & 1]
    return members, worst, gold_fail, ws


def classes(n, kind):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    base = 2 if kind == 'tournament' else 3
    perms = list(itertools.permutations(range(n)))
    pidx = {p: e for e, p in enumerate(pairs)}
    seen = set()
    for code in range(base ** len(pairs)):
        if code in seen:
            continue
        digits = []
        c = code
        for _ in pairs:
            digits.append(c % base)
            c //= base
        arc = [[False] * n for _ in range(n)]
        for e, (i, j) in enumerate(pairs):
            dgt = digits[e]
            if dgt in (0, 2):
                arc[i][j] = True
            if dgt in (1, 2):
                arc[j][i] = True
        for p in perms:
            nd = [0] * len(pairs)
            for e, (i, j) in enumerate(pairs):
                a, b = p[i], p[j]
                dgt = digits[e]
                if a < b:
                    nd[pidx[(a, b)]] = dgt
                else:
                    nd[pidx[(b, a)]] = {0: 1, 1: 0, 2: 2}[dgt]
            seen.add(sum(dd * base ** e for e, dd in enumerate(nd)))
        yield code, arc


@nb.njit(cache=True)
def paley_dp(q, adj, constrained):
    N = 1 << q
    f = np.zeros((N, q), dtype=np.int64)
    for v in range(q):
        if constrained and v == 1:
            continue
        f[1 << v, v] = 1
    for S in range(1, N):
        for v in range(q):
            c = f[S, v]
            if c == 0:
                continue
            for w in range(q):
                if (S >> w) & 1 or not (adj[v] >> w) & 1:
                    continue
                if constrained and w == 1 and not (S & 1):
                    continue
                f[S | (1 << w), w] += c
    tot = 0
    for v in range(q):
        tot += f[N - 1, v]
    return tot


def main():
    t0 = time.time()
    # 1. the stated n = 6 witness (session code tuple: 0 -> i->j, 1 -> j->i)
    n = 6
    code = (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    arc = [[False] * n for _ in range(n)]
    for e, (i, j) in enumerate(pairs):
        if code[e] == 0:
            arc[i][j] = True
        else:
            arc[j][i] = True
    HP = ham_paths(n, arc)
    W = [(0, 1, 5, 4, 2, 3), (4, 2, 5, 0, 1, 3), (5, 0, 1, 3, 4, 2), (5, 0, 1, 4, 2, 3)]
    rel = [(x, y) for x in range(n) for y in range(n) if x != y and all(p.index(x) < p.index(y) for p in W)]
    cons = [p for p in HP if all(p.index(x) < p.index(y) for x, y in rel)]
    splits = sorted(set(min(sum(p.index(x) < p.index(y) for p in W), 4 - sum(p.index(x) < p.index(y) for p in W)) for x, y in pairs))
    print(f'witness: H(T) = {len(HP)}; the 4 listed sequences are Hamiltonian paths: {all(w in HP for w in W)}; '
          f'common precedences {rel}; HP(T) cap Ext(Q_min) = listed set: {sorted(cons) == sorted(W)}; pair splits min(c, 4-c) in {splits} -> delta = 1/4')
    # 2. full families
    for kind, nmax in (('tournament', 6), ('semicomplete', 4)):
        for n in range(3, nmax + 1):
            t1 = time.time()
            ncls = mem = gf = 0
            worst = Fraction(1, 2)
            for code, arc in classes(n, kind):
                ncls += 1
                m, w, g, ws = family_stats(n, arc)
                mem += m
                gf += g
                if w < worst:
                    worst = w
            print(f'{kind} n = {n}: {ncls} classes, {mem} family members with |S| >= 2, min delta {worst}, gold failures {gf}  [{time.time() - t1:.0f} s]', flush=True)
    # 3. Paley P_19
    for q in (7, 11, 19):
        sq = {(x * x) % q for x in range(1, q)}
        adj = np.array([sum(1 << y for y in range(q) if (y - x) % q in sq) for x in range(q)], dtype=np.int64)
        H = int(paley_dp(q, adj, False))
        H01 = int(paley_dp(q, adj, True))
        print(f'Paley P_{q} (subset DP): H = {H}, #(0 before 1) = {H01}, p = {Fraction(H01, H)}', flush=True)
    print(f'[{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
