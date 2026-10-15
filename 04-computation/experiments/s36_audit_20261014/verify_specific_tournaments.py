#!/usr/bin/env python3
"""S36 referee audit, part B: brute force over ALL n! vertex orders (pure Python, exact integers) for specific tournaments.

(1) The two n = 9 counterexamples to the stochastic symmetry (SS) claimed in HYP-9280 (vertex w = 8, in-degree 7).
(2) The HYP-9277 Claim A counterexamples (n = 8, 9, 10) quoted in HYP-9280 ('a single inverted arc of probability
    13/45 and 17/81, and M(T) has 10 and 12 three-cycles').
(3) The extremal tournament T_n (transitive + reversed extreme arc) for n = 3..9: closed forms H = 2^(n-2) + 1,
    best 3-cycle certificate, vertex lemma value, min D_in, max Pr[w first].
Encoding: adj[v] = bitmask of out-neighbours of v.
"""
import itertools
from fractions import Fraction as Fr


def ham_paths(adj):
    n = len(adj)
    out = []
    for perm in itertools.permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if not (adj[perm[i]] >> perm[i + 1]) & 1:
                ok = False
                break
        if ok:
            out.append(perm)
    return out


def is_tournament(adj):
    n = len(adj)
    for x in range(n):
        if (adj[x] >> x) & 1:
            return False
        for y in range(x + 1, n):
            if ((adj[x] >> y) & 1) + ((adj[y] >> x) & 1) != 1:
                return False
    return True


def strongly_connected(n, arcs):
    """arcs: set of (x, y); Kosaraju-free check: reach from 0 forward and backward (BFS on explicit lists)"""
    fw = {v: [y for (x, y) in arcs if x == v] for v in range(n)}
    bw = {v: [x for (x, y) in arcs if y == v] for v in range(n)}
    for nb in (fw, bw):
        seen = {0}
        st = [0]
        while st:
            v = st.pop()
            for u in nb[v]:
                if u not in seen:
                    seen.add(u)
                    st.append(u)
        if len(seen) != n:
            return False
    return True


def analyse(adj, w_ss=None, label=''):
    n = len(adj)
    assert is_tournament(adj), 'not a tournament'
    P = ham_paths(adj)
    H = len(P)
    pos = [dict((v, i) for i, v in enumerate(p)) for p in P]
    C = [[sum(1 for q in pos if q[x] < q[y]) for y in range(n)] for x in range(n)]
    arcsT = {(x, y) for x in range(n) for y in range(n) if x != y and (adj[x] >> y) & 1}
    strongT = strongly_connected(n, arcsT)
    M = {(x, y) for x in range(n) for y in range(n) if x != y and 2 * C[x][y] > H}
    Fw = arcsT & M
    inv = sorted((x, y, Fr(C[x][y], H)) for (x, y) in arcsT if 2 * C[x][y] < H)
    ms = [sum(1 for y in range(n) if (x, y) in M) for x in range(n)]
    m3 = n * (n - 1) * (n - 2) // 6 - sum(s * (s - 1) // 2 for s in ms)
    fwd3 = any((x, y) in Fw and (y, z) in Fw and (z, x) in Fw for x in range(n) for y in range(n) for z in range(n)
               if len({x, y, z}) == 3)
    print(f'--- {label} adj = {adj}')
    print(f'    n = {n}, H = {H} (odd: {H % 2 == 1}), T strong: {strongT}, F strong: {strongly_connected(n, Fw)}, '
          f'M transitive: {m3 == 0}, #3-cycles of M = {m3}, forward 3-cycle of T: {fwd3}')
    print(f'    inverted arcs (tail, head, p(tail first)): {[(x, y, str(p)) for (x, y, p) in inv]}')
    if w_ss is not None:
        w = w_ss
        I = [u for u in range(n) if u != w and (adj[u] >> w) & 1]
        d = len(I)
        dist = [0] * (d + 1)
        for q in pos:
            dist[sum(1 for u in I if q[u] < q[w])] += 1
        EN = Fr(sum(k * dist[k] for k in range(d + 1)), H)
        viol = [k for k in range(d + 1) if sum(dist[d - k:]) < sum(dist[:k + 1])]
        best = max(Fr(C[u][w], H) for u in I)
        print(f'    vertex w = {w}: in-degree d = {d}, distribution of N = {dist}, E[N] = {EN} (d/2 = {Fr(d, 2)}), '
              f'best in-neighbour p = {best}, SS violated at k = {viol}')
        if viol:
            k = viol[0]
            print(f'      Pr[N >= {d - k}] = {sum(dist[d - k:])}/{H}  vs  Pr[N <= {k}] = {sum(dist[:k + 1])}/{H}')
    return H, C


def Tn(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    adj[0] &= ~(1 << (n - 1))
    adj[n - 1] |= 1
    return adj


if __name__ == '__main__':
    print('(1) SS counterexamples at n = 9 (HYP-9280)')
    analyse([502, 508, 504, 241, 480, 448, 384, 256, 8], w_ss=8, label='SS-1')
    analyse([250, 508, 497, 500, 480, 448, 384, 256, 1], w_ss=8, label='SS-2')
    print('\n(2) HYP-9277 Claim A counterexamples')
    analyse([190, 176, 250, 242, 64, 144, 163, 16], label='ClaimA n=8')
    analyse([460, 493, 328, 256, 367, 461, 264, 92, 128], label='ClaimA n=9')
    print('\n(3) T_n closed forms')
    for n in range(3, 10):
        adj = Tn(n)
        P = ham_paths(adj)
        H = len(P)
        pos = [dict((v, i) for i, v in enumerate(p)) for p in P]
        C = [[sum(1 for q in pos if q[x] < q[y]) for y in range(n)] for x in range(n)]
        cyc = [min(C[x][y], C[y][z], C[z][x]) for x in range(n) for y in range(n) for z in range(n)
               if len({x, y, z}) == 3 and (adj[x] >> y) & 1 and (adj[y] >> z) & 1 and (adj[z] >> x) & 1]
        vl = min(max(C[u][w] for u in range(n) if u != w and (adj[u] >> w) & 1) for w in range(n))
        din = min(sum(2 * C[u][w] - H for u in range(n) if u != w and (adj[u] >> w) & 1) for w in range(n))
        fm = max(sum(1 for p in P if p[0] == w) for w in range(n))
        print(f'    n = {n}: H = {H} (2^(n-2)+1 = {2 ** (n - 2) + 1}); cert = {Fr(max(cyc), H)}; vertex lemma = {Fr(vl, H)}; '
              f'min D_in = {Fr(din, H)}; max Pr[first] = {Fr(fm, H)}; formula {Fr(2 ** (n - 3) + 1, 2 ** (n - 2) + 1)}, '
              f'{Fr(1, 2 ** (n - 2) + 1)}, {Fr(2 ** (n - 3), 2 ** (n - 2) + 1)}; p(last->first arc) = {Fr(C[n - 1][0], H)}')
