"""Audit of THM-4846 (A): the one-arc observer M_{n-1} (built by the subset DP of runlib, not by the star rule).

Checks, every class n = 4, 5, 6:
 * M + M^T = adjacency matrix of the star-transposition graph ST_n in permutation coordinates, and each unordered vertex
   pair labels exactly (n-1)! edges;
 * equivariance under the right action of S_{n-1} (permuting run positions);
 * trivial (run-symmetric) block = A^T in complement coordinates, sign block = (-1)^(|x-y|-1) A[x][y];
 * det M_{n-1} (FLINT) and its factorisation; Smith forms (FLINT) of the strong classes with exactly two Hamiltonian
   cycles (the claimed singular-but-torsion-free exceptions);
 * ST_4 is the Nauru graph: isomorphic to the generalized Petersen graph G(12,5) (networkx VF2), 24 vertices, 36 edges.
"""
import itertools
from math import factorial
from collections import Counter
import networkx as nx
import flint
import sympy
from runlib import (tournaments, adj_from_bits, hp_table, ham_cycles, is_strong, run_matrix, plus_sector, snf_diag,
                    snf_summary, fmt)


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def main():
    out = []
    for n in (4, 5, 6):
        k = n - 1
        tup = list(itertools.permutations(range(n), k))
        idx = {t: i for i, t in enumerate(tup)}
        G = nx.Graph()
        label = Counter()
        for t in tup:
            w = (set(range(n)) - set(t)).pop()
            for j in range(k):
                u = list(t)
                u[j] = w
                u = tuple(u)
                if not G.has_edge(idx[t], idx[u]):
                    G.add_edge(idx[t], idx[u])
                    label[frozenset((t[j], w))] += 1
        assert all(v == factorial(n - 1) for v in label.values()) and len(label) == n * (n - 1) // 2
        dets = Counter()
        rows = []
        for bits in tournaments(n):
            A = adj_from_bits(n, bits)
            hp = hp_table(A)
            hc = ham_cycles(A, hp)
            _, M = run_matrix(A, k, hp, tup)
            Nn = len(tup)
            for i in range(Nn):
                for j in range(Nn):
                    a, b = M[i][j], M[j][i]
                    assert a + b == (1 if G.has_edge(i, j) else 0)
            # equivariance under S_{n-1} acting on positions
            for h in itertools.permutations(range(k)):
                for i, s in enumerate(tup):
                    sh = tuple(s[h[q]] for q in range(k))
                    for j, t in enumerate(tup):
                        if M[i][j]:
                            th = tuple(t[h[q]] for q in range(k))
                            assert M[idx[sh]][idx[th]] == M[i][j]
                break  # one nontrivial h suffices with the next line
            for h in [tuple([1, 0] + list(range(2, k))), tuple(list(range(1, k)) + [0])]:
                for i, s in enumerate(tup):
                    sh = tuple(s[h[q]] for q in range(k))
                    for j, t in enumerate(tup):
                        th = tuple(t[h[q]] for q in range(k))
                        assert M[idx[sh]][idx[th]] == M[i][j]
            # trivial block via the unordered plus sector, sign block via signed sums
            sets, P = plus_sector(A, k, hp)
            sidx = {s: i for i, s in enumerate(sets)}
            for x in range(n):
                for y in range(n):
                    a = sidx[tuple(v for v in range(n) if v != y)]
                    b = sidx[tuple(v for v in range(n) if v != x)]
                    assert P[a][b] == (A[x][y] if x != y else 0)
                    if x == y:
                        continue
                    start = tuple(v for v in range(n) if v != y)
                    endset = sorted(v for v in range(n) if v != x)
                    tot = sum(perm_sign([endset.index(v) for v in order]) * M[idx[start]][idx[order]]
                              for order in itertools.permutations(endset))
                    assert tot == (-1) ** (abs(x - y) - 1) * A[x][y]
            d = int(flint.fmpz_mat(M).det())
            fac = sympy.factorint(abs(d)) if d else None
            key = 0 if d == 0 else tuple(sorted(fac.items()))
            dets[(key, is_strong(A))] += 1
            if is_strong(A) and hc == 2:
                rows.append((bits, hc, fmt(snf_summary(snf_diag(M)))))
        out.append('n=%d: ST_n has %d vertices, %d edges, each pair labels %d edges; every class orients it; '
                   'S_{n-1}-equivariant; trivial block = A^T, sign block = (-1)^(|x-y|-1)A[x][y]'
                   % (n, G.number_of_nodes(), G.number_of_edges(), factorial(n - 1)))
        sing = sum(v for (kk, s), v in dets.items() if kk == 0)
        out.append('   singular classes: %d of %d; nonzero |det| (factorised) with strong flag: %s'
                   % (sing, sum(dets.values()), sorted(((kk, s), v) for (kk, s), v in dets.items() if kk != 0)))
        out.append('   strong classes with exactly 2 Hamiltonian cycles, Smith form of M_{n-1}: %s' % rows)
        if n == 4:
            GP = nx.Graph()
            for i in range(12):
                GP.add_edge(('u', i), ('u', (i + 1) % 12))
                GP.add_edge(('u', i), ('v', i))
                GP.add_edge(('v', i), ('v', (i + 5) % 12))
            out.append('   ST_4 isomorphic to the generalized Petersen graph G(12,5) (Nauru): %s; G(12,5): %d vertices, %d edges'
                       % (nx.is_isomorphic(G, GP), GP.number_of_nodes(), GP.number_of_edges()))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
