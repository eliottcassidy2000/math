#!/usr/bin/env python3
"""Referee checks of assorted exact claims: Paley clique numbers (own bitset branch-and-bound, cross-checked with
networkx where cheap), primary pseudoperfect numbers below 60000, Fermat product, Sylvester, 1807, 139, (2,3,k),
tau(n) = 4 abundancy, 2^32 - 1, collector denominators, mod-21 residues of 1093 and 3511, C(14,3)."""
import math, itertools
from fractions import Fraction as Fr


def paley_adj(q):
    sq = {(x * x) % q for x in range(1, q)}
    return [sum(1 << j for j in range(q) if j != i and (j - i) % q in sq) for i in range(q)], sq


def max_clique(adj, cand_mask):
    """maximum clique size within cand_mask (bitsets), greedy-colouring bound branch and bound"""
    best = [0]

    def colour_bound(P):
        # greedy sequential colouring; returns list of (vertex, colour) in colour order
        order = []
        col = 0
        U = P
        while U:
            col += 1
            Q = U
            while Q:
                v = (Q & -Q).bit_length() - 1
                Q &= ~(1 << v)
                Q &= ~adj[v]
                U &= ~(1 << v)
                order.append((v, col))
        return order

    def expand(size, P):
        order = colour_bound(P)
        for v, c in reversed(order):
            if size + c <= best[0]:
                return
            newP = P & adj[v]
            if newP == 0:
                if size + 1 > best[0]:
                    best[0] = size + 1
            else:
                expand(size + 1, newP)
            P &= ~(1 << v)

    expand(0, cand_mask)
    return best[0]


def paley_clique_number(q):
    adj, sq = paley_adj(q)
    # vertex-transitive: a maximum clique can be taken to contain 0; search inside N(0)
    return 1 + max_clique(adj, adj[0])


out = []
for q in (5, 13, 17, 29, 37, 41, 257):
    w = paley_clique_number(q)
    out.append(f'P_{q}: clique number {w}')
try:
    import networkx as nx
    for q in (13, 17, 29, 37, 41):
        adj, sq = paley_adj(q)
        G = nx.Graph()
        G.add_nodes_from(range(q))
        for i in range(q):
            for j in range(i + 1, q):
                if (j - i) % q in sq:
                    G.add_edge(i, j)
        w = max(len(c) for c in nx.find_cliques(G))
        out.append(f'   networkx cross-check P_{q}: {w}')
except ImportError:
    out.append('networkx unavailable')

# primary pseudoperfect: sum_{p | n} 1/p + 1/n = 1, n squarefree (by definition) -- own factorisation
def factor(n):
    f = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                m //= d
                e += 1
            f.append((d, e))
        d += 1
    if m > 1:
        f.append((m, 1))
    return f


pp = []
for n in range(2, 60000):
    f = factor(n)
    if sum(Fr(1, p) for p, e in f) + Fr(1, n) == 1:
        pp.append((n, all(e == 1 for p, e in f)))
out.append(f'primary pseudoperfect (sum 1/p + 1/n = 1) below 60000: {pp}')
out.append(f'3*5*17*257*65537 = {3*5*17*257*65537}, 2^32 - 1 = {2**32 - 1}')
for k in range(5):
    Fk = 2 ** (2 ** k) + 1
    d = next(e for e in range(1, Fk) if pow(2, e, Fk) == 1)
    prod = math.prod(2 ** (2 ** j) + 1 for j in range(k + 1))
    out.append(f'F_{k} = {Fk}: ord(2) = {d} = 2^{k+1}: {d == 2 ** (k + 1)}; 2^d - 1 = F_0...F_k: {2 ** d - 1 == prod}')
s = [2]
for _ in range(4):
    s.append(s[-1] ** 2 - s[-1] + 1)
out.append(f'Sylvester: {s}; 1807 = {factor(1807)}; 3^7 - 2^11 = {3**7 - 2**11}; Heegner numbers among them: '
           f'{[x for x in s if x in (1, 2, 3, 7, 11, 19, 43, 67, 163)]}; sum of reciprocals 1/2+1/3+1/7+1/43+1/1806 = {Fr(1,2)+Fr(1,3)+Fr(1,7)+Fr(1,43)+Fr(1,1806)}')
for k in range(3, 9):
    chi = Fr(1, 2) + Fr(1, 3) + Fr(1, k) - 1
    out.append(f'(2,3,{k}): chi = {chi}' + (f', 2/chi = {2 / chi}' if chi > 0 else '') + (f', (2g-2)/|chi| = {2 / (-chi)}(g-1)' if chi < 0 else ''))
# tau(n) = 4: max abundancy
best = max(((Fr(sum(d for d in range(1, n + 1) if n % d == 0), n), n) for n in range(2, 3000)
            if sum(1 for d in range(1, n + 1) if n % d == 0) == 4))
out.append(f'tau(n) = 4, n < 3000: max sigma(n)/n = {best[0]} at n = {best[1]}')
out.append(f'1093 mod 3 = {1093 % 3}, mod 7 = {1093 % 7}, mod 21 = {1093 % 21}; 3511 mod 3 = {3511 % 3}, mod 7 = {3511 % 7}, mod 21 = {3511 % 21};'
           f' (p-1) mod 13: {1092 % 13}, {3510 % 13}; C(14,3) = {math.comb(14, 3)}; (13-1)! + 1 mod 169 = {(math.factorial(12) + 1) % 169}')
out.append(f'1093 mod 4 = {1093 % 4}')
print('\n'.join(out))
