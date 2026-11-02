#!/usr/bin/env python3
"""paley_circulant_counts.py (klein-2026-11-02) -- pure-Python Hamiltonian path counts (subset DP), standard library:
(1) Paley tournaments P_q, q = 3, 7, 11: independent recount of scripts/hampath.c, and the quotient H / (q(q-1)/2);
(2) the rotational tournaments R_q = Circ(Z_q, {1..(q-1)/2}) for q = 3..13 and ALL circulant tournaments of orders 7, 9,
    11: H is a multiple of q (the rotation group acts freely on Hamiltonian paths) and the quotient is odd (Redei);
(3) the witnesses: the two order-6 tournaments with H = 17 and the first order-8 one with H = 257 (from the census),
    and the order-11 and order-12 hill-climb witnesses with H = 65537, recounted here with their score sequences."""
import itertools
def hcount(adj):
    n = len(adj); full = (1 << n) - 1
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            av = adj[v]
            for w in range(n):
                if av[w] and not (S >> w & 1): f[S | 1 << w][w] += c
    return sum(f[full])
def circ(q, S): return [[1 if ((j - i) % q) in S else 0 for j in range(q)] for i in range(q)]
def paley(q): return circ(q, {x * x % q for x in range(1, q)})
def decode(line, n):
    adj = [[0] * n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[idx] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj
def strong(adj):
    n = len(adj)
    def reach(A):
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in range(n):
                if A[u][v] and v not in seen: seen.add(v); st.append(v)
        return len(seen) == n
    return reach(adj) and reach([[adj[j][i] for j in range(n)] for i in range(n)])
print('(1) Paley tournaments')
for q in (3, 7, 11):
    h = hcount(paley(q)); aut = q * (q - 1) // 2
    print('    P_%d: H = %d = %d * %d (|Aut| = q(q-1)/2 = %d, quotient odd %s)' % (q, h, aut, h // aut, aut, h // aut % 2 == 1))
print('(2) circulant tournaments')
for q in (3, 5, 7, 9, 11, 13):
    h = hcount(circ(q, set(range(1, (q - 1) // 2 + 1))))
    print('    R_%d: H = %d = %d * %d, quotient odd %s' % (q, h, q, h // q, (h // q) % 2 == 1))
for q in (7, 9, 11):
    half = list(range(1, (q - 1) // 2 + 1)); vals = {}
    for signs in itertools.product([1, -1], repeat=len(half)):
        S = frozenset((s * a) % q for s, a in zip(signs, half))
        if 1 in S: vals.setdefault(hcount(circ(q, S)), []).append(tuple(sorted(S)))
    print('    all circulants of order %d (connection sets with 1): H values %s; all multiples of %d: %s; quotients %s; maximiser(s) %s' % (
        q, sorted(vals), q, all(h % q == 0 for h in vals), [h // q for h in sorted(vals)], vals[max(vals)]))
print('(3) witnesses')
for n, line, want in ((6, '111011111111110', 17), (6, '101011111111101', 17), (8, '1111101111001111111111110110', 257),
                      (11, '0010010010100011000001001100011011000110101001101110111', 65537),
                      (11, '1011101001100010111010111010100010001100101010100011000', 65537),
                      (11, '1000001010011010110001010001010000111000011110110011011', 65537),
                      (12, '000010001001001011101100010011010001011101111101110111110111100110', 65537)):
    adj = decode(line, n); h = hcount(adj)
    print('    order %2d %s: H = %d (expected %d) %s; scores %s; strong %s' % (n, line, h, want, 'OK' if h == want else 'MISMATCH', sorted(sum(r) for r in adj), strong(adj)))
