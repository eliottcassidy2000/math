#!/usr/bin/env python3
"""halved.py (Q-12): rooted pairs (h, B) of all tournaments of orders <= 8; which strong D of orders <= 7 have a halved
partner (h(T), B_t(T)) = (h(D), P_D(2)/2)?  Also lists the partners of the (15, 216) source."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def Bvec(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T))
    return H[N-1], [sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) for t in range(n)]
def classes1(n): return [[[0]]] if n == 1 else ([[[0, 1], [0, 0]]] if n == 2 else classes(n))
pairs = {}
for n in range(1, 9):
    for T in classes1(n):
        h, B = Bvec(T)
        for b in B:
            if (h, b) not in pairs: pairs[(h, b)] = (n, T)
hp = [tuple(map(int, l.split())) for l in open('hp8.txt')]
strongD = [(n, h, P) for n, s, h, P in hp if s and n <= 7]
have = [(n, h, P) for n, h, P in strongD if (h, P // 2) in pairs]
print(f"rooted pairs (h, B) of orders <= 8: {len(pairs)}; strong D of orders 1-7: {len(strongD)}; with a halved partner of order <= 8: {len(have)}")
print("C_3 (3, 12) partner (3, 6)?", (3, 6) in pairs, "; strong4 (5, 32) partner (5, 16)?", (5, 16) in pairs, "; (15, 216) partner (15, 108)?", (15, 108) in pairs, "-> order", pairs.get((15, 108), (None,))[0])
by_order = {}
for n, h, P in strongD: by_order.setdefault(n, [0, 0]); by_order[n][0] += 1; by_order[n][1] += ((h, P // 2) in pairs)
print("by order of D (total, with partner):", by_order)
print("partners found (n_D, h, P, order of partner):", sorted((n, h, P, pairs[(h, P // 2)][0]) for n, h, P in have)[:30])
