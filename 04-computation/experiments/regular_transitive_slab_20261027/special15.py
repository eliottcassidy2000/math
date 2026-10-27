"""Rows (hex out-masks, vertex 0 = bit 0) of the special regular 15-tournaments of earlier work, for
paley_geodesic_search mode 3: R5[C3] and C3[R5] (THM-4812, the vertex-transitive non-Cayley classes) and the doubly
regular doubling D(P_7) (THM-4557, Aut = F_21)."""
def lex(outer_n, outer_S, inner_n, inner_S):
    V = [(i, a) for i in range(outer_n) for a in range(inner_n)]; idx = {x: k for k, x in enumerate(V)}
    rows = []
    for (i, a) in V:
        r = 0
        for (j, b) in V:
            if (i, a) == (j, b): continue
            if (i != j and (j - i) % outer_n in outer_S) or (i == j and (b - a) % inner_n in inner_S): r |= 1 << idx[(j, b)]
        rows.append(r)
    return rows
def doubling_P7():
    n = 7; Q = {1, 2, 4}; A = lambda i, j: (j - i) % n in Q
    # vertices: i (0..6) -> k = i; 0' -> 7; i' -> 8 + i
    N = 15; rows = [0] * N
    def arc(u, w): rows[u] |= 1 << w
    for i in range(n):
        for j in range(n):
            if i != j and A(i, j): arc(i, j); arc(8 + j, 8 + i)        # T, and T' converse: i' -> j' iff j -> i
            if A(i, j) or i == j: arc(i, 8 + j)                         # i -> j' iff i -> j or i = j
            if i != j and A(i, j): arc(8 + i, j)                        # i' -> j iff i -> j
        arc(7, i); arc(8 + i, 7)                                        # 0' -> T, T' -> 0'
    return rows
def check(rows, N):
    for u in range(N):
        assert bin(rows[u]).count('1') == (N - 1) // 2, (u, rows)
        for w in range(N):
            if u != w: assert ((rows[u] >> w) & 1) + ((rows[w] >> u) & 1) == 1
    return rows
for name, rows in (("R5[C3]", lex(5, {1, 2}, 3, {1})), ("C3[R5]", lex(3, {1}, 5, {1, 2})), ("D(P7)", doubling_P7())):
    check(rows, 15); print(name, ' '.join(f"{r:x}" for r in rows))
