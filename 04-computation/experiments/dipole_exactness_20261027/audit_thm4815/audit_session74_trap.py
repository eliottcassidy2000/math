"""audit_session74_trap.py -- the Session74 'L=1 trap' (supplement section 4, rows as out-neighbour bitmasks) checked
against THM-4815 B5/D: strong?, cut vertices, endpoint relation, L-descent.  Run: python3 -I audit_session74_trap.py"""
rows = [4, 2037, 512, 1015, 709, 725, 517, 1093, 757, 129, 893]
n = len(rows)
A = [[(rows[u] >> v) & 1 for v in range(n)] for u in range(n)]
assert all(A[u][v] + A[v][u] == 1 for u in range(n) for v in range(n) if u != v) and all(A[u][u] == 0 for u in range(n))
s = [sum(r) for r in A]
def reach(vs, v):
    seen = {v}; st = [v]
    while st:
        u = st.pop()
        for w in vs:
            if A[u][w] and w not in seen: seen.add(w); st.append(w)
    return seen
def strong(vs):
    vs = set(vs)
    if len(vs) <= 1: return True
    if len(vs) == 2: return False
    return all(reach(vs, v) == vs for v in vs)
full = (1 << n) - 1
E = [[0] * n for _ in range(1 << n)]
for v in range(n): E[1 << v][v] = 1 << v
for S in range(1, 1 << n):
    for v in range(n):
        st = E[S][v]
        if not st: continue
        for w in range(n):
            if not (S >> w) & 1 and A[v][w]: E[S | 1 << w][w] |= st
Pi = {(a, b) for b in range(n) for a in range(n) if (E[full][b] >> a) & 1}
def L2(t): t = sorted(t); return sum(abs(t[i] - i) for i in range(n))
desc = [(a, b) for (a, b) in Pi if L2([n - 1 - s[v] + (v == a) - (v == b) for v in range(n)]) == L2(s) - 2]
print('scores', s, 'sorted', sorted(s), 'L =', L2(s) // 2, 'strong', strong(range(n)))
print('cut vertices', [v for v in range(n) if not strong([w for w in range(n) if w != v])])
print('Hamiltonian-path endpoint pairs', len(Pi), '; L-descending pairs in Pi:', desc)
