# mac-mini-2026-11-06-microcosm. Gate basics (THM-4900 context; the gate itself is Session210 graded-spectrum-gates): G(T) strong, gate unique for n >= 2, H(G(T)) = (n+2) H(T) + a_1(T); G(TT_n) has 2^n + 1 paths.
import itertools, random
def descents_vector(n, arc):
    """a_k = number of permutations of range(n) with exactly k non-arcs between consecutive entries."""
    a = [0] * n
    for p in itertools.permutations(range(n)):
        k = sum(1 for i in range(n - 1) if not arc[p[i]][p[i + 1]])
        a[k] += 1
    return a
def gate(n, arc):
    m = n + 2; s, t = n, n + 1
    A = [[False] * m for _ in range(m)]
    for i in range(n):
        for j in range(n):
            A[i][j] = arc[i][j]
        A[s][i] = True; A[i][t] = True
    A[t][s] = True
    return m, A
def is_strong(m, A):
    def reach(u, adjf):
        seen = {u}; st = [u]
        while st:
            x = st.pop()
            for y in range(m):
                if adjf(x, y) and y not in seen: seen.add(y); st.append(y)
        return seen
    return len(reach(0, lambda x, y: A[x][y])) == m and len(reach(0, lambda x, y: A[y][x])) == m
def gates(m, A):
    out = []
    for s in range(m):
        ins = [x for x in range(m) if x != s and A[x][s]]
        if len(ins) == 1:
            t = ins[0]
            outs = [y for y in range(m) if y != t and A[t][y]]
            if outs == [s]: out.append((s, t))
    return out
random.seed(1)
ok = True
for n in range(1, 7):
    for trial in range(40):
        arc = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: arc[i][j] = True
                else: arc[j][i] = True
        a = descents_vector(n, arc)
        m, A = gate(n, arc)
        b = descents_vector(m, A)
        H, a1 = a[0], (a[1] if n > 1 else 0)
        ok &= b[0] == (n + 2) * H + a1
        ok &= is_strong(m, A)
        g = gates(m, A)
        ok &= (len(g) == 1) if n >= 2 else (len(g) == 3)
print('H(G(T)) = (n+2) H(T) + a_1(T), G(T) strong, gate unique for n >= 2 (3 gates when n = 1):', ok)
for n in range(1, 8):
    arc = [[i < j for j in range(n)] for i in range(n)]
    m, A = gate(n, arc)
    print('transitive TT_%d: H(G) =' % n, descents_vector(m, A)[0], ' 2^n + 1 =', 2 ** n + 1)
