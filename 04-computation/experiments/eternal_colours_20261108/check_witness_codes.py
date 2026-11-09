# independent check of audit K1's order-11 witnesses (census convention: bit b <-> b-th pair i<j lexicographic, set iff i->j)
import sys
def decode(code, n):
    A = [[0]*n for _ in range(n)]; b = 0
    for i in range(n):
        for j in range(i+1, n):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A
def H(A):
    n = len(A); f = [[0]*n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            x = f[S][v]
            if not x: continue
            for u in range(n):
                if not (S >> u & 1) and A[v][u]: f[S | 1 << u][u] += x
    return sum(f[(1 << n) - 1])
def strong(A):
    n = len(A)
    def reach(adj):
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for u in range(n):
                if adj(v, u) and u not in seen: seen.add(u); st.append(u)
        return len(seen) == n
    return reach(lambda v, u: A[v][u]) and reach(lambda v, u: A[u][v])
def anti_iso(A):
    n = len(A); sc = [sum(r) for r in A]; phi = [-1]*n; used = [False]*n
    def bt(k):
        if k == n: return True
        for w in range(n):
            if used[w] or sc[w] != n - 1 - sc[k]: continue
            ok = all((A[k][j] == A[phi[j]][w]) and (A[j][k] == A[w][phi[j]]) for j in range(k))
            if ok:
                phi[k] = w; used[w] = True
                if bt(k + 1): return True
                used[w] = False
        phi[k] = -1; return False
    return list(phi) if bt(0) else None
if __name__ == "__main__":
  for code, exp in [(26516080775662589, 1551), (26520624255470589, 1627)]:
    A = decode(code, 11); p = anti_iso(A)
    print(code, "H =", H(A), "(expected", exp, ") strong:", strong(A), "anti-isomorphism:", p)
