"""audit_n6_strong_lmin.py -- list the strong L-local minima of order 6 and the Hamiltonian paths from the min-score
vertices of the named class 111011111111101 (brute force; reads gentourng output on stdin).
Run: gentourng -q 6 | python3 -I audit_n6_strong_lmin.py"""
import sys
from itertools import permutations
n = 6
def decode(s):
    A = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A
def hps(A):
    res = []
    def rec(p, used):
        if len(p) == n: res.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]: p.append(v); rec(p, used | 1 << v); p.pop()
    for s0 in range(n): rec([s0], 1 << s0)
    return res
def iso(A, B):
    return any(all(A[i][j] == B[p[i]][p[j]] for i in range(n) for j in range(n) if i != j) for p in permutations(range(n)))
def naut(A):
    return sum(1 for p in permutations(range(n)) if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j))
def strong(A):
    for s0 in range(n):
        seen = {s0}; st = [s0]
        while st:
            u = st.pop()
            for v in range(n):
                if A[u][v] and v not in seen: seen.add(v); st.append(v)
        if len(seen) < n: return False
    return True
L2 = lambda t: sum(abs(x - i) for i, x in enumerate(sorted(t)))
for line in sys.stdin:
    s = line.strip()
    if not s or s[0] not in '01': continue
    A = decode(s); sc = [sum(r) for r in A]
    if not strong(A) or L2(sc) == 0: continue
    P = hps(A)
    desc = any(L2([n-1-sc[v] + (v == p[0]) - (v == p[-1]) for v in range(n)]) == L2(sc) - 2 for p in P)
    if not desc:
        conv = [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]
        print(f"strong L-local minimum {s}: scores {sorted(sc)} L={L2(sc)//2} |Aut|={naut(A)} self-converse={iso(A, conv)}")
    if s == '111011111111101':
        print('   named class: Hamiltonian paths starting at 3:', [p for p in P if p[0] == 3], ' starting at 5:', [p for p in P if p[0] == 5])
