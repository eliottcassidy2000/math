# (1) exact identity used for THM-4756:  sum_(X subset V) H(T[X]) H(T[V\X]) = sum_R 4^|R| 2^(n - |V(R)|),
#     R over sets of pairwise vertex-disjoint directed odd cycles (H(empty) = 1);  all classes n <= 6.
# (2) Eulerian polynomials: A_n(x) = (x-1)^(n-1) mod 3 (n odd), (x-1)^(n-2)(x+1) mod 3 (n even), n <= 40.
# mac-mini-2026-10-22-sierp.
import itertools, subprocess
from functools import lru_cache
def tourns(n):
    for code in subprocess.run(['gentourng', str(n)], capture_output=True, text=True).stdout.split():
        if code.startswith('>'): continue
        adj = [[0]*n for _ in range(n)]; p = 0
        for i in range(n):
            for j in range(i+1, n):
                if code[p] == '1': adj[i][j] = 1
                else: adj[j][i] = 1
                p += 1
        yield adj
def H_of(adj, verts):
    verts = list(verts)
    if not verts: return 1
    return sum(1 for perm in itertools.permutations(verts) if all(adj[perm[t]][perm[t+1]] for t in range(len(perm)-1)))
def odd_cycles(adj, n):
    cyc = []
    for k in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), k):
            first = sub[0]
            for perm in itertools.permutations(sub[1:]):
                seq = (first,) + perm
                if all(adj[seq[t]][seq[(t+1) % k]] for t in range(k)): cyc.append(frozenset(sub))
    return cyc
bad = 0; tot = 0
for n in range(1, 7):
    for adj in tourns(n) if n >= 2 else [[[0]]]:
        V = range(n)
        lhs = sum(H_of(adj, X) * H_of(adj, [v for v in V if v not in X]) for r in range(n + 1) for X in itertools.combinations(V, r))
        cyc = odd_cycles(adj, n)
        rhs = 0
        def rec(i, used, cnt, size):
            global rhs
            if i == len(cyc):
                rhs += 4 ** cnt * 2 ** (n - size); return
            rec(i + 1, used, cnt, size)
            if not (cyc[i] & used): rec(i + 1, used | cyc[i], cnt + 1, size + len(cyc[i]))
        rec(0, frozenset(), 0, 0)
        tot += 1
        if lhs != rhs: bad += 1
print(f"(1) exact identity sum_X H(X)H(X^c) = sum_R 4^|R| 2^(n-|V(R)|): {tot} classes n <= 6, {bad} mismatches")
# Eulerian numbers
def eulerian(n):
    A = [1]
    for m in range(2, n + 1):
        B = [0] * m
        for k in range(m):
            B[k] = (k + 1) * (A[k] if k < len(A) else 0) + (m - k) * (A[k - 1] if k - 1 >= 0 else 0)
        A = B
    return A
def polymod3_power(a, e):  # (x + a)^e mod 3 as list low->high
    r = [1]
    for _ in range(e):
        r2 = [0] * (len(r) + 1)
        for i, c in enumerate(r): r2[i] = (r2[i] + a * c) % 3; r2[i + 1] = (r2[i + 1] + c) % 3
        r = r2
    return r
ok = True
for n in range(2, 41):
    A = [c % 3 for c in eulerian(n)]
    if n % 2: T = polymod3_power(-1 % 3, n - 1)
    else:
        P = polymod3_power(2, n - 2); T = [0] * (len(P) + 1)
        for i, c in enumerate(P): T[i] = (T[i] + c) % 3; T[i + 1] = (T[i + 1] + c) % 3
    if A != T[:len(A)] or any(T[len(A):]): ok = False; print("fail at", n)
print("(2) A_n(x) = (x-1)^(n-1) [odd] / (x-1)^(n-2)(x+1) [even] mod 3 for 2 <= n <= 40:", ok)
