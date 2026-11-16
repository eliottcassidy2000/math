#!/usr/bin/env python3
"""closure_check.py: (P-11) the near-join family is closed under left ordinal multiplication by counts:
H(R -> (T -> M_k -> z)) = h(R) H(T -> M_k -> z) for R, T of orders <= 4, every root, k = 1, 2;  (P-12) no strong tournament of
orders 4-8 has two digit-0 vertices;  (P-14) the single square-sum edge among the non-supersingular primes below 100."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
import math
def classes1(n): return [[[0]]] if n == 1 else classes(n)
def ordinal(A, B):
    a, b = len(A), len(B); n = a + b; U = [[0]*n for _ in range(n)]
    for i in range(a):
        for j in range(a): U[i][j] = A[i][j]
        for j in range(b): U[i][a+j] = 1
    for i in range(b):
        for j in range(b): U[a+i][a+j] = B[i][j]
    return U
def nearjoin(T, t0, k):
    M = [[1 if i < j else 0 for j in range(k)] for i in range(k)]; base = ordinal(T, M); n = len(base); U = [r[:] + [1] for r in base] + [[0]*(n+1)]
    U[n][t0] = 1; U[t0][n] = 0; return U
def H(U): return ham_all(U)[1][-1]
bad = tot = 0
for nr in range(1, 5):
    for R in classes1(nr):
        hR = H(R)
        for nt in range(1, 5):
            for T in classes1(nt):
                for t0 in range(nt):
                    for k in (1, 2):
                        NJ = nearjoin(T, t0, k); tot += 1
                        if H(ordinal(R, NJ)) != hR * H(NJ): bad += 1
print(f"(P-11) closure H(R -> near-join) = h(R) H(near-join): {tot} cases, {bad} failures")
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
multi = 0; tot = 0
for n in range(4, 9):
    for S in classes(n):
        if not strong(S): continue
        tot += 1; N = 1 << n; dp, Hm = ham_all(S); dpr, _ = ham_all(rev(S)); h = Hm[N-1]
        d0 = sum(1 for t in range(n) if sum(dpr[m][t] * Hm[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) == h)
        if d0 >= 2: multi += 1
print(f"(P-12) strong tournaments of orders 4-8: {tot}; with two or more digit-0 vertices: {multi}")
non = [37, 43, 53, 61, 67, 73, 79, 83, 89, 97]
print("(P-14) square-sum edges among the non-supersingular primes below 100:", [(p, q) for p in non for q in non if p < q and math.isqrt(p + q)**2 == p + q])
