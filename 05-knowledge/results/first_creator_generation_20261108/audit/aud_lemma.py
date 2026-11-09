"""aud_lemma.py - auditor's check of PROOF.md section 3.2 (lambda >= 3 iff min semi-degree >= 3 for n <= 13) and of its
sharpness at n = 14, plus random checks of the general form lambda = delta0 whenever n <= 4*delta0 + 1 (THM-4811 (B)).
Brute-force lambda over all 2^n - 2 cuts."""
import random, itertools
def lam(n, A):
    best = 10**9; full = (1 << n) - 1
    for X in range(1, full):
        b = 0
        for x in range(n):
            if X >> x & 1:
                for y in range(n):
                    if not (X >> y & 1) and A[y][x]: b += 1
        best = min(best, b)
    return best
def mindeg(n, A):
    return min(min(sum(A[v]), n - 1 - sum(A[v])) for v in range(n))
def circ(m, S):  # circulant tournament on Z_m with connection set S
    return [[1 if (j - i) % m in S else 0 for j in range(m)] for i in range(m)]
def osum(A, B):  # ordinal sum A => B
    a, b = len(A), len(B); n = a + b
    C = [[0] * n for _ in range(n)]
    for i in range(a):
        for j in range(a): C[i][j] = A[i][j]
        for j in range(b): C[i][a + j] = 1
    for i in range(b):
        for j in range(b): C[a + i][a + j] = B[i][j]
    return C
R7 = circ(7, {1, 2, 3})
T = osum(R7, R7)
print("n=14, R7 => R7: mindeg %d, lambda %d" % (mindeg(14, T), lam(14, T)))
for r in (1, 2):
    U = [row[:] for row in T]
    for k in range(r):   # reverse arc (k) -> (7 + k): a back arc from the second block into the first
        U[k][7 + k] = 0; U[7 + k][k] = 1
    print("n=14, R7 => R7 with %d back arc(s): mindeg %d, lambda %d" % (r, mindeg(14, U), lam(14, U)))
# random checks of lambda = delta0 for n <= 4 delta0 + 1 (biased toward near-regular tournaments)
random.seed(20261108); bad = 0; tested = 0
for n in range(5, 14):
    for t in range(60 if n <= 11 else 20):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: A[i][j] = 1
                else: A[j][i] = 1
        d = mindeg(n, A)
        if n <= 4 * d + 1:
            tested += 1
            if lam(n, A) != d: bad += 1
print("random tournaments n=5..13 with n <= 4*delta0+1: tested %d, lambda != delta0: %d" % (tested, bad))
