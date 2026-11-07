#!/usr/bin/env python3
# audit J1 (HYP-9500): switching construction (Paley and non-Paley doubly regular cores), and checks (6 <= n <= 96) of the construction that proves (a) at every even n >= 6: substitution R[T2] (R regular odd order >= 3, T2 equicentric) plus skew-Walsh towers at n = 2^e (proof in the J1 report).
import numpy as np

def strong(A):
    n = len(A); R = (A.astype(bool) | np.eye(n, dtype=bool))
    for k in range(n): R = R | (R[:, [k]] & R[[k], :])
    return bool(R.all())
def equic(A):
    M = A - A.T; y = M @ (M @ np.ones(len(A), dtype=np.int64)); return bool((y == y[0]).all()), int(-y[0])
def skewhad(A):
    n = len(A); H = np.eye(n, dtype=np.int64) + A - A.T; return bool((H.T @ H == n * np.eye(n, dtype=np.int64)).all())
def from_M(M): return (M == 1).astype(np.int64)
def circulant_regular(m):                      # i -> i+1, ..., i+(m-1)/2  (m odd)
    return np.array([[1 if 1 <= (j - i) % m <= (m - 1) // 2 else 0 for j in range(m)] for i in range(m)], dtype=np.int64)
def paley(q):
    QR = {(i * i) % q for i in range(1, q)}
    return np.array([[1 if (j - i) % q in QR else 0 for j in range(q)] for i in range(q)], dtype=np.int64)
def ordered_sum(A1, A2):                       # A1 => A2
    a, b = len(A1), len(A2); A = np.zeros((a + b, a + b), dtype=np.int64)
    A[:a, :a] = A1; A[a:, a:] = A2; A[:a, a:] = 1; return A
def lex(A1, A2):                               # A1[A2]: replace each vertex of A1 by a copy of A2
    n1, n2 = len(A1), len(A2); return np.kron(A1, np.ones((n2, n2), dtype=np.int64)) + np.kron(np.eye(n1, dtype=np.int64), A2)
def skew_double(A):                            # THM-447: M' = [[M, M+I], [M-I, -M]]
    n = len(A); M = A - A.T; I = np.eye(n, dtype=np.int64)
    return from_M(np.block([[M, M + I], [M - I, -M]]))
def doubly_regular_core(A):                    # Reid-Brown: switch a skew-Hadamard tournament so vertex 0 is a source; the rest is doubly regular
    n = len(A); M = A - A.T; d = np.where(M[0] == -1, -1, 1); d[0] = 1; D = np.diag(d)
    M2 = D @ M @ D; assert (M2[0, 1:] == 1).all(); return from_M(M2)[1:, 1:]
def switch(A, w):
    n = len(A); M = A - A.T; D = np.eye(n, dtype=np.int64); D[w, w] = -1; return from_M(D @ M @ D)
def bordered_switched(R, w=0):                 # HYP-9500: T = R + source s, then switch at w in R
    m = len(R); T = np.zeros((m + 1, m + 1), dtype=np.int64); T[:m, :m] = R; T[m, :m] = 1
    return T, switch(T, w)
def is_doubly_regular(R):
    m = len(R); AA = R @ R.T; off = AA[~np.eye(m, dtype=bool)]
    return bool((R.sum(1) == (m - 1) // 2).all() and (off == off[0]).all()), int(off[0])

print("1. Switching construction (HYP-9500 proof), every switch vertex w checked:")
cores = [(f"Paley {q}", paley(q)) for q in (7, 11, 19, 23, 31, 43, 47)]
H16 = skew_double(skew_double(skew_double(skew_double(np.zeros((1, 1), dtype=np.int64)))))
cores.append(("core of skew-Walsh 16 (non-Paley order 15)", doubly_regular_core(H16)))
for name, R in cores:
    dr, lam_c = is_doubly_regular(R); m = len(R)
    T, _ = bordered_switched(R, 0)
    res = [bordered_switched(R, w)[1] for w in range(m)]
    print(f"   {name}: doubly regular {dr} (common out-degree {lam_c} = (m-3)/4: {lam_c == (m - 3) // 4}); R+source skew-Hadamard {skewhad(T)}, strong {strong(T)};"
          f" switched at w: all skew-Hadamard {all(skewhad(A) for A in res)}, all strong {all(strong(A) for A in res)}, all equicentric {all(equic(A)[0] for A in res)}, lambda {sorted({equic(A)[1] for A in res})}")

print("2. Lexicographic product with a regular factor: R1[T2] is equicentric with lambda(T2), and strong (R1 strong, |R1| >= 3):")
T4 = ordered_sum(np.zeros((1, 1), dtype=np.int64), circulant_regular(3))          # R_1 => R_3, equicentric lambda 3, non-strong
for n1, T2 in ((3, T4), (5, T4), (3, ordered_sum(circulant_regular(3), circulant_regular(5))), (3, from_M(np.array([[0, 1], [-1, 0]])))):
    T = lex(circulant_regular(n1), T2); e, lam = equic(T)
    print(f"   C_{n1}[T2], |T2| = {len(T2)}, lambda(T2) = {equic(T2)[1]}: n = {len(T)}, equicentric {e}, lambda {lam}, strong {strong(T)}")

print("3. A strong equicentric tournament at EVERY even n, 6 <= n <= 96 (construction below):")
def construct(n):
    e = (n & -n).bit_length() - 1; odd = n >> e                      # n = 2^e * odd
    if odd >= 3:
        if e == 1: T2 = from_M(np.array([[0, 1], [-1, 0]]))         # single arc, lambda 1
        else: T2 = ordered_sum(np.zeros((1, 1), dtype=np.int64), circulant_regular(2 ** e - 1))   # R_1 => R_(2^e - 1)
        return lex(circulant_regular(odd), T2), f"C_{odd}[R_1=>R_{2**e - 1}]"
    H = np.zeros((1, 1), dtype=np.int64)
    for _ in range(e): H = skew_double(H)                             # skew-Walsh tower, order 2^e (THM-447(2))
    R = doubly_regular_core(H)
    return bordered_switched(R, 0)[1], f"switch(core(D^{e}(K1)) + source)"
fails = []
for n in range(6, 97, 2):
    T, how = construct(n); e, lam = equic(T); s = strong(T)
    if not (e and s and len(T) == n): fails.append(n)
    if n <= 20 or n in (32, 64, 96): print(f"   n={n}: {how}: equicentric {e} (lambda {lam}), strong {s}")
print(f"   failures in 6..96: {fails}")

print("4. HYP-9500's claimed strong n = 8, lambda = 3 class with x = (-3,-1,-1,-1,1,1,1,3):")
import subprocess
out = subprocess.run(["/opt/homebrew/bin/gentourng", "-q", "8"], capture_output=True).stdout
a = np.frombuffer(out, dtype=np.uint8).reshape(-1, 29)[:, :28] - 48; iu = np.triu_indices(8, 1)
for row in a:
    A = np.zeros((8, 8), dtype=np.int64); A[iu] = row; A[(iu[1], iu[0])] = 1 - row
    e, lam = equic(A)
    if e and strong(A): print(f"   strong equicentric: lambda {lam}, sorted x {sorted((A - A.T).sum(1).tolist())}, skew-Hadamard {skewhad(A)}")

print("5. THM-4911 (3c) construction as written: regular S, regular S-bar, s_a -> t_b iff (b - a) mod k in {0..(k-1)/2}:")
for k in (3, 5, 7, 9, 11, 13):
    R = circulant_regular(k); n = 2 * k; A = np.zeros((n, n), dtype=np.int64)
    A[:k, :k] = R; A[k:, k:] = R
    for a_ in range(k):
        for b_ in range(k):
            if (b_ - a_) % k <= (k - 1) // 2: A[a_, k + b_] = 1
            else: A[k + b_, a_] = 1
    e, lam = equic(A); sc = sorted(set(A.sum(1).tolist()))
    print(f"   k={k} (n={n}): scores {sc}, equicentric {e}, lambda {lam}, strong {strong(A)}")
