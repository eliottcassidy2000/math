"""C7: openai/math #197's 'Fano gadget' (complements of the 7 Fano lines, pairwise intersections 2)
is the closed in-neighbourhood design {N^-[i]} of the Paley tournament P_7 = T_3.  Check that every
doubly regular tournament of order n = 7 (mod 8) gives the same parity gadget: blocks of even size
(n+1)/2, pairwise intersections (n+1)/4 (even), each point in (n+1)/2 blocks, each pair in (n+1)/4 blocks,
and the full point set meets each block evenly but itself oddly.  Also the F_2-rank of the blocks.
Tournaments: the repo's Sierpinski skew-Hadamard tower T_k (THM-4557 doubling from T_3 = P_7) and Paley P_q.
"""
import itertools
import numpy as np

def paley(q):
    QR = {(i * i) % q for i in range(1, q)}
    return np.array([[1 if (j - i) % q in QR else 0 for j in range(q)] for i in range(q)], dtype=np.int64)

def double(A):
    """THM-4557 doubling D(T): vertices T, 0', T'; arcs: inside T as T; T' converse; 0' -> T; T' -> 0';
    i -> j' iff i -> j or i = j; i' -> j iff i -> j."""
    n = A.shape[0]
    N = 2 * n + 1
    B = np.zeros((N, N), dtype=np.int64)
    B[:n, :n] = A                       # T
    for i in range(n):
        B[n, i] = 1                     # 0' -> T
        B[n + 1 + i, n] = 1             # T' -> 0'
        for j in range(n):
            if A[j, i]:
                B[n + 1 + i, n + 1 + j] = 1   # i' -> j' iff j -> i
            if A[i, j] or i == j:
                B[i, n + 1 + j] = 1           # i -> j'
            if A[i, j]:
                B[n + 1 + i, j] = 1           # i' -> j
    return B

def is_tournament(A):
    n = A.shape[0]
    return all(A[i, i] == 0 for i in range(n)) and all(A[i, j] + A[j, i] == 1 for i in range(n) for j in range(n) if i != j)

def is_drt(A):
    n = A.shape[0]
    C = A @ A.T   # common out-neighbours
    t = (n - 3) // 4
    return all(C[i, j] == t for i in range(n) for j in range(n) if i != j)

def f2_rank(M):
    M = M.copy() % 2
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r

def gadget(A, name):
    n = A.shape[0]
    D = (A.T + np.eye(n, dtype=np.int64))  # row i = indicator of N^-[i] = {i} u {j : j -> i}
    sizes = set(D.sum(1))
    G = D @ D.T
    inter = {G[i, j] for i in range(n) for j in range(n) if i != j}
    pts = set(D.sum(0))
    pair = {int(sum(D[:, a] * D[:, b])) for a, b in itertools.combinations(range(n), 2)}
    even = all(s % 2 == 0 for s in sizes) and all(v % 2 == 0 for v in inter)
    print("  %-10s n=%3d  n mod 8=%d  tournament=%s DRT=%s | block sizes %s, pairwise intersections %s, "
          "point degrees %s, pair degrees %s, even gadget=%s, F2-rank=%d"
          % (name, n, n % 8, is_tournament(A), is_drt(A), sorted(sizes), sorted(inter), sorted(pts), sorted(pair),
             even, f2_rank(D)))
    return D

print("Fano check: P_7 closed in-neighbourhoods vs complements of the lines {i+1, i+2, i+4}")
A7 = paley(7)
D7 = A7.T + np.eye(7, dtype=np.int64)
lines = [{(i + d) % 7 for d in (1, 2, 4)} for i in range(7)]
comps = [set(range(7)) - L for L in lines]
blocks = [set(np.nonzero(D7[i])[0]) for i in range(7)]
print("  closed in-neighbourhoods == Fano line complements (as set systems):",
      sorted(map(sorted, blocks)) == sorted(map(sorted, comps)))
print("\nGadget table")
T = A7
for k in range(3, 8):
    gadget(T, "T_%d" % k)
    T = double(T)
for q in (11, 19, 23, 31, 43, 47, 71, 79):
    gadget(paley(q), "Paley P_%d" % q)
