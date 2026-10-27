"""Emit tournaments in gentourng's upper-triangle format (one line per tournament; char for pair (i,j), i<j, row order,
'1' iff i -> j), for feeding lemmaB_all_cycles / construct_DE beyond the orders gentourng can list.
  python3 gen_tournaments.py circ n               every circulant Z_n{S}, S ∩ -S = ∅, |S| = (n-1)/2 (2^((n-1)/2) sets)
  python3 gen_tournaments.py regular n count seed rotational circulant + 60 n^2 random 3-cycle reversals (scores kept)
  python3 gen_tournaments.py balanced n count seed  (n even) a random regular (n+1)-tournament minus vertex n, then
                                                     60 n^2 random 3-cycle reversals (scores m-1, m kept)
"""
import random, sys
def emit(A, n):
    print("".join("1" if A[i][j] else "0" for i in range(n) for j in range(i + 1, n)))
def circulant(n, S): return [[1 if i != j and (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
def shuffle3(A, n, rng):
    for _ in range(60 * n * n):
        a, b, c = rng.sample(range(n), 3)
        if A[a][b] and A[b][c] and A[c][a]:
            A[a][b] = A[b][c] = A[c][a] = 0; A[b][a] = A[c][b] = A[a][c] = 1
    return A
mode = sys.argv[1]; n = int(sys.argv[2])
if mode == "circ":
    m = (n - 1) // 2
    for mask in range(1 << m):
        S = {d if mask >> (d - 1) & 1 else n - d for d in range(1, m + 1)}
        emit(circulant(n, S), n)
elif mode == "regular":
    cnt, rng = int(sys.argv[3]), random.Random(int(sys.argv[4])); m = (n - 1) // 2
    for _ in range(cnt): emit(shuffle3(circulant(n, set(range(1, m + 1))), n, rng), n)
elif mode == "balanced":
    cnt, rng = int(sys.argv[3]), random.Random(int(sys.argv[4])); N = n + 1; m = n // 2
    for _ in range(cnt):
        A = shuffle3(circulant(N, set(range(1, m + 1))), N, rng)
        B = [row[:n] for row in A[:n]]
        assert all(sum(r) in (m - 1, m) for r in B)
        emit(shuffle3(B, n, rng), n)
