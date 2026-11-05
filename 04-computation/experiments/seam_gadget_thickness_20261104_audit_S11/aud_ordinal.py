"""aud_ordinal.py -- audit S11 of the note's section 4 (d): the four ordinal sums of a rooted tournament (X, a) with C3,
by permutation brute force over every labelled rooted tournament X of order 1..4.
  (1) X => C3 (every X-vertex beats every C3-vertex), exit a in X        claimed: v -> 3 v
  (2) X => C3, exit c in C3                                              R = 3 * C2(X), C2(X) = sum_U H(U) H(X - U)
  (3) C3 => X, exit a in X                                               v -> (3H, 0, 12R)   (non-scalar, inside the triple)
  (4) C3 => X, exit c in C3                                              v -> (3H, H, 3H)
Also tests whether the case-(2) state is a function of v(X, a) (i.e. stays in the triple calculus)."""
import itertools
from collections import defaultdict

def tours(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def is_path(A, seq): return all(A[seq[t]][seq[t + 1]] for t in range(len(seq) - 1))

def state(A, n, a):
    H = S = R = 0
    for p in itertools.permutations(range(n)):
        if is_path(A, p): H += 1; S += p[0] == a
        for c in range(1, n + 1):
            if p[c - 1] == a and is_path(A, p[:c]) and is_path(A, p[c:]): R += 1
    return (H, S, R)

def Hcount(A, verts):
    return sum(1 for p in itertools.permutations(verts) if is_path(A, p)) if verts else 1

def C2(A, n):
    tot = 0
    for mask in range(1 << n):
        U = [v for v in range(n) if mask >> v & 1]; W = [v for v in range(n) if not mask >> v & 1]
        tot += Hcount(A, U) * Hcount(A, W)
    return tot

def ordinal(A, n, x_first):
    """X on 0..n-1, C3 on n..n+2 (n -> n+1 -> n+2 -> n); x_first: every X vertex beats every C3 vertex"""
    N = n + 3; B = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n): B[i][j] = A[i][j]
    B[n][n + 1] = B[n + 1][n + 2] = B[n + 2][n] = 1
    for i in range(n):
        for c in range(n, N):
            if x_first: B[i][c] = 1
            else: B[c][i] = 1
    return B, N

bad = defaultdict(int); tested = 0; case2 = defaultdict(set)
for n in range(1, 5):
    for A in tours(n):
        c2 = C2(A, n)
        for a in range(n):
            H, S, R = v = state(A, n, a); tested += 1
            B1, N = ordinal(A, n, True); B3, _ = ordinal(A, n, False)
            if state(B1, N, a) != (3 * H, 3 * S, 3 * R): bad[1] += 1
            w2 = state(B1, N, n)
            if w2 != (3 * H, 0, 3 * c2): bad[2] += 1
            case2[v].add(w2)
            if state(B3, N, a) != (3 * H, 0, 12 * R): bad[3] += 1
            if state(B3, N, n) != (3 * H, H, 3 * H): bad[4] += 1
print(f"rooted tournaments X of order 1-4 tested: {tested}; formula failures per case: {dict(bad) or 'none'}")
amb = {v: sorted(s) for v, s in case2.items() if len(s) > 1}
print(f"case (2): input states v(X,a) whose output state is NOT determined by v(X,a): {len(amb)}; examples: {list(amb.items())[:3]}")
