#!/usr/bin/env python3
"""Independent brute verification of a two-level .sets witness for
R(B_{n-1}, B_n) >= 4n-1, plus Epoch adjacency-string emission."""
import numpy as np
import sys

def load(fn):
    with open(fn) as f:
        head = f.readline().split()
        n, m = int(head[1]), int(head[3])
        U0 = [int(c) for c in f.readline().split()[1]]
        U1 = [int(c) for c in f.readline().split()[1]]
        D = [int(c) for c in f.readline().split()[1]]
    return n, m, np.array(U0), np.array(U1), np.array(D)

def build(n, m, U0, U1, D):
    N = 4 * n - 2
    A = np.zeros((N, N), dtype=np.int8)
    for i in range(m):
        for j in range(m):
            dd = (j - i) % m
            if i != j:
                if U0[dd]: A[i, j] = 1
                if U1[dd]: A[m + i, m + j] = 1
            if D[dd]:
                A[i, m + j] = 1
                A[m + j, i] = 1
    assert (A == A.T).all() and (np.diag(A) == 0).all()
    return A

def verify(A, n):
    N = A.shape[0]
    M = A.astype(np.int64)
    C = M @ M
    comp = 1 - M
    np.fill_diagonal(comp, 0)
    Cc = comp @ comp
    won = -1; woff = -1; bad = 0
    for i in range(N):
        for j in range(i + 1, N):
            if M[i, j]:
                won = max(won, C[i, j])
                if C[i, j] > n - 2: bad += 1
            else:
                woff = max(woff, Cc[i, j])
                if Cc[i, j] > n - 1: bad += 1
    return bad == 0, won, woff, bad

if __name__ == "__main__":
    fn = sys.argv[1]
    n, m, U0, U1, D = load(fn)
    A = build(n, m, U0, U1, D)
    ok, won, woff, bad = verify(A, n)
    N = A.shape[0]
    print(f"{fn}: n={n} N={N} |U0|={U0.sum()} |U1|={U1.sum()} |D|={D.sum()}")
    print(f"VERIFY: ok={ok} worst_edge_codeg={won} (need <= {n-2}) "
          f"worst_nonedge_cononbrs={woff} (need <= {n-1}) violations={bad}")
    if ok:
        out = fn.replace(".sets", "_adjstring.txt")
        with open(out, "w") as f:
            f.write("".join(str(A[i, j]) for j in range(N) for i in range(N)) + "\n")
        print(f"witness OK — Epoch adjacency string written to {out}")
