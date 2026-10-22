# Doubly regular tournaments (S S^T = nI - J) from skew-Hadamard matrices by doubling:
# A skew-Hadamard of order m (A + A^T = 2I, A A^T = mI)  =>  K = [[A, A], [-A^T, A^T]] skew-Hadamard of order 2m.
# Normalising K at row r (conjugate by the diagonal sign matrix making row r of K - I equal to +1 off the diagonal)
# and deleting r gives a doubly regular tournament of order 2m - 1.  The loop expansion (THM-4757) predicts that
# all DRTs of one order agree in H through three loops, i.e. to relative O(n^-3).  Usage: drt_test.py k (order 2^k).
import sys, itertools, math
import numpy as np

def skew_hadamard(k):
    A = np.array([[1, 1, 1, 1], [-1, 1, 1, -1], [-1, -1, 1, 1], [-1, 1, -1, 1]])   # order 4: A + A^T = 2I
    assert (A + A.T == 2 * np.eye(4)).all() and (A @ A.T == 4 * np.eye(4)).all()
    while A.shape[0] < 2 ** k:
        A = np.block([[A, A], [-A.T, A.T]])
    m = A.shape[0]
    assert (A + A.T == 2 * np.eye(m)).all() and (A @ A.T == m * np.eye(m)).all()
    return A

def drt_from(A, r):
    m = A.shape[0]; S = A - np.eye(m, dtype=int)
    d = np.array([S[r, j] if j != r else 1 for j in range(m)])
    S2 = (d[:, None] * S) * d[None, :]          # row r of S2 is +1 off diagonal
    keep = [j for j in range(m) if j != r]
    T = S2[np.ix_(keep, keep)]
    n = m - 1
    assert (T @ T.T == n * np.eye(n) - np.ones((n, n))).all()
    return T

def held_karp(S):
    n = S.shape[0]
    out = [sum(1 << w for w in range(n) if S[v, w] == 1) for v in range(n)]
    dp = np.zeros((1 << n, n), dtype=object)
    for v in range(n): dp[1 << v, v] = 1
    for msk in range(1, 1 << n):
        row = dp[msk]
        for v in range(n):
            c = row[v]
            if not c: continue
            nxt = out[v] & ~msk
            while nxt:
                w = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                dp[msk | (1 << w), w] += c
    return int(sum(dp[(1 << n) - 1]))

if __name__ == '__main__':
    k = int(sys.argv[1]); A = skew_hadamard(k); m = A.shape[0]; n = m - 1
    seen = {}
    for r in range(m):
        T = drt_from(A, r)
        # cheap invariant: sorted list over vertices of the number of 3-cycles in the out-neighbourhood is constant for DRT;
        # use the multiset of 4-cycle-ish counts of S^3 entries and H itself
        key = tuple(sorted(np.round(np.abs(np.linalg.eigvals((T == 1).astype(float)) ), 6)))
        seen.setdefault(r, T)
    print(f"skew-Hadamard order {m}; DRTs of order {n} from the {m} row normalisations")
    if n <= 15:
        Hs = {}
        for r in range(m):
            H = held_karp(seen[r]); Hs[r] = H
            print(f"  row {r:2d}: H = {H}  kappa = {H * 2 ** (n - 1) / math.factorial(n):.8f}", flush=True)
        vals = sorted(set(Hs.values()))
        print(f"distinct H values: {vals}")
        if len(vals) > 1: print(f"relative spread (max-min)/mean = {(vals[-1] - vals[0]) / np.mean(vals):.3e}; n^-3 = {n ** -3:.3e}")
    np.save(f'drt_{n}.npy', np.array([seen[r] for r in range(m)]))
