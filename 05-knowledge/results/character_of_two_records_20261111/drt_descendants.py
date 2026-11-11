#!/usr/bin/env python3
"""Doubly regular tournaments obtained from one skew Hadamard matrix: for each index c of the Cayley-Dickson sign
matrix G of order 2^k (G + G^T = 2I after the sign normalisation of the diagonal), normalise row and column c and take
the tournament on the other 2^k - 1 indices.  Prints, for each descendant: regular / doubly regular check, H,
H mod 4, the numbers of 3-, 5-cycles' parity data, whether it is self-converse, |Aut|, and an isomorphism label.
Usage: drt_descendants.py k   (k = 3 or 4)"""
import sys
from cayley_dickson import gamma_table
def ham(out, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            m = out[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1; dp[S | (1 << w)][w] += c; m &= m - 1
    return sum(dp[full])
def count_maps(A, B, n, anti=False, limit=None):
    """number of isomorphisms A -> B (anti=False) or anti-isomorphisms (a->b iff f(b)->f(a))."""
    arcA = [[(A[a] >> b) & 1 for b in range(n)] for a in range(n)]
    arcB = [[(B[a] >> b) & 1 for b in range(n)] for a in range(n)]
    perm = [-1] * n; used = [False] * n; cnt = 0
    def rec(i):
        nonlocal cnt
        if i == n: cnt += 1; return limit is not None and cnt >= limit
        for c in range(n):
            if used[c]: continue
            ok = True
            for j in range(i):
                if anti:
                    if arcA[i][j] != arcB[perm[j]][c]: ok = False; break
                else:
                    if arcA[i][j] != arcB[c][perm[j]]: ok = False; break
            if not ok: continue
            perm[i] = c; used[c] = True
            if rec(i + 1): return True
            perm[i] = -1; used[c] = False
        return False
    rec(0); return cnt
if __name__ == "__main__":
    k = int(sys.argv[1]); N = 1 << k; g = gamma_table(k)
    # skew Hadamard form: M = (D G)^T with D = diag(1, -1, ..., -1).  Then M + M^T = 2I, M is Hadamard, and for
    # a != b both non-zero M[a][b] = -G[b][a] = G[a][b], so the descendant at c = 0 is CD_k itself.
    M = [[(g[b][a] if b == 0 else -g[b][a]) for b in range(N)] for a in range(N)]
    assert all(M[a][a] == 1 for a in range(N))
    assert all(M[a][b] == g[a][b] for a in range(1, N) for b in range(1, N) if a != b)
    assert all(M[a][b] == -M[b][a] for a in range(N) for b in range(N) if a != b)
    assert all(sum(M[a][c] * M[b][c] for c in range(N)) == 0 for a in range(N) for b in range(a + 1, N))
    tours = []
    for c in range(N):
        d = [M[c][i] if i != c else 1 for i in range(N)]
        Mp = [[d[a] * M[a][b] * d[b] for b in range(N)] for a in range(N)]
        assert all(Mp[c][i] == 1 for i in range(N) if i != c) and all(Mp[i][c] == -1 for i in range(N) if i != c)
        idx = [i for i in range(N) if i != c]; n = N - 1
        out = [0] * n
        for x, a in enumerate(idx):
            for y, b in enumerate(idx):
                if a != b and Mp[a][b] == 1: out[x] |= 1 << y
        scores = set(bin(o).count('1') for o in out)
        common = set(bin(out[a] & out[b]).count('1') for a in range(n) for b in range(a + 1, n))
        tours.append((c, out, scores, common))
    reps = []
    for c, out, scores, common in tours:
        n = N - 1; H = ham(out, n)
        aut = count_maps(out, out, n); anti = count_maps(out, out, n, anti=True, limit=1)
        label = None
        for r, (lab, o2) in enumerate(reps):
            if count_maps(out, o2, n, limit=1): label = lab; break
        if label is None:
            label = chr(ord('A') + len(reps)); reps.append((label, out))
        print("index %2d: scores %s, common out-neighbours %s, H = %d, H mod 4 = %d, |Aut| = %d, self-converse: %s, class %s" % (
            c, sorted(scores), sorted(common), H, H % 4, aut, "yes" if anti else "no", label)); sys.stdout.flush()
    # converse relations between classes
    for i, (la, oa) in enumerate(reps):
        for lb, ob in reps[i:]:
            if count_maps(oa, ob, N - 1, anti=True, limit=1): print("classes %s and %s are converse to each other" % (la, lb))
