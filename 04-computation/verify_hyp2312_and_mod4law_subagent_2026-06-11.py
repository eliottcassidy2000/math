#!/usr/bin/env python3
"""
Independent verification (fresh code, no imports from repo scripts).

TASK A: n=6 exhaustive (2^15 tournaments):
  - H(T) = number of directed Hamiltonian paths (DP over subsets)
  - find max H, ALL labeled tournaments attaining it, group into iso classes
    (canonical form = min bit-encoding over all 720 relabelings)
  - per class: |Pf| = sqrt(det(S)), S = A - A^T (skew adjacency)
  Tests HYP-2312: "max-H => |Pf| = 1" (universal vs existential reading).

TASK B: MOD-4 SCORE LAW, S = A - A^T, s_i = out-degree:
  odd n:  adj(S)_ij == (-1)^(s_i+s_j)  (mod 4)   for all i,j
  even n: u = adj(S)@1,  u_i == +(-1)^(s_i) (mod 4)  -- sign convention probed
  Exhaustive n=5; random n=7,9 (odd) and n=6,8 (even), 200 samples each.
  adj computed EXACTLY via integer cofactors (Bareiss determinants).

BONUS: n=8 stochastic search for high-H tournaments; |Pf| of best found.
"""
import random
from math import isqrt
from itertools import permutations, combinations

random.seed(20260611)

# ---------- exact integer determinant (Bareiss, fraction-free) ----------
def det_bareiss(Min):
    M = [row[:] for row in Min]
    n = len(M)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def adjugate(M):
    """adj(M)_ij = (-1)^(i+j) * minor(M, row j deleted, col i deleted). Exact."""
    n = len(M)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            sub = [[M[r][c] for c in range(n) if c != i] for r in range(n) if r != j]
            adj[i][j] = (-1) ** (i + j) * det_bareiss(sub)
    return adj

# ---------- tournament helpers ----------
def pairs(n):
    return list(combinations(range(n), 2))

def bits_to_adj(bits, n, P):
    """bit b=1 means i->j for pair (i,j), else j->i. Returns 0/1 matrix A."""
    A = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(P):
        if (bits >> b) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def skew(A):
    n = len(A)
    return [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

def scores(A):
    return [sum(row) for row in A]

def ham_paths(A):
    """count directed Hamiltonian paths via DP over subsets."""
    n = len(A)
    out = [0] * n
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                out[i] |= 1 << j
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c and (mask >> v) & 1:
                nxt = out[v] & ~mask
                w = 0
                m2 = nxt
                while m2:
                    w = (m2 & -m2).bit_length() - 1
                    dp[mask | (1 << w)][w] += c
                    m2 &= m2 - 1
    return sum(dp[full])

def canon(bits, n, P):
    """canonical form: min encoding over all relabelings."""
    A = bits_to_adj(bits, n, P)
    best = None
    for perm in permutations(range(n)):
        # relabeled tournament A'[x][y] = A[perm[x]][perm[y]]
        enc = 0
        for b, (i, j) in enumerate(P):
            if A[perm[i]][perm[j]]:
                enc |= 1 << b
        if best is None or enc < best:
            best = enc
    return best

# =====================================================================
print("=" * 72)
print("TASK A: n=6 exhaustive — all H-maximizing classes and their |Pf|")
print("=" * 72)
n = 6
P = pairs(n)
m = len(P)
maxH = -1
winners = []
for bits in range(1 << m):
    A = bits_to_adj(bits, n, P)
    H = ham_paths(A)
    if H > maxH:
        maxH = H
        winners = [bits]
    elif H == maxH:
        winners.append(bits)
print(f"n=6: max H over all 2^{m} tournaments = {maxH}")
print(f"labeled tournaments attaining max H: {len(winners)}")

classes = {}
for bits in winners:
    c = canon(bits, n, P)
    classes.setdefault(c, []).append(bits)
print(f"isomorphism classes attaining H={maxH}: {len(classes)}")
for c, members in sorted(classes.items()):
    A = bits_to_adj(c, n, P)
    S = skew(A)
    d = det_bareiss(S)
    pf = isqrt(d)
    assert pf * pf == d, f"det(S)={d} not a perfect square!"
    H = ham_paths(A)
    sc = sorted(scores(A))
    print(f"  class canon={c:5d}: labeled count={len(members):4d}, H={H}, "
          f"det(S)={d}, |Pf|={pf}, scores={sc}")

# =====================================================================
print()
print("=" * 72)
print("TASK B: MOD-4 SCORE LAW  (S = A - A^T, exact integer adjugate)")
print("=" * 72)

def check_odd(A):
    """returns list of violations (i,j,adj_ij mod4, expected)"""
    S = skew(A)
    adj = adjugate(S)
    s = scores(A)
    bad = []
    for i in range(len(A)):
        for j in range(len(A)):
            expect = 1 if (s[i] + s[j]) % 2 == 0 else 3   # (-1)^(s_i+s_j) mod 4
            if adj[i][j] % 4 != expect:
                bad.append((i, j, adj[i][j] % 4, expect))
    return bad

def check_even(A):
    """returns (plus_ok, minus_ok, u mod 4 list, expected list)"""
    S = skew(A)
    adj = adjugate(S)
    s = scores(A)
    nn = len(A)
    u = [sum(adj[i][j] for j in range(nn)) for i in range(nn)]
    plus_ok = all(u[i] % 4 == (1 if s[i] % 2 == 0 else 3) for i in range(nn))
    minus_ok = all(u[i] % 4 == (3 if s[i] % 2 == 0 else 1) for i in range(nn))
    return plus_ok, minus_ok

# ---- odd n=5 exhaustive ----
n5, P5 = 5, pairs(5)
viol = 0
for bits in range(1 << len(P5)):
    if check_odd(bits_to_adj(bits, n5, P5)):
        viol += 1
print(f"odd n=5 EXHAUSTIVE (1024 tournaments): violations = {viol}")

# ---- odd n=7, 9 random ----
for nn in (7, 9):
    Pn = pairs(nn)
    viol = 0
    for _ in range(200):
        bits = random.getrandbits(len(Pn))
        if check_odd(bits_to_adj(bits, nn, Pn)):
            viol += 1
    print(f"odd n={nn} RANDOM (200 samples): violations = {viol}")

# ---- even n=6 exhaustive (sign census) ----
n6, P6 = 6, pairs(6)
cnt_plus = cnt_minus = cnt_neither = 0
for bits in range(1 << len(P6)):
    p, mn = check_even(bits_to_adj(bits, n6, P6))
    if p:
        cnt_plus += 1
    elif mn:
        cnt_minus += 1
    else:
        cnt_neither += 1
print(f"even n=6 EXHAUSTIVE (32768): +sign holds {cnt_plus}, "
      f"global -sign {cnt_minus}, neither {cnt_neither}")

# ---- even n=8 random ----
n8, P8 = 8, pairs(8)
cnt_plus = cnt_minus = cnt_neither = 0
for _ in range(200):
    bits = random.getrandbits(len(P8))
    p, mn = check_even(bits_to_adj(bits, n8, P8))
    if p:
        cnt_plus += 1
    elif mn:
        cnt_minus += 1
    else:
        cnt_neither += 1
print(f"even n=8 RANDOM (200): +sign holds {cnt_plus}, "
      f"global -sign {cnt_minus}, neither {cnt_neither}")

# =====================================================================
print()
print("=" * 72)
print("BONUS: n=8 hill-climb for high-H tournaments; |Pf| of maximizers found")
print("=" * 72)
nn, Pn = 8, pairs(8)
mm = len(Pn)
bestH, best_examples = -1, {}
for restart in range(60):
    bits = random.getrandbits(mm)
    A = bits_to_adj(bits, nn, Pn)
    H = ham_paths(A)
    improved = True
    while improved:
        improved = False
        for b in range(mm):
            nb = bits ^ (1 << b)
            H2 = ham_paths(bits_to_adj(nb, nn, Pn))
            if H2 > H:
                bits, H = nb, H2
                improved = True
    if H > bestH:
        bestH = H
        best_examples = {bits}
    elif H == bestH:
        best_examples.add(bits)
print(f"n=8 hill-climb best H found = {bestH} "
      f"(session/exhaustive claim: 661); examples: {len(best_examples)}")
pfs = set()
for bits in best_examples:
    A = bits_to_adj(bits, nn, Pn)
    S = skew(A)
    d = det_bareiss(S)
    pf = isqrt(d)
    assert pf * pf == d
    pfs.add(pf)
    print(f"  H={ham_paths(A)}, det(S)={d}, |Pf|={pf}, scores={sorted(scores(A))}")
print(f"|Pf| values among best-found n=8: {sorted(pfs)}")
