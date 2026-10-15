# THM-4655 1(f): in a strong transitive-plus-one tournament (m >= 3), P(u_i before z) lies in [1/3, 2/3] for every
# 2 <= i <= m-1.  Proof = THM-4655 1(c) for i = 2, a crude geometric bound for m >= 8 and 3 <= i <= (m+1)/2, reversal
# symmetry, and this exact check for m <= 14.  The script
#   (a) cross-checks the class decomposition against explicit DFS enumeration of Hamiltonian paths (m <= 9);
#   (b) checks the statement exactly for every strong word with m <= 14;
#   (c) checks the crude-bound premises X <= W_M/2 and Y <= W_M/2 for 3 <= i <= (m+1)/2, m = 8..14 (where the proof uses them),
#       with X = K_(>=i) + W_(alpha=i) + W_(beta>i), Y = K_(<i) + W_(alpha<i) + W_(beta=i), W_M = W_(beta<i<alpha).
from fractions import Fraction as Fr
from itertools import product
from tpo_cycle_arcs import classes, probs, bal
def tournament(w):
    m = len(w); n = m + 1; z = 0                       # vertex 0 = z, vertex i = u_i
    A = [[False] * n for _ in range(n)]
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1): A[i][j] = True
        if w[i - 1] == 1: A[0][i] = True
        else: A[i][0] = True
    return A
def brute(w):
    A = tournament(w); n = len(A); cnt = {}; H = 0
    def dfs(p, used):
        nonlocal H
        if len(p) == n:
            H += 1; pos = {v: k for k, v in enumerate(p)}
            for i in range(1, n): cnt[i] = cnt.get(i, 0) + (pos[i] < pos[0])
            return
        for u in range(n):
            if not used >> u & 1 and A[p[-1]][u]: dfs(p + [u], used | 1 << u)
    for v in range(n): dfs([v], 1 << v)
    return H, {i: Fr(cnt[i], H) for i in cnt}
mism = 0
for m in range(2, 10):
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        H, pA, _ = probs(w); Hb, pb = brute(w)
        if H != Hb or any(pA[i] != pb[i] for i in range(1, m + 1)): mism += 1
print("(a) class decomposition vs DFS enumeration, m <= 9: mismatches", mism)
for m in range(3, 15):
    bad = 0
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        H, pA, _ = probs(w)
        bad += sum(1 for i in range(2, m) if not bal(pA[i]))
    print(f"(b) m={m}: (word, i) with P(u_i before z) outside [1/3,2/3], 2 <= i <= m-1: {bad}")
prem_bad = 0; checked = 0
for m in range(8, 15):
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        for i in range(3, (m + 1) // 2 + 1):
            K = sorted(set([0, m] + [k for k in range(1, m) if w[k - 1] == 0 and w[k] == 1]))
            X = sum(1 for k in K if k >= i); Y = sum(1 for k in K if k < i); WM = 0
            for b in range(1, m + 1):
                if w[b - 1] != 1: continue
                for a in range(b + 1, m + 1):
                    if w[a - 1] != 0: continue
                    wt = 2 ** (a - b - 1)
                    if a < i: Y += wt
                    elif a == i: X += wt
                    elif b < i: WM += wt
                    elif b == i: Y += wt
                    else: X += wt
            checked += 1
            if not (2 * X <= WM and 2 * Y <= WM): prem_bad += 1
print(f"(c) crude-bound premises 2X <= W_M and 2Y <= W_M for m = 8..14, 3 <= i <= (m+1)/2: {checked} cases, failures {prem_bad}")
