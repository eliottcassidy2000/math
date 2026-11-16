#!/usr/bin/env python3
"""padding_check.py: the ordinal-sum laws for the root response B_t, the end count e_t and the charge a_t = B_t - e_t - h + 1.
For U = A -> B (every vertex of A beats every vertex of B):
   t in B: B_t(U) = h(A) B_t(B), e_t(U) = h(A) e_t(B), so a_t(U) = h(A)(a_t(B) - 1) + 1 and sigma_t(U) = h(A) B_t(B) - h(A)h(B);
   t in A: B_t(U) = B_t(A) P_B(2), e_t(U) = 0, a_t(U) = B_t(A) P_B(2) - h(A) h(B) + 1, where P_B(2) = sum_S h(B[S]) h(B[V-S]).
Corollaries: transitive prefix TT_j -> S leaves B_t, e_t, h unchanged; transitive suffix S -> TT_k gives B_t = 2^k B_t(S),
e_t = 0: digit sigma_t = 2^k B_t(S) - h, charge a_t = 2^k B_t(S) - h + 1.  Checked on random pairs."""
import random, itertools
random.seed(7)
def ham_all(adj, n):
    N = 1 << n; dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    for m in range(1, N):
        for l in range(n):
            a = dp[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dp[m|(1<<v)][v] += a
    H = [sum(dp[m]) for m in range(N)]; H[0] = 1
    return dp, H
def Bt_et(adj, n, t):
    dp, H = ham_all(adj, n); N = 1 << n
    dps = [[0]*n for _ in range(N)]; dps[1<<t][t] = 1
    for m in range(1, N):
        if not (m >> t) & 1: continue
        for l in range(n):
            a = dps[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dps[m|(1<<v)][v] += a
    B = sum(sum(dps[m]) * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1)
    e = dp[N-1][t]; h = H[N-1]
    P2 = sum(H[S] * H[(N-1) & ~S] for S in range(N))
    return h, B, e, P2
def random_tournament(n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: adj[i][j] = 1
            else: adj[j][i] = 1
    return adj
def ordinal(A, B):
    nA, nB = len(A), len(B); n = nA + nB; adj = [[0]*n for _ in range(n)]
    for i in range(nA):
        for j in range(nA): adj[i][j] = A[i][j]
    for i in range(nB):
        for j in range(nB): adj[nA+i][nA+j] = B[i][j]
    for i in range(nA):
        for j in range(nB): adj[i][nA+j] = 1
    return adj
bad = 0; checks = 0
for trial in range(60):
    nA, nB = random.randint(1, 4), random.randint(1, 4); A, B = random_tournament(nA), random_tournament(nB); U = ordinal(A, B)
    hA = ham_all(A, nA)[1][-1]; hB = ham_all(B, nB)[1][-1]
    for t in range(nB):
        h, Bu, eu, _ = Bt_et(U, nA+nB, nA+t); hb, Bb, eb, _ = Bt_et(B, nB, t); checks += 1
        if not (Bu == hA*Bb and eu == hA*eb and h == hA*hB): bad += 1
    for t in range(nA):
        h, Bu, eu, _ = Bt_et(U, nA+nB, t); ha, Ba, ea, _ = Bt_et(A, nA, t); _, _, _, P2B = Bt_et(B, nB, 0); checks += 1
        if not (Bu == Ba*P2B and eu == 0 and h == hA*hB): bad += 1
print(f"ordinal-sum laws: {checks} checks, {bad} failures")
# the transitive-suffix law for digits and charges
bad2 = 0
for trial in range(40):
    nS = random.randint(2, 5); S = random_tournament(nS)
    for k in (1, 2, 3):
        TT = [[1 if i < j else 0 for j in range(k)] for i in range(k)]; U = ordinal(S, TT)
        for t in range(nS):
            h, Bu, eu, _ = Bt_et(U, nS+k, t); hs, Bs, es, _ = Bt_et(S, nS, t)
            if not (Bu == (2**k)*Bs and eu == 0 and h == hs): bad2 += 1
print("transitive-suffix law B_t -> 2^k B_t, e_t -> 0:", "holds" if bad2 == 0 else f"{bad2} failures")
