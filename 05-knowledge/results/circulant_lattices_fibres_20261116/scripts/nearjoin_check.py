#!/usr/bin/env python3
"""nearjoin_check.py: (1) near-join T' = T -> M_k -> z with z -> t_0: verify H(T') = h + 2^k B_{t_0}(T) and B_z(T') = H(T'),
and whether z's strong component carries the whole count; (2) the identities sigma_u = #backward-junction two-path covers
and sigma_z = e_{t_0}(S - z) for sole-exit z; (3) the converse test: strong S with a digit-0 vertex u -- is H(S) a near-join
value h(S-u) + 2^k B_t(S-u)?  Orders <= 6 for (1), <= 7 for (2), (3)."""
import subprocess, itertools
def parse(s, n):
    adj = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj
def classes(n): return [parse(s, n) for s in subprocess.run(['gentourng','-q',str(n)],capture_output=True,text=True).stdout.split()]
def ham_all(adj):
    n = len(adj); N = 1 << n; dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    for m in range(1, N):
        for l in range(n):
            a = dp[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dp[m|(1<<v)][v] += a
    H = [sum(dp[m]) for m in range(N)]; H[0] = 1
    return dp, H
def starts(adj, t):
    n = len(adj); N = 1 << n; dps = [[0]*n for _ in range(N)]; dps[1<<t][t] = 1
    for m in range(1, N):
        if not (m >> t) & 1: continue
        for l in range(n):
            a = dps[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dps[m|(1<<v)][v] += a
    return [sum(dps[m]) for m in range(N)], dps
def B_of(adj, t):
    n = len(adj); N = 1 << n; dp, H = ham_all(adj); st, _ = starts(adj, t)
    return sum(st[m] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1), H[N-1], dp[N-1][t]
def strong(adj):
    n = len(adj)
    def reach(s, A):
        r = {s}; f = [s]
        while f:
            x = f.pop()
            for y in range(n):
                if A[x][y] and y not in r: r.add(y); f.append(y)
        return r
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, adj)) == n and len(reach(0, rev)) == n
def sub(adj, S):
    S = sorted(S); return [[adj[i][j] for j in S] for i in S]
def scc_of(adj, v):
    n = len(adj)
    def reach(s, A):
        r = {s}; f = [s]
        while f:
            x = f.pop()
            for y in range(n):
                if A[x][y] and y not in r: r.add(y); f.append(y)
        return r
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    return reach(v, adj) & reach(v, rev)
# (1)
bad_formula = bad_Bz = 0; checks = 0; z_full = 0; z_not_full = 0; z_not_full_nonstrongT = 0
for n in (2, 3, 4, 5, 6):
    for T in classes(n):
        for t0 in range(n):
            B0, h, _ = B_of(T, t0)
            for k in (1, 2, 3):
                m = n + k + 1; adj = [[0]*m for _ in range(m)]
                for i in range(n):
                    for j in range(n): adj[i][j] = T[i][j]
                for i in range(n):
                    for j in range(n, m): adj[i][j] = 1       # T -> M_k, T -> z
                for i in range(n, m):
                    for j in range(i+1, m): adj[i][j] = 1     # chain, M_k -> z
                z = m - 1; adj[t0][z] = 0; adj[z][t0] = 1       # reverse z -> t0
                Bz, H, _ = B_of(adj, z); checks += 1
                if H != h + (2**k) * B0: bad_formula += 1
                if Bz != H: bad_Bz += 1
                comp = scc_of(adj, z); Hc = ham_all(sub(adj, comp))[1][-1]
                if Hc == H: z_full += 1
                else:
                    z_not_full += 1
                    if not strong(T): z_not_full_nonstrongT += 1
print(f"(1) near-joins checked {checks}: formula failures {bad_formula}, B_z = H failures {bad_Bz}; z's component carries the full count in {z_full}, not in {z_not_full} (of which base T non-strong: {z_not_full_nonstrongT})")
# (2) identities
bad_cov = bad_sole = 0; tested = 0; tested_sole = 0
for n in (3, 4, 5, 6, 7):
    for S in classes(n):
        dp, H = ham_all(S); N = 1 << n
        for u in range(n):
            B, h, e = B_of(S, u); sigma = B - h
            # backward-junction covers: pairs (P1 Hamiltonian path of S[A], A nonempty subset of V-u, P2 from u on V-A) with u beating last(P1)
            st, dps = starts(S, u)
            cnt = 0
            for A in range(1, N):
                if (A >> u) & 1: continue
                # paths of S[A] ending at x with u -> x
                for x in range(n):
                    if (A >> x) & 1 and S[u][x] and dp[A][x]: cnt += dp[A][x] * st[(N-1) & ~A]
            tested += 1
            if cnt != sigma: bad_cov += 1
            outs = [x for x in range(n) if S[u][x]]
            if len(outs) == 1:
                t0 = outs[0]; rest = [v for v in range(n) if v != u]; Sr = sub(S, rest); t0r = rest.index(t0)
                dpr, Hr = ham_all(Sr); tested_sole += 1
                if dpr[(1 << (n-1)) - 1][t0r] != sigma: bad_sole += 1
print(f"(2) sigma_u = backward-junction covers: {tested} checks, {bad_cov} failures; sole-exit identity sigma_z = e_t0(S - z): {tested_sole} checks, {bad_sole} failures")
# (3) converse test
withd0 = 0; reproduced = 0
for n in (4, 5, 6, 7):
    for S in classes(n):
        if not strong(S): continue
        for u in range(n):
            B, h, e = B_of(S, u)
            if B != h: continue
            withd0 += 1
            rest = [v for v in range(n) if v != u]; Sr = sub(S, rest); hr = ham_all(Sr)[1][-1]
            ok = False
            for t in range(n-1):
                Bt, _, _ = B_of(Sr, t)
                for k in range(1, 8):
                    if hr + (2**k) * Bt == h: ok = True
            if ok: reproduced += 1
print(f"(3) strong (S, u) with a digit-0 vertex u (orders 4-7): {withd0}; H(S) = h(S-u) + 2^k B_t(S-u) for some t, k: {reproduced}")
