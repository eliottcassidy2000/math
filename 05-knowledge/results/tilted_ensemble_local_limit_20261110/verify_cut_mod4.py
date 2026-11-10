#!/usr/bin/env python3
"""Independent brute-force check of the statements of THM-5127 (C) on small tournaments.

For every tournament R on m <= 5 labelled vertices (all 2^C(m,2) of them) and 400 seeded random R on 6 and 7 vertices:
  (1) cut quadratic: H(R + v_X) = H(R) + sum_{i in X} a_i - sum_{{i,j} in X} mu(i,j)   for all 2^m out-sets X,
      where a_i = H(R + v_{i}) - H(R) and mu(i,j) = a_i + a_j - (H(R + v_{i,j}) - H(R));
  (2) mu(i,j) = M[i][j] + M[j][i] with M[j][i] = #ordered two-path covers (P1 ending at j, P2 starting at i);
  (3) mu(i,j) is even and >= 2;
  (4) #{X : q_R(X) = 0} - #{X : q_R(X) = 1} is 0 or +-2^(m - rank/2), q_R(X) = (H(R+v_X) - H(R))/2 mod 2,
      rank = F_2-rank of beta = (mu/2 mod 2).
H is counted by plain permutation-free subset DP written independently of the C programs.
"""
import itertools, random, sys

def ham(n, out):
    """number of Hamiltonian paths; out[v] = bitmask of out-neighbours"""
    full = (1 << n) - 1
    f = [[0]*n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = f[S][v]
            if not c: continue
            nxt = out[v] & ~S
            while nxt:
                w = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                f[S | (1 << w)][w] += c
    return sum(f[full])

def two_cover_matrix(m, out):
    """M[j][i] by direct enumeration over subsets with independent DP tables"""
    full = (1 << m) - 1
    # e[S][v]: paths covering S ending at v ; s[S][v]: paths covering S starting at v
    e = [[0]*m for _ in range(1 << m)]; s = [[0]*m for _ in range(1 << m)]
    inn = [sum(1 << u for u in range(m) if out[u] >> v & 1) for v in range(m)]
    for S in range(1, 1 << m):
        for v in range(m):
            if not S >> v & 1: continue
            R = S & ~(1 << v)
            if R == 0: e[S][v] = s[S][v] = 1; continue
            e[S][v] = sum(e[R][u] for u in range(m) if R >> u & 1 and inn[v] >> u & 1)
            s[S][v] = sum(s[R][u] for u in range(m) if R >> u & 1 and out[v] >> u & 1)
    M = [[0]*m for _ in range(m)]
    for j in range(m):
        for i in range(m):
            if i == j: continue
            for S in range(1, full):
                if S >> j & 1 and not S >> i & 1:
                    M[j][i] += e[S][j] * s[full ^ S][i]
    return M

def rank_f2(rows):
    rows = rows[:]; r = 0
    for c in range(len(rows)):
        p = next((k for k in range(r, len(rows)) if rows[k] >> c & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for k in range(len(rows)):
            if k != r and rows[k] >> c & 1: rows[k] ^= rows[r]
        r += 1
    return r

def check(m, out):
    H0 = ham(m, out)
    def ext(X):
        o = out[:] + [X]
        for u in range(m):
            if not X >> u & 1: o[u] |= 1 << m
        return ham(m + 1, o)
    vals = [ext(X) for X in range(1 << m)]
    assert vals[0] == H0 and vals[(1 << m) - 1] == H0
    a = [vals[1 << i] - H0 for i in range(m)]
    mu = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j: mu[i][j] = a[i] + a[j] - (vals[(1 << i) | (1 << j)] - H0)
    M = two_cover_matrix(m, out)
    for i in range(m):
        for j in range(m):
            if i != j:
                assert mu[i][j] == M[i][j] + M[j][i], "pair term is not the two-cover count"
                assert mu[i][j] % 2 == 0 and mu[i][j] >= 2, "mu not even >= 2"
    for X in range(1 << m):
        pred = H0 + sum(a[i] for i in range(m) if X >> i & 1) - sum(mu[i][j] for i in range(m) for j in range(i) if X >> i & 1 and X >> j & 1)
        assert pred == vals[X], "cut quadratic fails"
    rows = [sum(((mu[i][j] // 2) & 1) << j for j in range(m) if j != i) for i in range(m)]
    rk = rank_f2(rows)
    diff = sum(1 if ((v - H0) // 2) % 2 == 0 else -1 for v in vals)
    assert rk % 2 == 0 and abs(diff) in (0, 1 << (m - rk // 2)), "mod-4 Gauss sum prediction fails"
    return rk, diff

def all_tournaments(m):
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    for code in range(1 << len(pairs)):
        out = [0]*m
        for k, (i, j) in enumerate(pairs):
            if code >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out

if __name__ == "__main__":
    total = 0
    for m in range(1, 6):
        cnt = 0
        for out in all_tournaments(m):
            check(m, out); cnt += 1
        print(f"m={m}: all {cnt} labelled tournaments pass (cut quadratic, two-cover pair term, mu even >= 2, mod-4 Gauss sum)")
        total += cnt
    rng = random.Random(20261110)
    for m, reps in ((6, 300), (7, 100)):
        for _ in range(reps):
            out = [0]*m
            for i in range(m):
                for j in range(i + 1, m):
                    if rng.random() < 0.5: out[i] |= 1 << j
                    else: out[j] |= 1 << i
            check(m, out)
        print(f"m={m}: {reps} seeded random tournaments pass")
        total += reps
    print("TOTAL tournaments checked:", total, " failures: 0")
