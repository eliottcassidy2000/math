"""aud_single_apex.py - audit check of single_apex_spectrum.out (opus-2026-11-08-S86), written independently.

One-apex tournament: transitive chain v_1 -> ... -> v_m plus an apex x; bit w_i = 1 iff x -> v_i.  It is strong iff
w_1 = 1 and w_m = 0 ("reduced word").  Every odd cycle uses x; their number is computed here by brute-force cycle
enumeration for m <= 9 and by a chain-path DP for all m <= 16 (the two are compared), and the minima per reduced length
and the non-attained alpha_1 values below 8192 are recomputed.
usage: python aud_single_apex.py
"""
from itertools import permutations


def alpha1_bruteforce(w):
    m = len(w)
    n = m + 1                      # vertex 0 = apex, vertex i = v_i
    adj = [[False] * n for _ in range(n)]
    for i in range(1, n):
        for j in range(i + 1, n):
            adj[i][j] = True       # v_i -> v_j for i < j
    for i in range(1, n):
        if w[i - 1]:
            adj[0][i] = True
        else:
            adj[i][0] = True
    # every cycle contains the apex (the chain is transitive); enumerate simple paths from the apex
    cnt = 0

    def dfs(v, seen, length):
        nonlocal cnt
        for u in range(n):
            if adj[v][u]:
                if u == 0 and length >= 3 and length % 2 == 1:
                    cnt += 1
                elif u != 0 and not seen[u]:
                    seen[u] = True
                    dfs(u, seen, length + 1)
                    seen[u] = False
    seen = [False] * n
    seen[0] = True
    dfs(0, seen, 1)
    return cnt


def chain_even_counts(dmax):
    E = [0] * (dmax + 1)
    for d in range(1, dmax + 1):
        even = [0] * (d + 1)
        odd = [0] * (d + 1)
        odd[0] = 1
        for j in range(1, d + 1):
            even[j] = sum(odd[:j])
            odd[j] = sum(even[:j])
        E[d] = even[d]
    return E


def main():
    E = chain_even_counts(16)

    def a1(w):
        return sum(E[j - i] for i in range(len(w)) if w[i] for j in range(i + 1, len(w)) if not w[j])
    bad = 0
    for m in range(1, 10):
        for mask in range(1 << m):
            w = [(mask >> i) & 1 for i in range(m)]
            if alpha1_bruteforce(w) != a1(w):
                bad += 1
    print(f"brute force vs chain-path DP on all words of length <= 9: {bad} disagreements")
    vals = set()
    for m in range(2, 17):
        mn = None
        cnt = 0
        for mask in range(1 << (m - 2)):
            w = [1] + [(mask >> i) & 1 for i in range(m - 2)] + [0]
            v = a1(w)
            cnt += 1
            mn = v if mn is None else min(mn, v)
            if v < 8192:
                vals.add(v)
        print(f"reduced length m = {m:2d} (order {m + 1:2d}): {cnt} words, min alpha_1 = {mn}, 2^(m-3) = {2 ** (m - 3) if m >= 3 else 0.5}")
    vals.add(0)
    holes = [x for x in range(8192) if x not in vals]
    print(f"alpha_1 values in [0, 8192) not attained by one-apex tournaments: {len(holes)}; first: {holes[:20]}")


if __name__ == "__main__":
    main()
