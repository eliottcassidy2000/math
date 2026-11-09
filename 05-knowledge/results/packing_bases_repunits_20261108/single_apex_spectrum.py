"""Single-apex tournaments: a transitive tournament v_1 -> v_2 -> ... -> v_m (v_i -> v_j for i < j) plus one apex a,
with a -> v_i iff w_i = 1.  Every cycle passes through a (the rest is acyclic), so no two odd cycles are disjoint:
alpha_2 = 0 and, by OCF, H = 1 + 2 alpha_1 and I(Omega, b) = 1 + b alpha_1 for every b.

alpha_1(w) = number of odd cycles = sum over i < j with w_i = 1, w_j = 0 of E(j - i - 1), where E(0) = 1 and
E(r) = 2^(r-1) (number of even-size subsets of an r-set): the cycle is a -> v_i -> (even number of intermediate
vertices, increasing) -> v_j -> a.

Leading 0s and trailing 1s contribute nothing, so only reduced words 1 u 0 matter, and a reduced word of length m has
alpha_1 >= E(m - 2) >= 2^(m - 3).  Hence all values alpha_1 <= X come from reduced words of length <= log2(X) + 3.

This script enumerates every reduced word of length <= MMAX, records the attained alpha_1 values below 2^(MMAX-3)
(complete in that range), cross-checks alpha_1 against a brute-force cycle count and H against a Hamiltonian-path
DP for small m, and reports the values NOT attained.
usage: python single_apex_spectrum.py MMAX
"""
import sys
import itertools
import numpy as np


def alpha1_word(w):
    """Exact alpha_1 by the pair formula."""
    m = len(w)
    s = 0
    for i in range(m):
        if w[i] != 1:
            continue
        for j in range(i + 1, m):
            if w[j] == 0:
                r = j - i - 1
                s += 1 if r == 0 else 1 << (r - 1)
    return s


def tournament_of_word(w):
    """Adjacency (out-neighbour bit masks) of the single-apex tournament; apex = vertex m."""
    m = len(w)
    n = m + 1
    out = [0] * n
    for i in range(m):
        for j in range(i + 1, m):
            out[i] |= 1 << j
    for i in range(m):
        if w[i] == 1:
            out[m] |= 1 << i
        else:
            out[i] |= 1 << m
    return n, out


def ham_paths(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            nb = out[v] & ~S
            while nb:
                u = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                dp[S | (1 << u)][u] += c
    return sum(dp[full])


def odd_cycles_brute(n, out):
    """Count directed odd cycles (each cycle once) by DP over subsets with a fixed least vertex."""
    total = 0
    for s in range(n):
        # paths starting at s using only vertices > s
        dp = {(1 << s, s): 1}
        frontier = [(1 << s, s)]
        counts = {}
        for _ in range(n):
            nxt = {}
            for (S, v), c in dp.items():
                nb = out[v] & ~S
                while nb:
                    u = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    if u <= s:
                        continue
                    key = (S | (1 << u), u)
                    nxt[key] = nxt.get(key, 0) + c
            for (S, v), c in dp.items():
                if (out[v] >> s) & 1:
                    k = bin(S).count("1")
                    if k >= 3 and k % 2 == 1:
                        total += c
            dp = nxt
            if not dp:
                break
    return total


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    # cross-check on all words of length <= 9
    checked = 0
    for m in range(0, 10):
        for w in itertools.product((0, 1), repeat=m):
            a1 = alpha1_word(w)
            n, out = tournament_of_word(w)
            assert odd_cycles_brute(n, out) == a1, (w,)
            assert ham_paths(n, out) == 1 + 2 * a1, (w,)
            checked += 1
    print(f"cross-check: alpha_1 formula = brute-force odd-cycle count and H = 1 + 2 alpha_1 on all {checked} words of length <= 9")

    bound = 1 << (mmax - 3)          # complete for alpha_1 < bound
    seen = np.zeros(bound, dtype=bool)
    # reduced words 1 u 0 of length m >= 2; vectorised over the middle u
    for m in range(2, mmax + 1):
        k = m - 2                      # free middle letters
        U = np.arange(1 << k, dtype=np.int64)
        bits = ((U[:, None] >> np.arange(k - 1, -1, -1)) & 1).astype(np.int64) if k > 0 else np.zeros((1, 0), np.int64)
        W = np.concatenate([np.ones((U.size, 1), np.int64), bits, np.zeros((U.size, 1), np.int64)], axis=1)
        # alpha_1 = sum_{i<j} [w_i=1][w_j=0] E(j-i-1)
        E = np.array([1] + [1 << (r - 1) for r in range(1, m)], dtype=np.int64)
        val = np.zeros(U.size, dtype=np.int64)
        for i in range(m):
            wi = W[:, i]
            for j in range(i + 1, m):
                val += wi * (1 - W[:, j]) * E[j - i - 1]
        small = val[val < bound]
        seen[small] = True
        print(f"m = {m:2d}: {U.size} reduced words, min alpha_1 = {val.min()}", flush=True)
    seen[0] = True                      # the transitive tournament (empty reduced word)
    missing = np.flatnonzero(~seen)
    print(f"alpha_1 values in [0, {bound}) NOT attained by single-apex tournaments: {missing.tolist()[:60]}"
          f"{' ...' if missing.size > 60 else ''} (count {missing.size})")
    print(f"so H = 1 + 2 alpha_1: odd H in [1, {2 * bound}) not attained: {(1 + 2 * missing).tolist()[:60]}")


if __name__ == "__main__":
    main()
