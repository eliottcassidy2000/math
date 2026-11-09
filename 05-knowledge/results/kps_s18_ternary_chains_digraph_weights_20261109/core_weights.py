#!/usr/bin/env python3
"""Weighted cycle decompositions of the strong semicomplete cores with H = 7 and H = 21 (THM-5071 (4)).

For a loopless semicomplete digraph D (pairs: one arc = T-pair, both arcs = B-pair), THM-5071 (2) gives
    H(D) = sum over packings of vertex-disjoint admissible cycles of prod w(C),
admissible = directed cycles of D; w(C) = 2 if C uses only T-pairs and has odd length, 0 if it uses only
T-pairs and has even length, 1 if it uses at least one B-pair (any length >= 2).

This script enumerates every semicomplete digraph of order n <= 5, keeps the strong ones with H in {7, 21},
reduces them to isomorphism classes, and prints for each class: b (reciprocal pairs), the numbers of
admissible cycles of each kind, and the packing decomposition of H.
"""
import itertools
import sys


def ham(n, A):
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, N):
        row = dp[s]
        for v in range(n):
            c = row[v]
            if c:
                for w in range(n):
                    if not s >> w & 1 and A[v][w]:
                        dp[s | 1 << w][w] += c
    return sum(dp[N - 1])


def strong(n, A):
    for s in range(n):
        seen = {s}
        stack = [s]
        while stack:
            v = stack.pop()
            for w in range(n):
                if A[v][w] and w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) < n:
            return False
    return True


def canon(n, A):
    best = None
    for p in itertools.permutations(range(n)):
        code = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or code < best:
            best = code
    return best


def cycles(n, A):
    """All directed cycles (as tuples starting at their minimum vertex), length >= 2."""
    out = []
    for k in range(2, n + 1):
        for sub in itertools.combinations(range(n), k):
            s = sub[0]
            for perm in itertools.permutations(sub[1:]):
                seq = (s,) + perm
                if all(A[seq[i]][seq[(i + 1) % k]] for i in range(k)):
                    out.append(seq)
    return out


def weight(A, cyc):
    k = len(cyc)
    usesB = any(A[cyc[i]][cyc[(i + 1) % k]] and A[cyc[(i + 1) % k]][cyc[i]] for i in range(k))
    if usesB:
        return 1, "B"
    return (2, "Todd") if k % 2 else (0, "Teven")


def decomposition(n, A):
    cyc = cycles(n, A)
    wts = [weight(A, c) for c in cyc]
    masks = [sum(1 << v for v in c) for c in cyc]
    adm = [(m, w, kind, c) for m, (w, kind), c in zip(masks, wts, cyc) if w != 0]
    total = 0
    terms = {}

    def rec(i, used, prod, sig):
        nonlocal total
        if i == len(adm):
            total += prod
            terms[sig] = terms.get(sig, 0) + prod
            return
        rec(i + 1, used, prod, sig)
        m, w, kind, _ = adm[i]
        if not m & used:
            rec(i + 1, used | m, prod * w, tuple(sorted(sig + (kind,))))

    rec(0, 0, 1, ())
    counts = {}
    for (w, kind) in wts:
        counts[kind] = counts.get(kind, 0) + 1
    return total, counts, terms


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for n in range(2, nmax + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        classes = {}
        for code in range(3 ** len(pairs)):
            A = [[0] * n for _ in range(n)]
            c = code
            for (i, j) in pairs:
                st = c % 3 + 1
                c //= 3
                A[i][j] = st & 1
                A[j][i] = st >> 1 & 1
            h = ham(n, A)
            if h in (7, 21) and strong(n, A):
                key = canon(n, A)
                if key not in classes:
                    classes[key] = (h, [row[:] for row in A])
        for key, (h, A) in sorted(classes.items(), key=lambda t: t[1][0]):
            total, counts, terms = decomposition(n, A)
            assert total == h, (total, h)
            b = sum(1 for i in range(n) for j in range(i + 1, n) if A[i][j] and A[j][i])
            rows = "/".join("".join(map(str, r)) for r in A)
            parts = " + ".join(f"{v}[{'+'.join(k) if k else 'empty'}]" for k, v in sorted(terms.items(),
                              key=lambda t: (len(t[0]), t[0])))
            print(f"n={n} H={h} b={b} rows={rows} cycles={counts} H = {parts}")


if __name__ == "__main__":
    main()
