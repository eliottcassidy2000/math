#!/usr/bin/env python3
"""Lane thinorbit (S93), finishing pass.  Check of the PROOF of the parity lemma (PROOFS.txt, Lemma P), not only of
its statement.

Lemma P (base 2): for every marked block (B, w, u2),
    N(B, w, u2) = [[1, 0, 0], [e, A0, A1], [o, A1, A0]]  (mod 2),
where e / o = number of directed paths of B starting at u2 with an even / odd number of vertices, and
A0 / A1 = number of pairs (Q, R) with Q a directed path of B from u2 ending at w, R a directed path of
B - V(Q) or the empty path, and |V(Q)| + |V(R)| even / odd.
(Proof: Redei gives H = 1 mod 2 for every tournament, so (h, E, O) = (1, #even paths from the port, #odd paths from
the port) mod 2; a path of the glued tournament from u2 either stays in B or is Q, then a path of T from the old
port, then R; the three base states are independent mod 2.)

This script recomputes (e, o, A0, A1) by brute-force path enumeration for every representative block of a library
file (MATFILE.reps: nine integers | k rows w u2) and compares with the stored matrix mod 2.
Usage: finish_parity_lemma_check.py MATFILE
"""
import sys


def paths_from(adj, k, start, allowed):
    """yield vertex masks (and last vertex) of all directed paths from start inside the vertex set allowed (mask)."""
    stack = [(1 << start, start)]
    while stack:
        mask, v = stack.pop()
        yield mask, v
        for t in range(k):
            if adj[v][t] and (allowed >> t) & 1 and not (mask >> t) & 1:
                stack.append((mask | (1 << t), t))


def main():
    path = sys.argv[1]
    nb = 0
    bad = 0
    pat = {}
    for line in open(path + ".reps"):
        a, b = line.split("|")
        N = list(map(int, a.split()))
        ks, rows, ws, us = b.split()
        k, w, u2 = int(ks), int(ws), int(us)
        adj = [[int(rows[i * k + j]) for j in range(k)] for i in range(k)]
        full = (1 << k) - 1
        e = o = 0
        A = [0, 0]
        # number of directed paths (nonempty) by parity of order inside a vertex set, memoised
        memo = {}

        def tails(allowed):
            if allowed in memo:
                return memo[allowed]
            c = [1, 0]                       # the empty path
            for s in range(k):
                if (allowed >> s) & 1:
                    for mask, _ in paths_from(adj, k, s, allowed):
                        c[bin(mask).count("1") % 2] += 1
            memo[allowed] = c
            return c

        for mask, v in paths_from(adj, k, u2, full):
            sz = bin(mask).count("1")
            if sz % 2 == 0:
                e += 1
            else:
                o += 1
            if v == w:
                c = tails(full & ~mask)
                A[sz % 2] += c[0]
                A[(sz + 1) % 2] += c[1]
        pred = [1, 0, 0, e % 2, A[0] % 2, A[1] % 2, o % 2, A[1] % 2, A[0] % 2]
        nb += 1
        if [x % 2 for x in N] != pred:
            bad += 1
            if bad <= 5:
                print("MISMATCH", N, b.strip(), pred)
        key = (A[0] % 2, A[1] % 2)
        pat[key] = pat.get(key, 0) + 1
    print(f"{path}: {nb} representative blocks; matrix mod 2 = [[1,0,0],[e,A0,A1],[o,A1,A0]] with the path counts of the "
          f"proof: {nb - bad} agree, {bad} disagree; (A0, A1) mod 2 counts: {dict(sorted(pat.items()))}")


if __name__ == "__main__":
    main()
