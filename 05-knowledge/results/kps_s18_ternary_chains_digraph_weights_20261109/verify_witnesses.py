#!/usr/bin/env python3
"""Independent check of the base-b witnesses for HYP-9530 (B) (THM-5070 (6)).

Input lines (from the hill-climb `hunt`):  base B m=M target=V : FOUND n=N rows ROW1 ROW2 ...
For each witness this script checks that the rows form a tournament and recomputes from scratch:
- c[S] = the number of directed Hamiltonian cycles on each odd vertex set S (|S| >= 3), by a DP from min(S);
- I(Omega, x) = sum over set partitions into packed odd-cycle blocks (independence polynomial of Omega);
- the Hamiltonian path count H, checked against OCF (H = I(Omega, 2)).
It reports any witness whose I(Omega, B) differs from V = 1 + B*M.

Usage: verify_witnesses.py FILE [FILE ...]
"""
import sys


def check(n, A, base):
    N = 1 << n
    cyc = [0] * N
    for s in range(n):
        low = (1 << s) - 1
        dp = {}
        dp[(1 << s, s)] = 1
        # process subsets containing s with min s in increasing order of size
        order = sorted((S for S in range(N) if S >> s & 1 and not S & low), key=lambda S: bin(S).count("1"))
        for S in order:
            for v in range(n):
                c = dp.get((S, v), 0)
                if not c:
                    continue
                if A[v][s] and bin(S).count("1") >= 3:
                    cyc[S] += c
                for w in range(s + 1, n):
                    if not S >> w & 1 and A[v][w]:
                        key = (S | 1 << w, w)
                        dp[key] = dp.get(key, 0) + c

    def I(x):
        F = [0] * N
        F[0] = 1
        for S in range(1, N):
            m = (S & -S).bit_length() - 1
            rest = S & ~(1 << m)
            val = F[rest]
            T = rest
            while True:
                U = T | (1 << m)
                k = bin(U).count("1")
                if k % 2 == 1 and k >= 3 and cyc[U]:
                    val += x * cyc[U] * F[S & ~U]
                if T == 0:
                    break
                T = (T - 1) & rest
            F[S] = val
        return F[N - 1]

    # Hamiltonian paths
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if not S >> w & 1 and A[v][w]:
                        dp[S | 1 << w][w] += c
    H = sum(dp[N - 1])
    return I(base), I(2), H


def main():
    total = bad = 0
    for path in sys.argv[1:]:
        for line in open(path):
            if "FOUND" not in line:
                continue
            parts = line.split()
            base = int(parts[1])
            target = int(parts[3].split("=")[1])
            n = int(line.split("n=")[1].split()[0])
            rows = line.split("rows")[1].split()
            A = [[int(c) for c in r] for r in rows]
            assert len(A) == n and all(len(r) == n for r in A)
            assert all(A[i][i] == 0 for i in range(n))
            assert all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n) if i != j)
            ib, i2, H = check(n, A, base)
            total += 1
            ok = ib == target and i2 == H
            if not ok:
                bad += 1
                print("MISMATCH", line.strip()[:80], ib, i2, H)
            else:
                print(f"base {base} alpha-sum {(target - 1) // base}: verified (order {n}, I(Omega,{base}) = {ib}, H = {H})")
    print(f"{total} witnesses checked, {bad} mismatches")


if __name__ == "__main__":
    main()
