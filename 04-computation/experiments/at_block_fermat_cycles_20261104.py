"""(AT_m, A = {m-1, m}) as a block in opus-S75's core theory (Lemma B4).
Checks for 4 <= m <= MMAX:
  (i)  the block conditions: AT_m strong, |A| >= 2, every vertex of W - A has out-degree >= 2 in W, and no non-empty Z
       inside W has (#arcs Z -> W - Z) + |Z & A| = 1 (all 2^m - 1 subsets Z);
  (ii) the type: end_(m-1)(AT_m) + end_m(AT_m) = 2^(m-3) + 1;
  (iii) the smallest tournament built on it: D = W + u (u beaten by A, u beats W - A), T = D + l (only u -> l leaves D;
       l beats W). T has no sole exit and hc(T) = 2^(m-3) + 1 (Hamiltonian cycles counted by DP, independent of (ii)).
kind-pasteur-2026-11-04-S10."""
import sys

def at_rows(m):  # vertices 0..m-1 stand for 1..m
    rows = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if (i, j) == (0, m - 1): rows[m - 1] |= 1 << 0
            else: rows[i] |= 1 << j
    return rows

def ham_paths_end(rows):
    n = len(rows); full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nx = rows[v] & ~S
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[S | 1 << w][w] += c
    return dp[full]

def hcycles(rows):
    n = len(rows); M = 1 << n
    dp = [[0] * n for _ in range(M)]
    dp[1][0] = 1
    for S in range(1, M, 2):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nx = rows[v] & ~S
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[S | 1 << w][w] += c
    return sum(dp[M - 1][v] for v in range(1, n) if rows[v] & 1)

def strong(rows):
    n = len(rows)
    for s in range(n):
        seen = 1 << s; stack = [s]
        while stack:
            v = stack.pop(); nx = rows[v] & ~seen
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1; seen |= 1 << w; stack.append(w)
        if seen != (1 << n) - 1: return False
    return True

def main(mmax):
    for m in range(4, mmax + 1):
        W = at_rows(m); A = (1 << (m - 2)) | (1 << (m - 1))
        cond_out = all(bin(W[v]).count("1") >= 2 for v in range(m) if not (A >> v & 1))
        bad_Z = 0
        for Z in range(1, 1 << m):
            leave = sum(bin(W[v] & ~Z & ((1 << m) - 1)).count("1") for v in range(m) if Z >> v & 1)
            if leave + bin(Z & A).count("1") == 1: bad_Z += 1
        ends = ham_paths_end(W)
        typ = ends[m - 2] + ends[m - 1]
        # T: W on 0..m-1, u = m, l = m+1
        n = m + 2; u, l = m, m + 1
        T = [r for r in W] + [0, 0]
        for v in range(m):
            if A >> v & 1: T[v] |= 1 << u          # A beats u
            else: T[u] |= 1 << v                   # u beats W - A
        T[u] |= 1 << l                             # the one arc leaving D
        for v in range(m): T[l] |= 1 << v          # l beats all of W
        assert all(bin(T[v]).count("1") + sum(1 for w in range(n) if T[w] >> v & 1) == n - 1 for v in range(n))
        min_out = min(bin(r).count("1") for r in T)
        hc = hcycles(T)
        ok = strong(W) and cond_out and bad_Z == 0 and typ == 2 ** (m - 3) + 1 and min_out >= 2 and hc == typ
        print(f"m={m:2d} block conditions: strong={strong(W)} outdeg>=2 off A={cond_out} bad Z={bad_Z}; "
              f"type end_(m-1)+end_m={typ} (2^(m-3)+1={2 ** (m - 3) + 1}); T order {n}: min out-degree {min_out}, hc={hc}  ok={ok}")
        assert ok

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
