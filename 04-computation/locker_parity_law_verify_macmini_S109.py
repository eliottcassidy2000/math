#!/usr/bin/env python3
"""locker_parity_law_verify_macmini_S109.py -- mac-mini-2026-07-15-S109.
INDEPENDENT verification of the locker-parity refutation:
(A) t_m (odd cycles through top vertex m of D_m) via subset-DP (independent of the S109 DFS),
    with exit-cell breakdown; m = 3..17.
(B) H(D_n) by Held-Karp path DP, n = 3..15; check H mod 4 == (1 + 2*alpha_1) mod 4 with
    alpha_1(D_n) = sum_{m<=n} t_m.
(C) the divisor-pairing parity table: for composite m, does exit-cell parity match d <-> m/d?
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

def build_out(n):
    out = [[False]*(n+1) for _ in range(n+1)]
    for u in range(1, n+1):
        for v in range(u+1, n+1):
            if v % u == 0 or v == u+1: out[v][u] = True
            else: out[u][v] = True
    return out

def t_m_subsetDP(m):
    """odd cycles through m, by exit b: paths b -> u inside {2..m-1}, u -> m an arc,
    cycle edges = |S|+1 where S = path vertex set; odd cycle <=> |S| even."""
    out = build_out(m)
    verts = list(range(2, m))                    # vertex 1 never on a cycle (outdeg 0)
    idx = {v: i for i, v in enumerate(verts)}
    V = len(verts)
    exits = [b for b in range(2, m) if out[m][b]]        # m -> b
    entries = [u for u in range(2, m) if out[u][m]]      # u -> m
    entry_set = set(entries)
    total_by_exit = {}
    for b in exits:
        dp = [dict() for _ in range(1 << V)]     # dp[S][v] = #paths b->v with vertex set S
        b_i = idx[b]
        dp[1 << b_i][b] = 1
        odd_count = 0
        for S in range(1 << V):
            if not dp[S]: continue
            popc = bin(S).count("1")
            for v, c in dp[S].items():
                if v in entry_set and popc % 2 == 0:     # |S| even -> odd cycle
                    odd_count += c
                for w in verts:
                    wi = idx[w]
                    if not (S >> wi) & 1 and out[v][w]:
                        dp[S | (1 << wi)][w] = dp[S | (1 << wi)].get(w, 0) + c
        # NOTE: odd_count added per (S,v) BEFORE extension -- but every dp cell is visited
        # exactly once when S reached, so each path counted once. Correct.
        total_by_exit[b] = odd_count
    return total_by_exit

def ham_paths(n):
    out = build_out(n)
    N = n
    adj = [0]*N
    for v in range(1, n+1):
        for w in range(1, n+1):
            if v != w and out[v][w]: adj[v-1] |= (1 << (w-1))
    dp = [[0]*N for _ in range(1 << N)]
    for v in range(N): dp[1 << v][v] = 1
    for S in range(1 << N):
        row = dp[S]
        for v in range(N):
            c = row[v]
            if not c or not (S >> v) & 1: continue
            av = adj[v]
            rest = av & ~S
            while rest:
                lb = rest & (-rest); u = lb.bit_length() - 1
                dp[S | lb][u] += c
                rest ^= lb
    return sum(dp[(1 << N) - 1][v] for v in range(N))

def main():
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 17
    HMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    print("(A) t_m via subset-DP, with exit cells:")
    t = {}
    for m in range(3, MMAX+1):
        be = t_m_subsetDP(m)
        t[m] = sum(be.values())
        divs = sorted(d for d in range(2, m) if m % d == 0)
        pairing = []
        for d in divs:
            e = m // d
            if d <= e and e in be:
                pairing.append(f"{d}<->{e}: {be[d]%2}/{be[e]%2} {'MATCH' if be[d]%2==be[e]%2 else 'MISMATCH'}")
        print(f"m={m:2d}: t_m={t[m]:10d} (par {t[m]%2})  by_exit={ {k: be[k] for k in sorted(be)} }")
        if pairing: print(f"       divisor-pairing parity: {pairing}   m-1 cell: {be.get(m-1,0)%2}")
    print()
    print("(B) H(D_n) Held-Karp vs OCF digit (H mod 4 = 1 + 2*alpha1 mod 4):")
    a1 = 0
    for n in range(3, HMAX+1):
        a1 += t.get(n, 0)
        if n <= MMAX:
            H = ham_paths(n)
            pred = (1 + 2*a1) % 4
            print(f"n={n:2d}: H={H:14d}  H mod 4 = {H%4}  alpha1={a1:8d} (par {a1%2})  "
                  f"pred {pred}  {'OK' if H%4==pred else 'FAIL'}"
                  + ("   <-- LAW BREAKS (H = 3 mod 4)" if H % 4 == 3 else ""))
    print("\nDONE")

if __name__ == "__main__":
    main()
