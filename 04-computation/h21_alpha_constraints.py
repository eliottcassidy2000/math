#!/usr/bin/env python3
"""
H=21 gap analysis: why is H=21 impossible at n=7?

H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ...
H=21 requires 2*alpha_1 + 4*alpha_2 + ... = 20

Exhaustive check of (alpha_1, alpha_2) constraints at n=7 shows:
  alpha_1=10, alpha_2 needed=0: IMPOSSIBLE (always alpha_2=2)
  alpha_1=8, alpha_2 needed=1: IMPOSSIBLE (always alpha_2=0)
  alpha_1=6, alpha_2 needed=2: IMPOSSIBLE (alpha_2 only 0 or 1)
  alpha_1=4, alpha_2 needed=3: IMPOSSIBLE (always alpha_2=0)

Conclusion: NO (alpha_1, alpha_2) pair gives H=21 at n=7.

Hypothesis: H=21 is a PERMANENT gap for all n.
Evidence: absent at n<=6 (exhaustive), n=7,8,9 (200k-500k samples each).

kind-pasteur-2026-03-07-S28
"""
import random
from itertools import combinations
from collections import defaultdict

def count_cycles_dp(A, n, cl):
    """Count directed cl-cycles using DP on vertex subsets."""
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for mask in range(1, 1 << cl):
            for v in range(cl):
                if not (mask & (1<<v)) or dp[mask][v] == 0: continue
                for u in range(cl):
                    if mask & (1<<u): continue
                    if sub[v][u]: dp[mask|(1<<u)][u] += dp[mask][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def get_all_cycles_as_sets(A, n):
    """Get all directed odd cycles as vertex sets (with multiplicity)."""
    cycles = []
    for cl in range(3, n+1, 2):
        for verts in combinations(range(n), cl):
            sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
            dp = [[0]*cl for _ in range(1 << cl)]
            dp[1][0] = 1
            for mask in range(1, 1 << cl):
                for v in range(cl):
                    if not (mask & (1<<v)) or dp[mask][v] == 0: continue
                    for u in range(cl):
                        if mask & (1<<u): continue
                        if sub[v][u]: dp[mask|(1<<u)][u] += dp[mask][v]
            full = (1 << cl) - 1
            count = sum(dp[full][v] for v in range(1, cl) if sub[v][0])
            for _ in range(count):
                cycles.append(set(verts))
    return cycles

def compute_alphas(cycles):
    """Compute alpha_0, alpha_1, alpha_2 from cycle list."""
    nc = len(cycles)
    alpha = defaultdict(int)
    alpha[0] = 1  # empty set

    if nc == 0:
        return dict(alpha)

    # Build adjacency
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True

    # Enumerate independent sets (limit to small nc)
    if nc > 25:
        # Just count alpha_1 and alpha_2
        alpha[1] = nc
        pairs = 0
        for i in range(nc):
            for j in range(i+1, nc):
                if not adj[i][j]:
                    pairs += 1
        alpha[2] = pairs
        return dict(alpha)

    for mask in range(1, 1 << nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok: break
        if ok:
            alpha[len(nodes)] = alpha.get(len(nodes), 0) + 1

    return dict(alpha)


def main():
    n = 7
    m = n*(n-1)//2
    rng = random.Random(42)

    print("=" * 70)
    print(f"H=21 GAP ANALYSIS AT n={n}")
    print("=" * 70)

    # Collect (alpha_1, alpha_2) -> count and H values
    a1_a2_to_H = {}
    num_samples = 5000

    print(f"\nSampling {num_samples} random tournaments...")
    for trial in range(num_samples):
        bits = rng.getrandbits(m)
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        cycles = get_all_cycles_as_sets(A, n)
        alphas = compute_alphas(cycles)
        alpha1 = alphas.get(1, 0)
        alpha2 = alphas.get(2, 0)
        H = sum(2**k * v for k, v in alphas.items())

        key = (alpha1, alpha2)
        if key not in a1_a2_to_H:
            a1_a2_to_H[key] = {'count': 0, 'H': H}
        a1_a2_to_H[key]['count'] += 1

    print(f"\n(alpha_1, alpha_2) -> H mapping:")
    print(f"{'alpha_1':>8} {'alpha_2':>8} {'H':>6} {'count':>8}")
    for (a1, a2) in sorted(a1_a2_to_H.keys()):
        d = a1_a2_to_H[(a1, a2)]
        needed_a2 = (20 - 2*a1) / 4
        marker = " *** H=21 would need this" if needed_a2 == a2 else ""
        print(f"{a1:8d} {a2:8d} {d['H']:6d} {d['count']:8d}{marker}")

    # Check which (alpha_1, alpha_2) pairs would give H=21
    print(f"\n\nH=21 requires 2*a1 + 4*a2 = 20:")
    for a1 in range(0, 15):
        if (20 - 2*a1) % 4 == 0:
            a2_needed = (20 - 2*a1) // 4
            if a2_needed >= 0:
                observed = a1_a2_to_H.get((a1, a2_needed), None)
                if observed:
                    print(f"  a1={a1}, a2={a2_needed}: FOUND with H={observed['H']} ({observed['count']} times)")
                else:
                    # Check what a2 values exist for this a1
                    a2_vals = {a2 for (aa1, a2) in a1_a2_to_H if aa1 == a1}
                    print(f"  a1={a1}, a2={a2_needed}: NOT FOUND (observed a2 values: {sorted(a2_vals)})")

if __name__ == "__main__":
    main()
