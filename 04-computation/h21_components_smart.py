#!/usr/bin/env python3
"""
Smart component analysis: don't compute I(C,2) directly.
Use H = product of I(C_i, 2) to infer component values.

For single-component: I(C,2) = H.
For multi-component: need to determine which factor goes where.
  -> Actually, for 2 components: H = I(C1,2) * I(C2,2).
     We can compute I(C1,2) by restricting to the cycles in C1.
     But this is still exponential.

SIMPLER APPROACH: just count components and check if H is divisible by 3,7.
If single-component and H != 21: OK.
If multi-component: check if any factorization allows I=7 or I=21.

Actually, the simplest approach: just find which H values come from
multi-component Omega, and what the component sizes are.

For the H=21 question: since H=21 is NEVER achieved at n<=7,
and the only way to get I=21 is H=21 (single comp) or H=21*k (multi-comp with I=21),
we just need to show H=21 never occurs (which we know).

But the PROOF needs to show WHY. The component approach says:
- 21 = 21 (single component impossible)
- 21 = 3*7 (impossible because I=7 is impossible for a component)
So proving I(C,2)=7 impossible would complete the proof!

Let me focus on: can I(C,2)=7 occur as a component value?

At n=5: H values = {1,3,5,9,11,13,15}. Single-component gives I = H.
  Multi-component: need to find tournaments with disconnected Omega.
  If multi-comp and one has I=7, then H >= 7*3 = 21 > 15. Impossible.

At n=6: max H = 45. I=7 in multi-comp gives H >= 21.
  H=21 not achievable, H=35 (=7*5) not achievable at n=6.
  So I=7 not achievable at n=6.

At n=7: max H = 189. H=35 (=7*5) IS achievable. Is it from a (7,5) factorization?
  If H=35 comes from a SINGLE connected component, then I=35, no I=7 needed.
  If from 2 components: 35 = 5*7 or 7*5. Both need I=7.

So the key question: does H=35 at n=7 come from a connected or disconnected Omega?

Let me check by computing Omega connectivity at n=7 for tournaments with H=35.
For efficiency: use only 3-cycles for Omega (much faster) and check connectivity.
If 3-cycle-only Omega is connected, full Omega is also connected (adding 5/7-cycles
only adds more edges).

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict
import time


def held_karp(n, adj_bits):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj_bits[v] & (1 << u):
                    dp[S | (1 << u)][u] += c
    return sum(dp[(1 << n) - 1])


def get_three_cycles(n, adj_bits):
    """Get list of 3-cycle vertex sets."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] & (1 << j)):
                A[i][j] = 1

    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append(frozenset({a, b, c}))
    return cycles


def count_components_3cycle(cycles):
    """Count connected components of 3-cycle conflict graph."""
    if not cycles:
        return 0, []

    nc = len(cycles)
    # Build adjacency
    adj = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)

    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        queue = [start]
        visited[start] = True
        while queue:
            v = queue.pop()
            comp.append(v)
            for u in adj[v]:
                if not visited[u]:
                    visited[u] = True
                    queue.append(u)
        components.append(comp)

    return len(components), components


def analyze_n7_h35():
    """Check if H=35 at n=7 comes from connected or disconnected Omega."""
    n = 7
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m  # 2,097,152

    h35_connected = 0
    h35_disconnected = 0
    h35_examples = []

    # Also track multi-component tournaments generally
    multi_comp_H = defaultdict(int)

    t0 = time.time()

    for bits in range(total):
        if bits % 200000 == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  {bits}/{total} ({elapsed:.0f}s, ETA {eta:.0f}s)")

        adj_bits = [0]*n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj_bits[j] |= (1 << i)
            else:
                adj_bits[i] |= (1 << j)

        # Quick filter: only compute H for potentially interesting cases
        cycles = get_three_cycles(n, adj_bits)

        if not cycles:
            continue

        ncomp, components = count_components_3cycle(cycles)

        if ncomp > 1:
            # Multi-component! Compute H.
            H = held_karp(n, adj_bits)
            multi_comp_H[H] += 1

            if H == 35:
                h35_disconnected += 1
                if len(h35_examples) < 5:
                    comp_sizes = [len(c) for c in components]
                    comp_verts = [[set(cycles[i]) for i in comp] for comp in components]
                    h35_examples.append((bits, comp_sizes, comp_verts))

        else:
            # Single component in 3-cycle graph.
            # Full Omega is also connected (5/7-cycles only add edges).
            # H = I(full Omega, 2), which equals the single component value.
            # Check if H=35
            H = held_karp(n, adj_bits)
            if H == 35:
                h35_connected += 1

    elapsed = time.time() - t0
    print(f"\nn=7 analysis in {elapsed:.0f}s")
    print(f"\nH=35 tournaments:")
    print(f"  Connected Omega: {h35_connected}")
    print(f"  Disconnected Omega (3-cycle graph): {h35_disconnected}")

    if h35_examples:
        print(f"\n  Disconnected H=35 examples:")
        for bits, sizes, verts in h35_examples[:3]:
            print(f"    bits={bits}: component sizes {sizes}")
            for i, cv in enumerate(verts):
                print(f"      comp {i}: {cv}")

    print(f"\nAll multi-component H values: {sorted(multi_comp_H.keys())}")
    print(f"Multi-component H frequency:")
    for h in sorted(multi_comp_H.keys()):
        print(f"  H={h}: {multi_comp_H[h]} tournaments")


if __name__ == '__main__':
    analyze_n7_h35()
