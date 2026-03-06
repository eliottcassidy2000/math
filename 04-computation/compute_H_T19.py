#!/usr/bin/env python3
"""
Compute H(T_19): number of Hamiltonian paths in the Paley tournament on 19 vertices.

The Paley tournament T_p (p = 3 mod 4 prime) has vertex set {0,...,p-1}
with edge i->j iff (j-i) is a quadratic residue mod p.

For p=19, QR = {1, 4, 5, 6, 7, 9, 11, 16, 17}.

Uses bitmask DP: dp[mask][last] = # of Hamiltonian paths ending at `last`
using exactly the vertices in `mask`.

Since pure Python is too slow for 2^19 * 19 ~ 10M states, we:
1. Generate and compile a C program for the heavy DP
2. Parse its output
3. Also compute odd cycle counts c_3, c_5
"""

import subprocess
import os
import sys
import tempfile
import time
from math import comb

P = 19

def quadratic_residues(p):
    """Return set of quadratic residues mod p (nonzero squares)."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    return qr

def build_adjacency(p):
    """Build adjacency matrix for Paley tournament T_p.
    adj[i][j] = 1 iff i->j (i.e., (j-i) mod p is a QR)."""
    qr = quadratic_residues(p)
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                adj[i][j] = 1
    return adj

def count_directed_3cycles(adj, n):
    """Count directed 3-cycles (i->j->k->i)."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 2 orientations of 3-cycle on {i,j,k}
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_directed_5cycles(adj, n):
    """Count directed 5-cycles. A directed 5-cycle visits 5 distinct vertices
    in order v0->v1->v2->v3->v4->v0. We count ordered cycles and divide by 5
    (for rotational equivalence) but not by 2 (since the tournament fixes orientation)."""
    count = 0
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or not adj[v0][v1]:
                continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1 or not adj[v1][v2]:
                    continue
                for v3 in range(n):
                    if v3 == v0 or v3 == v1 or v3 == v2 or not adj[v2][v3]:
                        continue
                    for v4 in range(n):
                        if v4 == v0 or v4 == v1 or v4 == v2 or v4 == v3:
                            continue
                        if adj[v3][v4] and adj[v4][v0]:
                            count += 1
    # Each 5-cycle is counted 5 times (once per rotation)
    return count // 5

def generate_c_code(adj, n):
    """Generate C code for bitmask DP to count Hamiltonian paths."""
    # Build adjacency bitmask: for each vertex v, adj_mask[v] = bitmask of vertices v points to
    adj_out = []
    for i in range(n):
        mask = 0
        for j in range(n):
            if adj[i][j]:
                mask |= (1 << j)
        adj_out.append(mask)

    code = f"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N {n}
#define FULL ((1 << N) - 1)

// dp[mask][last] = number of Hamiltonian paths using vertices in mask, ending at last
// We use 64-bit integers since the count can be large
// Memory: 2^19 * 19 * 8 bytes = ~80 MB -- feasible

static long long dp[1 << N][N];

// adj_out[v] = bitmask of vertices that v has an edge TO
static const int adj_out[N] = {{
"""
    for i in range(n):
        code += f"    0x{adj_out[i]:05x}"
        if i < n - 1:
            code += ","
        code += f"  /* vertex {i} */\n"
    code += "};\n\n"

    # Also build adj_in for the DP (predecessors)
    adj_in = []
    for j in range(n):
        mask = 0
        for i in range(n):
            if adj[i][j]:
                mask |= (1 << i)
        adj_in.append(mask)

    code += "// adj_in[v] = bitmask of vertices that have an edge TO v\n"
    code += "static const int adj_in[N] = {\n"
    for j in range(n):
        code += f"    0x{adj_in[j]:05x}"
        if j < n - 1:
            code += ","
        code += f"  /* vertex {j} */\n"
    code += "};\n\n"

    code += """
int main() {
    memset(dp, 0, sizeof(dp));

    // Base case: paths of length 1 (single vertex)
    for (int v = 0; v < N; v++) {
        dp[1 << v][v] = 1;
    }

    // Fill DP
    for (int mask = 1; mask <= FULL; mask++) {
        int popcount = __builtin_popcount(mask);
        if (popcount < 2) continue;  // need at least 2 vertices for a transition

        for (int last = 0; last < N; last++) {
            if (!(mask & (1 << last))) continue;  // last must be in mask

            // Predecessors of last that are in mask
            int prev_mask = mask ^ (1 << last);  // mask without last
            int candidates = adj_in[last] & prev_mask;  // vertices in prev_mask that point to last

            while (candidates) {
                int prev = __builtin_ctz(candidates);  // lowest set bit
                dp[mask][last] += dp[prev_mask][prev];
                candidates &= candidates - 1;  // clear lowest bit
            }
        }
    }

    // Sum over all ending vertices for the full mask
    long long total = 0;
    for (int v = 0; v < N; v++) {
        total += dp[FULL][v];
    }

    printf("H(T_%d) = %lld\\n", N, total);

    // Also print per-endpoint counts for verification
    printf("\\nPer-endpoint counts:\\n");
    for (int v = 0; v < N; v++) {
        printf("  endpoint %d: %lld\\n", v, dp[FULL][v]);
    }

    return 0;
}
"""
    return code

def main():
    print(f"=== Computing H(T_{P}) ===")
    print(f"Paley tournament on {P} vertices")
    print()

    # Build adjacency
    qr = quadratic_residues(P)
    print(f"Quadratic residues mod {P}: {sorted(qr)}")
    adj = build_adjacency(P)

    # Print adjacency for verification
    print(f"\nAdjacency (i->j means entry is 1):")
    print("   ", " ".join(f"{j:2d}" for j in range(P)))
    for i in range(P):
        row = " ".join(f"{adj[i][j]:2d}" for j in range(P))
        print(f"{i:2d}: {row}")

    # Verify: each vertex should have out-degree (p-1)/2 = 9
    for i in range(P):
        out_deg = sum(adj[i])
        assert out_deg == (P - 1) // 2, f"Vertex {i} has out-degree {out_deg}, expected {(P-1)//2}"
    print(f"\nAll vertices have out-degree {(P-1)//2} (verified)")

    # Compute cycle counts
    print(f"\n=== Odd cycle counts ===")

    c3 = count_directed_3cycles(adj, P)
    print(f"c_3(T_{P}) = {c3}")
    # For Paley tournament: c_3 = p(p-1)(p-5)/24 when p >= 7
    expected_c3 = P * (P - 1) * (P - 5) // 24
    print(f"  Expected (formula): {expected_c3}")

    print(f"\nComputing c_5 (this may take a moment for n={P})...")
    t0 = time.time()
    c5 = count_directed_5cycles(adj, P)
    t1 = time.time()
    print(f"c_5(T_{P}) = {c5}  (computed in {t1-t0:.1f}s)")

    # Automorphism group
    aut_size = P * (P - 1) // 2
    print(f"\n=== Automorphism group ===")
    print(f"|Aut(T_{P})| = {P} * {(P-1)//2} = {aut_size}")

    # Generate and compile C code for Hamiltonian path count
    print(f"\n=== Hamiltonian path computation ===")
    print(f"State space: 2^{P} * {P} = {(1 << P) * P:,} entries")
    print(f"Memory: ~{(1 << P) * P * 8 / 1024 / 1024:.0f} MB")

    c_code = generate_c_code(adj, P)
    c_file = "/home/e/Documents/claude/math/04-computation/compute_H_T19.c"
    exe_file = "/home/e/Documents/claude/math/04-computation/compute_H_T19"

    with open(c_file, 'w') as f:
        f.write(c_code)
    print(f"C code written to {c_file}")

    # Compile
    print("Compiling...")
    result = subprocess.run(
        ["gcc", "-O2", "-o", exe_file, c_file],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        sys.exit(1)
    print("Compilation successful")

    # Run
    print(f"\nRunning DP (this may take a few minutes)...")
    t0 = time.time()
    result = subprocess.run(
        [exe_file],
        capture_output=True, text=True,
        timeout=600  # 10 minute timeout
    )
    t1 = time.time()

    if result.returncode != 0:
        print(f"Execution failed: {result.stderr}")
        sys.exit(1)

    print(result.stdout)
    print(f"Computation time: {t1-t0:.1f}s")

    # Parse H value
    for line in result.stdout.split('\n'):
        if line.startswith('H('):
            h_val = int(line.split('=')[1].strip())
            print(f"\n=== Summary ===")
            print(f"H(T_{P}) = {h_val}")
            print(f"|Aut(T_{P})| = {aut_size}")
            print(f"H(T_{P}) / |Aut(T_{P})| = {h_val} / {aut_size} = {h_val / aut_size}")
            if h_val % aut_size == 0:
                print(f"  = {h_val // aut_size} (exact integer)")
            else:
                print(f"  (not an integer)")
            print(f"c_3 = {c3}")
            print(f"c_5 = {c5}")
            # Check parity
            print(f"\nH(T_{P}) mod 2 = {h_val % 2}")
            print(f"H(T_{P}) is {'odd' if h_val % 2 == 1 else 'even'}")
            break

if __name__ == "__main__":
    main()
