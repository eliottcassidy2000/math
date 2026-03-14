"""
h21_n8_sumgap.py — Investigate whether H=21 occurs at n=8.

H(T) = number of Hamiltonian paths in tournament T.
Uses DP on bitmasks with numpy vectorization for speed.

Generates 500k random tournaments on n=8 and checks:
1. Does H=21 ever occur?
2. Distribution of H values near 21
"""

import numpy as np
from collections import Counter
import time
import ctypes
import os
import tempfile
import subprocess
import sys

N = 8
FULL = (1 << N) - 1
NUM_SAMPLES = 500_000

# Write and compile a C extension for the DP
C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define N 8
#define FULL ((1<<N)-1)

/* Count Hamiltonian paths in a tournament given by adj[i][j] */
int count_hp(const int adj[N][N]) {
    /* dp[mask][v] = # of Hamiltonian paths using vertices in mask, ending at v */
    int dp[1<<N][N];
    memset(dp, 0, sizeof(dp));

    for (int v = 0; v < N; v++)
        dp[1<<v][v] = 1;

    for (int mask = 3; mask <= FULL; mask++) {
        int pc = __builtin_popcount(mask);
        if (pc < 2) continue;
        for (int v = 0; v < N; v++) {
            if (!(mask & (1<<v))) continue;
            int prev = mask ^ (1<<v);
            int total = 0;
            for (int u = 0; u < N; u++) {
                if ((prev & (1<<u)) && adj[u][v])
                    total += dp[prev][u];
            }
            dp[mask][v] = total;
        }
    }

    int h = 0;
    for (int v = 0; v < N; v++)
        h += dp[FULL][v];
    return h;
}

/* Process batch of tournaments.
   arcs: flat array of num_tournaments * N*(N-1)/2 bytes (0 or 1).
   For pair (i,j) with i<j, arcs[t*28 + idx] = 1 means i->j, 0 means j->i.
   results: output array of num_tournaments ints.
*/
void process_batch(const unsigned char *arcs, int *results, int num_tournaments) {
    int adj[N][N];
    int pair_idx;

    for (int t = 0; t < num_tournaments; t++) {
        memset(adj, 0, sizeof(adj));
        pair_idx = 0;
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {
                if (arcs[t * 28 + pair_idx]) {
                    adj[i][j] = 1;
                } else {
                    adj[j][i] = 1;
                }
                pair_idx++;
            }
        }
        results[t] = count_hp(adj);
    }
}
"""

def compile_c_extension():
    """Compile the C code into a shared library."""
    src_path = "/tmp/_hp_count.c"
    lib_path = "/tmp/_hp_count.so"

    with open(src_path, "w") as f:
        f.write(C_CODE)

    result = subprocess.run(
        ["gcc", "-O3", "-march=native", "-shared", "-fPIC", "-o", lib_path, src_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("C compilation failed:", result.stderr)
        return None

    return ctypes.CDLL(lib_path)


def main():
    print(f"Investigating H=21 at n={N}")
    print(f"Generating {NUM_SAMPLES} random tournaments")
    print()

    # Compile C extension
    print("Compiling C extension...", flush=True)
    lib = compile_c_extension()
    if lib is None:
        print("FATAL: Could not compile C extension")
        sys.exit(1)
    print("C extension compiled successfully.")

    lib.process_batch.argtypes = [
        ctypes.POINTER(ctypes.c_ubyte),  # arcs
        ctypes.POINTER(ctypes.c_int),    # results
        ctypes.c_int                     # num_tournaments
    ]
    lib.process_batch.restype = None

    rng = np.random.default_rng(42)
    num_pairs = N * (N - 1) // 2  # 28

    h_counter = Counter()
    t0 = time.time()
    batch_size = 50000

    for batch_start in range(0, NUM_SAMPLES, batch_size):
        batch_end = min(batch_start + batch_size, NUM_SAMPLES)
        actual_batch = batch_end - batch_start

        # Generate random arcs: for each pair (i<j), 1 means i->j
        arcs = rng.integers(0, 2, size=(actual_batch, num_pairs), dtype=np.uint8)
        arcs_flat = np.ascontiguousarray(arcs.reshape(-1))

        results = np.zeros(actual_batch, dtype=np.int32)

        lib.process_batch(
            arcs_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
            results.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.c_int(actual_batch)
        )

        for h_val in results:
            h_counter[h_val] += 1

        elapsed = time.time() - t0
        done = batch_end
        rate = done / elapsed if elapsed > 0 else 0
        eta = (NUM_SAMPLES - done) / rate if rate > 0 else 0
        print(f"  {done}/{NUM_SAMPLES} done  ({rate:.0f}/s, ETA {eta:.0f}s)", flush=True)

    total_time = time.time() - t0
    print(f"\nCompleted in {total_time:.1f}s ({NUM_SAMPLES/total_time:.0f} tournaments/s)")

    # Results
    print("\n" + "="*60)
    print("FULL H DISTRIBUTION (sorted by H value)")
    print("="*60)
    for h_val in sorted(h_counter.keys()):
        count = h_counter[h_val]
        pct = 100.0 * count / NUM_SAMPLES
        print(f"  H={h_val:6d}: {count:8d} ({pct:6.3f}%)")

    print("\n" + "="*60)
    print("H VALUES NEAR 21")
    print("="*60)
    for h_val in range(10, 36):
        count = h_counter.get(h_val, 0)
        marker = " <--- TARGET" if h_val == 21 else ""
        print(f"  H={h_val:4d}: {count:8d}{marker}")

    print("\n" + "="*60)
    print("KEY QUESTION: Does H=21 occur?")
    print("="*60)
    h21_count = h_counter.get(21, 0)
    if h21_count > 0:
        print(f"  YES! H=21 occurred {h21_count} times out of {NUM_SAMPLES}")
    else:
        print(f"  NO. H=21 never occurred in {NUM_SAMPLES} samples.")

    # Check parity
    print("\n" + "="*60)
    print("H PARITY ANALYSIS (n=8 is even)")
    print("="*60)
    odd_vals = sorted([h for h in h_counter if h % 2 == 1])
    even_vals = sorted([h for h in h_counter if h % 2 == 0])
    odd_count = sum(h_counter[h] for h in odd_vals)
    even_count = sum(h_counter[h] for h in even_vals)
    print(f"  Odd H values:  {len(odd_vals)} distinct, {odd_count} total ({100*odd_count/NUM_SAMPLES:.1f}%)")
    print(f"  Even H values: {len(even_vals)} distinct, {even_count} total ({100*even_count/NUM_SAMPLES:.1f}%)")
    if odd_vals:
        print(f"  Odd H list: {odd_vals[:40]}{'...' if len(odd_vals)>40 else ''}")
    if even_vals:
        print(f"  Even H list: {even_vals[:40]}{'...' if len(even_vals)>40 else ''}")

    # Summary statistics
    all_h = []
    for h_val, count in h_counter.items():
        all_h.extend([h_val] * count)
    all_h = np.array(all_h)
    print(f"\n  min(H)    = {all_h.min()}")
    print(f"  max(H)    = {all_h.max()}")
    print(f"  mean(H)   = {all_h.mean():.2f}")
    print(f"  median(H) = {np.median(all_h):.1f}")
    print(f"  std(H)    = {all_h.std():.2f}")

    # Missing values in range
    h_min, h_max = int(all_h.min()), int(all_h.max())
    missing = [h for h in range(h_min, h_max+1) if h not in h_counter]
    if missing:
        print(f"\n  Missing H values in [{h_min}, {h_max}]:")
        # Show in chunks
        for i in range(0, len(missing), 20):
            chunk = missing[i:i+20]
            print(f"    {chunk}")

    # H mod 4 analysis
    print("\n" + "="*60)
    print("H MOD 4 ANALYSIS")
    print("="*60)
    for r in range(4):
        vals = sorted([h for h in h_counter if h % 4 == r])
        total = sum(h_counter[h] for h in vals)
        print(f"  H ≡ {r} (mod 4): {len(vals)} distinct, {total} total ({100*total/NUM_SAMPLES:.1f}%)")


if __name__ == "__main__":
    main()
