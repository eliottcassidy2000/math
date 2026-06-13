"""
paley_betti_direct.py — Compute Paley tournament Betti numbers via direct method

Use the standard full_chain_complex_modp for T_p. For exact results,
we need max_deg = p-1, but for large p this is infeasible.

Key insight from MISTAKE-020: beta at max_deg is an UPPER BOUND.
For degrees 0 through max_deg-1, results are exact.

Strategy: compute at increasing max_deg to get as many exact values as possible.

Known results:
  T_3: beta = (1, 1, 0)
  T_7: beta = (1, 0, 0, 0, 6, 0, 0)
  T_11: beta = (1, 0, 0, 0, 0, ?, ..., ?)

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import full_chain_complex_modp


def paley_adj(p):
    """Construct Paley tournament for prime p = 3 mod 4."""
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


def main():
    print("PALEY TOURNAMENT BETTI NUMBERS (DIRECT METHOD)")
    print("=" * 60)

    for p in [3, 7, 11]:
        print(f"\n--- T_{p} ---")
        A = paley_adj(p)

        # For each max_deg, compute and show which bettis are exact
        for md in range(3, min(p, 10)):
            t0 = time.time()
            cc = full_chain_complex_modp(A, p, max_p=md)
            elapsed = time.time() - t0

            bettis = cc['bettis']
            omega_dims = cc['omega_dims']

            # bettis 0..md-1 are EXACT, bettis[md] is upper bound
            exact_str = ", ".join(f"b{k}={bettis[k]}" for k in range(md))
            ub_str = f"b{md}<={bettis[md]}"

            print(f"  max_deg={md} ({elapsed:.1f}s): {exact_str}, {ub_str}")
            print(f"    Omega dims: {[omega_dims.get(k,0) for k in range(md+1)]}")

            if elapsed > 300:  # 5 minutes timeout per increment
                print(f"    STOPPING: too slow for max_deg={md+1}")
                break


if __name__ == '__main__':
    main()
    print("\nDONE.")
