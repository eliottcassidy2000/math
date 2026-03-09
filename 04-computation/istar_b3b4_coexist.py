"""
istar_b3b4_coexist.py — Check i_*-injectivity for beta_3=beta_4=1 tournaments

These are RARE (~0.15% at n=8). The consecutive seesaw fails for them.
Key question: Does the LES still give beta_3 ≤ 1 even when beta_4(T) > 0?

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    full_chain_complex_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME

# Import the compute function from istar_n8_investigation
sys.path.insert(0, '.')
from istar_n8_investigation import compute_istar_all_degrees, compute_relative_bettis_from_les


def main():
    print("=" * 70)
    print("i_*-INJECTIVITY FOR beta_3=beta_4=1 TOURNAMENTS (n=8)")
    print("=" * 70)

    n = 8
    max_p = 7
    rng = np.random.RandomState(42)

    # Search for beta_3=beta_4=1 tournaments
    print(f"\n--- Searching for beta_3=beta_4=1 tournaments ---")
    found = 0
    target = 5
    t0 = time.time()
    trials = 0

    while found < target:
        A = random_tournament(n, rng)
        trials += 1
        cc = full_chain_complex_modp(A, n, max_p)
        b3 = cc['bettis'].get(3, 0)
        b4 = cc['bettis'].get(4, 0)

        if b3 > 0 and b4 > 0:
            found += 1
            profile = tuple(cc['bettis'].get(p, 0) for p in range(max_p + 1))
            scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"\n  FOUND #{found} (trial {trials}): bettis={profile}, scores={scores}")

            # Check ALL vertices
            for v in range(n):
                ri, bT, bTv = compute_istar_all_degrees(A, n, v, max_p)
                bTv_profile = tuple(bTv.get(p, 0) for p in range(max_p + 1))
                rel = compute_relative_bettis_from_les(ri, bT, bTv, max_p)
                rel_profile = tuple(rel.get(p, 0) for p in range(max_p + 1))

                ri_vals = tuple(ri.get(p, 0) for p in range(max_p + 1))
                print(f"    v={v}: b(T\\v)={bTv_profile}, rank(i_*)={ri_vals}, H_rel={rel_profile}")

            elapsed = time.time() - t0
            print(f"    ({elapsed:.1f}s total, {trials} trials)")

        if trials % 1000 == 0:
            elapsed = time.time() - t0
            print(f"  {trials} trials, {found}/{target} found, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  Total: {found} beta_3=beta_4>0 found in {trials} trials, {elapsed:.1f}s")
    print(f"  Frequency: {found/trials*100:.2f}%")


if __name__ == '__main__':
    main()
    print("\nDONE.")
