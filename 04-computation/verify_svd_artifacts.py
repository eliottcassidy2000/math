"""
verify_svd_artifacts.py — Check if S45 beta_3+beta_4 coexistence was an SVD artifact

S45 found 3 tournaments at n=8 with betti=[1,0,0,1,1,0,0,0] using SVD.
This script regenerates those exact tournaments and computes Betti numbers
using BOTH SVD and mod-p to compare.

If mod-p shows beta_3*beta_4=0 where SVD showed coexistence,
then the coexistence was a NUMERICAL ARTIFACT.

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp, full_chain_complex_svd
)


def main():
    print("=" * 70)
    print("SVD vs MOD-P VERIFICATION OF S45 COEXISTENCE")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(11111)  # Same seed as S45

    # The S45 output showed coexistence at trials 498, 851, 994
    target_trials = {498, 851, 994}

    print(f"\n  Regenerating n=8 tournaments with seed=11111...")
    print(f"  Checking trials: {sorted(target_trials)}")

    for trial in range(max(target_trials) + 1):
        A = random_tournament(n, rng)

        if trial in target_trials:
            print(f"\n  --- Trial {trial} ---")
            score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  Score: {score}")

            # Mod-p (exact)
            res_modp = full_chain_complex_modp(A, n, max_p=7)
            bettis_modp = tuple(res_modp['bettis'].get(p, 0) for p in range(8))
            print(f"  MOD-P bettis:  {bettis_modp}")

            # SVD (potentially inaccurate)
            try:
                res_svd = full_chain_complex_svd(A, n, max_p=7)
                bettis_svd = tuple(res_svd['bettis'].get(p, 0) for p in range(8))
                print(f"  SVD bettis:    {bettis_svd}")
            except Exception as e:
                print(f"  SVD failed: {e}")

            # Compare
            if bettis_modp != bettis_svd:
                print(f"  *** MISMATCH: SVD artifact detected! ***")
                b3_modp = bettis_modp[3]
                b4_modp = bettis_modp[4]
                b3_svd = bettis_svd[3] if 'bettis_svd' in dir() else '?'
                b4_svd = bettis_svd[4] if 'bettis_svd' in dir() else '?'
                if b3_modp * b4_modp == 0 and bettis_svd[3] * bettis_svd[4] > 0:
                    print(f"  SVD showed coexistence (b3={bettis_svd[3]}, b4={bettis_svd[4]})")
                    print(f"  Mod-p shows NO coexistence (b3={b3_modp}, b4={b4_modp})")
                    print(f"  => COEXISTENCE WAS AN SVD ARTIFACT")
            else:
                b3 = bettis_modp[3]
                b4 = bettis_modp[4]
                if b3 > 0 and b4 > 0:
                    print(f"  BOTH METHODS AGREE: coexistence is REAL (b3={b3}, b4={b4})")
                else:
                    print(f"  Both methods agree: no coexistence")

    # Also test 100 more random tournaments
    print(f"\n\n  --- Additional 200 random n=8 comparisons ---")
    rng2 = np.random.RandomState(777)
    mismatches = 0
    total_checked = 0
    for trial in range(200):
        A = random_tournament(n, rng2)
        total_checked += 1

        res_modp = full_chain_complex_modp(A, n, max_p=7)
        bettis_modp = tuple(res_modp['bettis'].get(p, 0) for p in range(8))

        try:
            res_svd = full_chain_complex_svd(A, n, max_p=7)
            bettis_svd = tuple(res_svd['bettis'].get(p, 0) for p in range(8))
        except:
            continue

        if bettis_modp != bettis_svd:
            mismatches += 1
            if mismatches <= 5:
                print(f"  MISMATCH trial {trial}: modp={bettis_modp} vs svd={bettis_svd}")

    print(f"\n  Mismatches: {mismatches}/{total_checked}")
    if mismatches > 0:
        print(f"  SVD has {mismatches/total_checked*100:.1f}% error rate at n=8!")
    else:
        print(f"  SVD agrees perfectly with mod-p at n=8")

    print("\nDONE.")


if __name__ == '__main__':
    main()
