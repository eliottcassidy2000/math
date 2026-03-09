"""
consecutive_seesaw_n8.py — Test HYP-394 (consecutive seesaw) at n=8 with mod-p

HYP-394 (opus-S54): beta_k * beta_{k+1} = 0 for ALL k>=1, ALL tournaments.
Verified at n<=7 (exhaustive n=6, sampled n=7).

S45 found beta_3=1, beta_4=1 at n=8 using SVD (3/2000 samples).
This script CONFIRMS or REFUTES using exact mod-p arithmetic.

Also tests consecutive seesaw more broadly: all pairs (beta_k, beta_{k+1}).

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp
)


def main():
    print("=" * 70)
    print("CONSECUTIVE SEESAW TEST AT n=8 (mod-p exact)")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(42)

    # Track consecutive pairs
    violations = []
    profiles_found = {}
    total = 0
    b3_count = 0
    b4_count = 0
    b3b4_count = 0

    t0 = time.time()
    N_TRIALS = 5000

    for trial in range(N_TRIALS):
        A = random_tournament(n, rng)
        total += 1

        res = full_chain_complex_modp(A, n, max_p=7)
        bettis = res['bettis']
        profile = tuple(bettis.get(p, 0) for p in range(8))

        b3 = bettis.get(3, 0)
        b4 = bettis.get(4, 0)

        if b3 > 0:
            b3_count += 1
        if b4 > 0:
            b4_count += 1
        if b3 > 0 and b4 > 0:
            b3b4_count += 1

        # Check ALL consecutive pairs
        for k in range(1, 7):
            bk = bettis.get(k, 0)
            bk1 = bettis.get(k+1, 0)
            if bk > 0 and bk1 > 0:
                score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                violations.append({
                    'trial': trial,
                    'k': k,
                    'bk': bk,
                    'bk1': bk1,
                    'profile': profile,
                    'score': score
                })
                key = (k, k+1)
                if key not in profiles_found:
                    profiles_found[key] = []
                profiles_found[key].append(profile)
                print(f"  VIOLATION at trial={trial}: k={k}, "
                      f"beta_{k}={bk}, beta_{k+1}={bk1}, "
                      f"profile={profile}, score={score}")

        # Track exotic profiles
        nonzero = [p for p in range(1, 8) if bettis.get(p, 0) > 0]
        if len(nonzero) >= 2:
            key = tuple(nonzero)
            if key not in profiles_found:
                profiles_found[key] = []
            if len(profiles_found[key]) < 3:
                profiles_found[key].append(profile)

        if (trial+1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (trial+1) / elapsed
            print(f"  {trial+1}/{N_TRIALS}, {elapsed:.1f}s ({rate:.0f}/s), "
                  f"b3={b3_count}, b4={b4_count}, b3&b4={b3b4_count}, "
                  f"violations={len(violations)}", flush=True)

    t1 = time.time()
    print(f"\n  Total: {total} tournaments in {t1-t0:.1f}s")
    print(f"  beta_3 > 0: {b3_count} ({b3_count*100/total:.2f}%)")
    print(f"  beta_4 > 0: {b4_count} ({b4_count*100/total:.2f}%)")
    print(f"  beta_3 AND beta_4 > 0: {b3b4_count} ({b3b4_count*100/total:.2f}%)")

    if violations:
        print(f"\n  CONSECUTIVE SEESAW VIOLATIONS: {len(violations)}")
        for v in violations[:20]:
            print(f"    trial={v['trial']}: beta_{v['k']}={v['bk']}, "
                  f"beta_{v['k']+1}={v['bk1']}, profile={v['profile']}")

        # Categorize by which pair is violated
        pairs = {}
        for v in violations:
            key = (v['k'], v['k']+1)
            pairs[key] = pairs.get(key, 0) + 1
        print(f"\n  Violation pairs:")
        for key, count in sorted(pairs.items()):
            print(f"    beta_{key[0]}*beta_{key[1]} > 0: {count} cases")
    else:
        print(f"\n  CONSECUTIVE SEESAW HOLDS at n=8 ({total} samples)")

    # Also check n=7 quickly for validation
    print(f"\n--- Quick n=7 check (1000 samples) ---")
    rng7 = np.random.RandomState(123)
    viol7 = 0
    for trial in range(1000):
        A = random_tournament(7, rng7)
        res = full_chain_complex_modp(A, 7, max_p=6)
        bettis = res['bettis']
        for k in range(1, 6):
            bk = bettis.get(k, 0)
            bk1 = bettis.get(k+1, 0)
            if bk > 0 and bk1 > 0:
                viol7 += 1
    print(f"  n=7: {viol7} violations in 1000 samples")

    print("\nDONE.")


if __name__ == '__main__':
    main()
