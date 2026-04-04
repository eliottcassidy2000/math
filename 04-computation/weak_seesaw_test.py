"""
weak_seesaw_test.py — opus-2026-04-04-S1

Test the WEAK SEESAW CONJECTURE: β₁ + β₃ ≤ 2 for all tournaments.
This replaces the strong seesaw β₁·β₃=0 which fails at n=8.

Key predictions:
  - β₁=1, β₃=2 is IMPOSSIBLE (would give sum=3)
  - Violations are ALWAYS β₁=β₃=1
  - Violations require SC score sequence

Sample at n=8 (5000) and n=9 (2000).
"""
import numpy as np
from collections import defaultdict
from seesaw_n8_analysis import (make_tournament, scores, count_3cycles,
                                 is_sc_score, compute_betti_1_3)

def test_weak_seesaw(n, num_samples, seed_offset=0):
    print(f"\n{'='*60}")
    print(f"WEAK SEESAW TEST: n={n}, {num_samples} samples")
    print(f"{'='*60}")

    stats = defaultdict(int)
    violations = []
    max_sum = 0

    for trial in range(num_samples):
        if trial % 500 == 0:
            print(f"  trial {trial}/{num_samples}...", flush=True)

        T = make_tournament(n, seed=trial + seed_offset)
        beta1, beta3 = compute_betti_1_3(T, n)
        s = beta1 + beta3

        pair = (beta1, beta3)
        stats[pair] += 1

        if s > max_sum:
            max_sum = s
            print(f"  NEW MAX β₁+β₃={s}: trial={trial}, β₁={beta1}, β₃={beta3}")

        if beta1 > 0 and beta3 > 0:
            sc = scores(T, n)
            c3 = count_3cycles(T, n)
            violations.append({
                'trial': trial,
                'beta1': beta1,
                'beta3': beta3,
                'scores': sc,
                'c3': c3,
                'sc_score': is_sc_score(sc),
            })

        if s > 2:
            print(f"  *** WEAK SEESAW VIOLATION! β₁+β₃={s} > 2 ***")
            sc = scores(T, n)
            print(f"  scores={sc}, c₃={count_3cycles(T, n)}")

    print(f"\n(β₁, β₃) pair distribution:")
    for pair in sorted(stats.keys()):
        pct = 100 * stats[pair] / num_samples
        print(f"  {pair}: {stats[pair]} ({pct:.2f}%)")

    print(f"\nMax β₁+β₃: {max_sum}")
    print(f"Weak seesaw β₁+β₃ ≤ 2: {'HOLDS' if max_sum <= 2 else 'FAILS'}")

    print(f"\nStrong seesaw violations (β₁·β₃>0): {len(violations)}/{num_samples}")
    if violations:
        all_sc = all(v['sc_score'] for v in violations)
        all_11 = all(v['beta1'] == 1 and v['beta3'] == 1 for v in violations)
        print(f"  All (1,1): {all_11}")
        print(f"  All SC score: {all_sc}")
        score_dist = defaultdict(int)
        for v in violations:
            score_dist[v['scores']] += 1
        print(f"  Score sequences:")
        for s, c in sorted(score_dist.items(), key=lambda x: -x[1]):
            print(f"    {s}: {c}")

    return max_sum, stats

if __name__ == "__main__":
    # n=8 with 5000 samples
    max8, stats8 = test_weak_seesaw(8, 5000, seed_offset=10000)

    # n=9 with 2000 samples
    max9, stats9 = test_weak_seesaw(9, 2000, seed_offset=20000)

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"n=8: max β₁+β₃ = {max8}, weak seesaw {'HOLDS' if max8 <= 2 else 'FAILS'}")
    print(f"n=9: max β₁+β₃ = {max9}, weak seesaw {'HOLDS' if max9 <= 2 else 'FAILS'}")

    if max8 <= 2 and max9 <= 2:
        print(f"\nCONJECTURE (Weak Seesaw): β₁(T) + β₃(T) ≤ 2 for all tournaments.")
        print(f"  Equivalent to: β₁ = 1 implies β₃ ≤ 1.")
        print(f"  Verified at n=8 (5000 samples) and n=9 (2000 samples).")
        print(f"  Previously proved: β₁·β₃=0 at n≤7 (exhaustive at n=7).")
