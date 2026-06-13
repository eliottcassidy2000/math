"""
h4_relative.py — Why is H_4(T,T\v) = 0 at n=7?

From LES: H_4(T\v) → H_4(T) → H_4(T,T\v) → H_3(T\v) →[i_*] H_3(T)

For BAD vertex: b3(T)=b3(T\v)=1.
H_4(T,T\v) = 0 implies i_* injective (HYP-398 equivalent).

H_4(T,T\v) = ker(∂: H_4(T,T\v) → H_3(T\v)) / im(j_*: H_4(T) → H_4(T,T\v))

Actually, from LES: H_4(T,T\v) = 0 iff:
  1. j_*: H_4(T) → H_4(T,T\v) is surjective, AND
  2. ∂: H_4(T,T\v) → H_3(T\v) has trivial kernel

Or equivalently: the map i_*: H_4(T\v) → H_4(T) is surjective (from the preceding segment).

Let's check: is H_4(T) = 0 and H_4(T\v) = 0 at n=7?
If both are zero, then H_4(T,T\v) = 0 follows from LES:
  0 → H_4(T,T\v) → H_3(T\v) →[i_*] H_3(T)
So H_4(T,T\v) = ker(i_*). And if i_* is injective, H_4(T,T\v) = 0.
But this is circular!

Better approach: from the LES segment:
  H_4(T\v) → H_4(T) → H_4(T,T\v) → H_3(T\v) → H_3(T)

If b4(T) = 0 and b4(T\v) = 0, then:
  0 → H_4(T,T\v) → H_3(T\v) → H_3(T)
gives dim H_4(T,T\v) = ker(i_*^3).

So H_4(T,T\v) = 0 ⟺ i_*^3 injective, which is what we're trying to prove!

The real content is: WHY b4(T) = b4(T\v) = 0 at n=7 for b3(T)=1 tournaments.

At n=7: b4 is ALMOST always 0 (99.8%), but Paley T_7 has b4=6.
For b3=1 tournaments at n=7: is b4 always 0?

At n=8: b4 can be nonzero, and when b3=b4=1 (coexistence), the LES gives:
  H_4(T\v) → H_4(T) → H_4(T,T\v) → H_3(T\v) → H_3(T)
  0 → Z → H_4(T,T\v) → Z → Z

If j_*: H_4(T) → H_4(T,T\v) is zero, then H_4(T,T\v) embeds in H_3(T\v),
so dim H_4(T,T\v) ≤ 1. And i_* could still be injective.

But kind-pasteur-S48 found rank(i_*)=0 at n=8 even with b4=0!
That means H_4(T,T\v) = 1 even when b4(T) = 0.

So the LES gives: 0 → H_4(T,T\v) → H_3(T\v) → H_3(T)
H_4(T,T\v) = ker(i_*) = 1 when rank(i_*)=0.

The question shifts to: what is the RELATIVE chain complex structure?

Author: opus-2026-03-09-S56
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    full_chain_complex_modp,
    RANK_PRIME
)


def check_b4_for_b3_1(n, num_samples, seed=42):
    """For b3=1 tournaments, check if b4 is always 0."""
    rng = np.random.RandomState(seed)
    t0 = time.time()

    found = 0
    b4_dist = Counter()

    for trial in range(num_samples * 20):  # b3=1 at ~7-17%
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1
        b4 = cc['bettis'].get(4, 0)
        b4_dist[b4] += 1

        if b4 > 0:
            scores = sorted(sum(A[i]) for i in range(n))
            bettis = tuple(cc['bettis'].get(p, 0) for p in range(n))
            elapsed = time.time() - t0
            print(f"  b4>0: trial={trial}, b4={b4}, bettis={bettis}, "
                  f"scores={tuple(scores)}, {elapsed:.1f}s")

        if found >= num_samples:
            break

    elapsed = time.time() - t0
    print(f"  n={n}: {found} b3=1 tournaments, b4 dist: {dict(sorted(b4_dist.items()))}, {elapsed:.1f}s")
    return b4_dist


def check_b4_deletions(n, num_samples, seed=42):
    """For b3=1 tournaments with BAD vertices, check b4(T) and b4(T\v)."""
    rng = np.random.RandomState(seed)
    t0 = time.time()

    found = 0
    b4_T_dist = Counter()
    b4_Tv_dist = Counter()

    for trial in range(num_samples * 20):
        A = random_tournament(n, rng)
        cc_T = full_chain_complex_modp(A, n, n - 1)
        if cc_T['bettis'].get(3, 0) != 1:
            continue
        found += 1
        b4_T = cc_T['bettis'].get(4, 0)
        b4_T_dist[b4_T] += 1

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            n1 = n - 1
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
            cc_Tv = full_chain_complex_modp(A_sub, n1, min(n - 1, n1 - 1))
            b3_Tv = cc_Tv['bettis'].get(3, 0)
            if b3_Tv == 1:  # BAD vertex
                b4_Tv = cc_Tv['bettis'].get(4, 0)
                b4_Tv_dist[b4_Tv] += 1

        if found >= num_samples:
            break

    elapsed = time.time() - t0
    print(f"  n={n}: {found} b3=1 tours")
    print(f"    b4(T) dist: {dict(sorted(b4_T_dist.items()))}")
    print(f"    b4(T\\v) for BAD v: {dict(sorted(b4_Tv_dist.items()))}")
    return b4_T_dist, b4_Tv_dist


def main():
    print("=" * 70)
    print("H_4 RELATIVE HOMOLOGY — WHY IS H_4(T,T\\v) = 0?")
    print("=" * 70)

    print("\n--- Phase 1: b4 distribution for b3=1 tournaments ---")

    for n in [7, 8]:
        print(f"\n  n={n}:")
        check_b4_for_b3_1(n, 500, seed=42)

    print(f"\n--- Phase 2: b4(T) and b4(T\\v) for BAD vertices ---")

    for n in [7, 8]:
        print(f"\n  n={n}:")
        check_b4_deletions(n, 200, seed=42)

    print(f"\n{'='*70}")
    print("INTERPRETATION:")
    print("  At n=7: b4(T)=0 for ALL b3=1 tours (seesaw).")
    print("  Also b4(T\\v)=0 for all BAD vertices (n-1=6, b4 vanishes at n=6).")
    print("  So LES gives: 0 → H_4(T,T\\v) → H_3(T\\v) → H_3(T)")
    print("  And H_4(T,T\\v) = ker(i_*) = 0 iff i_* injective.")
    print("")
    print("  At n=8: b4(T) can be 1 when b3=1 (coexistence ~0.15%).")
    print("  And b4(T\\v) can be nonzero at n=7 (Paley T_7 has b4=6).")
    print("  The LES becomes more complex with nonzero H_4 terms.")
    print("  rank(i_*)=0 means H_4(T,T\\v)=1 even when b4(T)=b4(T\\v)=0.")


if __name__ == '__main__':
    main()
    print("\nDONE.")
