"""
na_face_distribution.py — Why is the NA face distribution 25/50/25 for 3-paths?

For a 3-path (a,b,c,d) in a tournament:
  Face (b,c,d): always in A_2 (b->c, c->d)
  Face (a,c,d): in A_2 iff a->c (else c->a) AND c->d (always TRUE). So: NA iff c->a.
  Face (a,b,d): in A_2 iff a->b (TRUE) AND b->d. So: NA iff d->b.
  Face (a,b,c): always in A_2 (a->b, b->c)

So a 3-path (a,b,c,d) has:
  0 NA faces iff a->c AND b->d
  1 NA face iff exactly one of {c->a, d->b}
  2 NA faces iff c->a AND d->b

Each condition (a->c vs c->a, b->d vs d->b) is determined by the tournament.
Are these independent? If so, the distribution would be exactly 25/50/25.

But (a->c) and (b->d) involve DIFFERENT vertex pairs, so they ARE independent
in a RANDOM tournament. But for a FIXED tournament, the distributions over
all 3-paths might not be independent...

Let's check: is the 25/50/25 distribution EXACT or approximate?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths


def main():
    print("=" * 70)
    print("NA FACE DISTRIBUTION FOR 3-PATHS")
    print("=" * 70)

    # Part 1: Check if 25/50/25 is exact for every tournament
    print("\n--- Part 1: Per-tournament distribution of NA faces per 3-path ---")

    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        exact_25_50_25 = 0
        total = 0
        devs = []

        for trial in range(N):
            A = bits_to_adj(trial, n)
            a2 = enumerate_allowed_paths(A, n, 2)
            a3 = enumerate_allowed_paths(A, n, 3)
            a2_set = set(a2)

            if len(a3) == 0:
                continue

            counts = Counter()
            for path in a3:
                a, b, c, d = path
                na = 0
                # Face (a,c,d): NA iff c->a (since a->b->c->d, need a->c for (a,c,d) in A_2)
                if not A[a][c]:  # c->a
                    na += 1
                # Face (a,b,d): NA iff d->b (since a->b, need b->d for (a,b,d) in A_2)
                if not A[b][d]:  # d->b
                    na += 1
                counts[na] += 1

            total += 1
            p0 = counts[0] / len(a3)
            p1 = counts[1] / len(a3)
            p2 = counts[2] / len(a3)

            if abs(p0 - 0.25) < 1e-10 and abs(p1 - 0.5) < 1e-10 and abs(p2 - 0.25) < 1e-10:
                exact_25_50_25 += 1
            else:
                devs.append((trial, counts[0], counts[1], counts[2], len(a3)))

        print(f"  n={n}: exact 25/50/25 in {exact_25_50_25}/{total} tournaments")
        if devs:
            for trial, c0, c1, c2, tot_paths in devs[:5]:
                print(f"    trial {trial}: {c0}/{c1}/{c2} out of {tot_paths} "
                      f"({100*c0/tot_paths:.1f}/{100*c1/tot_paths:.1f}/{100*c2/tot_paths:.1f}%)")

    # Part 2: Algebraic explanation
    print("\n--- Part 2: Algebraic explanation ---")
    print("""
    For a 3-path (a,b,c,d) with a->b->c->d:
      NA iff c->a:  depends on edge (a,c)
      NA iff d->b:  depends on edge (b,d)

    These are two DIFFERENT edges. In a random tournament, they're independent.
    But for a FIXED tournament, the set of 3-paths is correlated with edges.

    However: summing over ALL 3-paths, each 3-path has two "skip" edges:
    (a,c) = "inner skip" and (b,d) = "outer skip".
    Each skip is independently forward or backward.

    The key question: does the existence of the 3-path (a->b->c->d) correlate
    with the direction of (a,c) or (b,d)?
    """)

    # Part 3: Direct combinatorial count
    # For a tournament T, count 3-paths (a,b,c,d) with a->c and b->d
    # vs a->c and d->b, etc.
    print("--- Part 3: Correlation between skip directions and path existence ---")
    n = 6
    N = 2**(n*(n-1)//2)

    # For each tournament, compute the conditional probability:
    # P(a->c | (a,b,c,d) is a 3-path) vs P(a->c | a,c distinct vertices)
    for trial in range(min(N, 5)):
        A = bits_to_adj(trial, n)
        a3 = enumerate_allowed_paths(A, n, 3)

        fwd_ac = sum(1 for a,b,c,d in a3 if A[a][c])
        fwd_bd = sum(1 for a,b,c,d in a3 if A[b][d])

        # Baseline: P(a->c) over all ordered pairs (a,c) with a != c
        total_fwd = sum(A[i][j] for i in range(n) for j in range(n) if i != j)
        total_pairs = n * (n-1)

        print(f"  trial {trial}: P(a->c|3-path)={fwd_ac/len(a3):.3f}, "
              f"P(b->d|3-path)={fwd_bd/len(a3):.3f}, "
              f"P(fwd|random)={total_fwd/total_pairs:.3f}, "
              f"|A_3|={len(a3)}")

    # Part 4: The 25/50/25 result OVER ALL TOURNAMENTS summed
    # (from the earlier computation) is clearly exact. But is it exact
    # for EACH tournament?
    print("\n--- Part 4: Checking per-tournament deviation from 25/50/25 ---")
    n = 5
    N = 2**(n*(n-1)//2)
    max_dev = 0
    for trial in range(N):
        A = bits_to_adj(trial, n)
        a3 = enumerate_allowed_paths(A, n, 3)
        if len(a3) == 0: continue
        counts = Counter()
        for path in a3:
            a, b, c, d = path
            na = (0 if A[a][c] else 1) + (0 if A[b][d] else 1)
            counts[na] += 1
        p0 = counts[0] / len(a3)
        dev = abs(p0 - 0.25)
        max_dev = max(max_dev, dev)

    print(f"  n=5: max deviation from 25% for '0 NA' class: {max_dev:.6f}")
    if max_dev < 1e-10:
        print("  EXACT 25/50/25 for EVERY tournament!")
    else:
        print(f"  NOT exact — deviates by up to {max_dev:.4f}")

    # Part 5: Check at larger n with sampling
    print("\n--- Part 5: Per-tournament exactness at n=7,8 ---")
    for n in [7, 8]:
        rng = np.random.RandomState(42 + n)
        max_dev = 0
        for trial in range(200):
            A = random_tournament(n, rng)
            a3 = enumerate_allowed_paths(A, n, 3)
            if len(a3) == 0: continue
            counts = Counter()
            for path in a3:
                a, b, c, d = path
                na = (0 if A[a][c] else 1) + (0 if A[b][d] else 1)
                counts[na] += 1
            p0 = counts[0] / len(a3)
            dev = abs(p0 - 0.25)
            max_dev = max(max_dev, dev)
        print(f"  n={n}: max deviation from 25%: {max_dev:.6f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
