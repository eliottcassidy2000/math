"""
h21_moon_bound_cycle_rich.py — Explore Moon's bound on t3 for cycle-rich tournaments.

For cycle-rich T (no source/sink, all in 3-cycle), what is the minimum t3?
Moon's formula: t3 = C(n,3) - sum_v C(s_v, 2)
where s_v = out-degree (score) of v.

Key question: does min t3 for cycle-rich grow fast enough with n to force H > 21?

Since H >= 1 + 2*alpha_1 and alpha_1 >= t3 (each 3-cycle vertex set contributes
at least one vertex to Omega), if t3 >= 11 then H >= 23 > 21.

Author: opus-2026-03-07-S43
"""
import random
from math import comb
from itertools import combinations, permutations

def random_tournament(n):
    """Random tournament as adjacency matrix."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def scores(A, n):
    return [sum(A[i]) for i in range(n)]

def is_cycle_rich(A, n):
    """No source/sink, every vertex in a 3-cycle."""
    sc = scores(A, n)
    for s in sc:
        if s == 0 or s == n-1:
            return False
    # Check every vertex is in a 3-cycle
    for v in range(n):
        found = False
        for a in range(n):
            if a == v:
                continue
            for b in range(n):
                if b == v or b == a:
                    continue
                # Check v->a->b->v
                if A[v][a] and A[a][b] and A[b][v]:
                    found = True
                    break
            if found:
                break
        if not found:
            return False
    return True

def count_3cycles(A, n):
    """Count number of directed 3-cycles (= number of 3-cycle vertex sets)."""
    count = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                # Check both orientations
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    count += 1
    return count

def count_3cycles_moon(sc, n):
    """Moon's formula: t3 = C(n,3) - sum C(s_v, 2)."""
    return comb(n, 3) - sum(comb(s, 2) for s in sc)

def min_t3_score_sequence(n):
    """
    Find score sequences in [1, n-2] that minimize t3.
    Must have sum = C(n,2) = n(n-1)/2.
    Minimize t3 = C(n,3) - sum C(s,2) => maximize sum C(s,2).
    By convexity, extreme scores (1 and n-2) maximize sum C(s,2).
    """
    target = n * (n-1) // 2
    best_t3 = float('inf')
    best_seq = None

    # Try sequences with k vertices at score 1 and k at score n-2
    # (they must be balanced because 1 + (n-2) = n-1 = average*2... no)
    # Average score = (n-1)/2.

    # For small n, enumerate:
    if n <= 8:
        # Generate all score sequences with sum = target, all in [1, n-2]
        from itertools import product as iprod
        # This is too slow for n=8, let's use a smarter approach
        pass

    # Analytical bound: put as many as possible at extremes 1 and n-2
    # k at score 1, k at score n-2: uses score k + k*(n-2) = k*(n-1)
    # Remaining n-2k vertices: need score n(n-1)/2 - k(n-1) = (n-1)(n/2 - k)
    # Average remaining score = (n-1)(n/2-k)/(n-2k)
    # Need this in [1, n-2].

    for k in range(1, n//2 + 1):
        remaining = n - 2*k
        if remaining < 0:
            break
        remaining_sum = target - k * (n-1)
        if remaining == 0:
            if remaining_sum == 0:
                seq = [1]*k + [n-2]*k
                t3 = comb(n,3) - k*comb(1,2) - k*comb(n-2,2)
                if t3 < best_t3:
                    best_t3 = t3
                    best_seq = sorted(seq)
            continue
        avg = remaining_sum / remaining
        if avg < 1 or avg > n-2:
            continue
        # Put remaining all at floor/ceil of avg
        low = int(avg)
        high = low + 1
        if low < 1:
            low = 1
        if high > n-2:
            high = n-2
        # How many at high vs low
        # n_high * high + n_low * low = remaining_sum
        # n_high + n_low = remaining
        if high == low:
            n_high = remaining
            n_low = 0
        else:
            n_high = remaining_sum - remaining * low
            n_low = remaining - n_high
        if n_high < 0 or n_low < 0:
            continue
        seq = [1]*k + [n-2]*k + [low]*n_low + [high]*n_high
        s_c2 = k*comb(1,2) + k*comb(n-2,2) + n_low*comb(low,2) + n_high*comb(high,2)
        t3 = comb(n,3) - s_c2
        if t3 < best_t3:
            best_t3 = t3
            best_seq = sorted(seq)

    return best_t3, best_seq

# Analytical lower bound on t3 for cycle-rich
print("=== Moon's bound: min t3 for score sequences in [1, n-2] ===")
print(f"{'n':>3} {'C(n,3)':>8} {'min_t3':>8} {'score_seq'}")
for n in range(3, 20):
    t3, seq = min_t3_score_sequence(n)
    print(f"{n:3d} {comb(n,3):8d} {t3:8d} {seq}")

print()
print("=== Key question: is min_t3 >= 11 for cycle-rich n >= ? ===")
print("(If t3 >= 11, then alpha_1 >= 11, so H >= 1+22 = 23 > 21)")
print()

# But wait: alpha_1 counts ALL odd cycles, not just 3-cycles.
# And alpha_1 >= t3 only if every 3-cycle vertex set has a directed 3-cycle.
# Actually, alpha_1 is the number of vertices in Omega(T), which is the number
# of directed odd cycles. t3 counts 3-cycle vertex sets, each giving exactly
# one directed 3-cycle (since both orientations count as 1 vertex set).
# Wait no: each 3-cycle vertex set gives exactly ONE directed 3-cycle (the
# cyclic one, not the transitive one). But the OTHER orientation is also a
# directed 3-cycle! No: {a,b,c} gives EITHER a->b->c->a OR a->c->b->a,
# not both (since the tournament has only one arc per pair). But we count
# both orientations as the same vertex set, so t3 = number of vertex sets.
# And each vertex set contributes exactly ONE vertex to Omega(T).
# So alpha_1 >= t3 (with equality when there are no 5-cycles, 7-cycles, etc.)

# Actually: if t3 >= 11, then Omega has >= 11 vertices from 3-cycles alone,
# so alpha_1 >= 11, and H >= 1 + 2*11 = 23. This is correct.

# But the analytical bound above is for ALL score sequences in [1, n-2],
# not specifically for cycle-rich (which also requires every vertex in 3-cycle).
# The "every vertex in 3-cycle" constraint may force MORE 3-cycles.

# Let's check computationally at n=8:
print("=== Computational check: min t3 for cycle-rich at small n ===")
for n in range(5, 9):
    min_t3_found = float('inf')
    count_cr = 0
    trials = 2000000 if n <= 7 else 500000
    for _ in range(trials):
        A = random_tournament(n)
        if not is_cycle_rich(A, n):
            continue
        count_cr += 1
        t3 = count_3cycles(A, n)
        if t3 < min_t3_found:
            min_t3_found = t3
            print(f"  n={n}: new min t3={t3}, scores={sorted(scores(A,n))}")
    print(f"n={n}: min_t3={min_t3_found} (found {count_cr} cycle-rich in {trials} trials)")

# Now check if t3 >= 11 at n >= 9
print("\n=== Sampling min t3 for cycle-rich at n=9,10,11 ===")
for n in [9, 10, 11, 12]:
    min_t3_found = float('inf')
    count_cr = 0
    trials = 500000
    for _ in range(trials):
        A = random_tournament(n)
        if not is_cycle_rich(A, n):
            continue
        count_cr += 1
        sc = scores(A, n)
        t3 = count_3cycles_moon(sc, n)
        if t3 < min_t3_found:
            min_t3_found = t3
            print(f"  n={n}: new min t3={t3}, scores={sorted(sc)}")
    print(f"n={n}: min_t3={min_t3_found} ({count_cr} cycle-rich in {trials} trials)")
    # Critical: if min_t3 >= 11 for n >= 9, then H >= 23 > 21 trivially!
