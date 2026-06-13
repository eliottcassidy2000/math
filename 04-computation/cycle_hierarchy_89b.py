#!/usr/bin/env python3
"""
CYCLE HIERARCHY — DO ODD CYCLE COUNTS DETERMINE H?
opus-2026-03-14-S89b

DISCOVERY: H = 1 + 2t3 for n≤4 (exact).
At n=5: t5 discriminates within each t3 class.
QUESTION: Does (t3, t5) determine H at n=5? At n=6?
Does (t3, t5, t_7, ...) always determine H?

This would give H as a function of ODD CYCLE COUNTS only —
a beautiful connection to the OCF (Odd Cycle Collection Formula).
"""

from math import factorial, comb
from collections import Counter, defaultdict
from itertools import permutations

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj.get((v, u), 0) == 1:
                    new_mask = mask | (1 << u)
                    dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

def tournament_from_bits(n, bits):
    adj = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
            idx += 1
    return adj

def count_directed_k_cycles(n, adj, k):
    """Count directed k-cycles in tournament."""
    from itertools import permutations
    count = 0
    vertices = list(range(n))
    for subset in _subsets_of_size(vertices, k):
        for perm in permutations(subset):
            # Check if perm[0]→perm[1]→...→perm[k-1]→perm[0]
            is_cycle = True
            for i in range(k):
                if adj.get((perm[i], perm[(i+1) % k]), 0) != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    # Each k-cycle is counted k times (k rotations)
    return count // k

def _subsets_of_size(vertices, k):
    """Generate all subsets of given size."""
    from itertools import combinations
    return combinations(vertices, k)

def count_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0)) or \
                   (adj.get((i,k),0) and adj.get((k,j),0) and adj.get((j,i),0)):
                    count += 1
    return count

def count_5cycles(n, adj):
    """Optimized 5-cycle counter."""
    from itertools import combinations
    count = 0
    for subset in combinations(range(n), 5):
        for perm in permutations(subset):
            is_cycle = True
            for i in range(5):
                if adj.get((perm[i], perm[(i+1) % 5]), 0) != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
        # Already counted rotations: 5 per cycle
    return count // 5

print("=" * 70)
print("CYCLE HIERARCHY — DO ODD CYCLE COUNTS DETERMINE H?")
print("opus-2026-03-14-S89b")
print("=" * 70)

# ======================================================================
# PART 1: n=5 — DOES (t3, t5) DETERMINE H?
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: n=5 — (t3, t5) → H")
print("=" * 70)

n = 5
m = 10

group = defaultdict(lambda: Counter())
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    t5 = count_5cycles(n, adj)
    group[(t3, t5)][h] += 1

print(f"\n  n={n}: Grouping by (t3, t5) → H:")
all_determined = True
for key in sorted(group.keys()):
    h_dist = group[key]
    determined = len(h_dist) == 1
    if not determined:
        all_determined = False
    mark = "✓" if determined else "✗"
    print(f"    (t3={key[0]}, t5={key[1]}): {dict(sorted(h_dist.items()))} {mark}")

print(f"\n  (t3, t5) determines H at n=5? {all_determined}")

# ======================================================================
# PART 2: n=6 — DOES (t3, t5) DETERMINE H?
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: n=6 — (t3, t5) → H")
print("=" * 70)

n = 6
m = 15

print(f"  Computing t3 and t5 for all {2**m} tournaments at n={n}...")

group6 = defaultdict(lambda: Counter())
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    t5 = count_5cycles(n, adj)
    group6[(t3, t5)][h] += 1

print(f"  Done. {len(group6)} distinct (t3, t5) pairs.\n")

non_det = 0
total = 0
for key in sorted(group6.keys()):
    h_dist = group6[key]
    total += 1
    determined = len(h_dist) == 1
    if not determined:
        non_det += 1
        print(f"    (t3={key[0]:2d}, t5={key[1]:3d}): {dict(sorted(h_dist.items()))}")

print(f"\n  Total groups: {total}")
print(f"  Non-determined: {non_det}")
print(f"  (t3, t5) determines H at n=6? {non_det == 0}")

if non_det > 0:
    print(f"\n  Need more invariants. Let's check (t3, t5, scores):")
    group6b = defaultdict(lambda: Counter())
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        h = compute_H(n, adj)
        t3 = count_3cycles(n, adj)
        t5 = count_5cycles(n, adj)
        scores = tuple(sorted(sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
        group6b[(t3, t5, scores)][h] += 1

    non_det_b = sum(1 for h_dist in group6b.values() if len(h_dist) > 1)
    total_b = len(group6b)
    print(f"  (t3, t5, scores): {total_b} groups, {non_det_b} non-determined")

    if non_det_b > 0:
        print(f"\n  Still non-determined groups (showing first 10):")
        count_shown = 0
        for key in sorted(group6b.keys()):
            h_dist = group6b[key]
            if len(h_dist) > 1:
                if count_shown < 10:
                    print(f"    {key}: {dict(sorted(h_dist.items()))}")
                    count_shown += 1

# ======================================================================
# PART 3: THE OCF CONNECTION
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: CONNECTION TO THE ODD CYCLE COLLECTION FORMULA")
print("=" * 70)

print("""
  The OCF states:
    H(T) ≡ 1 (mod 2)  always
    H(T) = sum over independent odd cycle sets of contributions

  If the OCF gives H as a function of CYCLE STRUCTURE,
  then the odd cycle counts (t3, t5, t7, ...) should suffice.

  But "cycle structure" is more than just COUNTS — it includes
  the INTERSECTION PATTERN of cycles (which cycles share vertices).

  The fact that t3 alone determines H for n≤4 but not n≥5
  shows that the intersection pattern matters starting at n=5
  (where the first non-trivial cycle intersections appear).

  The t5 discriminates at n=5 because the 5-cycle on all 5 vertices
  encodes the GLOBAL cycle structure (how the 3-cycles interact).

  CONJECTURE (Cycle Hierarchy):
  For each n, there exists a set of odd cycle counts
  (t3, t5, ..., t_{2k+1}) with 2k+1 ≤ n such that
  the tuple (t3, t5, ..., t_{2k+1}) determines H(T).

  At n=3: t3 suffices
  At n=4: t3 suffices
  At n=5: (t3, t5) suffices
  At n=6: (t3, t5, ?) — need to check

  But from our data: (t3, t5) does NOT determine H at n=6.
  The missing invariant is NOT just scores — it's something deeper.
  Likely the INTERSECTION GRAPH of 3-cycles (which pairs of
  3-cycles share an edge).
""")

# ======================================================================
# PART 4: EXACT H FORMULAS BY CYCLE TYPE
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: EXACT H FORMULAS")
print("=" * 70)

print("""
  For n ≤ 4:  H = 1 + 2·t3  (EXACT)

  For n = 5:  H = 1 + 2·t3 + 2·ε5(T)
  where ε5(T) depends on the 5-cycle count:
    If t3 ≤ 3: ε5 = 0 (t3 alone determines)
    If t3 = 4: ε5 = t5 - 1 (so H = 1 + 2·4 + 2·(t5-1) = 2t5 + 7)
    Actually: for t3=4, H = 11 when t5=1, H=13 when t5=2, H=15 when t5=3
    H = 2·t5 + 9? No: 2·1+9=11 ✓, 2·2+9=13 ✓, 2·3+9=15 ✓ YES!
    For t3=5: H=15, t5=?
""")

# Check t5 for t3=5 at n=5
n = 5
m = 10
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    t3 = count_3cycles(n, adj)
    if t3 == 5:
        h = compute_H(n, adj)
        t5 = count_5cycles(n, adj)
        print(f"  t3=5: H={h}, t5={t5}")
        break

# Comprehensive: all (t3, t5, H) at n=5
print("\n  Complete (t3, t5, H) table at n=5:")
for key in sorted(group.keys()):
    t3, t5 = key
    h_dist = group[key]
    for h, count in sorted(h_dist.items()):
        # Check: H = 2t5 + f(t3) ?
        if t5 > 0:
            offset = h - 2*t5
        else:
            offset = h
        print(f"    t3={t3}, t5={t5}: H={h} (count={count}), H-2t5={offset}")

# ======================================================================
# PART 5: THE UNIVERSAL FORMULA ATTEMPT
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: UNIVERSAL FORMULA — H = 1 + 2t3 + 2t5 + ... ?")
print("=" * 70)

# Check if H = 1 + 2*(t3 + t5) works at n=5
print(f"\n  Testing H = 1 + 2(t3 + t5) at n=5:")
n = 5
m = 10
correct = 0
total = 0
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    t5 = count_5cycles(n, adj)
    predicted = 1 + 2*(t3 + t5)
    if predicted == h:
        correct += 1
    total += 1

print(f"    Correct: {correct}/{total} ({correct/total:.4f})")

# What about H = 1 + 2*t3 + a*t5 for some other a?
print(f"\n  Testing H = 1 + 2t3 + a·t5 at n=5 for various a:")
for a_num in range(0, 10):
    a = a_num
    correct = 0
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        h = compute_H(n, adj)
        t3 = count_3cycles(n, adj)
        t5 = count_5cycles(n, adj)
        predicted = 1 + 2*t3 + a*t5
        if predicted == h:
            correct += 1
    print(f"    a={a}: {correct}/{total} correct ({correct/total:.4f})")

# Regression: H = c0 + c1*t3 + c2*t5
print(f"\n  Exact regression H = c0 + c1*t3 + c2*t5 at n=5:")
from fractions import Fraction

# Build the system
T3 = []
T5 = []
HH = []
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    t5 = count_5cycles(n, adj)
    T3.append(t3)
    T5.append(t5)
    HH.append(h)

N = len(HH)
# Normal equations:
# N*c0 + sum(t3)*c1 + sum(t5)*c2 = sum(H)
# sum(t3)*c0 + sum(t3^2)*c1 + sum(t3*t5)*c2 = sum(t3*H)
# sum(t5)*c0 + sum(t3*t5)*c1 + sum(t5^2)*c2 = sum(t5*H)

s0 = N
s_t3 = sum(T3)
s_t5 = sum(T5)
s_t3t3 = sum(t**2 for t in T3)
s_t5t5 = sum(t**2 for t in T5)
s_t3t5 = sum(T3[i]*T5[i] for i in range(N))
s_h = sum(HH)
s_t3h = sum(T3[i]*HH[i] for i in range(N))
s_t5h = sum(T5[i]*HH[i] for i in range(N))

# Solve 3x3 system using Cramer's rule with Fractions
A = [[Fraction(s0), Fraction(s_t3), Fraction(s_t5)],
     [Fraction(s_t3), Fraction(s_t3t3), Fraction(s_t3t5)],
     [Fraction(s_t5), Fraction(s_t3t5), Fraction(s_t5t5)]]
b = [Fraction(s_h), Fraction(s_t3h), Fraction(s_t5h)]

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
           -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
           +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

D = det3(A)
if D != 0:
    A0 = [[b[i] if j==0 else A[i][j] for j in range(3)] for i in range(3)]
    A1 = [[b[i] if j==1 else A[i][j] for j in range(3)] for i in range(3)]
    A2 = [[b[i] if j==2 else A[i][j] for j in range(3)] for i in range(3)]

    c0 = det3(A0) / D
    c1 = det3(A1) / D
    c2 = det3(A2) / D

    print(f"    c0 = {c0} = {float(c0):.6f}")
    print(f"    c1 = {c1} = {float(c1):.6f}")
    print(f"    c2 = {c2} = {float(c2):.6f}")
    print(f"    H ≈ {float(c0):.2f} + {float(c1):.2f}·t3 + {float(c2):.2f}·t5")

    # R² for this regression
    ss_res = sum((HH[i] - float(c0 + c1*T3[i] + c2*T5[i]))**2 for i in range(N))
    mean_h = sum(HH)/N
    ss_tot = sum((HH[i] - mean_h)**2 for i in range(N))
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    print(f"    R² = {r2:.6f}")

    # Check if residual is always 0
    max_resid = max(abs(HH[i] - float(c0 + c1*T3[i] + c2*T5[i])) for i in range(N))
    print(f"    Max residual: {max_resid:.6f}")
    if max_resid < 0.001:
        print(f"    *** H = {c0} + {c1}·t3 + {c2}·t5 EXACTLY! ***")

print("\n" + "=" * 70)
print("DONE — CYCLE HIERARCHY")
print("=" * 70)
