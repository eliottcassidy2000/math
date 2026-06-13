#!/usr/bin/env python3
"""
d5_parity_mod4.py — opus-2026-03-14-S71g

Deep investigation of the d₅ parity constraint.

KEY OBSERVATION from alpha1_parity_analysis.py:
  At n=5: t₃=3 → d₅ always odd; t₃ ∈ {0,1,2,5} → d₅ always even
  At n=6: t₃=3 → d₅ always odd; t₃=7 → d₅ always odd
          t₃ ∈ {0,1,2,5} → d₅ always even

HYPOTHESIS: t₃ ≡ 3 mod 4 forces d₅ odd.

This script:
1. Verifies the mod-4 pattern exhaustively at n=5,6
2. Searches for the underlying algebraic identity
3. Checks if the pattern extends to n=7 (sample)
4. Explores connection to the (z-2)(z-3) framework
"""

from itertools import permutations, combinations
from collections import defaultdict
import sys

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        elif A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def count_directed_5cycles(A, n):
    """Count directed 5-cycles with canonical start = min vertex."""
    d5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            ok = True
            for idx in range(5):
                if A[order[idx]][order[(idx+1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                d5 += 1
    return d5

def count_directed_7cycles(A, n):
    """Count directed 7-cycles with canonical start = min vertex."""
    d7 = 0
    for verts in combinations(range(n), 7):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            ok = True
            for idx in range(7):
                if A[order[idx]][order[(idx+1) % 7]] != 1:
                    ok = False
                    break
            if ok:
                d7 += 1
    return d7

def score_sequence(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

# ============================================================
# Analysis 1: d₅ mod 2 as function of t₃ mod 4
# ============================================================

print("=" * 70)
print("d₅ PARITY vs t₃ mod 4 — EXHAUSTIVE ANALYSIS")
print("=" * 70)

for n in [5, 6]:
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    print(f"\n--- n = {n} ({num_t} tournaments) ---")

    # Collect (t₃, d₅) for every tournament
    t3_mod4_to_d5_parities = defaultdict(lambda: {'even': 0, 'odd': 0, 'vals': set()})
    t3_to_d5_vals = defaultdict(lambda: defaultdict(int))

    for bits in range(num_t):
        A = make_tournament(bits, n)
        t3 = count_3cycles(A, n)
        d5 = count_directed_5cycles(A, n)

        key = t3 % 4
        if d5 % 2 == 0:
            t3_mod4_to_d5_parities[key]['even'] += 1
        else:
            t3_mod4_to_d5_parities[key]['odd'] += 1
        t3_mod4_to_d5_parities[key]['vals'].add(t3)

        t3_to_d5_vals[t3][d5] += 1

    print(f"\n  t₃ mod 4 → d₅ parity:")
    for r in range(4):
        if r in t3_mod4_to_d5_parities:
            g = t3_mod4_to_d5_parities[r]
            total = g['even'] + g['odd']
            det = "ALWAYS EVEN" if g['odd'] == 0 else ("ALWAYS ODD" if g['even'] == 0 else "BOTH")
            t3_vals = sorted(g['vals'])
            print(f"    t₃ ≡ {r} (mod 4): even={g['even']}, odd={g['odd']} [{det}]  (t₃ ∈ {t3_vals})")

    print(f"\n  Detailed (t₃, d₅) values:")
    for t3 in sorted(t3_to_d5_vals.keys()):
        d5_dict = t3_to_d5_vals[t3]
        d5_vals = sorted(d5_dict.keys())
        parities = set(v % 2 for v in d5_vals)
        par_str = "all even" if parities == {0} else ("all odd" if parities == {1} else "mixed")
        print(f"    t₃={t3} (≡{t3%4} mod 4): d₅ ∈ {d5_vals} [{par_str}]")

# ============================================================
# Analysis 2: Look for algebraic identity
# ============================================================

print(f"\n{'='*70}")
print("ALGEBRAIC IDENTITY SEARCH")
print("='*70")

# At n=5: for each 5-vertex set S, define:
#   hc(S) = # directed Hamiltonian cycles of T[S] (with canonical start)
#   t3(S) = # 3-cycles in T[S]
# Then d₅ = Σ_S hc(S) over all C(n,5) sets S.

# Key identity at n=5: d₅ = hc(T) since there's only one 5-set.
# And t₃ is the global count.

# For n=5, the relationship between t₃ and d₅ is:
# t₃ determines d₅ modulo 2 (but not the exact value) for t₃ ≤ 3 or t₃ = 5.

# Let's check: is d₅ ≡ something(t₃) mod 2 where something depends on ANOTHER invariant?

n = 5
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

print(f"\nn=5: Testing d₅ mod 2 formulas")

# Collect more invariants
formulas_tested = 0
formulas_working = 0

# Test: d₅ ≡ C(t₃, k) mod 2 for various k
from math import comb

for k in range(1, 6):
    matches = 0
    for bits in range(num_t):
        A = make_tournament(bits, n)
        t3 = count_3cycles(A, n)
        d5 = count_directed_5cycles(A, n)
        if d5 % 2 == comb(t3, k) % 2:
            matches += 1
    formulas_tested += 1
    if matches == num_t:
        print(f"  d₅ ≡ C(t₃, {k}) mod 2: EXACT MATCH ({matches}/{num_t})")
        formulas_working += 1
    elif matches > num_t * 0.9:
        print(f"  d₅ ≡ C(t₃, {k}) mod 2: {matches}/{num_t} ({100*matches/num_t:.1f}%)")

# Test: d₅ ≡ t₃*(t₃-1)*(t₃-2)/6 mod 2
# This is C(t₃, 3). Already tested above.

# Test: d₅ ≡ floor(t₃/3) mod 2
matches = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    d5 = count_directed_5cycles(A, n)
    if d5 % 2 == (t3 // 3) % 2:
        matches += 1
if matches == num_t:
    print(f"  d₅ ≡ floor(t₃/3) mod 2: EXACT MATCH")
else:
    print(f"  d₅ ≡ floor(t₃/3) mod 2: {matches}/{num_t}")

# Test: d₅ ≡ floor(t₃/2) mod 2  (this would detect the t₃=3 pattern)
matches = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    d5 = count_directed_5cycles(A, n)
    if d5 % 2 == (t3 // 2) % 2:
        matches += 1
if matches == num_t:
    print(f"  d₅ ≡ floor(t₃/2) mod 2: EXACT MATCH")
else:
    print(f"  d₅ ≡ floor(t₃/2) mod 2: {matches}/{num_t}")

# Test: d₅ ≡ (t₃ >> 1) & 1, i.e., bit 1 of t₃
matches = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    d5 = count_directed_5cycles(A, n)
    if d5 % 2 == (t3 >> 1) & 1:
        matches += 1
if matches == num_t:
    print(f"  d₅ ≡ bit_1(t₃) mod 2: EXACT MATCH *** ← t₃ ≡ 2,3 mod 4")
else:
    print(f"  d₅ ≡ bit_1(t₃) mod 2: {matches}/{num_t}")

# At n=5, this doesn't distinguish t₃=4 (undetermined).
# But "for determined cases" it would say: d₅ is odd when bit 1 of t₃ is 1, i.e., t₃ ∈ {2,3,6,7,...}
# But at n=5, t₃=2 gives d₅ even! So bit_1 alone doesn't work.

# Let me try involving more invariants.
# What about p₄ = number of directed 4-paths?
# Or transitive triples = C(n,3) - t₃?

# Try quadratic: d₅ ≡ t₃*(t₃-1)/2 mod 2
matches = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    d5 = count_directed_5cycles(A, n)
    if d5 % 2 == (t3*(t3-1)//2) % 2:
        matches += 1
if matches == num_t:
    print(f"  d₅ ≡ t₃(t₃-1)/2 mod 2: EXACT MATCH (= C(t₃,2))")
else:
    print(f"  d₅ ≡ t₃(t₃-1)/2 mod 2: {matches}/{num_t}")

# The key insight might be at the sub-tournament level.
# At n=6: d₅ = Σ_{|S|=5} hc(T[S]).
# t₃ of the sub-tournament T[S] is related to the global t₃ by:
# each global 3-cycle appears in (n-3) = 3 five-element subsets.

# ============================================================
# Analysis 3: Sub-tournament decomposition at n=6
# ============================================================

print(f"\n{'='*70}")
print("SUB-TOURNAMENT DECOMPOSITION at n=6")
print(f"{'='*70}")

n = 6
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

# For each tournament, compute:
# - t₃(T) = global 3-cycle count
# - For each 5-element subset S, compute t₃(T[S]) and hc(T[S])
# - d₅(T) = Σ hc(T[S])

# Check: d₅ mod 2 = Σ hc(T[S]) mod 2 = Σ (hc(T[S]) mod 2) mod 2

# Since at n=5, the parity of hc is determined by t₃ for most values:
# t₃(sub)=3 → hc odd (= 1 always at n=5)
# So: d₅ mod 2 = (# subtournaments with t₃(sub) = 3) mod 2?

# But at n=5, other t₃ values can also give odd hc.
# Actually looking at the n=5 data: t₃=4 can give d₅=1,2,3 (odd or even).

# Let me focus: at n=5, which t₃ values give odd hc?
# From the (t₃,d₅) joint:
# t₃=3: d₅=1 (odd). t₃=4: d₅ ∈ {1,2,3}. t₃=5: d₅=2 (even).
# So for t₃=4, hc can be 1 (odd) or 2 (even) or 3 (odd).

# At n=5, hc (directed Hamiltonian cycles, canonicalized) = d₅:
# 0 (even), 1 (odd), 2 (even), 3 (odd)

# Is hc mod 2 = t₃ mod 2... no, t₃=3→hc=1 (both odd), t₃=4→hc∈{1,2,3},
# t₃=5→hc=2. Doesn't match.

# Actually I think the cleanest approach is to just test ALL possible mod-2
# formulas involving a few common invariants.

# Let me try: the number of 4-paths P₄ as a complement to t₃.
# In a tournament on n=5: each ordered pair (i,j) with i→j contributes to
# some number of paths. Let me count P₄ = directed paths of length 4 (visiting 5 vertices).
# Actually H(T) at n=5 = P₅ (directed Hamiltonian paths).

# Let me try: use S₂ = Σ C(dᵢ,2) as the complementary invariant.
# t₃ = C(n,3) - S₂ (classical formula).

# Then d₅ mod 2 as function of (t₃, S₂)?
# But S₂ = C(n,3) - t₃, so they carry the same information!

# Try: number of "strongly connected components" at n=5?

# Actually, let me just try: d₅ mod 2 as function of score sequence at n=5.
print("\nn=5: d₅ mod 2 by (t₃, score_seq):")
n = 5
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

t3_ss_to_d5 = defaultdict(lambda: defaultdict(int))
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    ss = score_sequence(A, n)
    d5 = count_directed_5cycles(A, n)
    t3_ss_to_d5[(t3, ss)][d5 % 2] += 1

for key in sorted(t3_ss_to_d5.keys()):
    t3, ss = key
    parities = t3_ss_to_d5[key]
    det = "DET" if len(parities) == 1 else "MIX"
    print(f"  t₃={t3}, ss={ss}: {dict(parities)} [{det}]")

# Now n=6: is score sequence enough to determine d₅ mod 2?
print(f"\nn=6: d₅ mod 2 by (t₃, score_seq):")
n = 6
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

t3_ss_to_d5 = defaultdict(lambda: defaultdict(int))
count = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    ss = score_sequence(A, n)
    d5 = count_directed_5cycles(A, n)
    t3_ss_to_d5[(t3, ss)][d5 % 2] += 1
    count += 1
    if count % 10000 == 0:
        print(f"  ... processed {count}/{num_t}", file=sys.stderr)

mixed_count = 0
det_count = 0
for key in sorted(t3_ss_to_d5.keys()):
    t3, ss = key
    parities = t3_ss_to_d5[key]
    if len(parities) > 1:
        mixed_count += 1
        total_in_group = sum(parities.values())
        print(f"  t₃={t3}, ss={ss}: {dict(parities)} [MIXED, {total_in_group} tourn]")
    else:
        det_count += 1

print(f"\n  Determined by (t₃, ss): {det_count} groups")
print(f"  NOT determined: {mixed_count} groups")

# ============================================================
# Analysis 4: The sub-tournament HC parity sum
# ============================================================

print(f"\n{'='*70}")
print("SUB-TOURNAMENT HC PARITY SUM at n=6")
print(f"{'='*70}")

# d₅(T) = Σ_{|S|=5} hc(T[S])
# d₅ mod 2 = Σ hc(T[S]) mod 2

# Count: for how many of the 6 subtournaments is hc odd?
n = 6
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

# For each tournament, count the number of 5-vertex subtournaments with odd HC
hc_odd_count_dist = defaultdict(lambda: {'d5_even': 0, 'd5_odd': 0})

for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)

    # For each 5-vertex subset
    total_hc = 0
    num_odd_subs = 0
    for omit in range(n):
        verts = [v for v in range(n) if v != omit]
        # Build sub-tournament
        sub_A = [[0]*5 for _ in range(5)]
        for a in range(5):
            for b in range(5):
                if a != b:
                    sub_A[a][b] = A[verts[a]][verts[b]]

        # Count directed HCs of this 5-vertex tournament
        hc = 0
        for p in permutations(range(1, 5)):
            order = [0] + list(p)
            ok = True
            for idx in range(5):
                if sub_A[order[idx]][order[(idx+1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                hc += 1

        total_hc += hc
        if hc % 2 == 1:
            num_odd_subs += 1

    d5 = total_hc  # This IS d₅

    if d5 % 2 == 0:
        hc_odd_count_dist[num_odd_subs]['d5_even'] += 1
    else:
        hc_odd_count_dist[num_odd_subs]['d5_odd'] += 1

print("\n  # of 5-vertex subs with odd HC → d₅ parity:")
for k in sorted(hc_odd_count_dist.keys()):
    g = hc_odd_count_dist[k]
    total = g['d5_even'] + g['d5_odd']
    parity_of_k = "even" if k % 2 == 0 else "odd"
    d5_always_even = g['d5_odd'] == 0
    d5_always_odd = g['d5_even'] == 0
    consistent = (d5_always_even and k % 2 == 0) or (d5_always_odd and k % 2 == 1)
    marker = " ← d₅≡k mod 2!" if consistent and total > 0 else ""
    print(f"  k={k} ({parity_of_k}): d₅ even={g['d5_even']}, d₅ odd={g['d5_odd']} (total={total}){marker}")

# The point: d₅ = Σ hc(sub), so d₅ mod 2 = (# subs with odd hc) mod 2.
# This is tautological. But the INSIGHT is: which sub-tournament features determine hc parity?

# ============================================================
# Analysis 5: What invariant of T[S] determines hc(T[S]) mod 2?
# ============================================================

print(f"\n{'='*70}")
print("WHAT DETERMINES hc(T[S]) mod 2 at n=5?")
print(f"{'='*70}")

n = 5
total_edges = n * (n - 1) // 2
num_t = 2 ** total_edges

# For each 5-vertex tournament, we know t₃ and hc.
# t₃=3 → hc=1 (odd, always). But t₃=4 → hc ∈ {1,2,3}.
# What else separates hc=1 from hc=2 from hc=3 at t₃=4?

t3_4_data = []
for bits in range(num_t):
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    if t3 != 4:
        continue

    d5 = count_directed_5cycles(A, n)
    ss = score_sequence(A, n)

    # Compute: number of "kings" (vertices that beat everyone within 2 steps)
    # Actually, simpler: compute out-degree vector
    degs = tuple(sorted(sum(A[i]) for i in range(n)))

    # Compute: is the tournament "almost regular"?
    max_deg = max(sum(A[i]) for i in range(n))
    min_deg = min(sum(A[i]) for i in range(n))
    deg_range = max_deg - min_deg

    t3_4_data.append({
        'bits': bits, 'd5': d5, 'ss': ss, 'degs': degs,
        'deg_range': deg_range
    })

# Group by d5 value
from collections import Counter
d5_to_ss = defaultdict(list)
for d in t3_4_data:
    d5_to_ss[d['d5']].append(d['ss'])

print(f"\n  At t₃=4 (n=5), {len(t3_4_data)} tournaments:")
for d5 in sorted(d5_to_ss.keys()):
    ss_counts = Counter(d5_to_ss[d5])
    print(f"    d₅={d5} (hc mod 2 = {d5%2}): {len(d5_to_ss[d5])} tourn")
    for ss, cnt in sorted(ss_counts.items()):
        print(f"      ss={ss}: {cnt}")

# ============================================================
# Analysis 6: n=7 sample — d₅+d₇ parity
# ============================================================

print(f"\n{'='*70}")
print("n=7 SAMPLE: d₅+d₇ parity vs t₃ mod 4")
print(f"{'='*70}")

import random
random.seed(42)
n = 7
total_edges = n * (n - 1) // 2
sample_size = 200

t3_mod4_data = defaultdict(lambda: {'sum_even': 0, 'sum_odd': 0})

for trial in range(sample_size):
    bits = random.randint(0, 2**total_edges - 1)
    A = make_tournament(bits, n)
    t3 = count_3cycles(A, n)
    d5 = count_directed_5cycles(A, n)
    d7 = count_directed_7cycles(A, n)

    sum_d = d5 + d7
    key = t3 % 4
    if sum_d % 2 == 0:
        t3_mod4_data[key]['sum_even'] += 1
    else:
        t3_mod4_data[key]['sum_odd'] += 1

    if trial < 5:
        print(f"  trial {trial}: t₃={t3}, d₅={d5}, d₇={d7}, d₅+d₇={sum_d}, t₃%4={t3%4}")

print(f"\n  t₃ mod 4 → (d₅+d₇) parity ({sample_size} samples):")
for r in range(4):
    if r in t3_mod4_data:
        g = t3_mod4_data[r]
        total = g['sum_even'] + g['sum_odd']
        det = "ALWAYS EVEN" if g['sum_odd'] == 0 else ("ALWAYS ODD" if g['sum_even'] == 0 else "BOTH")
        print(f"    t₃ ≡ {r} (mod 4): even={g['sum_even']}, odd={g['sum_odd']} [{det}]")

# Also check: α₁ = t₃ + d₅ + d₇ mod 2
print(f"\n  α₁ mod 2 vs t₃ mod 4:")
for r in range(4):
    t3_mod4_alpha = {'even': 0, 'odd': 0}
    random.seed(42)
    for trial in range(sample_size):
        bits = random.randint(0, 2**total_edges - 1)
        A = make_tournament(bits, n)
        t3 = count_3cycles(A, n)
        if t3 % 4 != r:
            continue
        d5 = count_directed_5cycles(A, n)
        d7 = count_directed_7cycles(A, n)
        alpha1 = t3 + d5 + d7
        if alpha1 % 2 == 0:
            t3_mod4_alpha['even'] += 1
        else:
            t3_mod4_alpha['odd'] += 1
    total = t3_mod4_alpha['even'] + t3_mod4_alpha['odd']
    if total > 0:
        det = "ALWAYS EVEN" if t3_mod4_alpha['odd'] == 0 else ("ALWAYS ODD" if t3_mod4_alpha['even'] == 0 else "BOTH")
        print(f"    t₃ ≡ {r} (mod 4): α₁ even={t3_mod4_alpha['even']}, odd={t3_mod4_alpha['odd']} [{det}]")

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
