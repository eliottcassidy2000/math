#!/usr/bin/env python3
"""
Verify and prove: c3=3 forces all 3 triples to share a common vertex,
and span at most 5 vertices in their union.

For n=5,6: exhaustive.
For n=7: test all regular tournaments (2640) + sample non-regular with c3=3.

Also: prove the common-vertex property graph-theoretically.

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
import random

def get_cyclic_triples(T, n):
    """Return list of cyclic triples (as frozensets)."""
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

def scores_from_T(T, n):
    return tuple(sorted(sum(T[i]) for i in range(n)))

def analyze_c3_3(T, n):
    """Given T with c3=3, analyze the triple structure."""
    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        return None

    # Union of vertices
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)

    # Common vertex (in all 3 triples)
    common = triples[0] & triples[1] & triples[2]

    # Pairwise intersections
    pw = []
    for i in range(3):
        for j in range(i+1, 3):
            pw.append(len(triples[i] & triples[j]))

    return {
        'span': span,
        'common': common,
        'common_size': len(common),
        'pairwise_intersections': sorted(pw),
        'scores': scores_from_T(T, n)
    }

def generate_random_tournament(n):
    """Generate a random tournament on n vertices."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def c3_from_scores(scores, n):
    """Moon's formula: c3 = C(n,3) - sum C(s_i, 2)."""
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

# ============================================================
# PART 1: Exhaustive verification at n=5,6
# ============================================================
for n in [5, 6]:
    print(f"=== n={n}: Exhaustive verification ===")
    m = n*(n-1)//2
    count = 0
    all_share = True
    span_counts = {}
    common_size_counts = {}
    pw_patterns = {}
    score_set = set()

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        triples = get_cyclic_triples(T, n)
        if len(triples) != 3:
            continue
        count += 1
        info = analyze_c3_3(T, n)
        span_counts[info['span']] = span_counts.get(info['span'], 0) + 1
        cs = info['common_size']
        common_size_counts[cs] = common_size_counts.get(cs, 0) + 1
        if cs == 0:
            all_share = False
        pw = tuple(info['pairwise_intersections'])
        pw_patterns[pw] = pw_patterns.get(pw, 0) + 1
        score_set.add(info['scores'])

    print(f"  Total tournaments with c3=3: {count}")
    print(f"  Vertex span distribution: {span_counts}")
    print(f"  Common vertex sizes: {common_size_counts}")
    print(f"  All share a common vertex? {all_share}")
    print(f"  Pairwise intersection patterns: {pw_patterns}")
    print(f"  Score sequences: {score_set}")
    print()

# ============================================================
# PART 2: n=7 verification
# ============================================================
print("=== n=7: Verification ===")
n = 7

# First: which score sequences can give c3=3?
# c3 = C(7,3) - sum C(s_i,2) = 35 - sum s_i*(s_i-1)/2
# For c3=3: sum s_i*(s_i-1)/2 = 32
# sum s_i = C(7,2) = 21 always
# We need to find all compositions (s_0,...,s_6) with sum=21,
# 0 <= s_i <= 6, and sum s_i*(s_i-1)/2 = 32

print("Score sequences with c3=3 at n=7:")
print("  Need sum s_i(s_i-1)/2 = 32, sum s_i = 21")
valid_scores = []

def find_scores(n, target_sum, target_sq, prefix=[], min_val=0):
    """Find all sorted score sequences."""
    remaining = n - len(prefix)
    if remaining == 0:
        if target_sum == 0 and target_sq == 0:
            valid_scores.append(tuple(prefix))
        return
    current_sum = sum(prefix)
    for s in range(min_val, n):
        if target_sum - s < 0:
            break
        # Remaining values must be >= s and <= n-1
        min_remaining_sum = s * (remaining - 1)
        max_remaining_sum = (n-1) * (remaining - 1)
        if target_sum - s < min_remaining_sum or target_sum - s > max_remaining_sum:
            continue
        find_scores(n, target_sum - s, target_sq - s*(s-1)//2, prefix + [s], s)

find_scores(7, 21, 32)
for s in valid_scores:
    print(f"  {s}")

# Now enumerate tournaments with these score sequences at n=7
# For efficiency, use random sampling
print(f"\nSampling tournaments at n=7 with c3=3 (random)...")
random.seed(42)

found_c3_3 = 0
violations = 0
MAX_SAMPLES = 500000

for _ in range(MAX_SAMPLES):
    T = generate_random_tournament(n)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores, n)
    if c3 != 3:
        continue

    found_c3_3 += 1
    triples = get_cyclic_triples(T, n)
    assert len(triples) == 3, f"Moon's formula mismatch: {len(triples)} vs expected 3"

    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)

    common = triples[0] & triples[1] & triples[2]

    if len(common) == 0:
        print(f"  VIOLATION: No common vertex! triples={[sorted(t) for t in triples]}")
        violations += 1
    if span > 5:
        print(f"  VIOLATION: Span={span} > 5! triples={[sorted(t) for t in triples]}")
        violations += 1

print(f"  Sampled {MAX_SAMPLES} random tournaments, found {found_c3_3} with c3=3")
print(f"  Violations: {violations}")

# ============================================================
# PART 3: Exhaustive n=7 for score sequences that allow c3=3
# ============================================================
# Since full n=7 enumeration is 2^21 ~ 2M, that's feasible
print(f"\nn=7: Exhaustive verification (2^21 = {1<<21} tournaments)...")
n = 7
m = n*(n-1)//2  # 21
count = 0
all_share = True
span_counts = {}
common_size_counts = {}

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores, n)
    if c3 != 3:
        continue

    triples = get_cyclic_triples(T, n)
    assert len(triples) == 3
    count += 1

    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span = len(all_verts)
    span_counts[span] = span_counts.get(span, 0) + 1

    common = triples[0] & triples[1] & triples[2]
    cs = len(common)
    common_size_counts[cs] = common_size_counts.get(cs, 0) + 1
    if cs == 0:
        all_share = False
        print(f"  COUNTEREXAMPLE at n=7! bits={bits}, triples={[sorted(t) for t in triples]}")

print(f"  Total tournaments with c3=3: {count}")
print(f"  Span distribution: {span_counts}")
print(f"  Common vertex sizes: {common_size_counts}")
print(f"  All share a common vertex? {all_share}")

# ============================================================
# PART 4: Theoretical analysis - why must 3 triples share a vertex?
# ============================================================
print(f"\n=== THEORETICAL ANALYSIS ===")
print("""
Consider 3 cyclic triples T1, T2, T3 in a tournament.
Each triple is a set of 3 vertices. Their union has at most 9 vertices.

Claim: If exactly 3 cyclic triples exist, they must all share a common vertex.

Proof approach by contradiction:
Suppose no vertex is in all 3 triples.

Case 1: Some pair of triples is disjoint (share 0 vertices).
  Then T1={a,b,c}, T2={d,e,f}, T3 shares at most 2 vertices with each.
  Union spans >= 7 vertices.
  But with only 3 cyclic triples total, the remaining C(n,3)-3 triples
  must all be transitive.

Case 2: Every pair shares exactly 1 vertex, but no vertex is in all 3.
  T1={a,b,c}, T2={a,d,e} (share a), T3={b,d,f} (share b with T1, d with T2)
  Union = {a,b,c,d,e,f}, span = 6.

Case 3: Some pair shares 2 vertices.
  T1={a,b,c}, T2={a,b,d} (share a,b).
  T3 doesn't contain both a and b (else it'd share a common vertex with both).
  If T3 contains a but not b: T3={a,e,f}. Common vertex a is in T1,T2,T3 -> contradiction.
  If T3 contains b but not a: T3={b,e,f}. Similar.
  If T3 contains neither: T3={c,d,e} or similar.

Actually Case 3 is the key: if two triples share 2 vertices {a,b},
then {a,b,c} and {a,b,d} are both cyclic.
The edge a->b (WLOG) means:
  In T1: a->b->c->a (cyclic)
  In T2: a->b->d->a (cyclic)
So both c and d beat a and lose to b (in the cyclic orientation with a->b).
Wait: a->b->c->a means c->a, b->c. And a->b->d->a means d->a, b->d.
So we know: a->b, b->c, c->a, b->d, d->a.
What about c vs d?
If c->d: {c,d,a} has c->d->a->... wait, a->c? No: c->a. So c->d, d->a, a->c?
  That's NOT a cycle since we need c->d->a->c but a->c means the triple {a,c,d}
  has a->c, c->d, d->a which IS a cycle! That's a 4th cyclic triple.
If d->c: similarly {a,c,d} has d->c, c->a, a->d? We have d->a not a->d.
  Actually: d->c, c->a, and between a,d we have d->a.
  So the triple {a,c,d}: d->c->a and d->a. That's transitive (d beats both).

So if c->d, we get a 4th cyclic triple {a,c,d} -> contradiction with c3=3.
If d->c, {a,c,d} is transitive. OK.

Now T3 must not contain {a,b} (else the common vertex argument applies).
But T3 must create a cycle. With d->c forced, what constraints exist?
""")

# Let's verify the Case 2 analysis computationally
print("Verifying Case 2: Can 3 triples have pairwise intersection = {1,1,1}?")
print("(Each pair shares exactly 1 vertex, no common vertex to all 3)")

# At n=6, check if any c3=3 tournament has this pattern
n = 6
m = n*(n-1)//2
case2_found = False
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue
    common = triples[0] & triples[1] & triples[2]
    if len(common) == 0:
        pw = []
        for i in range(3):
            for j in range(i+1, 3):
                pw.append(len(triples[i] & triples[j]))
        print(f"  Case 2 example: bits={bits}, pw={sorted(pw)}, triples={[sorted(t) for t in triples]}")
        case2_found = True
        break

if not case2_found:
    print("  Case 2 NEVER occurs at n=6!")

# Also check: can two triples be disjoint?
print("\nCan two of the 3 triples be vertex-disjoint?")
disjoint_found = False
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue
    for i in range(3):
        for j in range(i+1, 3):
            if len(triples[i] & triples[j]) == 0:
                print(f"  Disjoint pair: {sorted(triples[i])}, {sorted(triples[j])}")
                disjoint_found = True
                break
        if disjoint_found:
            break
    if disjoint_found:
        break

if not disjoint_found:
    print("  No disjoint pairs at n=6!")

# Check what pairwise intersection patterns actually occur
print("\nPairwise intersection patterns at n=6 with c3=3:")
pw_counts = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue
    pw = []
    for i in range(3):
        for j in range(i+1, 3):
            pw.append(len(triples[i] & triples[j]))
    key = tuple(sorted(pw))
    pw_counts[key] = pw_counts.get(key, 0) + 1

for pw, cnt in sorted(pw_counts.items()):
    print(f"  {pw}: {cnt} tournaments")

print("\nDone.")
