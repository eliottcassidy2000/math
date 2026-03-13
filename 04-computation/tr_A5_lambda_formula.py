#!/usr/bin/env python3
"""
Derive the exact formula tr(A^5) in terms of lambda values.

Key identity: λ(i,j) = A_{ij}·(A²)_{ji} + A_{ji}·(A²)_{ij}
where (A²)_{ij} = #{k ≠ i,j : i→k→j}

Strategy: Express tr(A^5) = Σ_i Σ_{j,k,l,m} A_{ij}A_{jk}A_{kl}A_{lm}A_{mi}
by decomposing into cases based on which indices coincide.

For a tournament: A_{ij} + A_{ji} = 1 for i≠j, A_{ii} = 0.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
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

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

# First, let's understand the structure of tr(A^5) more carefully.
# tr(A^5) = Σ A_{i0,i1} A_{i1,i2} A_{i2,i3} A_{i3,i4} A_{i4,i0}
# where i0,...,i4 range over [n].
# Since A_{ii}=0, we need all consecutive pairs to be distinct.
# But non-consecutive indices CAN coincide.

# Decompose by the NUMBER of distinct vertices visited:
# 5 distinct: directed 5-cycles (counted 5 times each)
# 4 distinct: one repeated → "lollipop" type walks
# 3 distinct: two repeated → walks on 3 vertices
# 2 distinct: impossible (need 5 steps alternating between 2 vertices,
#             but A_{ij}A_{ji}=0 for tournaments)

# Let's verify this decomposition computationally at n=5.
print("=" * 60)
print("TRACE DECOMPOSITION OF tr(A^5)")
print("=" * 60)

n = 5
total_bits = n * (n - 1) // 2

# For each tournament, decompose tr(A^5) by number of distinct vertices
results = []
np.random.seed(42)

for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    A5 = np.linalg.matrix_power(A, 5)
    tr5 = int(np.trace(A5))
    L = lambda_graph(A, n)

    # Count 5-cycles directly
    c5_dir = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = all(A[path[i]][path[(i+1) % 5]] for i in range(5))
            if valid:
                c5_dir += 1

    # Decompose walks by distinct vertex count
    walk_by_dist = defaultdict(int)
    for i0 in range(n):
        for i1 in range(n):
            if i1 == i0: continue
            if A[i0][i1] == 0: continue
            for i2 in range(n):
                if i2 == i1: continue
                if A[i1][i2] == 0: continue
                for i3 in range(n):
                    if i3 == i2: continue
                    if A[i2][i3] == 0: continue
                    for i4 in range(n):
                        if i4 == i3: continue
                        if A[i3][i4] == 0: continue
                        if A[i4][i0] == 0: continue
                        if i4 == i0: continue  # A_{i0,i0}=0
                        ndist = len(set([i0,i1,i2,i3,i4]))
                        walk_by_dist[ndist] += 1

    total = sum(walk_by_dist.values())
    assert total == tr5, f"Mismatch: {total} vs {tr5}"

    # Lambda power sums
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)

    # Degree sequence
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
    D2 = sum(int(sum(A[i]))**2 for i in range(n))

    results.append({
        'bits': bits, 'tr5': tr5, 'c5': c5_dir,
        'w5': walk_by_dist.get(5, 0),
        'w4': walk_by_dist.get(4, 0),
        'w3': walk_by_dist.get(3, 0),
        'w2': walk_by_dist.get(2, 0),
        'S1': S1, 'S2': S2, 'S3': S3, 'D2': D2,
        'scores': scores
    })

# Verify: w5 = 5 * c5_dir always
for r in results:
    assert r['w5'] == 5 * r['c5'], f"Walk-5 != 5*c5: {r['w5']} vs {5*r['c5']}"
print("w5 = 5*c5_dir: VERIFIED for all 1024 tournaments")

# Check: is w4 (4-distinct-vertex walks) determined by lambda?
print(f"\n--- 4-vertex walks (w4) ---")
lam_groups_w4 = defaultdict(set)
for r in results:
    key = (r['S1'], r['S2'], r['S3'])
    lam_groups_w4[key].add(r['w4'])
ambig_w4 = sum(1 for v in lam_groups_w4.values() if len(v) > 1)
print(f"  w4 determined by (S1,S2,S3)? Ambiguous groups: {ambig_w4}")

# Check: is w3 determined by lambda?
lam_groups_w3 = defaultdict(set)
for r in results:
    key = (r['S1'], r['S2'], r['S3'])
    lam_groups_w3[key].add(r['w3'])
ambig_w3 = sum(1 for v in lam_groups_w3.values() if len(v) > 1)
print(f"  w3 determined by (S1,S2,S3)? Ambiguous groups: {ambig_w3}")

# Check: is tr(A^5) determined by lambda?
lam_groups_tr = defaultdict(set)
for r in results:
    key = (r['S1'], r['S2'], r['S3'])
    lam_groups_tr[key].add(r['tr5'])
ambig_tr = sum(1 for v in lam_groups_tr.values() if len(v) > 1)
print(f"  tr(A^5) determined by (S1,S2,S3)? Ambiguous groups: {ambig_tr}")

# Now let's try to find what DETERMINES w4.
# w4 counts closed walks of length 5 on exactly 4 distinct vertices.
# Pattern: one vertex is visited twice. The walk is:
# i0 → i1 → i2 → i3 → i4 → i0 where exactly one pair (ia, ib) with a<b are equal and (b-a) >= 2.
# Possible collapse patterns (which pair coincides):
# (i0=i2): walk = i0 → i1 → i0 → i3 → i4 → i0: but A_{i0,i1}A_{i1,i0}=0 for tournaments!
# Wait: A_{i1,i0} might be 0. The step i1 → i0 is A_{i1,i0}. For this walk to exist we need
# A_{i0,i1} * A_{i1,i0} = 1, which means BOTH edges exist — impossible in a tournament.

# Wait, that's wrong. i0=i2 means the walk is:
# i0 → i1 → i2=i0 → i3 → i4 → i0
# This requires A_{i0,i1}=1, A_{i1,i0}=1 (step i1→i2=i0), which is impossible.

# So i0=i2 is impossible! Similarly for any pair (ia, i_{a+1 mod 5}) — consecutive equals already excluded.
# And ia = i_{a+2 mod 5} requires a back-and-forth through the middle vertex — impossible in tournament.

# Let's enumerate which collapse patterns are possible:
# 5 vertices i0,i1,i2,i3,i4 with exactly 4 distinct. One pair equal, non-consecutive.
# Non-consecutive pairs: (i0,i2), (i0,i3), (i1,i3), (i1,i4), (i2,i4)
# But we also need no CONSECUTIVE pair equal (which is already excluded by A_{ii}=0).

# Check: (i0,i2) same → need A_{i1,i2=i0} = A_{i1,i0} and A_{i0,i1} both 1 → impossible.
# (i0,i3): walk = i0,i1,i2,i0,i4,i0 → need A_{i2,i0}=1 and A_{i0,i4}=1 and A_{i4,i0}=1 → last two contradict!
# Wait, A_{i0,i4}·A_{i4,i0} = 0. So i3=i0 means we need i3→i4→i0, so A_{i0,i4}=1 AND A_{i4,i0}=1? No:
# Step i3→i4: A_{i3,i4} = A_{i0,i4} (since i3=i0)
# Step i4→i0: A_{i4,i0}
# Need both = 1, so A_{i0,i4}=1 and A_{i4,i0}=1 → impossible in tournament!

# Hmm, so (i0,i3) is also impossible?

# Let me reconsider. Walk: i0→i1→i2→i3→i4→i0
# If i0=i3: walk is i0→i1→i2→i0→i4→i0
#   Need: A_{i2,i0}·A_{i0,i4}·A_{i4,i0} ≥ 1
#   But A_{i0,i4}·A_{i4,i0} = 0. So IMPOSSIBLE.

# If i1=i3: walk is i0→i1→i2→i1→i4→i0
#   Need: A_{i2,i1}·A_{i1,i4}: OK, both can be 1 independently
#   But also A_{i0,i1}·A_{i1,i2}·A_{i2,i1}: A_{i1,i2}·A_{i2,i1} = 0!

# If i1=i4: walk is i0→i1→i2→i3→i1→i0
#   Need: A_{i3,i1}·A_{i1,i0}: both can be 1
#   And A_{i0,i1}·A_{i1,i2}·A_{i2,i3}: all can be 1

# If i2=i4: walk is i0→i1→i2→i3→i2→i0
#   Need: A_{i3,i2}·A_{i2,i0}: both can be 1
#   And A_{i0,i1}·A_{i1,i2}·A_{i2,i3}·A_{i3,i2}: A_{i2,i3}·A_{i3,i2} = 0!

# So the ONLY possible 4-distinct-vertex pattern is i1=i4!
# Walk: i0→i1→i2→i3→i1→i0 with i0,i1,i2,i3 all distinct

print(f"\n--- Analyzing 4-vertex walk patterns ---")

# So w4 counts walks of the form: i0→i1→i2→i3→i1→i0
# with {i0,i1,i2,i3} = 4 distinct vertices
# Conditions: A_{01}=1, A_{12}=1, A_{23}=1, A_{31}=1, A_{10}=1
# Last two: A_{31}·A_{10}=1, meaning i3→i1 and i1→i0
# So: i0→i1, i1→i2, i2→i3, i3→i1, i1→i0
# The 3-cycle (i1,i2,i3) is a directed 3-cycle i1→i2→i3→i1
# Plus i0→i1→i0: but A_{01}·A_{10}=1 means both edges exist — IMPOSSIBLE!

# Wait, A_{01} = A[i0][i1] = 1 means i0→i1.
# A_{10} = A[i1][i0] = 1 means i1→i0.
# For a tournament, only one of these can be 1. CONTRADICTION.

# So w4 = 0 for all tournaments??

print(f"  w4 values: {set(r['w4'] for r in results)}")

# Check w3
print(f"  w3 values: {set(r['w3'] for r in results)}")

# And w2
print(f"  w2 values: {set(r['w2'] for r in results)}")

# So tr(A^5) = 5*c5_dir + w3 (if w4=w2=0)
# or tr(A^5) = 5*c5_dir + (3-vertex closed walks of length 5)

# 3-vertex closed walks of length 5:
# Walk on {a,b,c} of length 5: a→?→?→?→?→a
# Each step must be between different vertices of {a,b,c}.
# So the walk alternates, but on 3 vertices with 5 steps that form a cycle.
# Possible patterns: a→b→c→a→b→a (last step b→a)
# Actually, any 5-step walk on {a,b,c} with no consecutive repeats.
# Since we're on a tournament with 3 vertices, exactly one of the 2 directed 3-cycles exists.

# For a 3-element tournament:
# If a→b→c→a (directed 3-cycle): all arcs go "forward"
# Then the possible 5-walks that close are:
# a,b,c,a,b,c (6 steps = would need to end at a after 5 steps)
# Actually, let me think step by step. With arcs a→b, b→c, c→a:
# From a: can go to b only (since a→b, c→a means a doesn't go to c)
# Wait: a→b and c→a, but what about a→c? Since there's a 3-cycle a→b→c→a, we have a→b, b→c, c→a.
# But tournaments are complete: a→c or c→a? Since c→a in the cycle, we have c→a, so a does NOT go to c.
# So from a we can ONLY go to b. From b, only to c. From c, only to a.
# So the ONLY walk is a→b→c→a→b→c→a→... which is the cycle repeated.
# A closed walk of length 5 starting at a: a→b→c→a→b→c but this is 5 steps ending at c, not a. Not closed!
# a→b→c→a→b→c→a is 6 steps. So no 5-step closed walk on a 3-vertex tournament with a 3-cycle!

# What about 3 vertices WITHOUT a 3-cycle (transitive: a→b, a→c, b→c)?
# From a: can go to b or c.
# Possible walks: a→b→c→?  From c, c can go to... in transitive, c→nobody (c is the sink).
# Wait, c→a? No, a→c. c→b? No, b→c. So c is a sink — no outgoing edges. Dead end.
# From a: a→c→dead. a→b→c→dead.
# So NO closed walks of length 5 on 3 vertices!

# This means w3 = 0 too. So tr(A^5) = 5 * c5_dir?
print(f"\n  CONCLUSION: tr(A^5) = 5*c5_dir? {all(r['tr5'] == 5*r['c5'] for r in results)}")

# Now: is tr(A^5) lambda-determined?
print(f"\n--- Is tr(A^5) lambda-determined? ---")
# We need: tr(A^5) = some function of lambda values.
# Since A² involves lambda (λ(i,j) = A_{ij}(A²)_{ji} + A_{ji}(A²)_{ij}),
# and tr(A^5) = tr(A·A²·A²), this is plausible.

# tr(A^5) = tr(A · A^4) = tr(A · (A²)²)
# = Σ_i Σ_j A_{ij} (A^4)_{ji}
# = Σ_i Σ_j A_{ij} Σ_k (A²)_{jk} (A²)_{ki}

# Now (A²)_{jk} depends on the full adjacency, not just lambda.
# But we know: (A²)_{jk} for j≠k on a tournament satisfies
# (A²)_{jk} + (A²)_{kj} = n - 2 - λ(j,k) + λ(j,k) ... no that's not right.

# Let me compute (A²)_{jk} + (A²)_{kj} for j≠k:
# = #{m: j→m→k} + #{m: k→m→j}
# Every m ≠ j,k contributes to exactly one of these (since j→m or m→j, and m→k or k→m).
# Actually, m contributes to (A²)_{jk} iff j→m AND m→k, i.e., A_{jm}=1, A_{mk}=1.
# m contributes to (A²)_{kj} iff k→m AND m→j, i.e., A_{km}=1, A_{mj}=1.
# The four cases for m ≠ j,k (using A_{jm}+A_{mj}=1, A_{mk}+A_{km}=1):
# j→m, m→k: contributes to (A²)_{jk}
# j→m, k→m: neither
# m→j, m→k: neither
# m→j, k→m: contributes to (A²)_{kj}
# So (A²)_{jk} + (A²)_{kj} = #{m: j→m, m→k} + #{m: m→j, k→m}
# = #{m: j→m→k} + #{m: k→m→j}
# But there are also the other two cases:
# #{m: j→m, k→m} = #{common out-neighbors of j and k} (not counting each other)
# #{m: m→j, m→k} = #{common in-neighbors}

# Total: all four cases sum to n-2 (number of vertices besides j,k).
# So: (A²)_{jk} + (A²)_{kj} + #{common out-neighbors} + #{common in-neighbors} = n-2.

# Now, λ(j,k) = A_{jk}(A²)_{kj} + A_{kj}(A²)_{jk}
# If j→k: λ(j,k) = (A²)_{kj} = #{m: k→m→j}
# If k→j: λ(j,k) = (A²)_{jk} = #{m: j→m→k}
# In either case, λ(j,k) picks out one of the two terms.

# The OTHER term:
# If j→k: (A²)_{jk} = #{m: j→m→k} -- this is NOT λ(j,k)
# Hmm. (A²)_{jk} is always #{m: j→m→k}, regardless of the j-k edge direction.
# If j→k: λ(j,k) = (A²)_{kj}, so (A²)_{jk} = (n-2) - (A²)_{kj} - co - ci
#   where co = #{common out-neighbors}, ci = #{common in-neighbors}
#   But co and ci depend on more than just lambda...

# Actually, there IS a formula:
# If j→k: (A²)_{jk} = d_j - 1 - λ(j,k) ... let me check.
# (A²)_{jk} = #{m≠j,k: j→m, m→k} = #{m: j→m} - A_{jk} - #{m≠k: j→m, k→m}
# Wait this is getting complicated. Let me just verify computationally.

# Key question: can tr(A^5) = Σ f(lambda)?
# Since tr(A^5) = 5*c5_dir and c5 depends on (S1,S2,S3) at n=5,
# the answer is YES for n=5.

# But does this extend to n=6, n=7?
# At n=7, kind-pasteur found that c7 is NOT lambda-determined.
# But c5 IS lambda-determined at n=7 (50k samples, 0 ambiguities).
# The question: what makes c5 special?

# INSIGHT: tr(A^5) = 5*c5_dir (no degenerate walks!)
# So c5_dir = tr(A^5)/5.
# tr(A^5) = tr(A·(A²)²) = Σ_{i,j} A_{ij} ((A²)²)_{ji}
# = Σ_{i,j} A_{ij} Σ_k (A²)_{jk}(A²)_{ki}

# Now (A²)_{jk} for j≠k:
# Define P_{jk} = (A²)_{jk} = #{m: j→m→k}
# And note: P_{jk} depends on the adjacency matrix, NOT just lambda.
# But lambda IS a specific combination of P values.

# However, tr(A^5) = Σ A_{ij} Σ_k P_{jk} P_{ki}
# This is a 3-fold sum involving A and P.

# Can we express this in terms of λ? Let's try.
# Note: Σ_k P_{jk} P_{ki} = ((A²)²)_{ji} = (A^4)_{ji}
# So tr(A^5) = Σ_{i,j} A_{ij} (A^4)_{ji} = tr(A·A^4) = tr(A^5). Circular.

# Let's try a DIFFERENT decomposition.
# tr(A^5) = tr((A²)²·A) = Σ_{i,j} (A^4)_{ij} A_{ji}

# Key: A_{ji} = 1 - A_{ij} for i≠j. So:
# tr(A^5) = Σ_{i≠j} (A^4)_{ij}(1 - A_{ij}) = Σ_{i≠j} (A^4)_{ij} - Σ_{i≠j} (A^4)_{ij}A_{ij}
# = tr((A^4)·J') - tr(A^5) where J' = J - I (all-ones minus identity)
# Wait: Σ_{i≠j} (A^4)_{ij} = Σ_i (Σ_j (A^4)_{ij} - (A^4)_{ii}) = Σ_i (A^4 · 1)_i - tr(A^4)

# This is getting circular. Let me just verify the key claim at n=6,7.

print(f"\n{'='*60}")
print("VERIFICATION AT n=6: tr(A^5) = 5*c5_dir?")
print(f"{'='*60}")

n6 = 6
tb6 = n6*(n6-1)//2
np.random.seed(42)
count = 0
fail = 0
for trial in range(10000):
    bits = np.random.randint(0, 1 << tb6)
    A = bits_to_adj(bits, n6)
    A5 = np.linalg.matrix_power(A, 5)
    tr5 = int(np.trace(A5))

    # Count c5_dir directly
    c5_total = 0
    for combo in combinations(range(n6), 5):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = all(A[path[i]][path[(i+1) % 5]] for i in range(5))
            if valid:
                c5_total += 1

    if tr5 != 5 * c5_total:
        fail += 1
        if fail <= 5:
            # Decompose
            w_count = defaultdict(int)
            for i0 in range(n6):
                for i1 in range(n6):
                    if i1==i0 or A[i0][i1]==0: continue
                    for i2 in range(n6):
                        if i2==i1 or A[i1][i2]==0: continue
                        for i3 in range(n6):
                            if i3==i2 or A[i2][i3]==0: continue
                            for i4 in range(n6):
                                if i4==i3 or A[i3][i4]==0: continue
                                if i4==i0 or A[i4][i0]==0: continue
                                nd = len({i0,i1,i2,i3,i4})
                                w_count[nd] += 1
            print(f"  FAIL at trial {trial}: tr5={tr5}, 5*c5={5*c5_total}, walk decomp={dict(w_count)}")
    count += 1

print(f"  Tested {count} tournaments, failures: {fail}")
if fail == 0:
    print("  tr(A^5) = 5*c5_dir CONFIRMED at n=6!")
else:
    print(f"  tr(A^5) ≠ 5*c5_dir in {fail} cases — degenerate walks exist at n≥6")

# If tr(A^5) = 5*c5_dir at all n, then c5_dir = tr(A^5)/5.
# And tr(A^5) = tr(A · A^4) = tr(A² · A³) = tr(A² · A · A²)
# Since A² = f(lambda, scores), this gives a formula for c5 in terms of lambda.
# More precisely: (A²)_{ij} = P_{ij} depends on ALL edges, but tr(A^5) might still
# be expressible as a polynomial in lambda values through cancellation.

# Let's check: does tr(A^5) depend on lambda at n=6?
print(f"\n{'='*60}")
print("IS tr(A^5) LAMBDA-DETERMINED AT n=6?")
print(f"{'='*60}")

# This is expensive exhaustively. Sample.
np.random.seed(42)
lam_groups = defaultdict(set)
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb6)
    A = bits_to_adj(bits, n6)
    A5 = np.linalg.matrix_power(A, 5)
    tr5 = int(np.trace(A5))
    L = lambda_graph(A, n6)
    lam_key = tuple(sorted(L[i][j] for i in range(n6) for j in range(i+1, n6)))
    lam_groups[lam_key].add(tr5)

ambig = sum(1 for v in lam_groups.values() if len(v) > 1)
print(f"  Lambda multiset groups: {len(lam_groups)}, ambiguous: {ambig}")
if ambig > 0:
    for key, vals in sorted(lam_groups.items()):
        if len(vals) > 1:
            print(f"    multiset={key}: tr5 values = {sorted(vals)}")
            if ambig > 5:
                print("    ... (truncated)")
                break

# Now check at n=7
print(f"\n{'='*60}")
print("IS tr(A^5) LAMBDA-DETERMINED AT n=7?")
print(f"{'='*60}")

n7 = 7
tb7 = n7*(n7-1)//2
np.random.seed(42)
lam_groups_7 = defaultdict(set)
for trial in range(20000):
    bits = np.random.randint(0, 1 << tb7)
    A = bits_to_adj(bits, n7)
    A5 = np.linalg.matrix_power(A, 5)
    tr5 = int(np.trace(A5))
    L = lambda_graph(A, n7)
    lam_key = tuple(sorted(L[i][j] for i in range(n7) for j in range(i+1, n7)))
    lam_groups_7[lam_key].add(tr5)

ambig_7 = sum(1 for v in lam_groups_7.values() if len(v) > 1)
print(f"  Lambda multiset groups: {len(lam_groups_7)}, ambiguous: {ambig_7}")
if ambig_7 > 0:
    for key, vals in sorted(lam_groups_7.items()):
        if len(vals) > 1:
            print(f"    multiset={key}: tr5 values = {sorted(vals)}")
            if sum(1 for v in lam_groups_7.values() if len(v) > 1) > 5:
                break

# Now check tr(A^7) at n=7 — this should NOT be lambda-determined (since c7 isn't)
print(f"\n{'='*60}")
print("IS tr(A^7) LAMBDA-DETERMINED AT n=7?")
print(f"{'='*60}")

lam_groups_tr7 = defaultdict(set)
for trial in range(20000):
    bits = np.random.randint(0, 1 << tb7)
    A = bits_to_adj(bits, n7)
    A7 = np.linalg.matrix_power(A, 7)
    tr7 = int(np.trace(A7))
    L = lambda_graph(A, n7)
    lam_key = tuple(sorted(L[i][j] for i in range(n7) for j in range(i+1, n7)))
    lam_groups_tr7[lam_key].add(tr7)

ambig_tr7 = sum(1 for v in lam_groups_tr7.values() if len(v) > 1)
print(f"  Lambda multiset groups: {len(lam_groups_tr7)}, ambiguous: {ambig_tr7}")
if ambig_tr7 > 0:
    count_shown = 0
    for key, vals in sorted(lam_groups_tr7.items()):
        if len(vals) > 1:
            print(f"    multiset={key}: tr7 values = {sorted(vals)}")
            count_shown += 1
            if count_shown >= 5:
                break

print(f"\nDone.")
