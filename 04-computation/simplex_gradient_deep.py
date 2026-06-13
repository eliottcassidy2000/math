"""
simplex_gradient_deep.py -- kind-pasteur-2026-03-13-S61

Deep analysis of the 3D internal space gradient:
1. What distinguishes the 6 ambiguous simplex profiles?
2. The exact sigma:lambda:delta ratio as function of score sequence
3. Fiber bundle parallel transport: how Vitali atoms move on the simplex
4. The "hidden dimension" — what information BEYOND the simplex profile
   is needed to determine c7?
5. Connection to path homology: is the simplex position a homological invariant?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

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
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

def get_simplex_data(A, n):
    """Get full (sigma, lambda, delta) for all pairs."""
    L = lambda_graph(A, n)
    A2 = A @ A
    pairs = []
    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])
            delta = n - 2 - sig - lam
            pairs.append((sig, lam, delta))
    return pairs

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("DEEP SIMPLEX GRADIENT ANALYSIS at n=7")
print("=" * 60)

np.random.seed(42)

# Collect tournament data
all_data = []
for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    pairs = get_simplex_data(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    c3 = count_directed_k_cycles(A, n, 3)
    c5 = count_directed_k_cycles(A, n, 5)
    scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))

    # Simplex profile = sorted multiset of (sig, lam, delta)
    profile = tuple(sorted(pairs))
    # Totals
    ts = sum(s for s, l, d in pairs)
    tl = sum(l for s, l, d in pairs)
    td = sum(d for s, l, d in pairs)

    all_data.append({
        'bits': bits,
        'profile': profile,
        'pairs': pairs,
        'c7': c7,
        'c3': c3,
        'c5': c5,
        'scores': scores,
        'total_sigma': ts,
        'total_lambda': tl,
        'total_delta': td,
    })

# 1. Find the 6 ambiguous profiles
print("\n--- Ambiguous Simplex Profiles ---")
profile_c7 = defaultdict(set)
profile_examples = defaultdict(list)
for entry in all_data:
    profile_c7[entry['profile']].add(entry['c7'])
    profile_examples[entry['profile']].append(entry)

ambig_profiles = {p: c7s for p, c7s in profile_c7.items() if len(c7s) > 1}
print(f"Total distinct profiles: {len(profile_c7)}")
print(f"Ambiguous profiles: {len(ambig_profiles)}")

for i, (profile, c7s) in enumerate(sorted(ambig_profiles.items(), key=lambda x: -len(x[1]))):
    examples = profile_examples[profile]
    c7_counts = Counter(ex['c7'] for ex in examples)
    scores_set = set(ex['scores'] for ex in examples)
    c3_set = set(ex['c3'] for ex in examples)
    c5_set = set(ex['c5'] for ex in examples)

    ts = examples[0]['total_sigma']
    tl = examples[0]['total_lambda']
    td = examples[0]['total_delta']

    print(f"\n  Ambiguous profile {i+1}:")
    print(f"    c7 values: {sorted(c7s)}")
    print(f"    c7 distribution: {dict(sorted(c7_counts.items()))}")
    print(f"    Total (sig,lam,del) = ({ts}, {tl}, {td})")
    print(f"    Score sequences: {scores_set}")
    print(f"    c3 values: {sorted(c3_set)}")
    print(f"    c5 values: {sorted(c5_set)}")
    print(f"    n_examples: {len(examples)}")

    # What distinguishes them?
    # Check: is it the ORDERING of pairs (not just multiset)?
    pair_lists = defaultdict(set)
    for ex in examples:
        pair_lists[ex['c7']].add(tuple(ex['pairs']))  # ordered pairs
    print(f"    Distinct ordered pair-lists per c7:")
    for c7_val in sorted(c7s):
        print(f"      c7={c7_val}: {len(pair_lists[c7_val])} distinct ordered pair-lists")

# 2. The sigma:lambda:delta ratio by score sequence
print("\n\n--- Sigma:Lambda:Delta by Score Sequence ---")
score_simplex = defaultdict(list)
for entry in all_data:
    score_simplex[entry['scores']].append(entry)

print(f"Distinct score sequences: {len(score_simplex)}")
for scores in sorted(score_simplex.keys()):
    entries = score_simplex[scores]
    avg_ts = np.mean([e['total_sigma'] for e in entries])
    avg_tl = np.mean([e['total_lambda'] for e in entries])
    avg_td = np.mean([e['total_delta'] for e in entries])
    avg_c7 = np.mean([e['c7'] for e in entries])
    n_ex = len(entries)
    if n_ex >= 5:  # only show well-sampled scores
        print(f"  {scores}: n={n_ex:>4}, avg S:L:D = {avg_ts:.1f}:{avg_tl:.1f}:{avg_td:.1f}, "
              f"avg c7={avg_c7:.1f}")

# 3. Total_lambda = 3*c3 verification
print("\n\n--- Total Lambda = 3*c3 Verification ---")
match = 0
mismatch = 0
for entry in all_data[:1000]:
    if entry['total_lambda'] == 3 * entry['c3']:
        match += 1
    else:
        mismatch += 1
        print(f"  MISMATCH: total_lambda={entry['total_lambda']}, 3*c3={3*entry['c3']}")
print(f"  Match: {match}/{match+mismatch}")

# 4. Total sigma formula: sum_v [C(indeg(v),2) + C(outdeg(v),2)]
print("\n\n--- Total Sigma Formula Verification ---")
match = 0
for entry in all_data[:1000]:
    bits = entry['bits']
    A = bits_to_adj(bits, n)
    scores = [int(sum(A[i])) for i in range(n)]
    predicted_sigma = sum(s*(s-1)//2 + (n-1-s)*(n-1-s-1)//2 for s in scores)
    if predicted_sigma == entry['total_sigma']:
        match += 1
    else:
        print(f"  MISMATCH: predicted={predicted_sigma}, actual={entry['total_sigma']}")
print(f"  Total sigma = sum C(outdeg,2)+C(indeg,2): {match}/1000")

# 5. The exact formula: total_sigma = sum C(s_i, 2) + C(n-1-s_i, 2)
#    = sum [s_i^2 - s_i + (n-1-s_i)^2 - (n-1-s_i)] / 2
#    = sum [2*s_i^2 - 2*(n-1)*s_i + (n-1)^2 - (n-1)] / 2
#    = [2*sum(s_i^2) - 2*(n-1)*n*(n-1)/2 + n*((n-1)^2 - (n-1))] / 2
#    Since sum(s_i) = n*(n-1)/2:
#    = sum(s_i^2) - n*(n-1)^2/2 + n*(n-1)(n-2)/2
#    = sum(s_i^2) - n*(n-1)/2
# So total_sigma = sum(s_i^2) - n*(n-1)/2

print("\n--- Compact Total Sigma Formula ---")
match = 0
for entry in all_data[:1000]:
    bits = entry['bits']
    A = bits_to_adj(bits, n)
    scores = [int(sum(A[i])) for i in range(n)]
    ssq = sum(s*s for s in scores)
    predicted = ssq - n*(n-1)//2
    if predicted == entry['total_sigma']:
        match += 1
    else:
        print(f"  MISMATCH: ssq={ssq}, predicted={predicted}, actual={entry['total_sigma']}")
print(f"  total_sigma = sum(s_i^2) - n(n-1)/2: {match}/1000")

# Corollary: total_delta = 21*5 - total_sigma - total_lambda
#           = 105 - [sum(s^2) - 21] - 3*c3
#           = 126 - sum(s^2) - 3*c3
print("\n--- Total Delta Formula ---")
match = 0
for entry in all_data[:1000]:
    bits = entry['bits']
    A = bits_to_adj(bits, n)
    scores = [int(sum(A[i])) for i in range(n)]
    ssq = sum(s*s for s in scores)
    predicted_delta = n*(n-1)*(n-2)//2 + n*(n-1)//2 - ssq - 3*entry['c3']
    # = n(n-1)(n-2)/2 + n(n-1)/2 - ssq - 3c3
    # = n(n-1)(n-1)/2 - ssq - 3c3
    # Wait let me recompute:
    # total = 21*5 = 105
    # total_sigma = ssq - 21
    # total_lambda = 3*c3
    # total_delta = 105 - (ssq - 21) - 3*c3 = 126 - ssq - 3*c3
    predicted_delta = 126 - ssq - 3*entry['c3']
    if predicted_delta == entry['total_delta']:
        match += 1
    else:
        print(f"  MISMATCH: predicted={predicted_delta}, actual={entry['total_delta']}")
print(f"  total_delta = 126 - sum(s_i^2) - 3*c3: {match}/1000")

# 6. The gradient: how total_sigma, total_lambda, total_delta
#    relate to c7 as a FUNCTION
print("\n\n--- c7 as Function of (total_sigma, total_lambda) ---")
sl_to_c7 = defaultdict(list)
for entry in all_data:
    sl_to_c7[(entry['total_sigma'], entry['total_lambda'])].append(entry['c7'])

print("(total_sigma, total_lambda) -> c7 stats:")
for (ts, tl) in sorted(sl_to_c7.keys()):
    c7s = sl_to_c7[(ts, tl)]
    td = n*(n-1)*(n-2)//2 - ts - tl  # actually 21*5 - ts - tl
    # Fix: total = 21*(n-2) = 21*5 = 105
    td = 105 - ts - tl
    if len(c7s) >= 3:
        print(f"  ({ts:>2},{tl:>2},{td:>2}): n={len(c7s):>4}, c7 range=[{min(c7s)},{max(c7s)}], "
              f"mean={np.mean(c7s):.1f}, std={np.std(c7s):.1f}")

# 7. The KEY: does (total_sigma, total_lambda) = (scores, c3)?
# total_sigma = sum(s^2) - 21 depends ONLY on score sequence
# total_lambda = 3*c3 depends ONLY on c3 count
# So the "total simplex coordinates" are just (score_variance, c3)!
print("\n\n--- Total Simplex = (Score Variance, c3) ---")
print("total_sigma = sum(s_i^2) - n(n-1)/2 = n*Var(scores) + n*(n-1)^2/4 - n(n-1)/2")
print("  where Var = sum(s_i - mean)^2 / n, mean = (n-1)/2")
print("total_lambda = 3 * c3")
print("total_delta = n(n-1)(n-2)/2 + n(n-1)/2 - sum(s^2) - 3*c3")
print()
print("So the 'total simplex position' is completely determined by (score_variance, c3)!")
print("The INDIVIDUAL pair positions carry additional structure beyond these totals.")

# 8. Fiber bundle: for FIXED lambda profile, what does sigma vary?
print("\n\n--- Lambda-Fiber Analysis ---")
# Group by lambda profile (multiset of lambda values across all 21 pairs)
lambda_profile_data = defaultdict(list)
for entry in all_data:
    lam_profile = tuple(sorted(l for s, l, d in entry['pairs']))
    lambda_profile_data[lam_profile].append(entry)

print(f"Distinct lambda profiles: {len(lambda_profile_data)}")
n_ambig_lam = 0
for lam_profile, entries in sorted(lambda_profile_data.items(), key=lambda x: -len(x[1])):
    c7s = set(e['c7'] for e in entries)
    sigma_profiles = set(tuple(sorted(s for s, l, d in e['pairs'])) for e in entries)
    if len(c7s) > 1:
        n_ambig_lam += 1
    if len(entries) >= 5:
        print(f"  lam_profile={lam_profile}: n={len(entries)}, "
              f"c7 range=[{min(c7s)},{max(c7s)}], "
              f"distinct sigma profiles={len(sigma_profiles)}")
        if len(entries) <= 20 and len(c7s) <= 3:
            break

print(f"\n  Lambda profiles ambiguous for c7: {n_ambig_lam}/{len(lambda_profile_data)}")

# 9. The question: does sigma profile determine c7 WITHIN a fixed lambda profile?
print("\n--- Sigma Determines c7 Within Lambda Fiber? ---")
ambig_within_fiber = 0
total_fibers_multi = 0
for lam_profile, entries in lambda_profile_data.items():
    if len(entries) < 2:
        continue
    # Within this lambda fiber, group by sigma profile
    sig_c7 = defaultdict(set)
    for e in entries:
        sig_profile = tuple(sorted(s for s, l, d in e['pairs']))
        sig_c7[sig_profile].add(e['c7'])

    for sig_p, c7s in sig_c7.items():
        total_fibers_multi += 1
        if len(c7s) > 1:
            ambig_within_fiber += 1

print(f"  Total (lambda, sigma) pairs: {total_fibers_multi}")
print(f"  Ambiguous for c7: {ambig_within_fiber}")
if ambig_within_fiber == 0:
    print("  SIGMA PROFILE DETERMINES c7 WITHIN EACH LAMBDA FIBER!")
else:
    print(f"  {ambig_within_fiber} cases where same (lambda, sigma) profiles give different c7")

# 10. The converse: does lambda determine c7 within fixed sigma?
print("\n--- Lambda Determines c7 Within Sigma Fiber? ---")
sigma_profile_data = defaultdict(list)
for entry in all_data:
    sig_profile = tuple(sorted(s for s, l, d in entry['pairs']))
    sigma_profile_data[sig_profile].append(entry)

ambig_within_sigma = 0
total_sigma_fibers = 0
for sig_profile, entries in sigma_profile_data.items():
    if len(entries) < 2:
        continue
    lam_c7 = defaultdict(set)
    for e in entries:
        lam_profile = tuple(sorted(l for s, l, d in e['pairs']))
        lam_c7[lam_profile].add(e['c7'])
    for lam_p, c7s in lam_c7.items():
        total_sigma_fibers += 1
        if len(c7s) > 1:
            ambig_within_sigma += 1

print(f"  Total (sigma, lambda) pairs: {total_sigma_fibers}")
print(f"  Ambiguous for c7: {ambig_within_sigma}")

# 11. What IS the hidden dimension?
# If the simplex profile almost determines c7 (only 6 ambig out of 248),
# what extra bit of information resolves the ambiguity?
print("\n\n--- Resolving the Ambiguity ---")
for i, (profile, c7s) in enumerate(sorted(ambig_profiles.items(), key=lambda x: -len(x[1]))):
    examples = profile_examples[profile]
    if len(examples) < 2:
        continue

    print(f"\n  Ambiguous profile {i+1}: c7 in {sorted(c7s)}")

    # Check: do the ambiguous cases differ in WHICH pairs have which simplex point?
    # i.e., the labeled (not just multiset) simplex assignment
    for c7_val in sorted(c7s):
        exs = [e for e in examples if e['c7'] == c7_val]
        # Check c5
        c5s = set(e['c5'] for e in exs)
        scores = set(e['scores'] for e in exs)
        print(f"    c7={c7_val}: {len(exs)} examples, c5={sorted(c5s)}, scores={scores}")

    # Does c5 resolve it?
    c5_c7 = defaultdict(set)
    for e in examples:
        c5_c7[e['c5']].add(e['c7'])
    ambig_c5 = sum(1 for c7s in c5_c7.values() if len(c7s) > 1)
    print(f"    c5 resolves ambiguity? {'YES' if ambig_c5 == 0 else 'NO (' + str(ambig_c5) + ' still ambig)'}")

print("\n\nDone.")
