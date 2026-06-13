"""
overlap_weight_210_structure.py -- kind-pasteur-2026-03-13-S61

The "2, 1, 0" dimensional structure in tournaments:

Each vertex pair (u,v) in a tournament has lambda(u,v) in {0, 1, ..., n-2}.
The lambda values create a FILTERED structure:
  lambda=0: "orthogonal" pairs (no common 3-cycle)
  lambda=1: "tangent" pairs (exactly one common 3-cycle)
  lambda=2: "overlapping" pairs (two common 3-cycles)
  ...

The Vitali atom operates on the lambda=2 boundary: the (1,1,2,2) reversal
preserves lambda but shifts sigma by +-1. This is a "tangent space" operation.

This script explores:
1. The "overlap weight" graph: edge-weighted graph with lambda values as weights
2. How the lambda-filtration decomposes into "dimensional layers"
3. The connection between lambda=0,1,2 pairs and the witness matrix structure
4. The "hidden higher-dimensional structure" = the fiber over each lambda value
5. How Vitali atoms act as "parallel transport" in this fiber bundle
"""

import numpy as np
from itertools import combinations
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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

# =====================================================
# PART 1: The Lambda Filtration at n=7
# =====================================================
n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"OVERLAP WEIGHT (LAMBDA) FILTRATION AT n={n}")
print("=" * 60)

np.random.seed(42)

# Sample tournaments and analyze lambda distribution
lambda_distributions = []
lambda_sigma_pairs = []

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    # Lambda distribution
    lam_vals = []
    sig_vals = []
    for u in range(n):
        for v in range(u+1, n):
            lam_vals.append(int(L[u][v]))
            sigma = n - 2 - int(A2[u][v]) - int(A2[v][u])
            sig_vals.append(sigma)

    lambda_distributions.append(tuple(sorted(lam_vals)))

    # Lambda-sigma pairs
    for i, (u, v) in enumerate([(u, v) for u in range(n) for v in range(u+1, n)]):
        lambda_sigma_pairs.append((lam_vals[i], sig_vals[i]))

# Lambda value distribution
print("\n--- Lambda Value Distribution ---")
lam_flat = [l for dist in lambda_distributions[:1000] for l in dist]
lam_dist = Counter(lam_flat)
total_pairs = len(lam_flat)
print(f"Total pairs sampled: {total_pairs}")
for k in sorted(lam_dist.keys()):
    pct = 100 * lam_dist[k] / total_pairs
    print(f"  lambda={k}: {lam_dist[k]} ({pct:.1f}%)")

# Lambda-sigma joint distribution
print("\n--- Lambda-Sigma Joint Distribution ---")
ls_dist = Counter(lambda_sigma_pairs[:21000])  # 1000 tournaments * 21 pairs
for lam in sorted(set(l for l, s in ls_dist.keys())):
    entries = {s: cnt for (l, s), cnt in ls_dist.items() if l == lam}
    print(f"  lambda={lam}: sigma distribution = {dict(sorted(entries.items()))}")

# KEY IDENTITY: lambda + sigma = ?
print("\n--- Lambda + Sigma Relationship ---")
ls_sum = Counter()
for l, s in lambda_sigma_pairs[:21000]:
    ls_sum[(l, s, l + s)] += 1

# What is lambda + sigma?
sum_dist = Counter(l + s for l, s in lambda_sigma_pairs[:21000])
print(f"  lambda + sigma distribution: {dict(sorted(sum_dist.items()))}")

# Is lambda + sigma = constant? NO (different values exist)
# What about the relationship?
# sigma(u,v) = #{common successors} + #{common predecessors}
# lambda(u,v) = #{3-cycles containing u,v}
# For a witness w of (u,v): w is in a 3-cycle with u,v iff (u->v->w->u or u->w->v->u)
# Common successor: u->w AND v->w. Common predecessor: w->u AND w->v.
# Neither succ nor pred: u->w,w->v or w->u,v->w (one in, one out).
# These "tangent" witnesses don't contribute to sigma but DO contribute to lambda!
#
# Actually: lambda(u,v) counts witnesses w where {u,v,w} forms a 3-cycle.
# A 3-cycle on {u,v,w}: either u->v->w->u or v->u->w->v.
# For u->v->w->u: w is a common predecessor (w->u, but v->w so w is predecessor of u, successor of v... wait)
# Let me be more precise:
# u->v->w->u means A[u][v]=1, A[v][w]=1, A[w][u]=1
# Then: v->w (v is predecessor of w), w->u (w is predecessor of u)
# For pair (u,v): w witnesses a 3-cycle.
# Is w a common successor? u->w? No, w->u. v->w? Yes.
# Is w a common predecessor? w->u? Yes. w->v? No, v->w.
# So w is NEITHER common successor nor common predecessor of (u,v)!
# But w IS in a 3-cycle with u,v.
#
# For v->u->w->v: A[v][u]=1, A[u][w]=1, A[w][v]=1
# u->w? Yes. v->w? No (w->v). So w is NOT common successor.
# w->u? No (u->w). w->v? Yes. So w is NOT common predecessor.
# Again: w witnesses a 3-cycle but is neither common successor nor predecessor!

# So lambda and sigma count DIFFERENT things:
# lambda(u,v) = #{w: w in 3-cycle with u,v} = #{w: "tangential" to (u,v)}
# sigma(u,v) = #{w: common successor or predecessor of (u,v)}
# The remaining witnesses: #{w: one-in-one-out for (u,v)} = n-2 - lambda - sigma ? NO...

# Let me verify with the 4-category decomposition:
# For each w != u,v, exactly one of:
# (a) u->w, v->w (common successor): contributes to sigma
# (b) w->u, w->v (common predecessor): contributes to sigma
# (c) u->w, w->v (u->w->v: transitive through w): contributes to lambda if also forms 3-cycle
#     Wait: u->w->v AND (v->u or u->v). If u->v: path u->w->v agrees with u->v, so
#     {u,v,w} has arcs u->v, u->w, w->v. Is there a 3-cycle? u->w->v->u? Only if v->u.
#     But we said u->v. Contradiction. So u->w,w->v,u->v => transitive triple, no 3-cycle.
#     If v->u: arcs v->u, u->w, w->v. That's a 3-cycle v->u->w->v!
# (d) w->u, v->w (v->w->u: transitive through w): similar analysis.

# So for category (c) u->w, w->v:
#   If u->v: transitive triple, no 3-cycle. Contributes to NEITHER lambda nor sigma.
#   If v->u: 3-cycle v->u->w->v. Contributes to lambda.
# For category (d) w->u, v->w:
#   If v->u: transitive triple v->w->u, v->u. No 3-cycle. NEITHER.
#   If u->v: 3-cycle u->v->w->u. Contributes to lambda.

# Summary:
# sigma = (a) + (b) = common successor + common predecessor
# lambda = #{(c) with v->u} + #{(d) with u->v}
# neither = #{(c) with u->v} + #{(d) with v->u}
# Total: (a) + (b) + (c) + (d) = n-2

# So: sigma + lambda + delta = n-2 where delta = "transitive non-witnesses"
# delta = #{w: {u,v,w} is a transitive triple}

print("\n--- Verifying: sigma + lambda + delta = n-2 ---")
verify_count = 0
verify_pass = 0
for trial in range(500):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])

            # Count transitive triples
            delta = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                # (c) u->w, w->v, u->v: transitive
                if A[u][w] and A[w][v] and A[u][v]:
                    delta += 1
                # (d) w->u, v->w, v->u: transitive
                if A[w][u] and A[v][w] and A[v][u]:
                    delta += 1

            verify_count += 1
            if sig + lam + delta == n - 2:
                verify_pass += 1
            else:
                print(f"  FAIL: u={u},v={v}, sig={sig}, lam={lam}, delta={delta}, sum={sig+lam+delta}")
                break
    else:
        continue
    break

print(f"  Verified: {verify_pass}/{verify_count}")
if verify_pass == verify_count:
    print("  CONFIRMED: sigma + lambda + delta = n-2 for all pairs!")
    print("  This is the FUNDAMENTAL DECOMPOSITION:")
    print("    n-2 = sigma + lambda + delta")
    print("    n-2 = (common succ/pred) + (3-cycles) + (transitive triples)")
    print()
    print("  This means the 'overlap weight' lambda and sigma are COMPLEMENTARY")
    print("  views of the same n-2 dimensional space of witnesses.")

# =====================================================
# PART 2: What Vitali Atoms Do to the Decomposition
# =====================================================
print(f"\n{'='*60}")
print("VITALI ATOM EFFECT ON (sigma, lambda, delta) DECOMPOSITION")
print(f"{'='*60}")

np.random.seed(42)
decomp_changes = []

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        A2_A = A @ A
        A2_B = B @ B

        for u in range(n):
            for v in range(u+1, n):
                sig_A = n - 2 - int(A2_A[u][v]) - int(A2_A[v][u])
                sig_B = n - 2 - int(A2_B[u][v]) - int(A2_B[v][u])
                lam = int(L[u][v])  # preserved!

                # Compute delta for A and B
                delta_A = 0
                delta_B = 0
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if A[u][w] and A[w][v] and A[u][v]:
                        delta_A += 1
                    if A[w][u] and A[v][w] and A[v][u]:
                        delta_A += 1
                    if B[u][w] and B[w][v] and B[u][v]:
                        delta_B += 1
                    if B[w][u] and B[v][w] and B[v][u]:
                        delta_B += 1

                u_in = u in S
                v_in = v in S
                ptype = "SS" if (u_in and v_in) else ("SX" if (u_in or v_in) else "XX")

                dsig = sig_B - sig_A
                ddelta = delta_B - delta_A

                if dsig != 0 or ddelta != 0:
                    decomp_changes.append({
                        'ptype': ptype,
                        'dsig': dsig,
                        'ddelta': ddelta,
                        'lam': lam,
                    })

        break

    if len(decomp_changes) > 3000:
        break

print(f"Collected {len(decomp_changes)} changes")

# Since lambda is preserved: dsig + ddelta = 0 !
# (sigma + lambda + delta = n-2, dlambda=0 => dsig = -ddelta)
print("\n--- Verifying dsig = -ddelta ---")
balance_check = all(c['dsig'] == -c['ddelta'] for c in decomp_changes)
print(f"  dsig = -ddelta always? {balance_check}")

if balance_check:
    print("  CONFIRMED: Vitali atoms TRANSFER weight between sigma and delta!")
    print("  When sigma goes up by 1, delta goes down by 1 (and vice versa).")
    print("  Lambda is the CONSERVED quantity; sigma and delta are conjugate variables.")

# Distribution by pair type
print("\n--- Change Distribution by Pair Type ---")
for ptype in ['XX', 'SX', 'SS']:
    entries = [c for c in decomp_changes if c['ptype'] == ptype]
    if entries:
        dsig_dist = Counter(c['dsig'] for c in entries)
        print(f"  {ptype}: dsig dist = {dict(sorted(dsig_dist.items()))}")
    else:
        print(f"  {ptype}: no changes (confirmed)")

# What lambda values do the changing pairs have?
print("\n--- Lambda Values of Changed Pairs ---")
for ptype in ['SX']:
    entries = [c for c in decomp_changes if c['ptype'] == ptype]
    lam_dist = Counter(c['lam'] for c in entries)
    print(f"  {ptype}: lambda distribution = {dict(sorted(lam_dist.items()))}")

    # dsig conditioned on lambda
    for lam in sorted(lam_dist.keys()):
        dsig_at_lam = Counter(c['dsig'] for c in entries if c['lam'] == lam)
        print(f"    lambda={lam}: dsig = {dict(sorted(dsig_at_lam.items()))}")

# =====================================================
# PART 3: The Fiber Bundle Interpretation
# =====================================================
print(f"\n{'='*60}")
print("FIBER BUNDLE INTERPRETATION")
print(f"{'='*60}")

print("""
The tournament's pair data has a FIBER BUNDLE structure:

  Base space: Lambda graph (edge-weighted graph with lambda values)
  Fiber over each pair (u,v): the split (sigma, delta) with sigma + delta = n-2-lambda(u,v)

  Vitali atoms act as PARALLEL TRANSPORT in this bundle:
    They move along the fiber (changing sigma <-> delta by +-1)
    while staying on the same base point (lambda preserved).

  The "hidden higher-dimensional structure" is precisely this fiber:
    - Lambda (the base) is 1-dimensional: just a non-negative integer per pair
    - The fiber is ALSO 1-dimensional: sigma ranges from 0 to n-2-lambda
    - Together: (lambda, sigma) is a 2D lattice point per pair
    - The Vitali atom shifts sigma by +-1 for SX pairs only

  The topology of the fiber bundle:
    - The fiber is discrete (integer sigma values)
    - The structure group is Z (translations by +-1)
    - Vitali atoms generate Z-actions on the fiber
    - The "parallel transport" around a closed loop of Vitali atoms
      may give a non-trivial holonomy (net sigma shift)
""")

# =====================================================
# PART 4: What does delta = transitive triples actually count?
# =====================================================
print(f"\n{'='*60}")
print("DELTA = TRANSITIVE TRIPLES: THE 'HIDDEN' COMPONENT")
print(f"{'='*60}")

# delta(u,v) = #{w: {u,v,w} is transitive}
# These are the witnesses that are "invisible" to both sigma and lambda.
# They represent the "dark matter" of tournament pair structure.

# In the sigma-algebra hierarchy:
#   Score -> Lambda -> (Lambda, Sigma) -> (Lambda, Sigma, Delta) = A^2
# But sigma + delta = n-2-lambda, so delta is determined by (lambda, sigma).
# Therefore (Lambda, Sigma) determines Delta, and we don't get new information.
# The hierarchy is: Score -> Lambda -> (Lambda, Sigma) -> A^2

# But delta has a nice interpretation: it counts "ordered transitive triples"
# through the pair (u,v). A transitive triple means u, v, and w have
# a total order consistent with the arc directions.

# Let's check: does the specific SET of transitive witnesses carry more info
# than just their count?
print("Does the SET of delta-witnesses (not just count) carry extra info?")

np.random.seed(42)
ls_with_delta_set = defaultdict(set)

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])

            # Identify each w's category
            categories = {}
            for w in range(n):
                if w == u or w == v:
                    continue
                uw = A[u][w]
                wv = A[w][v]
                wu = A[w][u]
                vw = A[v][w]
                uv = A[u][v]
                vu = A[v][u]

                if uw and vw:  # common successor
                    categories[w] = 'S'
                elif wu and wv:  # common predecessor
                    categories[w] = 'P'
                elif (uw and wv and uv) or (wu and vw and vu):
                    categories[w] = 'D'  # delta (transitive)
                else:
                    categories[w] = 'L'  # lambda (3-cycle witness)

            key = (lam, sig)
            cat_pattern = tuple(sorted(categories.values()))
            ls_with_delta_set[key].add(cat_pattern)

# How many distinct category patterns per (lambda, sigma)?
print("\nCategory patterns per (lambda, sigma):")
for key in sorted(ls_with_delta_set.keys()):
    patterns = ls_with_delta_set[key]
    print(f"  (lam={key[0]}, sig={key[1]}): {len(patterns)} distinct witness patterns")
    if len(patterns) <= 3:
        for p in patterns:
            print(f"    {p}")

# =====================================================
# PART 5: The 2-1-0 Structure
# =====================================================
print(f"\n{'='*60}")
print("THE 2-1-0 STRUCTURE: WINNER-DRAW-LOSER LAYERS")
print(f"{'='*60}")

# For a vertex v in a tournament on n vertices:
#   score(v) ranges from 0 to n-1
#   "regular" tournaments have all scores = (n-1)/2

# The 2-1-0 pattern appears in the (1,1,2,2) Vitali atom:
#   score 2: "winners" (dominate 2 out of 3 others)
#   score 1: "balanced" (beat 1, lose to 2)
#   score 0: doesn't appear (minimum score in (1,1,2,2) is 1)

# But more generally, for ANY pair (u,v):
#   Their shared witnesses w fall into categories:
#     w beats both u and v: w is a "2-winner" (dominates the pair)
#     w beats one, loses to other: w is a "1-mixer"
#     w loses to both: w is a "0-loser" (dominated by the pair)

# Let's count these categories
print("\nWitness w categories relative to pair (u,v):")
print("  2-winner: w beats both u and v (common predecessor of u,v)")
print("  0-loser: w loses to both u and v (common successor of u,v)")
print("  1-mixer: w beats one, loses to other")
print()

# For the 1-mixers: there are two subtypes
# 1a: w beats u, loses to v (w->u, v->w) => need to check u->v or v->u
# 1b: w beats v, loses to u (w->v, u->w) => need to check

# How do these relate to lambda and delta?
# If u->v:
#   1a (w->u, v->w): then v->w->u is a 3-cycle with v->u? No, we have u->v not v->u.
#                     So arcs: u->v, w->u, v->w. This is the 3-cycle u->v->w->u?
#                     v->w->u->v: A[v][w]=1, A[w][u]=1, A[u][v]=1. Yes! 3-cycle.
#                     So this contributes to lambda.
#   1b (u->w, w->v): arcs u->v, u->w, w->v. Transitive triple! Contributes to delta.
# If v->u:
#   1a (w->u, v->w): arcs v->u, w->u, v->w. Transitive (v->w->u and v->u). Delta.
#   1b (u->w, w->v): arcs v->u, u->w, w->v. 3-cycle v->u->w->v. Lambda.

# So: lambda = #{1-mixers where the "right" direction creates a 3-cycle}
#     delta = #{1-mixers where the direction creates a transitive triple}
#     sigma = #{2-winners} + #{0-losers}

# The sigma and (lambda+delta) are COMPLEMENTARY:
#   sigma counts 2-winners + 0-losers (extremal witnesses)
#   lambda + delta counts 1-mixers (balanced witnesses)
#   Total = n-2

# And among 1-mixers:
#   lambda vs delta depends on whether the mixer creates a CYCLE or TRANSITIVITY
#   with the pair's own arc direction.

print("FUNDAMENTAL THEOREM:")
print("  For pair (u,v) with arc u->v:")
print("  sigma(u,v) = #{w: w beats both u,v} + #{w: w loses to both u,v}")
print("  lambda(u,v) = #{w: (w->u, v->w) forming 3-cycle v->w->u->v}")
print("  delta(u,v) = #{w: (u->w, w->v) forming transitive u->w->v}")
print()
print("  The 3 categories (sigma, lambda, delta) = (extremal, cyclic, transitive)")
print("  partition all n-2 witnesses of each pair.")
print()
print("  Vitali atoms shift weight from 'transitive' to 'cyclic' (or vice versa)")
print("  for SX pairs, while preserving the extremal count... wait, no!")
print("  Since sigma changes by +-1 and lambda is preserved,")
print("  delta changes by -/+1. So the TRANSFER is between sigma and delta:")
print("  extremal <-> transitive, NOT cyclic!")
print()
print("  This means Vitali atoms convert TRANSITIVE witnesses into EXTREMAL witnesses")
print("  (or vice versa) for SX pairs, leaving CYCLIC witnesses unchanged.")

# Verify this interpretation with actual data
print("\n--- Verification: Which Category Shifts? ---")
np.random.seed(42)
category_shifts = []

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)

        for u in range(n):
            for v in range(u+1, n):
                u_in = u in S
                v_in = v in S
                if not (u_in != v_in):  # only SX pairs change
                    continue

                # Categorize each witness w in A and B
                for w in range(n):
                    if w == u or w == v:
                        continue

                    # In A
                    if A[u][w] and A[v][w]:
                        cat_A = 'S'  # common successor
                    elif A[w][u] and A[w][v]:
                        cat_A = 'P'  # common predecessor
                    elif (A[u][w] and A[w][v] and A[u][v]) or (A[w][u] and A[v][w] and A[v][u]):
                        cat_A = 'D'  # delta
                    else:
                        cat_A = 'L'  # lambda

                    # In B
                    if B[u][w] and B[v][w]:
                        cat_B = 'S'
                    elif B[w][u] and B[w][v]:
                        cat_B = 'P'
                    elif (B[u][w] and B[w][v] and B[u][v]) or (B[w][u] and B[v][w] and B[v][u]):
                        cat_B = 'D'
                    else:
                        cat_B = 'L'

                    if cat_A != cat_B:
                        category_shifts.append({
                            'from': cat_A,
                            'to': cat_B,
                            'w_in_S': w in S,
                            'u_in_S': u_in,
                        })

        break

    if len(category_shifts) > 5000:
        break

shift_dist = Counter((c['from'], c['to']) for c in category_shifts)
print(f"\nCategory transitions (from -> to):")
for (f, t), cnt in sorted(shift_dist.items()):
    print(f"  {f} -> {t}: {cnt}")

# Check if witness is in S or outside
print(f"\nCategory transitions by witness location:")
for w_in in [True, False]:
    label = "w in S" if w_in else "w outside S"
    entries = [c for c in category_shifts if c['w_in_S'] == w_in]
    if entries:
        sd = Counter((c['from'], c['to']) for c in entries)
        print(f"  {label}: {dict(sd)}")

print("\nDone.")
