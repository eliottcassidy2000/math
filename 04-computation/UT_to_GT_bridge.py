#!/usr/bin/env python3
"""
Bridge between U_T (Grinberg-Stanley symmetric function) and G_T(t,x) (inflated independence polynomial).

KEY QUESTION: Is G_T(t,x) a specialization of U_T? If so, what specialization?

G_T(t,x) = A_n(t) + sum_I x^|I| * A_{f+1}(t) * (t-1)^S
where f = n-1-S, S = sum(l_i - 1)

U_T = sum_{sigma in S(T,T^op)} (-1)^phi(sigma) * p_{type(sigma)}

The OCF specialization: p_1->1, p_{odd>=3}->2, p_even->0 gives H(T).
G_T(1,x) = n! for all x (the corrections vanish at t=1).
G_T(0,x) = I(Omega(T), x) [up to sign/scaling at t=0 which means u=infinity].

Actually: G_T(t,x)/t^m at t=1 means u=2, so P(2,x) = n! (corrections vanish via (u-2)^S/2).
And G_T(t,0) = A_n(t) for ALL tournaments (no cycle dependence).

So the specialization should:
- At x=0: give A_n(t) regardless of T (i.e., kill all nontrivial-cycle terms)
- At t=1, any x: give n! (universal)
- At t=1, x=2: give n! (consistent with H(T) coming from a different evaluation)

Wait. G_T(t,x) at t=0 is NOT I(Omega,x). Let me reclarify.

G_T(t,x) = t^m * P(u,x) where u = t + 1/t, m = (n-1)/2.
P(u,x) = P_n(u,0) + sum_I c_I * x^|I| * P_{n-S}(u,0) * (u-2)^{S/2}

At u=2 (t=1): P(2,x) = n! (all corrections vanish).
At x=0: P(u,0) = P_n(u,0) = A_n(t)/t^m (Eulerian polynomial in u-space).

The CORRECT evaluation for H(T) is G_T(0,2):
  G_T(0,x) = lim_{t->0} G_T(t,x) = lim_{t->0} [A_n(t) + sum_I x^|I| A_{f+1}(t) (t-1)^S]
  A_n(0) = A(n,0) = 1 (Eulerian number at k=0)
  A_{f+1}(0) = 1
  (0-1)^S = (-1)^S = 1 (since S is always even!)
  So G_T(0,x) = 1 + sum_I x^|I| * 1 * 1 = I(Omega(T), x)

YES! G_T(0,x) = I(Omega(T), x) exactly! And G_T(0,2) = I(Omega(T),2) = H(T).

So G_T gives a 2-parameter interpolation between:
- G_T(0,2) = H(T) [Hamiltonian path count]
- G_T(1,x) = n! [universal constant]
- G_T(t,0) = A_n(t) [Eulerian polynomial, tournament-independent]

Now: U_T under the specialization p_1->1, p_{odd}->2 gives H(T).
Can we find a 2-parameter specialization p_k -> f(k; t, x) such that
applying it to U_T gives G_T(t,x)?

Attempt: Since G_T(t,0) = A_n(t) = sum_sigma t^{des(sigma)} (for ALL sigma in S_n),
and U_T at x=0 should only include the identity permutation contribution...

Actually, no. U_T involves permutations with ALL nontrivial cycles being T or T^op cycles.
At "x=0" (killing nontrivial cycles), only the identity permutation contributes.
The identity has type (1^n) and phi=0, contributing p_1^n.
Under p_1 -> t (?), this would give t^n, not A_n(t).

That doesn't work. The relationship is more subtle.

Let me try a DIFFERENT approach: compute U_T's power-sum coefficients for small
tournaments, then check if G_T(t,x) can be recovered by a specific specialization.

kind-pasteur-2026-03-07-S30
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import factorial, comb

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def redei_berge_typed_coeffs(n, A):
    """Compute U_T's power-sum expansion: coefficients indexed by partition type.

    Returns dict: partition (tuple, sorted desc) -> integer coefficient.
    """
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    coeffs = defaultdict(int)

    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        # Check: each nontrivial cycle must be a T-cycle or T^op-cycle
        valid = True
        phi = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                valid = False
                break
            if is_T:
                phi += len(cyc) - 1

        if valid:
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            coeffs[cycle_type] += (-1)**phi

    return dict(coeffs)

def compute_GT(A, n, t_val, x_val):
    """Compute G_T(t, x) directly from formula."""
    # Find all directed odd cycles
    cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    cycles.append((canon, length))

    seen = set()
    unique_cycles = []
    for c, l in cycles:
        if c not in seen:
            seen.add(c)
            unique_cycles.append((set(c), l))

    # Build conflict graph and enumerate independent sets
    nc = len(unique_cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if unique_cycles[i][0] & unique_cycles[j][0]:
                adj[i][j] = adj[j][i] = True

    indep_sets = []
    for mask in range(2**nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            indep_sets.append(nodes)

    d = n - 1
    result = sum(eulerian_number(n, k) * t_val**k for k in range(n))

    for s in indep_sets:
        if not s:
            continue
        parts = len(s)
        sum_len = sum(unique_cycles[i][1] for i in s)
        S_val = sum_len - parts  # sum(l_i - 1)
        f = d - S_val

        A_f1 = sum(eulerian_number(f + 1, k) * t_val**k for k in range(f + 1))
        correction = x_val**parts * A_f1 * (t_val - 1)**S_val
        result += correction

    return result

def specialize_UT(coeffs, spec_func):
    """Apply a specialization function to U_T coefficients.

    spec_func(partition) -> value
    """
    total = 0
    for lam, coeff in coeffs.items():
        total += coeff * spec_func(lam)
    return total


# === Main computation ===
print("=" * 70)
print("BRIDGE: U_T <-> G_T(t,x)")
print("=" * 70)

# Test on C_5
n = 5
A = [[0]*n for _ in range(n)]
for i in range(n):
    A[i][(i+1) % n] = 1
    A[i][(i+2) % n] = 1

H = count_ham_paths(A, n)
print(f"\nC_5: H = {H}")

# Compute U_T power-sum coefficients
ut_coeffs = redei_berge_typed_coeffs(n, A)
print(f"U_T coefficients (partition -> coeff):")
for lam in sorted(ut_coeffs.keys()):
    if ut_coeffs[lam] != 0:
        print(f"  p_{lam} : {ut_coeffs[lam]}")

# Verify OCF specialization
def ocf_spec(lam):
    result = 1
    for part in lam:
        if part == 1:
            result *= 1
        elif part % 2 == 1:
            result *= 2
        else:
            result *= 0  # even parts -> 0
    return result

H_from_UT = specialize_UT(ut_coeffs, ocf_spec)
print(f"\nOCF specialization of U_T: {H_from_UT}")
print(f"H(T) = {H}, match: {H_from_UT == H}")

# Compute G_T at several (t, x) values
print(f"\nG_T values:")
for t_val in [0, 1, 2]:
    for x_val in [0, 1, 2]:
        gt = compute_GT(A, n, t_val, x_val)
        print(f"  G_T({t_val}, {x_val}) = {gt}")

# Now the key question: what specialization p_lambda -> f(lambda; t, x)
# gives G_T(t, x)?
#
# We need: for each partition lambda appearing in U_T,
# sum_lambda ut_coeffs[lambda] * f(lambda; t, x) = G_T(t, x)
#
# The partitions appearing are of n=5: (5,), (3,1,1), (1,1,1,1,1)
# Let's denote: a = f((5,)), b = f((3,1,1)), c = f((1,1,1,1,1))
#
# Then: ut[(5,)]*a + ut[(3,1,1)]*b + ut[(1,1,1,1,1)]*c = G_T(t,x)
#
# We have 3 unknowns {a, b, c} that are functions of (t,x).
# But we also need these to be "symmetric-function-natural" —
# i.e., f(lambda) = prod_i g(lambda_i) for some function g on parts.

print(f"\n{'='*70}")
print(f"SPECIALIZATION SEARCH")
print(f"{'='*70}")

# For C_5, the nonzero partitions are (5,), (3,1,1), (1,1,1,1,1)
# And possibly (3,3) if T^op cycles contribute at n=5... wait, n=5 so max 5.
# Actually: permutations can have cycle type (3,1,1) if cycle length 3.

# Let's see what partition types appear:
print(f"\nNonzero partition types in U_T(C_5):")
for lam, c in sorted(ut_coeffs.items()):
    if c != 0:
        # Check: is this all-odd parts?
        all_odd = all(p % 2 == 1 for p in lam)
        print(f"  {lam}: coeff = {c}, all-odd: {all_odd}")

# For a multiplicative specialization p_lambda -> prod p_{lambda_i}:
# We need p_1 and p_3 and p_5 values (as functions of t, x).
# G_T(t,x) for C_5:
# = A_5(t) + 2*x*(t-1)^2*A_3(t) + 2*x*(t-1)^4
# where 2 counts the cycle types (5 3-cycles or 2 5-cycles)
# Actually let me be precise about the cycle counts

from sympy import Symbol, expand, collect

t = Symbol('t')
x = Symbol('x')

# C_5 has: 5 directed 3-cycles, 2 directed 5-cycles (both directions of the full cycle)
# Wait, C_5 as cyclic tournament: edges are i->(i+1) and i->(i+2).
# 3-cycles: any 3 consecutive vertices form a cycle. There are 5 vertex-sets,
# each supporting exactly 1 directed 3-cycle in C_5 (the one going with the flow).
# Actually let me count properly.

# For C_5 with edges i->(i+1), i->(i+2) mod 5:
# 3-cycle on {0,1,2}: need 0->1->2->0 or 0->2->1->0
#   0->1 (yes), 1->2 (yes), 2->0 (yes, diff=3=2 mod 5? No, diff=-2 mod 5 = 3, and QR={1,2})
#   Actually edge set: 0->1, 0->2, 1->2, 1->3, 2->3, 2->4, 3->4, 3->0, 4->0, 4->1
#   So 0->1->2->0: 0->1 yes, 1->2 yes, 2->0 yes? We need A[2][0]=1.
#   A[2][0]: diff = 0-2 = -2 mod 5 = 3. Is 3 in {1,2}? No. So 0->2 exists but 2->0 doesn't.
#   0->2->1->0: A[0][2]=1, A[2][1]? diff=1-2=-1=4, not in {1,2}. No.
#   Hmm, so {0,1,2} may not support a 3-cycle in C_5.

# Let me just count directly
cycles_c5 = []
for length in [3, 5]:
    for subset in combinations(range(n), length):
        for perm in permutations(subset):
            is_cycle = True
            for i in range(length):
                if A[perm[i]][perm[(i+1)%length]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                min_idx = list(perm).index(min(perm))
                canon = perm[min_idx:] + perm[:min_idx]
                if canon not in [c[0] for c in cycles_c5]:
                    cycles_c5.append((canon, length))

c3 = sum(1 for c, l in cycles_c5 if l == 3)
c5 = sum(1 for c, l in cycles_c5 if l == 5)
print(f"\nC_5 cycle counts: c3={c3}, c5={c5}")

# Now compute G_T(t,x) symbolically
A_5t = sum(eulerian_number(5, k) * t**k for k in range(5))  # = 1 + 26t + 66t^2 + 26t^3 + t^4
A_3t = sum(eulerian_number(3, k) * t**k for k in range(3))  # = 1 + 4t + t^2
A_1t = 1  # A_1(t) = 1

GT_sym = A_5t
# 3-cycles: each has S = l-1 = 2, f = 4-2 = 2
# 5-cycles: each has S = l-1 = 4, f = 4-4 = 0
# Need to add corrections for ALL independent sets
# alpha_0 = 1 (counted in A_5t)
# alpha_1 contributions: c3 three-cycles with S=2, c5 five-cycles with S=4
# alpha_2 contributions: disjoint pairs (only 3-cycles can be disjoint at n=5)

# Actually, I should use the independence polynomial structure
# For C_5: alpha_0=1, alpha_1 = c3+c5 = ?, alpha_2 = #{disjoint pairs}

print(f"C_5 cycles:")
for c, l in cycles_c5:
    print(f"  length {l}: {c}")

# Check which pairs are disjoint
disjoint_pairs = 0
for i in range(len(cycles_c5)):
    for j in range(i+1, len(cycles_c5)):
        if not (set(cycles_c5[i][0]) & set(cycles_c5[j][0])):
            disjoint_pairs += 1
print(f"Disjoint pairs: {disjoint_pairs}")

# Build G_T(t,x) for C_5 symbolically
print(f"\n=== Symbolic G_T(t,x) for C_5 ===")
GT_sym = expand(A_5t)
for s_nodes in [[i] for i in range(len(cycles_c5))]:
    parts = 1
    sum_len = cycles_c5[s_nodes[0]][1]
    S_val = sum_len - parts
    f = 4 - S_val
    A_f1 = sum(eulerian_number(f + 1, k) * t**k for k in range(f + 1))
    GT_sym += x * A_f1 * (t - 1)**S_val

# Add disjoint pair contributions
for i in range(len(cycles_c5)):
    for j in range(i+1, len(cycles_c5)):
        if not (set(cycles_c5[i][0]) & set(cycles_c5[j][0])):
            parts = 2
            sum_len = cycles_c5[i][1] + cycles_c5[j][1]
            S_val = sum_len - parts
            f = 4 - S_val
            A_f1 = sum(eulerian_number(f + 1, k) * t**k for k in range(max(0, f + 1)))
            GT_sym += x**2 * A_f1 * (t - 1)**S_val

GT_sym = expand(GT_sym)
print(f"G_T(t,x) = {GT_sym}")

# Now check: can we express this as a specialization of U_T?
# U_T has coefficients on partitions. The product specialization
# p_lambda -> prod g(lambda_i) means:
#
# sum coeff[lambda] * prod g(lambda_i) = G_T(t,x)
#
# With partitions (1^5), (3,1^2), (5,):
# coeff[(1,1,1,1,1)] * g(1)^5 + coeff[(3,1,1)] * g(3)*g(1)^2 + coeff[(5,)] * g(5) = G_T
#
# Let me write this as a system and solve for g(1), g(3), g(5) in terms of t, x.

p_11111 = ut_coeffs.get((1,1,1,1,1), 0)
p_311 = ut_coeffs.get((3,1,1), 0)
p_5 = ut_coeffs.get((5,), 0)

print(f"\nU_T coefficients: p_(1^5) = {p_11111}, p_(3,1,1) = {p_311}, p_(5,) = {p_5}")

# Check: at the OCF specialization g(1)=1, g(3)=2, g(5)=2:
H_check = p_11111 * 1**5 + p_311 * 2 * 1**2 + p_5 * 2
print(f"OCF check: {p_11111}*1 + {p_311}*2 + {p_5}*2 = {H_check} (should be {H})")

# G_T(t,x) at (0,2): should be H
gt_02 = GT_sym.subs([(t, 0), (x, 2)])
print(f"G_T(0,2) = {gt_02} (should be {H})")

# G_T(1,x): should be 120 = 5!
gt_1x = GT_sym.subs(t, 1)
print(f"G_T(1,x) = {expand(gt_1x)} (should be 120)")

# For the multiplicative specialization, we need:
# p_11111 * g1^5 + p_311 * g3 * g1^2 + p_5 * g5 = G_T(t,x)
# This gives ONE equation in 3 unknowns. But for DIFFERENT tournaments,
# the coefficients change while g1, g3, g5 stay the same.
# So we need: for ALL tournaments T on 5 vertices,
#   ut_T[(1^5)] * g1^5 + ut_T[(3,1,1)] * g3*g1^2 + ut_T[(5,)] * g5 = GT_T(t,x)

# Let's compute for the transitive tournament too
print(f"\n{'='*70}")
print(f"TRANSITIVE T_5")
print(f"{'='*70}")

A_trans = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        A_trans[i][j] = 1

H_trans = count_ham_paths(A_trans, n)
print(f"H(T_trans) = {H_trans}")

ut_trans = redei_berge_typed_coeffs(n, A_trans)
print(f"U_T coefficients:")
for lam in sorted(ut_trans.keys()):
    if ut_trans[lam] != 0:
        print(f"  p_{lam} : {ut_trans[lam]}")

# Now we have TWO equations:
# For C_5: p_11111_c * g1^5 + p_311_c * g3*g1^2 + p_5_c * g5 = GT_c5(t,x)
# For Trans: p_11111_t * g1^5 + p_311_t * g3*g1^2 + p_5_t * g5 = GT_trans(t,x)

# If transitive has no cycles, then p_311_t = p_5_t = 0, and:
# p_11111_t * g1^5 = GT_trans(t, x) = A_5(t) (since no cycles)
# So g1 = (A_5(t) / p_11111_t)^{1/5}

# Check transitive cycle counts
cycles_trans = []
for length in [3, 5]:
    for subset in combinations(range(n), length):
        for perm in permutations(subset):
            is_cycle = True
            for i in range(length):
                if A_trans[perm[i]][perm[(i+1)%length]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                min_idx = list(perm).index(min(perm))
                canon = perm[min_idx:] + perm[:min_idx]
                if canon not in [c[0] for c in cycles_trans]:
                    cycles_trans.append((canon, length))

print(f"\nTransitive T_5 cycles: {len(cycles_trans)}")

# Transitive tournament has NO directed cycles at all!
# So U_T for transitive should only have (1^5) type.
# And GT_trans(t,x) = A_5(t) for all x.

# From transitive: coeff * g1^5 = A_5(t)
# So if coeff = 1 (identity permutation), g1^5 = A_5(t).
# A_5(t) = 1 + 26t + 66t^2 + 26t^3 + t^4
# g1 = A_5(t)^{1/5} which is NOT a polynomial!

# Hmm, that means the specialization is NOT multiplicative.
# p_lambda -> f(lambda, t, x) where f is NOT a product over parts.

# Alternative: maybe the relationship is more subtle.
# U_T is a symmetric function = sum c_lambda p_lambda.
# We can also write U_T in terms of other bases: e_lambda, h_lambda, s_lambda, etc.
# The quasisymmetric expansion gives:
# U_T = sum_{sigma} (-1)^phi * F_{Des(sigma), n}
# where F_{I,n} is the fundamental quasisymmetric function.
# Under the principal specialization ps(x_i = q^{i-1}):
#   ps(F_{I,n}) = q^{...} / (q-1)...(q^n-1)
# which gives the q-analog.

# At q=0: ps(F_{I,n})(0) = [I = emptyset] (only the empty descent set survives)
# So ps(U_T)(0) = sum_{sigma with Des(sigma)=empty} (-1)^phi
# Des(sigma) = empty means sigma is the identity (only perm with no descents)
# So ps(U_T)(0) = 1.

# This isn't giving us G_T(t,x) directly.

# KEY INSIGHT: G_T(t,x) is NOT a standard specialization of U_T.
# It's constructed differently. G_T uses:
# - The Eulerian structure of permutations (descents)
# - The cycle collection structure from the independence polynomial
# - The (t-1)^S factor which comes from expanding (t-1+1)^d in the correction term

# The connection is: BOTH U_T and G_T encode the same cycle information,
# but in different ways. U_T uses the symmetric function framework (power sums),
# while G_T uses the descent/Eulerian framework (inflated Eulerian polynomials).

# The SHARED data between them is the independence polynomial coefficients alpha_k.
# Both can recover I(Omega(T), x) and hence H(T) from their respective frameworks.

print(f"\n{'='*70}")
print(f"CONCLUSION: G_T is NOT a standard specialization of U_T")
print(f"{'='*70}")
print(f"U_T uses power-sum basis of Sym (tracks cycle types as symmetric functions)")
print(f"G_T uses Eulerian polynomial basis (tracks descents + cycle collections)")
print(f"Both encode the same cycle-disjoint collection data (independence polynomial)")
print(f"The bridge: alpha_S coefficients appear in BOTH frameworks with the same values")
print(f"U_T: alpha_S enters via p_lambda coefficients counting sigma with type lambda")
print(f"G_T: alpha_S enters via correction terms P_{{n-S}}(u,0) * (u-2)^{{S/2}}")
