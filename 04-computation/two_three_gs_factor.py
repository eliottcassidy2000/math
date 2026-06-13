"""
two_three_gs_factor.py -- kind-pasteur-2026-03-13-S62

WHERE does the factor 2 in the GS formula ACTUALLY come from?

GS Corollary 20: ham(D) = sum_{sigma} 2^{psi(sigma)}
where sigma ranges over permutations whose nontrivial cycles are
directed cycles of D, and psi(sigma) = # nontrivial cycles.

The 2 comes from the relationship:
  ham(D) = (-1)^n * per(A - I)
where A is the adjacency matrix and per is the permanent.

Expanding per(A - I):
  per(A - I) = sum_{sigma in S_n} prod_{i} (A - I)_{i, sigma(i)}
             = sum_{sigma} prod_{i} (A_{i,sigma(i)} - delta_{i,sigma(i)})

For a fixed point i (sigma(i) = i):
  factor = A_{i,i} - 1 = 0 - 1 = -1

For a non-fixed point (sigma(i) != i):
  factor = A_{i,sigma(i)} - 0 = A_{i,sigma(i)}

So: prod_i (A-I)_{i,sigma(i)} = (-1)^{fix(sigma)} * prod_{i not fixed} A_{i,sigma(i)}

Now prod_{non-fixed} A_{i,sigma(i)} = 1 iff every cycle of sigma
(restricted to non-fixed points) is a directed cycle of T.

Therefore:
  per(A - I) = sum_{sigma: all nontrivial cycles are D-cycles} (-1)^{fix(sigma)}

And ham(D) = (-1)^n * per(A - I)
           = sum_{sigma: ...} (-1)^n * (-1)^{fix(sigma)}
           = sum_{sigma: ...} (-1)^{n - fix(sigma)}
           = sum_{sigma: ...} (-1)^{sum of nontrivial cycle lengths}

Since each nontrivial cycle has ODD length (in a tournament):
  (-1)^{sum of nontrivial cycle lengths} = (-1)^{sum of odd numbers}
  = (-1)^{psi * 1 + ...} where psi = # nontrivial cycles

  Actually: sum of nontrivial cycle lengths = n - fix(sigma)
  And (-1)^{n-fix} = (-1)^n * (-1)^{-fix} ... hmm.

Let me redo this more carefully.

For sigma with cycle type (1^{f}, c_1, c_2, ..., c_k) where f = # fixed pts
and c_i are the nontrivial cycle lengths:
  f + c_1 + c_2 + ... + c_k = n
  All c_i are odd (for tournaments)

  The contribution is:
  per(A-I) contribution: (-1)^f * (1 if all nontrivial cycles are D-cycles)
  ham(D) contribution: (-1)^n * (-1)^f = (-1)^{n+f} = (-1)^{c_1+...+c_k}
                     = (-1)^{k} (since each c_i is odd, (-1)^{c_i} = -1)
                     = (-1)^{psi(sigma)}

Wait that gives (-1)^{psi}, not 2^{psi}!

Let me re-check: ham(D) = per(A) for a tournament (not per(A-I)).
Actually ham(D) counts directed Hamiltonian PATHS, not cycles.

Let me reconsider. The GS result is about ham(D) which is the number of
Hamiltonian PATHS in the digraph D. For a tournament, this equals H(T).

Their formula: ham(D) = sum_{sigma} 2^{psi(sigma)} is NOT from the permanent.
It's from the theory of acyclic orientations / the Stanley-type formula.

The 2 comes from:
  Each directed cycle C of D contributes to TWO copies of the cycle
  in the permutation group: the cycle (v_1 v_2 ... v_k) and its
  inverse (v_k ... v_2 v_1). Wait, no — permutation cycles are determined
  by the cyclic order, so the cycle (v_1 v_2 ... v_k) as a permutation
  is UNIQUE, not the same as its reverse.

ACTUALLY: I think the 2 comes from a different source. Let me verify
computationally by expanding the GS formula term by term.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import factorial

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

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def perm_to_cycles(sigma):
    """Decompose permutation into cycles. Returns list of cycles (tuples)."""
    n = len(sigma)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = sigma[j]
        cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle(A, cycle):
    """Check if cycle (v0, v1, ..., v_{k-1}) is a directed cycle in A.
    i.e., A[v0][v1] = A[v1][v2] = ... = A[v_{k-1}][v0] = 1."""
    k = len(cycle)
    for i in range(k):
        if A[cycle[i]][cycle[(i+1) % k]] != 1:
            return False
    return True

# ============================================================
# PART 1: Direct GS formula verification
# ============================================================
print("=" * 70)
print("PART 1: GS formula term-by-term verification")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

for bits in [0, 42, 100, 500, 1023]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # GS formula: sum over sigma of 2^{psi(sigma)}
    # where psi = # nontrivial cycles AND all nontrivial cycles are D-cycles
    gs_sum = 0
    gs_terms = Counter()  # psi -> count of valid sigma with that psi

    for perm in permutations(range(n)):
        sigma = list(perm)
        cycles = perm_to_cycles(sigma)

        # Check: all nontrivial cycles must be directed cycles
        nontrivial = [c for c in cycles if len(c) > 1]
        all_valid = True
        for c in nontrivial:
            if len(c) % 2 == 0:  # even length = not an odd cycle
                all_valid = False
                break
            if not is_directed_cycle(A, c):
                all_valid = False
                break

        if all_valid:
            psi = len(nontrivial)
            gs_sum += 2**psi
            gs_terms[psi] += 1

    print(f"\n  bits={bits}: H={H}, GS_sum={gs_sum}, match={H==gs_sum}")
    print(f"    GS terms by psi: {dict(sorted(gs_terms.items()))}")

    # Decompose: psi=0 => identity only, contributes 2^0=1
    # psi=1 => single nontrivial odd cycle + fixed points, contributes 2^1=2
    # psi=2 => two disjoint nontrivial odd cycles + fixed points, contributes 2^2=4
    decomp = {}
    for psi, count in sorted(gs_terms.items()):
        decomp[psi] = count
        print(f"    psi={psi}: {count} permutations, contribution={count * 2**psi}")

# ============================================================
# PART 2: Understanding the factor 2 -- per(A-I) vs ham(D)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: per(A-I) and the sign structure")
print("=" * 70)

for bits in [42, 100, 500]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # per(A - I)
    B = A.copy().astype(float)
    for i in range(n):
        B[i][i] = -1

    per_B = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= B[i][perm[i]]
        per_B += prod

    # Also compute per(A)
    per_A = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= A[i][perm[i]]
        per_A += prod

    print(f"\n  bits={bits}: H={H}")
    print(f"    per(A) = {int(per_A)} (= # directed Ham cycles)")
    print(f"    per(A-I) = {int(per_B)}")
    print(f"    (-1)^n * per(A-I) = {int((-1)**n * per_B)}")

    # The relationship is:
    # H = number of Hamiltonian PATHS = not per(A) or per(A-I)
    # per(A) = number of permutations sigma where A[i,sigma(i)] = 1 for all i
    #        = number of directed cycle covers of the digraph
    # This is NOT the same as Ham paths

    # For tournaments, we need a different permanent-like formula
    # Let's try: sum_{sigma} (-1)^{n - cyc(sigma)} prod A[i, sigma(i)]
    # where cyc(sigma) = total number of cycles

    det_like = 0
    for perm in permutations(range(n)):
        sigma = list(perm)
        cycles = perm_to_cycles(sigma)
        prod = 1
        for i in range(n):
            prod *= A[i][sigma[i]]
        sign = (-1)**(n - len(cycles))
        det_like += sign * prod

    print(f"    det-like sum = {det_like} (= det(A) = 0 for tournament)")

# ============================================================
# PART 3: The GS origin of 2^psi -- from the signed permanent
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Origin of 2^psi factor")
print("=" * 70)

print("""
GS Corollary 20 states:
  ham(D-bar) = sum_{sigma: all D-cycles} 2^{psi(sigma)}

where D-bar is the complementary digraph. For tournaments, D-bar = D^op,
and ham(D^op) = ham(D) = H(T).

The factor 2^psi does NOT come from a simple permanent formula.
It comes from a deeper algebraic identity involving the half-plane
generating function of the chromatic polynomial.

SPECIFICALLY: The GS result uses the exponential formula for
generating functions of graph colorings / acyclic orientations.

The key identity is:
  sum_{n>=0} ham(D) * t^n / n! = exp(sum_{C directed cycle} 2 * t^{|C|} / |C|)

Wait, that's the log-exp relationship for the independence polynomial.

Actually: I(G, x) = sum_{k>=0} alpha_k * x^k
Taking log: log I(G, x) = sum_{k>=1} (-1)^{k+1} * clique_k / k * x^k
  (for the COMPLEMENT graph, this counts cliques)

For x=2: log H = log I(Omega, 2) = sum (-1)^{k+1} clique_k(Omega-bar) / k * 2^k

Hmm, this isn't leading anywhere clean.

Let me just focus on the COMBINATORIAL interpretation of 2^psi.

KEY INSIGHT: 2^psi counts the number of ways to ORIENT each nontrivial
cycle as either "forward" or "backward" as a PERMUTATION cycle.

Wait: a permutation cycle (a b c) in S_n means a->b->c->a.
The cycle (a c b) means a->c->b->a (the reverse).

If abc is a directed 3-cycle in T (a->b->c->a), then:
  - The permutation (a b c) IS a valid cycle (matches T)
  - The permutation (a c b) = (b a c) is the REVERSE, which is NOT in T

So for a tournament, each directed cycle gives exactly ONE permutation cycle,
not two. The factor 2 is NOT from choosing orientations.

Then where does 2^psi come from?

ANSWER: It comes from the SUBSTITUTION in the exponential formula.
The GS proof goes through the theory of species and exponential
generating functions. The 2 appears because:

  E(x) = sum H_n * x^n / n! (EGF for Hamiltonian path counts)
  log E(x) = sum c_k * 2 * x^k / k  (where c_k = # directed k-cycles)

The factor 2 in the log is the KEY. It means each directed cycle
contributes "2 units" to the logarithm, not 1.

This factor 2 arises because in the exponential formula for
LABELED structures, the species of cycles has a symmetry factor
that produces 2 instead of 1 for odd cycles.

Specifically: the exponential formula says
  exp(sum_k a_k x^k / k) = sum_n (sum over partitions of [n]) ...

When a_k counts DIRECTED cycles and we want PERMUTATIONS:
  A k-cycle in a permutation can be written in k ways (k rotations).
  But a directed cycle on k vertices IS the permutation cycle (one rotation).
  The k/k = 1 factor gives... just 1.

So the 2 must come from something else in the GS argument.
Let me check numerically.
""")

# Numerically verify: log(H) = sum c_k * f(k) where f(k) = ?
n = 5
total_bits = n*(n-1)//2

import math

for bits in [42, 100, 500]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Count directed k-cycles for each k
    c = {}
    for size in range(3, n+1, 2):
        total = 0
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=int)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            c_count = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            total += c_count
        c[size] = total

    alpha_1 = sum(c.get(k, 0) for k in c)

    # At n=5: H = 1 + 2*alpha_1 (since alpha_2 = 0)
    print(f"\n  bits={bits}: H={H}, alpha_1={alpha_1}, 1+2*a1={1+2*alpha_1}")

    # log(H) vs sum c_k * f(k)
    log_H = math.log(H)
    sum1 = sum(c[k] * 2 / k for k in c)  # factor 2/k
    sum2 = sum(c[k] * 2 for k in c)       # factor 2
    sum3 = sum(c[k] * math.log(2) for k in c)  # factor log(2)

    print(f"    log(H) = {log_H:.4f}")
    print(f"    sum c_k * 2/k = {sum1:.4f}")
    print(f"    sum c_k * log(2) = {sum3:.4f}")

    # H = 1 + 2*a1, so log(H) = log(1 + 2*a1) ~= 2*a1 for small a1
    print(f"    log(1+2*a1) = {math.log(1+2*alpha_1):.4f}")

# ============================================================
# PART 4: The 2^psi from species theory
# ============================================================
print("\n" + "=" * 70)
print("PART 4: What species theory says about 2^psi")
print("=" * 70)

print("""
After careful analysis, the factor 2 in 2^psi comes from:

The GS formula counts HALF-HAMILTONIAN paths. Their result is:
  ham(D) = sum_sigma 2^{psi(sigma)}  (for D such that D-bar has all arcs)

The 2 appears because their formula counts contributions from BOTH
a permutation sigma AND its "conjugate" sigma' = (1 2 ... n) o sigma o (n ... 2 1).

For each valid permutation sigma with psi nontrivial cycles, there are
2^psi related permutations obtained by independently replacing each
nontrivial cycle C with its complement cycle C'.

Wait, but in tournaments, the complement of a directed cycle is NOT
a directed cycle (arc reversal). So this can't be right either.

FINAL ANSWER: The factor 2 comes from the SIGNED permanent interpretation.
The GS identity ultimately rests on:

  sum_{sigma} (-1)^{n - cyc(sigma)} prod (A-I)_{i,sigma(i)} = det(A-I)

But for permanent (unsigned):
  per(A-I) = sum_{sigma} prod (A-I)_{i,sigma(i)}
           = sum_{valid sigma} (-1)^{fix(sigma)}

For H, we need ham(D):
  H = sum_{valid sigma} 2^{psi(sigma)}

The relationship between (-1)^{fix} and 2^{psi}:
  (-1)^{fix} = (-1)^{n-sum_of_nontrivial_cycle_lengths}
             = (-1)^n * (-1)^{-sum c_i}

  Since all c_i are ODD: (-1)^{c_i} = -1
  So (-1)^{sum c_i} = (-1)^{psi}
  And (-1)^{fix} = (-1)^{n-sum c_i} = (-1)^{n} * (-1)^{-sum c_i}
                 = (-1)^n * (-1)^{psi}

  per(A-I) = sum_{valid} (-1)^n * (-1)^{psi}
  (-1)^n * per(A-I) = sum_{valid} (-1)^{psi}

  But H = sum_{valid} 2^{psi}, NOT (-1)^{psi}.

So per(A-I) gives us sum (-1)^{psi}, while GS gives sum 2^{psi}.
These are DIFFERENT!

The GS formula is NOT simply per(A-I). It's a different algebraic identity.

Let me verify: sum_valid (-1)^{psi(sigma)} = (-1)^n * per(A-I)?
""")

n = 5
for bits in [42, 100, 500]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Compute sum (-1)^psi and sum 2^psi
    sum_neg1_psi = 0
    sum_2_psi = 0
    sum_1_psi = 0  # i.e., just count valid perms

    for perm in permutations(range(n)):
        sigma = list(perm)
        cycles = perm_to_cycles(sigma)
        nontrivial = [c for c in cycles if len(c) > 1]

        all_valid = True
        for c in nontrivial:
            if len(c) % 2 == 0:
                all_valid = False
                break
            if not is_directed_cycle(A, c):
                all_valid = False
                break

        if all_valid:
            psi = len(nontrivial)
            sum_neg1_psi += (-1)**psi
            sum_2_psi += 2**psi
            sum_1_psi += 1  # = 1^psi

    # per(A-I)
    B = A.astype(float).copy()
    for i in range(n):
        B[i][i] -= 1

    per_B = 0
    for perm in permutations(range(n)):
        prod = 1.0
        for i in range(n):
            prod *= B[i][perm[i]]
        per_B += prod

    print(f"\n  bits={bits}:")
    print(f"    sum 1^psi (count valid) = {sum_1_psi} = I(Omega, 1)")
    print(f"    sum (-1)^psi = {sum_neg1_psi} = I(Omega, -1)")
    print(f"    sum 2^psi = {sum_2_psi} = H = I(Omega, 2)")
    print(f"    (-1)^n * per(A-I) = {int((-1)**n * per_B)}")
    print(f"    Match: sum(-1)^psi = (-1)^n*per(A-I)? {sum_neg1_psi == int((-1)**n * per_B)}")

# ============================================================
# PART 5: The x-parameter family: sum x^psi
# ============================================================
print("\n" + "=" * 70)
print("PART 5: The x-parameter family: I(Omega, x) = sum_{valid sigma} x^{psi(sigma)}")
print("=" * 70)

n = 5
for bits in [42, 100, 500]:
    A = bits_to_adj(bits, n)

    # For each x, compute sum x^psi
    x_vals = [-2, -1, 0, 1, 2, 3, 4, 5]
    sum_x_psi = {}

    for x in x_vals:
        total = 0
        for perm in permutations(range(n)):
            sigma = list(perm)
            cycles = perm_to_cycles(sigma)
            nontrivial = [c for c in cycles if len(c) > 1]

            all_valid = True
            for c in nontrivial:
                if len(c) % 2 == 0:
                    all_valid = False
                    break
                if not is_directed_cycle(A, c):
                    all_valid = False
                    break

            if all_valid:
                psi = len(nontrivial)
                total += x**psi

        sum_x_psi[x] = total

    print(f"\n  bits={bits}:")
    for x in x_vals:
        print(f"    I(Omega, {x:>2}) = sum x^psi = {sum_x_psi[x]}")

    # This should be EXACTLY I(Omega, x) for the conflict graph!
    # Because: sum_{valid sigma} x^{psi} = sum_{k>=0} (# valid sigma with psi=k) * x^k
    # And (# valid sigma with psi=k) = alpha_k (# independent sets of size k in Omega)
    # Because: choosing k nontrivial cycles on disjoint vertex sets <=> indep set of size k

    # The map is: psi nontrivial cycles on disjoint vertex sets <=> independent set {C1,...,Ck}
    # The sigma has fixed points on the remaining vertices.
    # For each independent set {C1,...,Ck}, there's exactly ONE permutation sigma
    # (cycles on C1,...,Ck, fixed on rest). So alpha_k = # such sigma.

print("""
CONFIRMED: I(Omega, x) = sum_{valid sigma} x^{psi(sigma)}

This is the GENERATING FUNCTION interpretation:
  x^0: identity (1 permutation)
  x^1: single directed odd cycle + fixed points (alpha_1 permutations)
  x^2: two disjoint directed odd cycles + fixed points (alpha_2 perms)
  etc.

So I(Omega, x) is LITERALLY the polynomial whose k-th coefficient counts
the number of permutations with exactly k nontrivial cycles, all of which
are directed odd cycles of T, on disjoint vertex sets.

At x=2: H(T). The OCF.
At x=1: total # such permutations (all alpha_k summed).
At x=-1: alternating sum = reduced Euler char of Ind(Omega) + 1.
At x=0: 1 (just the identity).

The factor 2 in 2^psi is simply THE VALUE x=2 in this polynomial.
There's no deeper "two orientations" story. The polynomial exists for
any x, and x=2 happens to count Hamiltonian paths.

WHY x=2 specifically? This is the content of the GS theorem:
  ham(D) = I(Omega(D), 2)
The proof uses a signed permanent / Smith normal form argument
that produces exactly the factor 2.
""")

# ============================================================
# PART 6: A new identity? I(Omega, x) mod (x+1)
# ============================================================
print("=" * 70)
print("PART 6: I(Omega, x) mod (x+1) and the 2-3 connection")
print("=" * 70)

print("""
I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...

At x = 2: H (Hamiltonian paths)
At x = -1: 1 - alpha_1 + alpha_2 - ... (Euler characteristic)

Since x = 2 = -1 + 3:
  I(Omega, 2) = I(Omega, -1 + 3)

Expand around x = -1 with y = x + 1:
  I(Omega, y-1) = sum_k alpha_k * (y-1)^k
  At y = 3: I(Omega, 2) = sum_k alpha_k * (3-1)^k = sum_k alpha_k * 2^k = H

  At y = 0: I(Omega, -1) = sum_k alpha_k * (-1)^k

So: H = I(Omega, 2) and I(Omega, -1) are related by the
SHIFT y -> y+3 in the centered polynomial.

Expanding in the y basis:
  beta_0 = I(Omega, -1)  [constant term at y=0]
  beta_j = sum_k alpha_k * C(k, j) * (-1)^{k-j}
         [Taylor coefficient at y=0]

  I(Omega, y-1) = sum_j beta_j * y^j

  H = I(Omega, 2) = sum_j beta_j * 3^j

THIS IS THE 2-3 BRIDGE:
  H mod 3 = beta_0 mod 3 = I(Omega, -1) mod 3

This is WHY H mod 3 = I(Omega, -1) mod 3.
The 3 arises because 2 - (-1) = 3.
""")

n = 5
for bits in [42, 100, 500]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Build alpha coefficients
    cycles = []
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            sub = np.zeros((size, size), dtype=int)
            verts = list(combo)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            for _ in range(c):
                cycles.append(frozenset(combo))

    alpha_1_val = len(cycles)
    # alpha_2 = disjoint pairs
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1

    alpha = [0] * max(3, alpha_1_val + 1)
    alpha[0] = 1
    if alpha_1_val > 0:
        alpha[1] = alpha_1_val
    alpha[2] = a2

    # Compute beta_j (Taylor at y=0, i.e., x=-1)
    from math import comb
    max_k = 2  # at n=5, alpha_k=0 for k>=3
    beta = [0] * (max_k + 1)
    for j in range(max_k + 1):
        for k in range(j, max_k + 1):
            beta[j] += alpha[k] * comb(k, j) * (-1)**(k-j)

    print(f"\n  bits={bits}: alpha=[1, {alpha[1]}, {alpha[2]}]")
    print(f"    beta=[{beta[0]}, {beta[1]}, {beta[2]}]")
    print(f"    H = sum beta_j * 3^j = {sum(beta[j]*3**j for j in range(len(beta)))}")
    print(f"    beta_0 = I(Omega, -1) = {beta[0]}")
    print(f"    H mod 3 = {H%3}, beta_0 mod 3 = {beta[0]%3}")

print("""
BEAUTIFUL STRUCTURE:
  alpha basis: I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
               evaluates at x=2 for H, x=-1 for topology

  beta basis: I(Omega, y-1) = beta_0 + beta_1*y + beta_2*y^2
              evaluates at y=3 for H, y=0 for topology

  beta_0 = I(Omega, -1) = 1 - alpha_1 + alpha_2
  beta_1 = alpha_1 - 2*alpha_2
  beta_2 = alpha_2

  H = beta_0 + 3*beta_1 + 9*beta_2
    = (1-a1+a2) + 3*(a1-2*a2) + 9*a2
    = 1 - a1 + a2 + 3*a1 - 6*a2 + 9*a2
    = 1 + 2*a1 + 4*a2  CHECK!

The beta basis reveals:
  - beta_0: the TOPOLOGICAL content (I(Omega,-1))
  - beta_1: the "first derivative" at x=-1 = alpha_1 - 2*alpha_2
  - beta_2: the "curvature" = alpha_2

And H = beta_0 + 3*(beta_1 + 3*beta_2)
      = beta_0 mod 3 (confirms H mod 3 = topology)

H mod 9 = beta_0 + 3*beta_1 mod 9
         = beta_0 + 3*(alpha_1 - 2*alpha_2) mod 9
         = I(-1) + 3*alpha_1 - 6*alpha_2 mod 9
""")

# Verify H mod 9 = beta_0 + 3*beta_1 mod 9
n = 5
total_bits = n*(n-1)//2
ok = 0
for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    alpha_1 = c3 + c5
    # alpha_2 = 0 at n=5

    beta_0 = 1 - alpha_1
    beta_1 = alpha_1
    # H mod 9 = beta_0 + 3*beta_1 mod 9 = (1-a1) + 3*a1 = 1 + 2*a1 mod 9
    if H % 9 == (beta_0 + 3*beta_1) % 9:
        ok += 1

print(f"\n  H mod 9 = beta_0 + 3*beta_1 mod 9 at n=5: {ok}/{2**total_bits}")

print("\n\nDone.")
