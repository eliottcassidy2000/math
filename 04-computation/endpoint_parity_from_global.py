#!/usr/bin/env python3
"""
BRIDGING THE GAP: From Global Even Cycle Vanishing to Endpoint Even-r-Powers

The central question: Can we derive M[a,b] has only even r-powers from
the GLOBAL fact that U_T uses only odd-part cycle types?

Key insight: The Grinberg-Stanley formula gives
  ham(D) = sum_S det(A_bar[S]) * per(A[S^c])

Our transfer matrix is the ENDPOINT-CONDITIONED version:
  M[a,b] = sum_S (-1)^|S| E_a(S+a) * B_b(R+b)

where E_a(S+a) = sum of Ham paths through S+a ending at a
      B_b(R+b) = sum of Ham paths through R+b starting at b

APPROACH: Express M[a,b] directly in terms of the CYCLE STRUCTURE
of permutations sigma in S_V(T, T^op), then check if the
even-cycle-vanishing involution works at the endpoint level.

kind-pasteur-2026-03-06-S23
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, Rational
from collections import defaultdict

def setup(n):
    """Setup symbolic framework with r = c/2 and skew variables."""
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')

    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]

    def t(i, j):
        """Arc weight: r + s_ij."""
        if i == j: return 0
        return r + s(i, j)

    return r, sv, s, t

def hp(t_fn, vset, start=None, end=None):
    """Hamiltonian path sum through vertex set."""
    vl = sorted(vset)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for p in permutations(vl):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= t_fn(p[i], p[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, n, a, b):
    """Transfer matrix entry M[a,b]."""
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        ea = hp(t_fn, set(S)|{a}, end=a)
        bb = hp(t_fn, set(R)|{b}, start=b)
        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("ENDPOINT PARITY FROM GLOBAL EVEN CYCLE VANISHING")
print("=" * 70)

# ============================================================
# Part 1: Express M[a,b] as sum over TWO-PATH COVERS
# ============================================================
print("\n--- Part 1: M[a,b] as sum over 2-path covers ---")
print("""
A '2-path cover of (a,b)' is a partition of [n] into:
  Path_a: Hamiltonian path through S+a ending at a
  Path_b: Hamiltonian path through R+b starting at b
where S union R = U = [n] \ {a,b}.

M[a,b] = sum over all such covers of (-1)^|S| * product of arc weights

Each 2-path cover is a PERMUTATION sigma of [n] where:
  sigma restricted to S+a is a path ending at a (with sign (-1)^|S|)
  sigma restricted to R+b is a path starting at b
""")

# ============================================================
# Part 2: Verify that the r-power structure of M matches
# the global parity prediction
# ============================================================
print("\n--- Part 2: r-power structure of M[a,b] ---")

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1

    M_ab = transfer_M(t, n, a, b)
    p = Poly(M_ab, r)

    print(f"\n  n={n}: M[0,1] as polynomial in r=c/2:")
    for power in range(p.degree() + 1):
        coeff = expand(p.nth(power))
        parity = "EVEN" if power % 2 == 0 else "ODD"
        is_zero = " = 0!" if coeff == 0 else f" ({len(coeff.as_ordered_terms())} terms)"
        print(f"    r^{power} [{parity}]: {is_zero}")

# ============================================================
# Part 3: The KEY question — what is the s-degree structure?
# ============================================================
print("\n--- Part 3: s-degree structure of M[a,b] ---")
print("At each r-power, what are the s-degrees of the surviving terms?")

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1

    M_ab = transfer_M(t, n, a, b)
    p = Poly(M_ab, r)

    print(f"\n  n={n}:")
    s_vars = list(sv.values())
    for power in range(p.degree() + 1):
        coeff = expand(p.nth(power))
        if coeff == 0:
            print(f"    r^{power}: ZERO")
            continue

        # Get degree in s-variables
        if s_vars:
            try:
                ps = Poly(coeff, *s_vars)
                terms_by_degree = defaultdict(int)
                for monom, c in ps.as_dict().items():
                    deg = sum(monom)
                    terms_by_degree[deg] += 1
                deg_str = ", ".join(f"deg{d}:{count}" for d, count in sorted(terms_by_degree.items()))
                print(f"    r^{power}: s-degrees [{deg_str}]")
            except Exception as e:
                print(f"    r^{power}: {len(coeff.as_ordered_terms())} terms (poly err)")
        else:
            print(f"    r^{power}: {coeff}")

# ============================================================
# Part 4: Check if the involution works ENDPOINT-SPECIFICALLY
# ============================================================
print("\n" + "=" * 70)
print("Part 4: Does cycle reversal involution preserve endpoints?")
print("=" * 70)
print("""
The global involution: reverse an even cycle c in sigma.
This changes phi by k-1 (odd for even k), flipping the sign.

For our endpoint M[a,b]: we have TWO paths, one ending at a, one starting at b.
If we reverse an even cycle, does this preserve the endpoint structure?

Key: the paths in a 2-path cover are NOT cycles.
The cycle reversal involution acts on the CYCLE TYPE of a permutation.
But our 2-path cover is not a permutation — it's two paths.

However, if we CONCATENATE the two paths (a-path and b-path) into a
single permutation, we get a permutation sigma of [n] with specific
cycle structure. The cycle reversal can then act on sigma.

Let's check: if sigma = (path ending at a) * (path starting at b),
what is the cycle structure of sigma?
""")

# At n=4, a=0, b=1: S can be empty, {2}, {3}, or {2,3}
n = 4; r_sym, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3]

print(f"\n  n=4, a=0, b=1:")
for mask in range(1 << len(U)):
    S_lst = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_lst = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_lst)

    print(f"\n    S={S_lst}, R={R_lst}, sign={sign:+d}")

    # List all 2-path covers for this S
    S_set = set(S_lst) | {a}
    R_set = set(R_lst) | {b}

    for p_a in permutations(sorted(S_set)):
        if p_a[-1] != a: continue
        for p_b in permutations(sorted(R_set)):
            if p_b[0] != b: continue

            # Build the permutation sigma from both paths
            # Path ending at a: p_a[0] -> p_a[1] -> ... -> a
            # Path starting at b: b -> p_b[1] -> ... -> p_b[-1]
            # As a permutation: sigma(p_a[i]) = p_a[i+1], sigma(a) = ???
            # sigma(p_b[i]) = p_b[i+1], sigma(p_b[-1]) = ???
            # The "missing" links are a -> ??? and ??? -> p_b[-1]
            # Natural completion: make it cyclic: sigma(a) = p_a[0], sigma(p_b[-1]) = b
            # No... this gives cycles, not a single permutation of interest

            # Actually, the PRODUCT of arc weights is what matters
            weight_a = 1
            for i in range(len(p_a)-1):
                weight_a *= t(p_a[i], p_a[i+1])
            weight_b = 1
            for i in range(len(p_b)-1):
                weight_b *= t(p_b[i], p_b[i+1])

            total_weight = expand(sign * weight_a * weight_b)

            # Expand in r
            pr = Poly(total_weight, r_sym)
            r1_coeff = expand(pr.nth(1)) if pr.degree() >= 1 else 0

            if r1_coeff != 0:
                print(f"      Path_a={p_a}, Path_b={p_b}: r^1 coeff = {r1_coeff}")

# ============================================================
# Part 5: Direct analysis — WHY does r^1 vanish?
# ============================================================
print("\n" + "=" * 70)
print("Part 5: Direct analysis of r^1 vanishing at n=4")
print("=" * 70)

n = 4; r_sym, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3]

# Collect ALL r^1 contributions
r1_terms = []
for mask in range(1 << len(U)):
    S_lst = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_lst = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_lst)

    S_set = set(S_lst) | {a}
    R_set = set(R_lst) | {b}

    for p_a in permutations(sorted(S_set)):
        if p_a[-1] != a: continue
        for p_b in permutations(sorted(R_set)):
            if p_b[0] != b: continue

            # Total number of arcs
            n_arcs = (len(p_a)-1) + (len(p_b)-1)

            # Each arc weight is r + s_ij. Product of n_arcs factors.
            # The r^1 coefficient: choose exactly 1 arc to contribute r,
            # the rest contribute s_ij.

            # r^1 = sum over each arc position of:
            #   r * product of s_ij for all other arcs

            arcs_a = [(p_a[i], p_a[i+1]) for i in range(len(p_a)-1)]
            arcs_b = [(p_b[i], p_b[i+1]) for i in range(len(p_b)-1)]
            all_arcs = arcs_a + arcs_b

            if len(all_arcs) == 0: continue

            # Product of all s values
            s_product = 1
            for arc in all_arcs:
                s_product *= s(arc[0], arc[1])
            s_product = expand(s_product)

            # r^1 coefficient: sum over positions, replacing s_ij by 1
            for k in range(len(all_arcs)):
                remaining = [s(arc[0], arc[1]) for j, arc in enumerate(all_arcs) if j != k]
                term = sign
                for factor in remaining:
                    term *= factor
                term = expand(term)
                r1_terms.append((S_lst, p_a, p_b, k, term))

# Now sum all r^1 terms
total_r1 = sum(t[-1] for t in r1_terms)
total_r1 = expand(total_r1)
print(f"\n  Total r^1 coefficient: {total_r1}")
print(f"  Number of individual terms: {len(r1_terms)}")

# Group by which arc is replaced by r
print("\n  Grouped by arc replaced:")
by_arc = defaultdict(lambda: 0)
for S_lst, p_a, p_b, k, term in r1_terms:
    arcs_a = [(p_a[i], p_a[i+1]) for i in range(len(p_a)-1)]
    arcs_b = [(p_b[i], p_b[i+1]) for i in range(len(p_b)-1)]
    all_arcs = arcs_a + arcs_b
    arc = all_arcs[k]
    by_arc[arc] += term

for arc in sorted(by_arc.keys()):
    val = expand(by_arc[arc])
    if val != 0:
        print(f"    Arc {arc}: {val}")

# ============================================================
# Part 6: The skew-symmetry cancellation mechanism
# ============================================================
print("\n" + "=" * 70)
print("Part 6: Skew-symmetry and the cancellation mechanism")
print("=" * 70)
print("""
At r^1: each contributing term has exactly one factor of r and the
rest are s_ij values. Since each s_ij = -s_ji, the product of
s-values has a PARITY (odd or even number of "reversed" arcs).

The key question: does the (-1)^|S| sign from the inclusion-exclusion
interact with the skew signs to cancel all r^1 terms?

If arc (i,j) contributes r, the remaining arcs contribute s-values.
The total s-degree of each such term is n-3 (= number of arcs minus 1).
For n=4: s-degree 1. For n=5: s-degree 2.

At n=4: r^1 has s-degree 1, and s-degree 1 monomials are just s_ij.
Total is 0, meaning the coefficient of each s_ij in the r^1 sum is 0.
""")

# Verify: decompose r^1 by s-monomial at n=4
print("\n  n=4: r^1 by s-monomial:")
for key in sorted(sv.keys()):
    var = sv[key]
    # Coefficient of this s-variable in total_r1
    if total_r1 == 0:
        print(f"    s{key}: (total is 0)")
        break

# Now do n=5 to see the pattern
print("\n  n=5: Computing r^1 coefficient...")
n = 5; r_sym, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3, 4]

M_ab = transfer_M(t, n, a, b)
p5 = Poly(M_ab, r_sym)

for power in range(p5.degree() + 1):
    coeff = expand(p5.nth(power))
    parity = "EVEN" if power % 2 == 0 else "ODD"
    is_zero = " = 0!" if coeff == 0 else ""
    print(f"    r^{power} [{parity}]: {'ZERO' if coeff == 0 else f'{len(coeff.as_ordered_terms())} terms'}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
THE BRIDGE between global and endpoint parity:

1. GLOBAL: U_T has only odd-part cycle types (Even Cycle Vanishing).
   Proved by cycle reversal involution.

2. ENDPOINT: M[a,b] has only even r-powers (our conjecture).
   This is an endpoint-conditioned version of the same parity.

3. CONNECTION: Both stem from the T <-> T^op involution.
   - Global: reversing an even cycle in sigma changes sign by (-1)^{k-1}
   - Endpoint: the (-1)^|S| signs filter out odd-parity s-terms

4. THE GAP: The global proof uses a clean involution on permutations.
   The endpoint version needs an involution on 2-PATH-COVERS.

   A 2-path-cover can be encoded as a permutation sigma of {1,...,n}
   by composing the two paths. The cycle structure of sigma encodes
   both the path structure AND the subset S.

   QUESTION: Does the cycle reversal involution on sigma preserve
   the endpoint structure (i.e., keep a as the end of one path
   and b as the start of the other)?

   If YES, the even-r-powers property follows immediately from
   the even cycle vanishing theorem.
""")
