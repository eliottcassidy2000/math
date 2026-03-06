"""
KEY DISCOVERY: M[a,b] = M[b,a] holds not just for tournaments (t_ij + t_ji = 1)
but for ALL "c-tournaments" where t_ij + t_ji = c for any constant c.

This script verifies this at n=3,4,5 symbolically and n=6,7 numerically.
"""

from itertools import permutations
from sympy import symbols, expand, Poly
import random

# ==========================================================
# Part 1: Symbolic verification at n=3
# ==========================================================
print("=" * 70)
print("c-TOURNAMENT SYMMETRY")
print("=" * 70)

def make_independent_vars(n):
    """Create independent arc variables for a digraph on n vertices."""
    t = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                t[(i,j)] = symbols(f't{i}{j}')
    return t

def hp_sym(t, vertex_set, start=None, end=None):
    """Symbolic Hamiltonian path count."""
    vl = sorted(vertex_set)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t[(perm[i], perm[i+1])]
        total += prod
    return expand(total)

def compute_M_entry(t, n, a, b):
    """Compute M[a,b] symbolically."""
    U = [v for v in range(n) if v != a and v != b]
    M_val = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1) ** len(S)
        ea = hp_sym(t, set(S)|{a}, end=a)
        bb = hp_sym(t, set(R)|{b}, start=b)
        M_val += sign * ea * bb
    return expand(M_val)

c = symbols('c')

for n in [3, 4]:
    print(f"\n--- n={n}: Symbolic c-tournament test ---")
    t = make_independent_vars(n)

    # Compute M[0,1] - M[1,0] with independent variables
    M01 = compute_M_entry(t, n, 0, 1)
    M10 = compute_M_entry(t, n, 1, 0)
    diff = expand(M01 - M10)
    print(f"  M[0,1] - M[1,0] (general digraph): {len(diff.as_ordered_terms()) if diff != 0 else 0} terms")

    # Substitute c-tournament constraint: t_ji = c - t_ij for i < j
    subs_c = {}
    for i in range(n):
        for j in range(i+1, n):
            subs_c[t[(j,i)]] = c - t[(i,j)]

    diff_c = expand(diff.subs(subs_c))
    print(f"  M[0,1] - M[1,0] (c-tournament): {diff_c}")
    print(f"  Zero for ALL c? {diff_c == 0}")

    if diff_c != 0:
        # Collect by powers of c
        forward_vars = [t[(i,j)] for i in range(n) for j in range(i+1, n)]
        p = Poly(diff_c, c, domain='ZZ[' + ','.join(str(v) for v in forward_vars) + ']')
        print(f"  Degree in c: {p.degree()}")
        for power in range(p.degree() + 1):
            coeff = expand(p.nth(power))
            print(f"    c^{power}: {coeff}")

    # Also check a different pair
    if n >= 4:
        M02 = compute_M_entry(t, n, 0, 2)
        M20 = compute_M_entry(t, n, 2, 0)
        diff2 = expand(M02 - M20)
        diff2_c = expand(diff2.subs(subs_c))
        print(f"  M[0,2] - M[2,0] (c-tournament): {diff2_c}")

# ==========================================================
# Part 2: Symbolic n=5 (may be slow but let's try)
# ==========================================================
print(f"\n--- n=5: Symbolic c-tournament test ---")
n = 5
t = make_independent_vars(n)
M01 = compute_M_entry(t, n, 0, 1)
M10 = compute_M_entry(t, n, 1, 0)
diff = expand(M01 - M10)
print(f"  M[0,1] - M[1,0] (general digraph): {len(diff.as_ordered_terms()) if diff != 0 else 0} terms")

subs_c = {}
for i in range(n):
    for j in range(i+1, n):
        subs_c[t[(j,i)]] = c - t[(i,j)]

diff_c = expand(diff.subs(subs_c))
print(f"  M[0,1] - M[1,0] (c-tournament): {diff_c}")
print(f"  Zero for ALL c? {diff_c == 0}")

# ==========================================================
# Part 3: Numerical verification of c-tournament symmetry
# ==========================================================
print(f"\n--- Part 3: Numerical c-tournament verification ---")
random.seed(42)

for n in [4, 5, 6, 7]:
    failures = 0
    for trial in range(200):
        c_val = random.uniform(-2, 5)  # Random c value (not just 0 or 1!)

        # Random c-tournament: t_ij random, t_ji = c - t_ij
        T = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                T[i][j] = random.uniform(-1, 2)
                T[j][i] = c_val - T[i][j]

        # Compute M[0,1] and M[1,0]
        def hp_num(vset, start=None, end=None):
            vl = sorted(vset)
            if len(vl) <= 1:
                if start is not None and (len(vl)==0 or vl[0] != start): return 0
                if end is not None and (len(vl)==0 or vl[0] != end): return 0
                return 1 if len(vl) == 1 else 0
            total = 0
            for p in permutations(vl):
                if start is not None and p[0] != start: continue
                if end is not None and p[-1] != end: continue
                prod = 1
                for i in range(len(p)-1):
                    prod *= T[p[i]][p[i+1]]
                total += prod
            return total

        a_t, b_t = 0, 1
        U_t = [v for v in range(n) if v != a_t and v != b_t]
        m_ab = 0
        m_ba = 0
        for mask in range(1 << len(U_t)):
            S = [U_t[i] for i in range(len(U_t)) if mask & (1 << i)]
            R = [U_t[i] for i in range(len(U_t)) if not (mask & (1 << i))]
            sign = (-1)**len(S)
            ea = hp_num(set(S)|{a_t}, end=a_t)
            bb = hp_num(set(R)|{b_t}, start=b_t)
            eb = hp_num(set(S)|{b_t}, end=b_t)
            ba = hp_num(set(R)|{a_t}, start=a_t)
            m_ab += sign * ea * bb
            m_ba += sign * eb * ba

        if abs(m_ab - m_ba) > 1e-8:
            failures += 1

    print(f"  n={n}: {failures}/200 random c-tournaments have M[0,1] != M[1,0]")

# ==========================================================
# Part 4: Does it hold when DIFFERENT pairs have DIFFERENT c?
# ==========================================================
print(f"\n--- Part 4: Non-uniform c (t_ij + t_ji = c_ij, different per pair) ---")
random.seed(42)

for n in [4, 5, 6]:
    failures = 0
    for trial in range(200):
        # Each pair {i,j} gets its OWN c_ij
        T = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                c_ij = random.uniform(-2, 5)  # DIFFERENT c for each pair
                T[i][j] = random.uniform(-1, 2)
                T[j][i] = c_ij - T[i][j]

        a_t, b_t = 0, 1
        U_t = [v for v in range(n) if v != a_t and v != b_t]
        m_ab = 0
        m_ba = 0
        for mask in range(1 << len(U_t)):
            S = [U_t[i] for i in range(len(U_t)) if mask & (1 << i)]
            R = [U_t[i] for i in range(len(U_t)) if not (mask & (1 << i))]
            sign = (-1)**len(S)
            ea = hp_num(set(S)|{a_t}, end=a_t)
            bb = hp_num(set(R)|{b_t}, start=b_t)
            eb = hp_num(set(S)|{b_t}, end=b_t)
            ba = hp_num(set(R)|{a_t}, start=a_t)
            m_ab += sign * ea * bb
            m_ba += sign * eb * ba

        if abs(m_ab - m_ba) > 1e-8:
            failures += 1

    print(f"  n={n}: {failures}/200 non-uniform c-digraphs have M[0,1] != M[1,0]")

# ==========================================================
# Part 5: What if c_ij only for pairs NOT involving a,b?
# ==========================================================
print(f"\n--- Part 5: Constraint only on U-pairs vs all pairs ---")
print("  Testing: what if we only constrain pairs {i,j} with both i,j in U (not a,b)?")
random.seed(42)

n = 4
a_t, b_t = 0, 1
U_t = [v for v in range(n) if v != a_t and v != b_t]
failures_u_only = 0
failures_ab_only = 0
failures_all_uniform = 0

for trial in range(500):
    c_val = random.uniform(-2, 5)
    T = [[0.0]*n for _ in range(n)]

    # Version 1: constrain ONLY U-pairs
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = random.uniform(-1, 2)
            if i in U_t and j in U_t:
                T[j][i] = c_val - T[i][j]
            else:
                T[j][i] = random.uniform(-1, 2)  # independent

    m_ab = 0
    m_ba = 0
    for mask in range(1 << len(U_t)):
        S = [U_t[ii] for ii in range(len(U_t)) if mask & (1 << ii)]
        R = [U_t[ii] for ii in range(len(U_t)) if not (mask & (1 << ii))]
        sign = (-1)**len(S)
        ea = hp_num(set(S)|{a_t}, end=a_t)
        bb = hp_num(set(R)|{b_t}, start=b_t)
        eb = hp_num(set(S)|{b_t}, end=b_t)
        ba = hp_num(set(R)|{a_t}, start=a_t)
        m_ab += sign * ea * bb
        m_ba += sign * eb * ba

    if abs(m_ab - m_ba) > 1e-8:
        failures_u_only += 1

print(f"  n=4, constrain only U-pairs: {failures_u_only}/500 failures")

# Version 2: constrain ONLY a-b adjacent pairs
for trial in range(500):
    c_val = random.uniform(-2, 5)
    T = [[0.0]*n for _ in range(n)]

    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = random.uniform(-1, 2)
            if i in [a_t, b_t] or j in [a_t, b_t]:
                T[j][i] = c_val - T[i][j]
            else:
                T[j][i] = random.uniform(-1, 2)

    m_ab = 0
    m_ba = 0
    for mask in range(1 << len(U_t)):
        S = [U_t[ii] for ii in range(len(U_t)) if mask & (1 << ii)]
        R = [U_t[ii] for ii in range(len(U_t)) if not (mask & (1 << ii))]
        sign = (-1)**len(S)
        ea = hp_num(set(S)|{a_t}, end=a_t)
        bb = hp_num(set(R)|{b_t}, start=b_t)
        eb = hp_num(set(S)|{b_t}, end=b_t)
        ba = hp_num(set(R)|{a_t}, start=a_t)
        m_ab += sign * ea * bb
        m_ba += sign * eb * ba

    if abs(m_ab - m_ba) > 1e-8:
        failures_ab_only += 1

print(f"  n=4, constrain only {a_t},{b_t}-adjacent pairs: {failures_ab_only}/500 failures")

# ==========================================================
# Part 6: The MINIMAL constraint for symmetry
# ==========================================================
print(f"\n--- Part 6: Which pairs need constraining? ---")
print("  At n=4, pairs are {0,2}, {0,3}, {1,2}, {1,3}, {2,3}")
print("  (arc {0,1} doesn't appear in M[0,1] at all)")
print()

# Test each subset of pairs
from itertools import combinations as combs
all_pairs = [(i,j) for i in range(n) for j in range(i+1, n) if not (i==a_t and j==b_t)]
print(f"  Relevant pairs: {all_pairs}")

random.seed(42)
for r in range(len(all_pairs)+1):
    for constrained_subset in combs(all_pairs, r):
        constrained_set = set(constrained_subset)
        n_fail = 0
        for trial in range(300):
            c_val = random.uniform(-2, 5)
            T = [[0.0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    T[i][j] = random.uniform(-1, 2)
                    if (i,j) in constrained_set:
                        T[j][i] = c_val - T[i][j]
                    else:
                        T[j][i] = random.uniform(-1, 2)

            m_ab = 0
            m_ba = 0
            for mask in range(1 << len(U_t)):
                S = [U_t[ii] for ii in range(len(U_t)) if mask & (1 << ii)]
                R2 = [U_t[ii] for ii in range(len(U_t)) if not (mask & (1 << ii))]
                sign = (-1)**len(S)
                ea = hp_num(set(S)|{a_t}, end=a_t)
                bb = hp_num(set(R2)|{b_t}, start=b_t)
                eb = hp_num(set(S)|{b_t}, end=b_t)
                ba = hp_num(set(R2)|{a_t}, start=a_t)
                m_ab += sign * ea * bb
                m_ba += sign * eb * ba

            if abs(m_ab - m_ba) > 1e-8:
                n_fail += 1

        if n_fail == 0:
            print(f"  WORKS with {r} constraints: {list(constrained_subset)}")

print(f"\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
The c-tournament symmetry theorem: M[a,b] = M[b,a] whenever
t_ij + t_ji = c for ALL pairs {i,j} (any constant c, same for all pairs).

This is a polynomial identity, not dependent on c=1 (tournaments).
The proof must use the UNIFORM constraint on arc pairs, but not its specific value.
""")
