#!/usr/bin/env python3
"""
Deeper isomorphism analysis of Cayley tournaments on F_21 = Z/7 x| Z/3.

From frobenius21_selfconverse.py:
- ALL 1024 connection sets give the same canonical hash
- ALL have exactly 1155 three-cycles
- ALL pass the hash match T vs T^op (self-converse test)

QUESTION: Are all 1024 Cayley tournaments on F_21 isomorphic?
If so, there's essentially ONE tournament on F_21.

We also check: is the hash match sufficient for self-converse?
We do this by finding an explicit anti-automorphism.

kind-pasteur-2026-03-06-S25e
"""

# Group operations for F_21
POW2 = [1, 2, 4]  # 2^j mod 7

def mul(a, b):
    return ((a[0] + POW2[a[1]] * b[0]) % 7, (a[1] + b[1]) % 3)

def inv(a):
    j_inv = (-a[1]) % 3
    i_inv = (-POW2[j_inv] * a[0]) % 7
    return (i_inv, j_inv)

elements = [(i, j) for j in range(3) for i in range(7)]
idx = {e: k for k, e in enumerate(elements)}
n = 21

# Enumerate all inverse pairs
nonid = [e for e in elements if e != (0, 0)]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        pairs.append((g, inv(g)))
        seen.add(g)
        seen.add(inv(g))

def adj_matrix(S_set):
    A = [[0]*n for _ in range(n)]
    for a in elements:
        for b in elements:
            if a == b:
                continue
            diff = mul(inv(a), b)
            if diff in S_set:
                A[idx[a]][idx[b]] = 1
    return A

def compute_invariant(A):
    """Compute a more detailed invariant vector for each vertex."""
    # For each vertex i, compute:
    # - out-degree
    # - number of 3-cycles through i
    # - for each out-neighbor j: out-degree of j restricted to N^+(i)
    n = len(A)
    vertex_profiles = []
    for i in range(n):
        out_deg = sum(A[i])
        in_deg = sum(A[j][i] for j in range(n))

        # 3-cycles through i: i->j->k->i
        tc = 0
        for j in range(n):
            if A[i][j] != 1:
                continue
            for k in range(n):
                if A[j][k] == 1 and A[k][i] == 1:
                    tc += 1

        # Edges within out-neighborhood
        out_set = [j for j in range(n) if A[i][j] == 1]
        edges_out = sum(1 for j in out_set for k in out_set if A[j][k] == 1)

        # Edges within in-neighborhood
        in_set = [j for j in range(n) if A[j][i] == 1]
        edges_in = sum(1 for j in in_set for k in in_set if A[j][k] == 1)

        vertex_profiles.append((out_deg, tc, edges_out, edges_in))

    return tuple(sorted(vertex_profiles))

# Test: are all 1024 really the same?
print("=" * 70)
print("ARE ALL F_21 CAYLEY TOURNAMENTS ISOMORPHIC?")
print("=" * 70)

invariants = set()
for bits in range(1024):
    S = set()
    for k, (g, gi) in enumerate(pairs):
        if (bits >> k) & 1:
            S.add(g)
        else:
            S.add(gi)
    A = adj_matrix(S)
    inv_vec = compute_invariant(A)
    invariants.add(inv_vec)

print(f"  Distinct invariant vectors across all 1024 sets: {len(invariants)}")
if len(invariants) == 1:
    print("  => All 1024 have IDENTICAL vertex profiles!")
    print("  => Extremely strong evidence that they're all isomorphic.")
else:
    print(f"  => {len(invariants)} distinct isomorphism classes (at least)")

# For a specific non-normal S, find an explicit anti-automorphism
print("\n" + "=" * 70)
print("FINDING EXPLICIT ANTI-AUTOMORPHISM FOR NON-NORMAL S")
print("=" * 70)

# Choose a specific non-normal S
S = set()
for k, (g, gi) in enumerate(pairs):
    if k < 5:
        S.add(g)
    else:
        S.add(gi)

# Check normality
def is_normal_set(S_set):
    conj_classes_map = {}
    for g in elements:
        cls = frozenset(mul(mul(h, g), inv(h)) for h in elements)
        if g in S_set:
            if cls in conj_classes_map:
                if conj_classes_map[cls] != 'in':
                    return False
            conj_classes_map[cls] = 'in'
        else:
            if g != (0, 0):
                if cls in conj_classes_map:
                    if conj_classes_map[cls] != 'out':
                        return False
                conj_classes_map[cls] = 'out'
    return True

normal = is_normal_set(S)
print(f"  S = {sorted(S)}")
print(f"  Normal (union of conjugacy classes)? {normal}")

A = adj_matrix(S)
AT = [[A[j][i] for j in range(n)] for i in range(n)]

# Try to find anti-automorphism by searching group-based maps
# An anti-auto sigma satisfies: A[sigma(i)][sigma(j)] = A[j][i]
# Since T is a Cayley tournament, try sigma = L_g (left translation by g)
# composed with inversion and/or group automorphisms.

# Group automorphisms of F_21 = Z/7 x| Z/3:
# Aut(F_21): the automorphism group of the Frobenius group.
# Inner automorphisms: conjugation by elements of F_21 (|Inn| = |F_21/Z(F_21)| = 21 since Z = {e})
# Outer automorphisms: automorphisms of Z/7 that commute with the Z/3 action.
# The Z/3 action is a -> 2a mod 7. Automorphisms of Z/7 are a -> ka for k in {1,...,6}.
# For alpha: a -> ka to commute with a -> 2a: alpha(2a) = 2*alpha(a), so 2ka = k*2a = 2ka. Always true!
# Wait: we need alpha to be an automorphism of the WHOLE group F_21, not just Z/7.
# If alpha fixes the Z/3 part: alpha(i, j) = (k*i, j), then:
# alpha((i1,j1)(i2,j2)) = alpha(i1+2^j1*i2, j1+j2) = (k(i1+2^j1*i2), j1+j2)
# alpha(i1,j1)*alpha(i2,j2) = (ki1, j1)(ki2, j2) = (ki1 + 2^j1*ki2, j1+j2) = (k(i1+2^j1*i2), j1+j2)
# These are equal! So alpha(i,j) = (ki, j) for k in {1,...,6} are all automorphisms.
# But also: alpha could change the Z/3 part. alpha(0,1) must have order 3.
# If alpha(0,1) = (c, 1) for some c, then alpha(i,0) = alpha((1,0)^i) = (ki, 0) as before.
# alpha(i,j) = alpha((i,0)*(0,j)) = (ki, 0) * alpha(0,1)^j.
# alpha(0,1) = (c, 1). So alpha(0,1)^2 = (c, 1)(c, 1) = (c+2c, 2) = (3c, 2).
# alpha(0,1)^3 = (3c, 2)(c, 1) = (3c + 4c, 0) = (0, 0) mod 7. Need 7c = 0 mod 7. Always true!
# So alpha(0,1) = (c, 1) for any c in Z/7.
# alpha(i, j) = (ki, 0) * (c, 1)^j.
# (c, 1)^0 = (0, 0), (c, 1)^1 = (c, 1), (c, 1)^2 = (3c, 2).
# So alpha(i, 0) = (ki, 0), alpha(i, 1) = (ki+c, 1), alpha(i, 2) = (ki+3c, 2).

# Total automorphisms: 6 choices for k * 7 choices for c = 42.
# But |Aut(F_21)| = 42 (since |Out(F_21)| = 2 and |Inn(F_21)| = 21? or different).
# Actually |Inn(F_21)| = |F_21/Z(F_21)| = 21 (Z = {e} since non-abelian).
# And Out: the k values give 6 outer automorphisms? Not exactly,
# since inner auts also act on Z/7 by conjugation (a -> 2^j * a).
# Conjugation by (0, j) maps (i, 0) -> (2^j * i, 0).
# So the inner automorphisms already include k = 1, 2, 4 (from conjugation by (0,0), (0,1), (0,2)).
# The remaining k values {3, 5, 6} give genuinely new (outer) automorphisms.

# For our purposes, let's try all combinations:
# sigma(g) = L_h(alpha(g^{-1})) = h * alpha(g^{-1}) for some h in F_21, alpha in Aut(F_21)

found = False
for k in range(1, 7):
    for c in range(7):
        # alpha(i, j): need the formula
        def auto(g, k=k, c=c):
            i, j = g
            if j == 0:
                return ((k * i) % 7, 0)
            elif j == 1:
                return ((k * i + c) % 7, 1)
            else:
                return ((k * i + 3 * c) % 7, 2)

        for h in elements:
            # sigma(g) = h * alpha(g^{-1})
            def sigma(g, h=h, auto=auto):
                return mul(h, auto(inv(g)))

            # Check: A[sigma(a)][sigma(b)] = A[b][a] for all a, b?
            ok = True
            for a in elements:
                sa = sigma(a)
                for b in elements:
                    if a == b:
                        continue
                    sb = sigma(b)
                    if A[idx[sa]][idx[sb]] != A[idx[b]][idx[a]]:
                        ok = False
                        break
                if not ok:
                    break

            if ok:
                print(f"  FOUND anti-automorphism: h={h}, k={k}, c={c}")
                print(f"  sigma(g) = {h} * alpha_{{k={k},c={c}}}(g^{{-1}})")
                found = True
                break
        if found:
            break
    if found:
        break

if not found:
    print("  No anti-automorphism of the form h*alpha(g^{-1}) found!")
    print("  The anti-automorphism might not be group-based.")

    # Try: just search for ANY anti-automorphism using the GROUP action.
    # Since G acts transitively, we can fix sigma(e_id) = some vertex h,
    # and then sigma is determined by its values on a generating set...
    # Actually no, sigma need not be a group element.

    # For n=21, searching all 21! permutations is impossible.
    # But we can use backtracking search.
    print("  Attempting backtracking search for anti-automorphism...")

    # Backtracking: assign sigma vertex by vertex
    def backtrack_anti_auto(A, AT, assignment, remaining):
        """Find anti-automorphism by backtracking.
        assignment[i] = sigma(i) for assigned vertices.
        We need A[assignment[i]][assignment[j]] = A[j][i] for all assigned i,j.
        """
        if not remaining:
            return assignment.copy()

        v = remaining[0]
        rest = remaining[1:]
        used = set(assignment.values())

        for target in range(n):
            if target in used:
                continue
            # Check consistency: for all already-assigned u,
            # A[assignment[u]][target] = A[v][u]
            # A[target][assignment[u]] = A[u][v]
            ok = True
            for u, su in assignment.items():
                if A[su][target] != A[v][u]:
                    ok = False
                    break
                if A[target][su] != A[u][v]:
                    ok = False
                    break
            if ok:
                assignment[v] = target
                result = backtrack_anti_auto(A, AT, assignment, rest)
                if result is not None:
                    return result
                del assignment[v]

        return None

    ordering = list(range(n))
    result = backtrack_anti_auto(A, AT, {}, ordering)
    if result is not None:
        print(f"  FOUND anti-automorphism via backtracking!")
        # Verify
        sigma = [result[i] for i in range(n)]
        ok = all(A[sigma[i]][sigma[j]] == A[j][i]
                 for i in range(n) for j in range(n) if i != j)
        print(f"  Verification: {ok}")
        print(f"  sigma = {sigma}")
    else:
        print(f"  NO ANTI-AUTOMORPHISM EXISTS!")
        print(f"  This F_21 Cayley tournament is NOT self-converse!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
