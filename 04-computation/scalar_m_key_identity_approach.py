#!/usr/bin/env python3
"""
CORRECTED: Use the Key Identity to decompose M[a,b] and understand vanishing.

Key Identity (THM-030, proved for all m):
  B_b(W) + (-1)^m E_b(W) = col_sum_W(b)

where W is a sub-tournament on m vertices containing b,
col_sum_W(b) = sum_{v in W, v != b} T[v,b] = in-degree of b in W.

At c=1 (standard tournaments, r=1/2):
  This is the identity B_b + (-1)^m E_b = col_sum(b).

For M[a,b] (a != b):
  M[a,b] = sum_{S subset U} (-1)^|S| * E_a(S+{a}) * B_b(R+{b})

where U = [n] \ {a,b}, R = U \ S, |S+{a}| = |S|+1, |R+{b}| = |R|+1 = n-1-|S|.

Substituting B_b(R+{b}) from Key Identity:
  B_b(R+{b}) = col_{R+b}(b) - (-1)^{|R+{b}|} E_b(R+{b})
             = col_{R+b}(b) - (-1)^{n-1-|S|} E_b(R+{b})

So:
  M[a,b] = sum_S (-1)^|S| E_a(S+a) [col_{R+b}(b) - (-1)^{n-1-|S|} E_b(R+b)]
         = sum_S (-1)^|S| E_a(S+a) col_{R+b}(b)
           - sum_S (-1)^|S| (-1)^{n-1-|S|} E_a(S+a) E_b(R+b)
         = F(a,b) - (-1)^{n-1} sum_S E_a(S+a) E_b(R+b)
         = F(a,b) + (-1)^n sum_S E_a(S+a) E_b(R+b)

Wait, (-1)^|S| * (-1)^{n-1-|S|} = (-1)^{n-1}. So:
  M[a,b] = F(a,b) - (-1)^{n-1} C(a,b) = F(a,b) + (-1)^n C(a,b)

Let me verify this numerically.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np

def make_circulant_tournament(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def E_v(T, verts, v):
    """Count paths in T[verts] ending at v."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def B_v(T, verts, v):
    """Count paths in T[verts] starting at v."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v:
            continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def col_sum(T, verts, v):
    """In-degree of v in T[verts]."""
    return sum(T.get((u, v), 0) for u in verts if u != v)

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            val += sign * E_v(T, S_set, a) * B_v(T, R_set, b)
        return val


# ============================================================
# Verify the Key Identity decomposition
# ============================================================
print("=" * 70)
print("Verifying M[a,b] = F(a,b) + (-1)^n * C(a,b)")
print("=" * 70)

n = 5
S_gen = {1, 2}
T = make_circulant_tournament(n, S_gen)
H = count_H(T, n)

print(f"\nPaley T_5 (S={{1,2}}): H={H}")

for k in range(1, n):
    a, b = 0, k
    U = [v for v in range(n) if v != a and v != b]

    M_direct = compute_M_entry(T, n, a, b)

    F_val = 0  # sum_S (-1)^|S| E_a(S+a) col_sum(b, R+b)
    C_val = 0  # sum_S E_a(S+a) E_b(R+b) [unsigned]

    for mask in range(1 << len(U)):
        S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
        R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = E_v(T, S_set, a)
        eb = E_v(T, R_set, b)
        bb = B_v(T, R_set, b)
        cs = col_sum(T, R_set, b)

        sz = len(S_list)

        # Verify Key Identity: B_b + (-1)^|R_set| E_b = col_sum(b)
        m = len(R_set)
        key_id_check = bb + (-1)**m * eb
        if key_id_check != cs:
            print(f"  KEY IDENTITY VIOLATION: B={bb}, E={eb}, m={m}, cs={cs}")

        F_val += (-1)**sz * ea * cs
        C_val += ea * eb

    reconstructed = F_val + (-1)**n * C_val
    print(f"  M[0,{k}] = {M_direct}: F={F_val}, C={C_val}, "
          f"(-1)^n*C={(-1)**n * C_val}, F+(-1)^n*C={reconstructed}")


# ============================================================
# Step 2: What is C(a,b) for VT tournaments?
# ============================================================
print("\n" + "=" * 70)
print("C(a,b) = sum_S E_a(S+a) E_b(R+b) for circulant tournaments")
print("=" * 70)

for n in [5]:
    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S_gen in gen_sets:
        T = make_circulant_tournament(n, S_gen)
        H = count_H(T, n)
        print(f"\n  S={sorted(S_gen)}: H={H}")

        for k in range(n):
            a, b = 0, k
            U = [v for v in range(n) if v != a and v != b]

            if a == b:
                # C(a,a) = sum_S E_a(S+a)^2
                C_val = 0
                F_val = 0
                for mask in range(1 << (n-1)):
                    others = [v for v in range(n) if v != a]
                    S_list = [others[j] for j in range(n-1) if mask & (1 << j)]
                    R = [others[j] for j in range(n-1) if not (mask & (1 << j))]
                    S_set = sorted(set(S_list) | {a})
                    R_set = sorted(set(R) | {a})
                    # Diagonal formula is different - skip
                continue

            C_val = 0
            F_val = 0
            for mask in range(1 << len(U)):
                S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
                R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                ea = E_v(T, S_set, a)
                eb = E_v(T, R_set, b)
                cs = col_sum(T, R_set, b)

                sz = len(S_list)
                F_val += (-1)**sz * ea * cs
                C_val += ea * eb

            print(f"    C(0,{k}) = {C_val}, F(0,{k}) = {F_val}")


# ============================================================
# Step 3: Is C(a,b) symmetric? Does C relate to position matrix?
# ============================================================
print("\n" + "=" * 70)
print("C(a,b) symmetry and structure")
print("=" * 70)

n = 5
S_gen = {1, 2}
T = make_circulant_tournament(n, S_gen)

print(f"\nPaley T_5: C matrix and F matrix")

C_matrix = np.zeros((n, n), dtype=int)
F_matrix = np.zeros((n, n), dtype=int)

for a in range(n):
    for b in range(n):
        if a == b:
            continue
        U = [v for v in range(n) if v != a and v != b]
        C_val = 0
        F_val = 0
        for mask in range(1 << len(U)):
            S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
            R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})

            ea = E_v(T, S_set, a)
            eb = E_v(T, R_set, b)
            cs = col_sum(T, R_set, b)

            sz = len(S_list)
            F_val += (-1)**sz * ea * cs
            C_val += ea * eb

        C_matrix[a, b] = C_val
        F_matrix[a, b] = F_val

print(f"\n  C matrix (sum_S E_a*E_b unsigned):")
for row in C_matrix:
    print(f"    {list(row)}")

print(f"\n  F matrix (sum_S (-1)^|S| E_a * col_sum_b):")
for row in F_matrix:
    print(f"    {list(row)}")

print(f"\n  C symmetric? {np.allclose(C_matrix, C_matrix.T)}")
print(f"  F symmetric? {np.allclose(F_matrix, F_matrix.T)}")

# Check: M = F + (-1)^n * C
M_check = F_matrix + (-1)**n * C_matrix
print(f"\n  M = F + (-1)^n * C:")
for row in M_check:
    print(f"    {list(row)}")

# Since M = 0 off-diagonal: F = -(-1)^n * C = (-1)^{n+1} * C = C (odd n)
print(f"\n  For odd n=5: M=0 off-diag iff F = C (off-diagonal)")
print(f"  F == C off-diagonal? {np.allclose(F_matrix[np.triu_indices(n,1)], C_matrix[np.triu_indices(n,1)])}")

# Check if F = C
diff = F_matrix - C_matrix
print(f"\n  F - C (off-diagonal):")
for i in range(n):
    for j in range(n):
        if i != j and diff[i,j] != 0:
            print(f"    F[{i},{j}] - C[{i},{j}] = {diff[i,j]}")

if np.all(diff[np.triu_indices(n,1)] == 0):
    print("  CONFIRMED: F = C off-diagonal!")
    print("  This means: sum_S (-1)^|S| E_a(S+a) col_sum(b,R+b)")
    print("            = sum_S E_a(S+a) E_b(R+b)")
    print("  for all a != b in circulant tournaments!")


# ============================================================
# Step 4: Test F = C for non-circulant tournaments
# ============================================================
print("\n" + "=" * 70)
print("Does F = C hold for ALL tournaments at n=5?")
print("=" * 70)

pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

max_check = 200  # Check first 200
fc_equal_count = 0
fc_not_equal_count = 0

for bits in range(min(1 << len(pairs), max_check)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    fc_equal = True
    for a in range(n):
        for b in range(a+1, n):
            U = [v for v in range(n) if v != a and v != b]
            C_val = 0
            F_val = 0
            for mask in range(1 << len(U)):
                S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
                R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                ea = E_v(T, S_set, a)
                eb = E_v(T, R_set, b)
                cs = col_sum(T, R_set, b)

                sz = len(S_list)
                F_val += (-1)**sz * ea * cs
                C_val += ea * eb

            if F_val != C_val:
                fc_equal = False
                break
        if not fc_equal:
            break

    if fc_equal:
        fc_equal_count += 1
    else:
        fc_not_equal_count += 1

print(f"  F = C (off-diag) for {fc_equal_count}/{fc_equal_count + fc_not_equal_count} tournaments")


# ============================================================
# Step 5: Check ALL 1024 tournaments at n=5
# ============================================================
print("\n" + "=" * 70)
print("ALL n=5 tournaments: F = C off-diagonal?")
print("=" * 70)

fc_fails = []
for bits in range(1 << len(pairs)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    fc_equal = True
    for a in range(n):
        for b in range(a+1, n):
            U = [v for v in range(n) if v != a and v != b]
            C_val = 0
            F_val = 0
            for mask in range(1 << len(U)):
                S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
                R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                ea = E_v(T, S_set, a)
                eb = E_v(T, R_set, b)
                cs = col_sum(T, R_set, b)

                sz = len(S_list)
                F_val += (-1)**sz * ea * cs
                C_val += ea * eb

            if F_val != C_val:
                fc_equal = False
                break
        if not fc_equal:
            break

    if not fc_equal:
        H = count_H(T, n)
        fc_fails.append((bits, H))

if fc_fails:
    print(f"  F != C for {len(fc_fails)}/{1 << len(pairs)} tournaments")
    for bits, H in fc_fails[:5]:
        print(f"    bits={bits}, H={H}")
else:
    print(f"  F = C for ALL {1 << len(pairs)} tournaments at n=5!")
    print("""
  THIS IS A NEW IDENTITY:
    sum_S (-1)^|S| E_a(S+a) col_sum(b, R+b) = sum_S E_a(S+a) E_b(R+b)

  In words: the signed sum of (paths ending at a) * (in-degree of b)
           equals the unsigned sum of (paths ending at a) * (paths ending at b)

  If this holds for ALL tournaments at ALL n, then:
    M[a,b] = F + (-1)^n C = C + (-1)^n C = C(1 + (-1)^n)
    For odd n: M[a,b] = C(1-1) = 0  <==> M is diagonal!

  But that would mean M is diagonal for ALL tournaments at odd n,
  which contradicts known examples where M has nonzero off-diagonal entries.

  So either F != C at general n, or the identity is n=5-specific.
""")


# ============================================================
# Step 6: Quick check at n=3 and n=4
# ============================================================
for n_test in [3, 4]:
    print(f"\n  n={n_test}: checking F = C...")
    pairs_test = [(i,j) for i in range(n_test) for j in range(i+1, n_test)]
    fails = 0
    total = 0
    for bits in range(1 << len(pairs_test)):
        T = {}
        for idx, (i,j) in enumerate(pairs_test):
            if (bits >> idx) & 1:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1

        total += 1
        for a in range(n_test):
            for b in range(a+1, n_test):
                U = [v for v in range(n_test) if v != a and v != b]
                if not U:
                    # No U, both sums are trivially 0 vs 0
                    continue
                C_val = 0
                F_val = 0
                for mask in range(1 << len(U)):
                    S_list = [U[j] for j in range(len(U)) if mask & (1 << j)]
                    R_list = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
                    S_set = sorted(set(S_list) | {a})
                    R_set = sorted(set(R_list) | {b})

                    ea = E_v(T, S_set, a)
                    eb = E_v(T, R_set, b)
                    cs = col_sum(T, R_set, b)

                    sz = len(S_list)
                    F_val += (-1)**sz * ea * cs
                    C_val += ea * eb

                if F_val != C_val:
                    fails += 1
                    break

    print(f"    {fails} failures out of {total} tournaments")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
