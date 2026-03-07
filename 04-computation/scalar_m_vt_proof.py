#!/usr/bin/env python3
"""
PROOF: M = (H/n)*I for vertex-transitive tournaments at odd n.

STRATEGY:
For a VT tournament T with automorphism group G acting transitively on [n]:
  (1) M[sigma(a), sigma(b)] = M[a,b] for all sigma in G, all a,b
      => M is a G-equivariant matrix
  (2) M is symmetric (THM-030)
  (3) At odd n: tr(M) = H(T), and by VT: M[a,a] = H/n for all a

CLAIM: For VT at odd n, M[a,b] = 0 whenever a != b.

PROOF APPROACH (representation theory):
  G acts transitively on [n]. The permutation representation C^n decomposes as:
    C^n = W_triv + W_1 + ... + W_k
  where W_triv = span(1,...,1) is the trivial representation.

  M is a G-equivariant symmetric operator on C^n.
  By Schur's lemma, M acts as a scalar on each irreducible component.

  But does M actually act as a SINGLE scalar on ALL components?
  Only if all eigenvalues are equal.

  For circulant (G = Z/nZ):
    C^n = chi_0 + chi_1 + ... + chi_{n-1} (all 1-dim irreps)
    M acts as lambda_j on chi_j.
    M = (H/n)*I iff lambda_j = H/n for all j.

  KEY: We need to show lambda_j = H/n for j >= 1.
  lambda_0 = H/n is guaranteed (M*1 = H*1 for... wait, does M*1 = H*1?)

Let me verify: Does M*[1,...,1] = H*[1,...,1]?

ACTUALLY: sum_b M[a,b] = M[a,a] + sum_{b!=a} M[a,b]
  = diag_a + off_sum_a

If M = (H/n)*I, then sum_b M[a,b] = H/n (just the diagonal).
But sum_b M[a,b] in general is NOT H.

Let me check what sum_b M[a,b] equals.

REVISED: M is the "transfer matrix" with M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b).
The row sum sum_b M[a,b] is:
  sum_b sum_S (-1)^|S| E_a(S+a) B_b(R+b)
  = sum_S (-1)^|S| E_a(S+a) * sum_b B_b(R+b)
  = sum_S (-1)^|S| E_a(S+a) * H(T[R+b])  ... no, b is in R+b.

Wait. For fixed a, the sum over b (b != a) involves different splittings.
For b = a: M[a,a] = sum_P (-1)^{pos(a,P)} over Ham paths P.

For b != a: M[a,b] = sum_{S subset U} (-1)^|S| E_a(S+a) B_b(R+b)
  where U = [n]\{a,b}, R = U\S.

The row sum is more complex. Let me just compute it.

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
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_entry(T, n, a, a)
        for b in range(a+1, n):
            val = compute_M_entry(T, n, a, b)
            M[a, b] = val
            M[b, a] = val
    return M

def position_matrix(T, n):
    """P[v, k] = number of Ham paths where vertex v is at position k."""
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P


# ============================================================
# Step 1: Verify M*1 vector and row sums
# ============================================================
print("=" * 70)
print("Step 1: Row sums of M and M*[1,...,1]")
print("=" * 70)

for n in [5]:
    print(f"\n  n={n}:")

    # Find all valid circulant generating sets
    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S in gen_sets:
        T = make_circulant_tournament(n, S)
        H = count_H(T, n)
        M = compute_M(T, n)

        row_sums = [sum(M[a, :]) for a in range(n)]
        M_times_ones = M @ np.ones(n)

        print(f"\n    S={sorted(S)}: H={H}")
        print(f"      M = {M.tolist()}")
        print(f"      Row sums: {row_sums}")
        print(f"      M*1 = {[int(x) for x in M_times_ones]}")

        P = position_matrix(T, n)
        print(f"      Position matrix P:")
        for v in range(n):
            print(f"        v={v}: {list(P[v, :])}")

        # Signed position count: sum_P (-1)^{pos(v,P)} for each v
        signed_pos = np.zeros(n, dtype=int)
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                for v in range(n):
                    pos = list(perm).index(v)
                    signed_pos[v] += (-1)**pos
        print(f"      M diagonal (signed pos): {list(signed_pos)}")


# ============================================================
# Step 2: The key identity — what makes off-diagonals vanish?
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Off-diagonal structure for circulant tournaments")
print("=" * 70)

n = 5
for S_tuple in combinations(range(1, n), (n-1)//2):
    S = set(S_tuple)
    complement = {(n - s) % n for s in S}
    if S & complement or (S | complement) != set(range(1, n)):
        continue

    T = make_circulant_tournament(n, S)
    H = count_H(T, n)

    print(f"\n  S={sorted(S)}: H={H}")

    # Compute E_a(S_set) and B_b(R_set) for specific (a,b) pair
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    print(f"    Off-diagonal M[{a},{b}] decomposition:")
    total = 0
    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = 0
        if len(S_set) == 1:
            ea = 1
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod

        bb2 = 0
        if len(R_set) == 1:
            bb2 = 1
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod

        contrib = sign * ea * bb2
        total += contrib
        if ea > 0 and bb2 > 0:
            print(f"      S={sorted(S_list)}, R={sorted(R)}: "
                  f"sign={sign:+d}, E_{a}({sorted(S_set)})={ea}, "
                  f"B_{b}({sorted(R_set)})={bb2}, "
                  f"contrib={contrib:+d}")

    print(f"      TOTAL M[{a},{b}] = {total}")


# ============================================================
# Step 3: Representation theory approach
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Eigenvalues of M for circulant tournaments (DFT analysis)")
print("=" * 70)

for n in [5, 7]:
    print(f"\n  n={n}:")

    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S in gen_sets[:3]:  # Limit at n=7 for speed
        T = make_circulant_tournament(n, S)
        H = count_H(T, n)
        M = compute_M(T, n)

        # Eigenvalues via DFT (M is circulant)
        row0 = M[0, :]
        evals_dft = np.fft.fft(row0)
        evals_sym = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]

        is_scalar = all(abs(M[i,j]) < 0.01 for i in range(n) for j in range(n) if i != j)

        print(f"\n    S={sorted(S)}: H={H}, scalar={is_scalar}")
        print(f"      row0 = {list(row0)}")
        print(f"      DFT evals = {[f'{e.real:.2f}' for e in evals_dft]}")
        print(f"      eigh evals = {[f'{e:.2f}' for e in evals_sym]}")

        if not is_scalar:
            print(f"      *** NOT SCALAR — this would be a counterexample! ***")


# ============================================================
# Step 4: The KEY algebraic identity to prove
# ============================================================
print("\n" + "=" * 70)
print("Step 4: Algebraic identity — why do off-diagonals vanish?")
print("=" * 70)

print("""
For circulant tournament with generating set S (i beats j iff (j-i) mod n in S):

  M[0,k] = sum_{subset W of U={all except 0,k}} (-1)^|W| E_0(W+{0}) B_k(U\\W+{k})

  where E_0(W+{0}) = # directed paths in T[W+{0}] ending at 0
        B_k(U\\W+{k}) = # directed paths in T[U\\W+{k}] starting at k

By cyclic symmetry sigma(x) = x+1 mod n:
  M[0,k] = M[1,k+1] = ... (circulant)

Also M[0,k] = M[0,n-k] (symmetry THM-030, since f(k) = f(n-k) for symmetric circulant)

CONJECTURE: For tournament circulants, the inclusion-exclusion in M[0,k] telescopes
to zero for every k != 0.

EVIDENCE: Let's look at the individual E_a and B_b values more carefully.
""")

n = 5
S_gen = {1, 2}  # Paley T_5
T = make_circulant_tournament(n, S_gen)
H = count_H(T, n)
print(f"  Paley T_5 (S={{1,2}}): H={H}")

for k in range(1, n):
    a, b = 0, k
    U = [v for v in range(n) if v != a and v != b]

    print(f"\n  M[0,{k}] decomposition:")

    # Organize by |W|
    by_size = {}
    for mask in range(1 << len(U)):
        W = [U[j] for j in range(len(U)) if mask & (1 << j)]
        R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]

        S_set = sorted(set(W) | {a})
        R_set = sorted(set(R) | {b})

        ea = 0
        if len(S_set) == 1:
            ea = 1
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1
                for pk in range(len(p)-1):
                    prod *= T.get((p[pk], p[pk+1]), 0)
                ea += prod

        bb2 = 0
        if len(R_set) == 1:
            bb2 = 1
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1
                for pk in range(len(p)-1):
                    prod *= T.get((p[pk], p[pk+1]), 0)
                bb2 += prod

        sz = len(W)
        if sz not in by_size:
            by_size[sz] = []
        by_size[sz].append((sorted(W), ea, bb2, (-1)**sz * ea * bb2))

    total = 0
    for sz in sorted(by_size.keys()):
        layer_sum = sum(x[3] for x in by_size[sz])
        total += layer_sum
        print(f"    |W|={sz}: layer_sum={layer_sum:+d}")
        for W, ea, bb2, contrib in by_size[sz]:
            if ea > 0 or bb2 > 0:
                print(f"      W={W}: E_0={ea}, B_{k}={bb2}, contrib={contrib:+d}")

    print(f"    TOTAL = {total}")


# ============================================================
# Step 5: Check if the alternating sum has a pattern
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Alternating layer sums for ALL circulant tournaments at n=5")
print("=" * 70)

n = 5
gen_sets = []
for S_tuple in combinations(range(1, n), (n-1)//2):
    S_set = set(S_tuple)
    complement = {(n - s) % n for s in S_set}
    if not (S_set & complement) and (S_set | complement) == set(range(1, n)):
        gen_sets.append(S_set)

for S_gen in gen_sets:
    T = make_circulant_tournament(n, S_gen)
    H = count_H(T, n)

    print(f"\n  S={sorted(S_gen)}: H={H}")

    for k in [1, 2]:  # Just check k=1,2 (others by symmetry)
        a, b = 0, k
        U = [v for v in range(n) if v != a and v != b]

        layer_sums = {}
        for mask in range(1 << len(U)):
            W = [U[j] for j in range(len(U)) if mask & (1 << j)]
            R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
            S_set_local = sorted(set(W) | {a})
            R_set = sorted(set(R) | {b})

            ea = 0
            if len(S_set_local) == 1:
                ea = 1
            else:
                for p in permutations(S_set_local):
                    if p[-1] != a: continue
                    prod = 1
                    for pk in range(len(p)-1):
                        prod *= T.get((p[pk], p[pk+1]), 0)
                    ea += prod

            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for pk in range(len(p)-1):
                        prod *= T.get((p[pk], p[pk+1]), 0)
                    bb2 += prod

            sz = len(W)
            layer_sums[sz] = layer_sums.get(sz, 0) + (-1)**sz * ea * bb2

        layers = [layer_sums.get(sz, 0) for sz in range(len(U)+1)]
        total = sum(layers)
        print(f"    M[0,{k}]: layers={layers}, total={total}")


# ============================================================
# Step 6: Position distribution vs M diagonal
# ============================================================
print("\n" + "=" * 70)
print("Step 6: Position distribution analysis")
print("=" * 70)

n = 5
# Also test non-circulant tournaments with M scalar
from itertools import permutations as perms

pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

scalar_count = 0
non_vt_scalar = 0

for bits in range(1 << len(pairs)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    H = count_H(T, n)
    if H != 15:  # Only maximizers
        continue

    M = compute_M(T, n)
    is_scalar = all(M[i,j] == 0 for i in range(n) for j in range(n) if i != j)

    if is_scalar:
        scalar_count += 1
        P = position_matrix(T, n)
        is_uniform = all(P[v,k] == H//n for v in range(n) for k in range(n))

        # Check VT
        is_vt = False
        for perm in perms(range(n)):
            if perm == tuple(range(n)):
                continue
            is_aut = True
            for i in range(n):
                for j in range(n):
                    if i != j:
                        if T.get((i,j), 0) != T.get((perm[i], perm[j]), 0):
                            is_aut = False
                            break
                if not is_aut:
                    break
            if is_aut:
                # Check if this automorphism + all others give transitivity
                pass

        scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

        if scalar_count <= 5:
            print(f"  bits={bits}: scores={scores}, uniform_pos={is_uniform}")

print(f"\n  Total H=15 with scalar M: {scalar_count}")


# ============================================================
# Step 7: PROOF IDEA — use the Key Identity (THM-030)
# ============================================================
print("\n" + "=" * 70)
print("Step 7: Using Key Identity for VT proof")
print("=" * 70)

print("""
KEY IDENTITY (THM-030, proved):
  B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)

  where W is any sub-tournament on m vertices containing b,
  r = 1/2 for standard tournaments,
  col_sum_W(b) = sum_{v in W, v != b} T[v,b].

At r = 1/2 (standard tournament):
  B_b(W) + (-1)^m E_b(W) = col_sum_W(b)

For odd m: B_b(W) - E_b(W) = col_sum_W(b)
For even m: B_b(W) + E_b(W) = col_sum_W(b)

This gives: B_b(W) = (col_sum_W(b) + (-1)^{m+1} * (-1)^m * ...

Actually, let me re-derive. The Key Identity at c=1 (r=1/2):
  B_b + (-1)^m E_b = col_sum(b)

So: B_b = col_sum(b) - (-1)^m E_b = col_sum(b) + (-1)^{m+1} E_b

For M[a,b] = sum_S (-1)^|S| E_a(S+{a}) B_b(R+{b}):
  Substitute B_b(R+{b}) = col_sum_{R+{b}}(b) + (-1)^{|R|} E_b(R+{b})
  where |R+{b}| = |R|+1 = n-1-|S|, so (-1)^{|R+{b}|} = (-1)^{n-1-|S|}

  => B_b = col_sum(b) + (-1)^{n-|S|} E_b

M[a,b] = sum_S (-1)^|S| E_a(S+a) * [col_sum_{R+b}(b) + (-1)^{n-|S|} E_b(R+b)]
        = sum_S (-1)^|S| E_a(S+a) * col_sum_{R+b}(b)
          + sum_S (-1)^|S| * (-1)^{n-|S|} E_a(S+a) E_b(R+b)
        = sum_S (-1)^|S| E_a(S+a) * col_sum_{R+b}(b)
          + (-1)^n sum_S E_a(S+a) E_b(R+b)

The second sum is UNSIGNED (no (-1)^|S| factor since (-1)^|S|*(-1)^{n-|S|} = (-1)^n).
By the S <-> R relabeling symmetry (swapping a,b):
  sum_S E_a(S+a) E_b(R+b) = sum_S E_b(S+b) E_a(R+a)

This is symmetric in (a,b)! So the second term is:
  (-1)^n * C(a,b) where C(a,b) = C(b,a)

For the FIRST term:
  F(a,b) = sum_S (-1)^|S| E_a(S+a) * col_sum_{R+b}(b)

col_sum_{R+b}(b) = sum_{v in R} T[v,b] = in-degree of b in T[R+{b}]

For VT: does this have enough symmetry to vanish?

Let me compute F(a,b) and C(a,b) separately for circulant T_5.
""")

n = 5
S_gen = {1, 2}  # Paley T_5
T = make_circulant_tournament(n, S_gen)
H = count_H(T, n)

for k in range(n):
    a, b = 0, k
    if a == b:
        continue
    U = [v for v in range(n) if v != a and v != b]

    F_sum = 0  # First term
    C_sum = 0  # Second term (unsigned)

    for mask in range(1 << len(U)):
        W = [U[j] for j in range(len(U)) if mask & (1 << j)]
        R = [U[j] for j in range(len(U)) if not (mask & (1 << j))]
        S_set = sorted(set(W) | {a})
        R_set = sorted(set(R) | {b})

        ea = 0
        if len(S_set) == 1:
            ea = 1
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1
                for pk in range(len(p)-1):
                    prod *= T.get((p[pk], p[pk+1]), 0)
                ea += prod

        eb = 0
        if len(R_set) == 1:
            eb = 1
        else:
            for p in permutations(R_set):
                if p[-1] != b: continue
                prod = 1
                for pk in range(len(p)-1):
                    prod *= T.get((p[pk], p[pk+1]), 0)
                eb += prod

        bb2 = 0
        if len(R_set) == 1:
            bb2 = 1
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1
                for pk in range(len(p)-1):
                    prod *= T.get((p[pk], p[pk+1]), 0)
                bb2 += prod

        # col_sum of b in R_set
        col_sum_b = sum(T.get((v, b), 0) for v in R)

        sz = len(W)
        F_sum += (-1)**sz * ea * col_sum_b
        C_sum += ea * eb  # unsigned

    M_val = compute_M_entry(T, n, a, b)
    reconstructed = F_sum + (-1)**n * C_sum

    print(f"  M[0,{k}] = {M_val}: F={F_sum}, (-1)^n*C={(-1)**n * C_sum}, F+(-1)^n*C={reconstructed}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
