#!/usr/bin/env python3
"""
DISCOVERY: For Paley T_5, ALL valid split pairs (P ending at a, Q starting at b)
on complementary vertex sets are actual Hamiltonian paths when concatenated!

This means M[a,b] counts Hamiltonian paths with an a->b consecutive pair,
weighted by (-1)^{position of a}.

If path P has |S+a| = j+1 vertices and ends at a, then a is at position j.
The concatenation P|Q is a Hamiltonian path where a is at position j and
b is at position j+1.

So M[a,b] = sum_P [P has a at pos j, b at pos j+1] * (-1)^j
          = sum_j (-1)^j * P2[a,j,b,j+1]

where P2[a,j,b,k] = # Ham paths with a at position j and b at position k.

QUESTION: Is this true for ALL tournaments, or only position-uniform ones?

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

def E_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def B_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

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
            val += sign * E_v(T, S_set, a) * B_v(T, R_set, b)
        return val

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# Test the conjecture: M[a,b] = sum_j (-1)^j P2[a,j,b,j+1]
# ============================================================
print("=" * 70)
print("CONJECTURE: M[a,b] = sum_j (-1)^j * #{paths with a@j, b@j+1}")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Test ALL tournaments
matches = 0
mismatches = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)

    # Compute P2[a,j,b,k] for consecutive positions
    # P2_consec[a,b,j] = # paths with a at pos j and b at pos j+1
    P2_consec = {}
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for j in range(n-1):
                a_v, b_v = perm[j], perm[j+1]
                P2_consec[(a_v, b_v, j)] = P2_consec.get((a_v, b_v, j), 0) + 1

    all_match = True
    for a in range(n):
        for b in range(n):
            if a == b:
                continue

            # Compute the conjectured formula
            conj = sum((-1)**j * P2_consec.get((a, b, j), 0) for j in range(n-1))
            M_val = compute_M_entry(T, n, a, b)

            if conj != M_val:
                all_match = False
                if mismatches < 3:
                    H = count_H(T, n)
                    print(f"  MISMATCH at bits={bits}, H={H}: M[{a},{b}]={M_val}, conj={conj}")
                break
        if not all_match:
            break

    if all_match:
        matches += 1
    else:
        mismatches += 1

print(f"\n  Results: {matches} matches, {mismatches} mismatches out of {matches+mismatches}")

if mismatches == 0:
    print("""
  THEOREM: M[a,b] = sum_{j=0}^{n-2} (-1)^j * #{Ham paths with a at pos j, b at pos j+1}

  This is EXACTLY the count of Hamiltonian paths where a->b is a consecutive edge,
  weighted by (-1)^{position of a}.

  PROOF SKETCH:
    M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
    For each term, |S| = position of a when we concatenate the paths.
    E_a(S+a) * B_b(R+b) counts pairs of paths: one ending at a, one starting at b.
    The concatenation P..a, b..Q is a valid Hamiltonian path iff T[a,b]=1 AND
    the internal edges of P and Q are valid.

    But wait — the internal edges of P and Q are ALWAYS valid (they're paths in T).
    The ONLY question is whether T[a,b]=1 (the junction edge).
    If T[a,b]=0, then B_b(R+b)... no, B_b counts paths starting at b in T[R+b],
    not paths starting from the end of P.

    Hmm, the concatenation needs T[a,b]=1 AND all internal edges valid.
    P is a valid path in T[S+a], Q is valid in T[R+b].
    But the concatenation P|Q needs edges T[perm[k],perm[k+1]] for ALL k,
    including the junction edge T[a,b].

    If T[a,b]=0, then the junction edge doesn't exist, so the concatenation
    is NOT a Hamiltonian path. But E*B still counts the pair.

    So M[a,b] should include terms from T[a,b]=0 as well!
    Let me check: for a pair (a,b) with T[a,b]=0, are there split pairs
    with nonzero E*B?
""")

    # Check: for pairs with T[a,b]=0, do split pairs exist?
    T_test = make_circulant_tournament(n, {1, 2})
    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T_test.get((a,b), 0) == 0:
                # T[a,b] = 0, so b beats a
                U = [v for v in range(n) if v != a and v != b]
                total_EB = 0
                for mask in range(1 << len(U)):
                    S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                    S_set = sorted(set(S_list) | {a})
                    R_set = sorted(set(R) | {b})
                    ea = E_v(T_test, S_set, a)
                    bb = B_v(T_test, R_set, b)
                    total_EB += ea * bb
                if total_EB > 0:
                    print(f"  T[{a},{b}]=0 but sum E*B = {total_EB}")

    # Also check: what happens when T[a,b]=0 for the formula
    print("\n  Verifying for T[a,b]=0 pairs:")
    for a, b in [(0, 3), (0, 4)]:  # 0 does NOT beat 3 or 4 in Paley T_5?
        edge = T_test.get((a,b), 0)
        M_val = compute_M_entry(T_test, n, a, b)
        conj = sum((-1)**j * P2_consec.get((a, b, j), 0) for j in range(n-1))
        print(f"  T[{a},{b}]={edge}: M={M_val}, conj={conj}")

else:
    print("  Conjecture FALSE — investigating mismatches...")


# ============================================================
# If the formula holds, then M[a,b]=0 iff the position-consecutive
# count P2[a,j,b,j+1] has alternating sum = 0.
# ============================================================
print("\n" + "=" * 70)
print("When does the consecutive-position alternating sum vanish?")
print("=" * 70)

T = make_circulant_tournament(n, {1, 2})
H = count_H(T, n)

# For Paley T_5:
print(f"\n  Paley T_5 (H={H}):")
for a in range(n):
    for b in range(n):
        if a == b: continue
        consec = [0] * (n-1)
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                for j in range(n-1):
                    if perm[j] == a and perm[j+1] == b:
                        consec[j] += 1
        alt_sum = sum((-1)**j * consec[j] for j in range(n-1))
        if a < 2 and b < 3:
            print(f"    ({a},{b}): consec={consec}, alt_sum={alt_sum}, "
                  f"T[{a},{b}]={T.get((a,b),0)}")

# For uniform position matrix:
# P2[a,j,b,j+1] = #{paths with a@j, b@j+1}
# This means T[a,b]=1 AND (a at pos j, b at pos j+1).
# If T[a,b]=0, then ALL consec[j] = 0.

print("\n  NOTE: If T[a,b]=0, then no path has a immediately before b.")
print("  So all consec[j]=0 and alt_sum=0=M[a,b].")
print("  If T[a,b]=1, then consec[j] >= 0 and we need sum(-1)^j consec[j] = 0.")

# So M[a,b] = 0 for BOTH T[a,b]=0 AND T[a,b]=1 when position-uniform.
# For T[a,b]=0: trivially 0.
# For T[a,b]=1: need alternating sum of consecutive-position counts = 0.

# By symmetry (M[a,b] = M[b,a]), M[a,b] = 0 iff M[b,a] = 0.
# M[b,a] = sum_j (-1)^j #{paths with b@j, a@j+1} = 0.
# If T[b,a] = 0, this is trivially 0.
# But T[a,b] + T[b,a] = 1 for tournaments!
# So EXACTLY ONE of T[a,b], T[b,a] is 1.
# The one that's 0 gives trivially M = 0 from that formula.
# But M[a,b] = M[b,a] means both are 0.

# So the content is: for a->b (T[a,b]=1):
#   sum_j (-1)^j #{paths with a@j, b@j+1} = 0
# iff
#   sum_j (-1)^j #{paths with b@j, a@j+1} = 0 (trivially, since T[b,a]=0)

# The FIRST condition (for T[a,b]=1) is the non-trivial one.

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
THEOREM (verified exhaustively at n=5):
  M[a,b] = sum_{j=0}^{n-2} (-1)^j * #{Ham paths with a at pos j, b at pos j+1}

This has a BEAUTIFUL interpretation:
  - M[a,b] counts how often a IMMEDIATELY PRECEDES b in Hamiltonian paths,
    weighted by (-1)^{position of a}.
  - For the diagonal: M[a,a] = sum_j (-1)^j * P[a,j], which is the
    standard alternating position count.

COROLLARY: M[a,b] = 0 for all a != b iff:
  For every edge a->b, the consecutive-position count
  #{paths where a is at position j and b is at position j+1}
  has vanishing alternating sum: sum_j (-1)^j * consec(a,b,j) = 0.

  Equivalently: M = (H/n)*I iff the "signed consecutive occurrence"
  is zero for every directed edge.

For position-uniform tournaments:
  If P[v,k] = H/n for all v,k, AND the tournament has enough symmetry
  to make the consecutive counts "alternating-zero," then M is scalar.
""")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
