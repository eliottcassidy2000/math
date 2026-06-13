#!/usr/bin/env python3
"""
WHY do the 40 non-regular n=5 H=15 tournaments have M=3*I
despite NON-UNIFORM E_a sums?

These have scores (1,2,2,2,3) and |Aut|=3 (Z/3Z).
The Z/3Z fixes two vertices (source-like and sink-like) and permutes three middle vertices.

The E_a sums are NOT uniform across all 5 vertices, but
M[a,b] = 0 for all a != b. So the inclusion-exclusion sum
cancels for reasons beyond pure Fourier uniformity.

HYPOTHESIS: The E*B product structure respects a hidden symmetry
that forces cancellation even when E alone is non-uniform.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np

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

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Find a non-regular H=15 tournament
for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    if H != 15:
        continue
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    if scores != (1, 2, 2, 2, 3):
        continue

    # Found one!
    print("=" * 70)
    print(f"Non-regular H=15 tournament (bits={bits})")
    print("=" * 70)

    # Identify vertices by score
    out_degrees = [sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)]
    print(f"\n  Out-degrees: {out_degrees}")

    # E_a sums by subset size for each vertex
    print(f"\n  E_a sums by subset size:")
    for a in range(n):
        others = [v for v in range(n) if v != a]
        e_sums = []
        for s in range(n):
            e_sum = 0
            for S_tuple in combinations(others, s):
                S_set = sorted(list(S_tuple) + [a])
                e_sum += E_v(T, S_set, a)
            e_sums.append(e_sum)
        alt_sum = sum((-1)**k * e_sums[k] for k in range(n))
        print(f"    a={a} (deg={out_degrees[a]}): e_s={e_sums}, alt_sum={alt_sum}")

    # B_a sums by subset size for each vertex
    print(f"\n  B_a sums by subset size:")
    for a in range(n):
        others = [v for v in range(n) if v != a]
        b_sums = []
        for s in range(n):
            b_sum = 0
            for S_tuple in combinations(others, s):
                S_set = sorted(list(S_tuple) + [a])
                b_sum += B_v(T, S_set, a)
            b_sums.append(b_sum)
        alt_sum = sum((-1)**k * b_sums[k] for k in range(n))
        print(f"    a={a} (deg={out_degrees[a]}): b_s={b_sums}, alt_sum={alt_sum}")

    # Now look at M[a,b] decomposition for specific pairs
    print(f"\n  M[a,b] decomposition for non-uniform pairs:")
    for a in range(n):
        for b in range(a+1, n):
            U = [v for v in range(n) if v != a and v != b]
            terms = []
            total = 0
            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})
                ea = E_v(T, S_set, a)
                bb = B_v(T, R_set, b)
                contrib = sign * ea * bb
                total += contrib
                if ea > 0 and bb > 0:
                    terms.append((sorted(S_list), ea, bb, contrib))

            if total != 0:
                print(f"    M[{a},{b}] = {total} *** NONZERO!")
            else:
                print(f"    M[{a},{b}] = 0: {len(terms)} nonzero terms, "
                      f"sum(+)={sum(c for _,_,_,c in terms if c>0)}, "
                      f"sum(-)={sum(c for _,_,_,c in terms if c<0)}")

    # Position matrix
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1

    print(f"\n  Position matrix P[v,k]:")
    for v in range(n):
        print(f"    v={v} (deg={out_degrees[v]}): {list(P[v,:])}")

    # Check: is position matrix uniform despite non-uniform E?
    is_uniform_pos = all(P[v,k] == 3 for v in range(n) for k in range(n))
    print(f"\n  Uniform position matrix: {is_uniform_pos}")

    # KEY: The position matrix IS uniform (P[v,k]=3 for all v,k)
    # but E_a(S+a) sums at intermediate levels are NOT uniform.
    # The cancellation in M happens at the FINAL inclusion-exclusion level,
    # not at the level of individual E values.

    # Let's look at the PAIRED structure: E_a * B_b for complementary subsets
    print(f"\n  E*B product analysis for a=0, b=1 (source-like, middle):")
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})
        ea = E_v(T, S_set, a)
        bb = B_v(T, R_set, b)
        eb = E_v(T, R_set, b)
        ba = B_v(T, S_set, a)
        sign = (-1)**len(S_list)
        print(f"    S={sorted(S_list)}: E_{a}(S+a)={ea}, B_{b}(R+b)={bb}, "
              f"E_{b}(R+b)={eb}, B_{a}(S+a)={ba}, "
              f"E*B={ea}*{bb}={ea*bb}, sign*E*B={sign*ea*bb:+d}")

    break  # Just analyze one representative


print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
For n=5 H=15 non-regular tournaments:
  - E_a(S+a) sums at intermediate subset sizes are NOT uniform
  - But the position matrix P[v,k] = 3 for all v,k (uniform positions)
  - And M[a,b] = 0 for all a != b (scalar transfer matrix)

This means: uniform POSITION distribution => scalar M,
even without uniform E at intermediate levels.

The cancellation in M[a,b] = sum_S (-1)^|S| E_a B_b uses the FULL
inclusion-exclusion, which telescopes to 0 despite individual
E_a(S+a) values being non-uniform.

This is consistent with M[a,a] = sum_P (-1)^pos(a,P) = 3 for all a.
The diagonal entries are uniform because every vertex appears 3 times
at each position, so the alternating sum is 3*(1-1+1-1+1) = 3.

Wait: 3*(1-1+1-1+1) = 3. Yes! Each vertex appears 3 times at each
of the 5 positions (0,1,2,3,4), so M[a,a] = 3*(1-1+1-1+1) = 3*1 = 3.

For the off-diagonal: M[a,b] = 0 must follow from some analogous
"uniform" property of the E*B product, but at the level of the
complete inclusion-exclusion, not at individual subset levels.

CONJECTURE: For ANY tournament at odd n with uniform position matrix
(P[v,k] = H/n for all v,k), the transfer matrix M = (H/n)*I.

This would unify all our observations:
  - VT => uniform positions (trivially)
  - Non-VT with Z/3Z aut at n=5 can also have uniform positions
  - Uniform positions => M scalar
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
