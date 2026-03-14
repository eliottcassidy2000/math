"""
cubic_transition_n9.py — The independence polynomial degree transition at n=9=3^2
and its connection to exponentiation, Cayley-Dickson, and forbidden H.
kind-pasteur-2026-03-14-S65

KEY FINDING from z_omega_tower.py Part 8:
  alpha(Omega) = 2 at n<=8 (quadratic I.P.)
  alpha(Omega) = 3 first possible at n=9 (cubic I.P.)
  Because 3+3+3 = 9 = minimum for 3 disjoint 3-cycles

This means:
  - The "quadratic era" (I.P. has at most degree 2) spans n=5 to n=8
  - The "cubic era" begins at n=9 = 3^2
  - The transition number 9 = 3^2 is the EXPONENTIATION of 3!

PART 1:  Verify alpha(Omega) values at n=6,7,8 (all <=2?)
PART 2:  Construct a tournament at n=9 with 3 disjoint 3-cycles
PART 3:  I(x) for cubic tournaments: the new coefficient a3
PART 4:  What H values become achievable with a3 > 0?
PART 5:  The forbidden H revisited: do some gaps CLOSE at n=9?
PART 6:  Connection: 9=3^2 and the degree-3 transition
PART 7:  The lattice of achievable (a1, a2, a3) at n=9
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
import math

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_directed_ham_cycles_sub(A_sub, m):
    if m < 3:
        return 0
    dp = {}
    start = 0
    dp[(1 << start, start)] = 1
    for mask_size in range(2, m + 1):
        for mask in range(1 << m):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(m):
                if not (mask & (1 << v)):
                    continue
                if v == start and mask_size < m:
                    continue
                prev_mask = mask ^ (1 << v)
                if not (prev_mask & 1):
                    continue
                total = 0
                for u in range(m):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = dp.get((mask, v), 0) + total
    full = (1 << m) - 1
    count = 0
    for v in range(m):
        if A_sub[v][start]:
            count += dp.get((full, v), 0)
    return count

def get_all_odd_cycles_detailed(A):
    """Return list of (vertex_set, length) for all directed odd cycles."""
    n = len(A)
    all_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            sub = [[A[verts[i]][verts[j]] for j in range(length)] for i in range(length)]
            hc = count_directed_ham_cycles_sub(sub, length)
            for _ in range(hc):
                all_cycles.append((frozenset(verts), length))
    return all_cycles

def get_independence_poly(A):
    """Compute full independence polynomial coefficients [1, a1, a2, a3, ...]."""
    cycles = get_all_odd_cycles_detailed(A)
    vertex_sets = [c[0] for c in cycles]
    a1 = len(vertex_sets)

    # Build conflict graph (Omega)
    n_cycles = len(vertex_sets)
    # Adjacency: two cycles conflict iff they share a vertex
    adj = [[False]*n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if not vertex_sets[i].isdisjoint(vertex_sets[j]):
                adj[i][j] = True
                adj[j][i] = True

    # Find independent sets by size (brute force for small n_cycles)
    max_indep = min(n_cycles, 10)  # cap for performance
    coeffs = [1]  # f_0 = 1 (empty set)

    for k in range(1, max_indep + 1):
        count = 0
        if k > n_cycles:
            break
        for combo in combinations(range(n_cycles), k):
            # Check independence
            is_indep = True
            for ii in range(len(combo)):
                for jj in range(ii+1, len(combo)):
                    if adj[combo[ii]][combo[jj]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                count += 1
        coeffs.append(count)
        if count == 0:
            break

    return coeffs

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
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
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


# ============================================================
# PART 1: Verify alpha(Omega) at n=6,7,8
# ============================================================
def part1():
    print("=" * 70)
    print("PART 1: VERIFY alpha(Omega) VALUES AT n=6,7,8")
    print("=" * 70)

    for n in [6, 7]:
        rng = np.random.default_rng(2026 + n)
        N_samples = 500 if n <= 7 else 100
        max_alpha = 0
        alpha_dist = defaultdict(int)

        for _ in range(N_samples):
            A = random_tournament(n, rng)
            coeffs = get_independence_poly(A)
            alpha = len(coeffs) - 1
            alpha_dist[alpha] = alpha_dist.get(alpha, 0) + 1
            max_alpha = max(max_alpha, alpha)

        print(f"\n  n={n} ({N_samples} samples):")
        print(f"    max alpha(Omega) = {max_alpha}")
        print(f"    alpha distribution: {dict(sorted(alpha_dist.items()))}")
        if max_alpha >= 3:
            print(f"    *** CUBIC I.P. FOUND AT n={n}! ***")

    # n=8 separately (smaller sample due to cost)
    n = 8
    rng = np.random.default_rng(2026 + n)
    N_samples = 100
    max_alpha = 0
    alpha_dist = defaultdict(int)
    cubic_example = None

    for trial in range(N_samples):
        A = random_tournament(n, rng)
        cycles = get_all_odd_cycles_detailed(A)
        # Quick check: can we find 3 mutually disjoint cycles?
        vertex_sets = [c[0] for c in cycles]
        n_cyc = len(vertex_sets)

        # Check max independent set size via greedy + small enumeration
        # Just check if alpha >= 3
        found_3 = False
        if n_cyc >= 3:
            for i in range(min(n_cyc, 50)):
                for j in range(i+1, min(n_cyc, 50)):
                    if not vertex_sets[i].isdisjoint(vertex_sets[j]):
                        continue
                    for k in range(j+1, min(n_cyc, 50)):
                        if (vertex_sets[i].isdisjoint(vertex_sets[k]) and
                            vertex_sets[j].isdisjoint(vertex_sets[k])):
                            found_3 = True
                            cubic_example = (A.copy(), cycles, i, j, k)
                            break
                    if found_3:
                        break
                if found_3:
                    break

        # Also compute alpha properly via independence poly
        coeffs = get_independence_poly(A)
        alpha = len(coeffs) - 1
        alpha_dist[alpha] += 1
        max_alpha = max(max_alpha, alpha)

        if (trial + 1) % 25 == 0:
            print(f"    n=8: {trial+1}/{N_samples} done, max alpha so far = {max_alpha}")

    print(f"\n  n=8 ({N_samples} samples):")
    print(f"    max alpha(Omega) = {max_alpha}")
    print(f"    alpha distribution: {dict(sorted(alpha_dist.items()))}")
    if max_alpha >= 3:
        print(f"    *** CUBIC I.P. EXISTS AT n=8! ***")
    else:
        print(f"    CONFIRMED: alpha(Omega) <= 2 at n=8 (quadratic era)")

# ============================================================
# PART 2: Construct tournament at n=9 with 3 disjoint 3-cycles
# ============================================================
def part2():
    print("\n" + "=" * 70)
    print("PART 2: CONSTRUCT n=9 WITH 3 DISJOINT 3-CYCLES")
    print("=" * 70)

    # Build a tournament where vertices {0,1,2}, {3,4,5}, {6,7,8}
    # each form a directed 3-cycle
    n = 9
    A = np.zeros((n, n), dtype=int)

    # 3-cycle on {0,1,2}: 0->1->2->0
    A[0][1] = 1; A[1][2] = 1; A[2][0] = 1

    # 3-cycle on {3,4,5}: 3->4->5->3
    A[3][4] = 1; A[4][5] = 1; A[5][3] = 1

    # 3-cycle on {6,7,8}: 6->7->8->6
    A[6][7] = 1; A[7][8] = 1; A[8][6] = 1

    # Fill inter-block edges: block 0 beats block 1, block 1 beats block 2, block 2 beats block 0
    # (circulant block structure)
    blocks = [[0,1,2], [3,4,5], [6,7,8]]
    for i in blocks[0]:
        for j in blocks[1]:
            A[i][j] = 1
    for i in blocks[1]:
        for j in blocks[2]:
            A[i][j] = 1
    for i in blocks[2]:
        for j in blocks[0]:
            A[i][j] = 1

    # Verify it's a valid tournament
    for i in range(n):
        for j in range(i+1, n):
            assert A[i][j] + A[j][i] == 1, f"Not a tournament at ({i},{j})"

    print(f"\nConstructed tournament with 3 disjoint 3-cycles on {{0,1,2}}, {{3,4,5}}, {{6,7,8}}")

    # Compute independence polynomial
    coeffs = get_independence_poly(A)
    print(f"  Independence polynomial: I(x) = {' + '.join(f'{c}*x^{k}' for k, c in enumerate(coeffs))}")
    print(f"  Degree = {len(coeffs)-1}")

    # Verify H
    H = count_ham_paths(A)
    I_at_2 = sum(c * 2**k for k, c in enumerate(coeffs))
    print(f"  H(T) = {H}")
    print(f"  I(2) = {I_at_2}")
    print(f"  H = I(2)? {H == I_at_2}")

    if len(coeffs) > 3:
        print(f"\n  *** CUBIC TERM: a3 = {coeffs[3]} ***")
        print(f"  This is the FIRST tournament with a3 > 0!")
    elif len(coeffs) > 2 and coeffs[-1] > 0:
        print(f"  Still quadratic: a2 = {coeffs[2] if len(coeffs) > 2 else 0}")

    # Check with random n=9 tournaments
    print(f"\nRandom sampling at n=9:")
    rng = np.random.default_rng(42)
    cubic_count = 0
    total = 50
    for trial in range(total):
        A = random_tournament(n, rng)
        coeffs = get_independence_poly(A)
        deg = len(coeffs) - 1
        if deg >= 3:
            cubic_count += 1
            if cubic_count <= 3:
                print(f"  Trial {trial}: I(x) = {' + '.join(f'{c}*x^{k}' for k, c in enumerate(coeffs))}")
    print(f"\n  Cubic I.P. frequency at n=9: {cubic_count}/{total} ({100*cubic_count/total:.1f}%)")

# ============================================================
# PART 3: New coefficient a3
# ============================================================
def part3():
    print("\n" + "=" * 70)
    print("PART 3: THE NEW COEFFICIENT a3 AT n=9")
    print("=" * 70)

    rng = np.random.default_rng(2026_0314)
    n = 9
    N_samples = 50

    print(f"\nCollecting (a1, a2, a3) statistics at n=9 ({N_samples} samples):")
    a3_vals = defaultdict(int)
    all_data = []

    for trial in range(N_samples):
        A = random_tournament(n, rng)
        coeffs = get_independence_poly(A)
        a1 = coeffs[1] if len(coeffs) > 1 else 0
        a2 = coeffs[2] if len(coeffs) > 2 else 0
        a3 = coeffs[3] if len(coeffs) > 3 else 0
        H = count_ham_paths(A)
        I_at_2 = sum(c * 2**k for k, c in enumerate(coeffs))
        a3_vals[a3] += 1
        all_data.append((a1, a2, a3, H, I_at_2))

        if trial < 5 or a3 > 0:
            print(f"  Trial {trial}: a1={a1:3d}, a2={a2:3d}, a3={a3:2d}, H={H:5d}, I(2)={I_at_2:5d}")
            if H != I_at_2:
                print(f"    *** MISMATCH: H={H} != I(2)={I_at_2} ***")

        if (trial+1) % 25 == 0:
            print(f"  ... {trial+1}/{N_samples}")

    print(f"\na3 distribution:")
    for a3 in sorted(a3_vals.keys()):
        print(f"  a3 = {a3}: {a3_vals[a3]} ({100*a3_vals[a3]/N_samples:.1f}%)")

    # H formula with a3
    print(f"\nH formula for cubic I.P.:")
    print(f"  H = I(2) = 1 + 2*a1 + 4*a2 + 8*a3")
    print(f"  The 8 = 2^3 coefficient of a3 connects to:")
    print(f"    - Cayley-Dickson: dim 8 = octonions")
    print(f"    - Exponentiation: 8 = 2^3")
    print(f"    - Each independent triple of cycles contributes 2^3 = 8 to H")

# ============================================================
# PART 4: What H values become achievable with a3 > 0?
# ============================================================
def part4():
    print("\n" + "=" * 70)
    print("PART 4: NEW H VALUES FROM a3 > 0")
    print("=" * 70)

    print(f"\nH = 1 + 2*a1 + 4*a2 + 8*a3")
    print(f"  With a3=0: H = 1 + 2*a1 + 4*a2 (quadratic)")
    print(f"  With a3=1: H = 9 + 2*a1 + 4*a2")
    print(f"  With a3=k: H = 1+8k + 2*a1 + 4*a2")

    print(f"\nCan a3 > 0 unlock previously forbidden H values?")
    print(f"  H=21: need 1 + 2*a1 + 4*a2 + 8*a3 = 21")
    print(f"         2*a1 + 4*a2 + 8*a3 = 20")
    print(f"         a1 + 2*a2 + 4*a3 = 10")
    print(f"    With a3=1: a1 + 2*a2 = 6")
    print(f"      Possible: (6,0), (4,1), (2,2), (0,3)")
    print(f"    With a3=2: a1 + 2*a2 = 2")
    print(f"      Possible: (2,0), (0,1)")

    print(f"\n  So H=21 COULD be achievable at n>=9 if (a1,a2,a3)=(6,0,1) etc.")
    print(f"  is realizable. This would mean:")
    print(f"    6 odd cycles, 0 disjoint pairs, but 1 independent TRIPLE")
    print(f"  Wait: if there are 0 disjoint pairs but 1 independent triple,")
    print(f"  that's impossible! An independent triple contains 3 mutually")
    print(f"  disjoint cycles, which means at least C(3,2)=3 disjoint pairs.")
    print(f"  So a3 > 0 implies a2 >= C(a3 choose 1, 2) = C(a3,2)... no,")
    print(f"  actually a3 independent triples means a2 >= 3*a3")
    print(f"  (each triple contributes 3 disjoint pairs).")

    print(f"\n  Constraint: a3 = 1 => a2 >= 3")
    print(f"  For H=21 with a3=1: a1 + 2*a2 = 6, a2 >= 3")
    print(f"    a2=3: a1=0. Need 0 cycles but 3 disjoint pairs and 1 triple.")
    print(f"    IMPOSSIBLE: 0 cycles means 0 pairs and 0 triples.")
    print(f"  So ALL a3=1 solutions for H=21 are ruled out!")

    print(f"\n  For a3=2: a1 + 2*a2 = 2, a2 >= C(2,1)*3 = 6. But 2*a2 <= 2, so a2 <= 1.")
    print(f"    IMPOSSIBLE: a2 >= 6 and a2 <= 1 contradicts.")

    print(f"\n  CONCLUSION: H=21 remains forbidden even with cubic I.P.!")
    print(f"  The constraint a3 > 0 => a2 >= 3*a3 makes it even harder.")

    # Same analysis for H=7
    print(f"\nH=7 with a3 > 0:")
    print(f"  1 + 2*a1 + 4*a2 + 8*a3 = 7")
    print(f"  2*a1 + 4*a2 + 8*a3 = 6")
    print(f"  a1 + 2*a2 + 4*a3 = 3")
    print(f"  a3 >= 1: 4*a3 >= 4 > 3. IMPOSSIBLE.")
    print(f"  H=7 remains forbidden for ALL n (no room for a3).")

    # H=63
    print(f"\nH=63 with a3 > 0:")
    print(f"  a1 + 2*a2 + 4*a3 = 31")
    print(f"  a3=1, a2>=3: a1 + 2*a2 = 27, a2 in [3, 13]")
    print(f"    Many solutions. Whether realizable depends on n=9 cycle structure.")
    print(f"  So H=63 MIGHT become achievable at some n>=9.")

# ============================================================
# PART 5: Forbidden H at n=9
# ============================================================
def part5():
    print("\n" + "=" * 70)
    print("PART 5: FORBIDDEN H AT n=9 — DO GAPS CLOSE?")
    print("=" * 70)

    rng = np.random.default_rng(2026_0314)
    n = 9
    N_samples = 200  # Limited due to n=9 computational cost

    h_values = defaultdict(int)
    print(f"\nSampling H values at n=9 ({N_samples} tournaments)...")

    for trial in range(N_samples):
        A = random_tournament(n, rng)
        H = count_ham_paths(A)
        h_values[H] += 1
        if (trial+1) % 50 == 0:
            print(f"  {trial+1}/{N_samples} done")

    h_list = sorted(h_values.keys())
    print(f"\n  H range: [{min(h_list)}, {max(h_list)}]")
    print(f"  Distinct H values: {len(h_list)}")

    # Check the n=7 forbidden values
    forbidden_n7 = [7, 21, 63, 107, 119, 149]
    print(f"\n  Status of n=7 forbidden values at n=9:")
    for h in forbidden_n7:
        count = h_values.get(h, 0)
        print(f"    H={h:3d}: {'FOUND' if count > 0 else 'still absent'} ({count} times)")

    # Gaps in small H
    print(f"\n  All odd H values up to 50 and their status:")
    for h in range(1, 51, 2):
        count = h_values.get(h, 0)
        marker = "" if count > 0 else "  <-- GAP"
        if h <= 50:
            print(f"    H={h:3d}: {count:4d}{marker}")

# ============================================================
# PART 6: 9 = 3^2 and degree-3 connection
# ============================================================
def part6():
    print("\n" + "=" * 70)
    print("PART 6: 9 = 3^2 AND THE DEGREE-3 TRANSITION")
    print("=" * 70)

    print(f"\nThe number 9 = 3^2 controls TWO transitions:")
    print(f"  1. alpha(Omega) can reach 3 (cubic I.P.)")
    print(f"  2. Real roots of I(Omega, x) first fail (THM-025)")
    print(f"")
    print(f"Both are consequences of n = 3^2:")
    print(f"  Transition 1: 3 disjoint 3-cycles need 3*3 = 9 vertices")
    print(f"  Transition 2: claw-free fails, complex roots emerge")
    print(f"")
    print(f"The 'square of 3' interpretation:")
    print(f"  3^1 = 3: the 3-cycle (fundamental building block)")
    print(f"  3^2 = 9: three independent copies of the building block")
    print(f"  The exponentiation 3^2 literally means 'square the building block count'")

    print(f"\nCayley-Dickson alignment:")
    print(f"  dim = evaluation points needed = degree + 1")
    print(f"  Quadratic I.P. (n<=8): degree 2, need 3 points = 3^1")
    print(f"  Cubic I.P. (n=9+): degree 3, need 4 points = 2^2")
    print(f"  The quadratic -> cubic transition happens at n = 3^2 = 9")
    print(f"  Next structural transition: quartic I.P. needs 5 points")
    print(f"    4 disjoint 3-cycles need 12 vertices, so n >= 12")
    print(f"    Or 3 disjoint 5-cycles need 15 vertices")
    print(f"    Or mixed: 3-cycle + 3-cycle + 3-cycle + 3-cycle = 12")
    print(f"    Actually: 4 disjoint odd cycles, smallest is 3+3+3+3=12")
    print(f"    But also: 3+3+3+5=14, 3+3+5+5=16, etc.")
    print(f"    Minimum for alpha >= 4: n = 12 = 3*4 = 3*2^2")

    print(f"\n  Minimum n for alpha(Omega) >= k:")
    for k in range(1, 8):
        min_n = 3 * k  # k disjoint 3-cycles need 3k vertices
        print(f"    alpha >= {k}: n >= {min_n} (= 3*{k})")

    print(f"\n  So the 'Cayley-Dickson' thresholds for I.P. degree are:")
    print(f"    degree >= 1: n >= 3 = 3*1")
    print(f"    degree >= 2: n >= 6 = 3*2 = 2*3")
    print(f"    degree >= 3: n >= 9 = 3*3 = 3^2")
    print(f"    degree >= 4: n >= 12 = 3*4 = 2^2*3")
    print(f"    degree >= 5: n >= 15 = 3*5 = 3*5")
    print(f"    These are the MULTIPLES OF 3!")
    print(f"    Degree d requires n >= 3*d (or 3*(d+1) for independence number >= d+1)")

    print(f"\n  CORRECTION: alpha >= d requires d disjoint cycles,")
    print(f"  smallest is d copies of 3-cycles needing 3*d vertices.")
    print(f"  But I.P. has degree alpha, so degree >= d requires alpha >= d.")
    print(f"  So degree d I.P. first appears at n = 3*d.")
    print(f"  The THRESHOLD SEQUENCE is 3, 6, 9, 12, 15, ... = 3*N.")
    print(f"  This is an ARITHMETIC PROGRESSION with common difference 3!")

# ============================================================
# PART 7: The (a1, a2, a3) lattice at n=9
# ============================================================
def part7():
    print("\n" + "=" * 70)
    print("PART 7: THE (a1, a2, a3) LATTICE AT n=9")
    print("=" * 70)

    rng = np.random.default_rng(2026_0314)
    n = 9
    N_samples = 100

    triplets = defaultdict(int)
    print(f"\nCollecting (a1, a2, a3) at n=9 ({N_samples} samples)...")

    for trial in range(N_samples):
        A = random_tournament(n, rng)
        coeffs = get_independence_poly(A)
        a1 = coeffs[1] if len(coeffs) > 1 else 0
        a2 = coeffs[2] if len(coeffs) > 2 else 0
        a3 = coeffs[3] if len(coeffs) > 3 else 0
        triplets[(a1, a2, a3)] += 1

        if (trial+1) % 50 == 0:
            print(f"  {trial+1}/{N_samples}")

    print(f"\n  Distinct (a1, a2, a3) triplets: {len(triplets)}")

    # Show some with a3 > 0
    cubic_trips = {k: v for k, v in triplets.items() if k[2] > 0}
    print(f"  Triplets with a3 > 0: {len(cubic_trips)}")
    if cubic_trips:
        for (a1, a2, a3), count in sorted(cubic_trips.items())[:10]:
            H = 1 + 2*a1 + 4*a2 + 8*a3
            print(f"    (a1={a1:3d}, a2={a2:3d}, a3={a3:2d}): H={H:5d}, count={count}")

    # The constraint a3 > 0 => a2 >= 3*a3
    print(f"\n  Checking constraint a3 > 0 => a2 >= 3*a3:")
    violations = 0
    for (a1, a2, a3), count in triplets.items():
        if a3 > 0 and a2 < 3 * a3:
            violations += 1
            print(f"    VIOLATION: (a1={a1}, a2={a2}, a3={a3}): a2={a2} < 3*a3={3*a3}")
    if violations == 0:
        print(f"    No violations! a2 >= 3*a3 holds for all observed triplets.")

    # Show the H formula verification
    print(f"\n  H = 1 + 2*a1 + 4*a2 + 8*a3 verification:")
    rng2 = np.random.default_rng(99)
    all_match = True
    for trial in range(20):
        A = random_tournament(n, rng2)
        H = count_ham_paths(A)
        coeffs = get_independence_poly(A)
        I_at_2 = sum(c * 2**k for k, c in enumerate(coeffs))
        if H != I_at_2:
            print(f"    MISMATCH trial {trial}: H={H}, I(2)={I_at_2}")
            all_match = False
    if all_match:
        print(f"    All 20 trials match: H = I(2) confirmed at n=9")

    print(f"\n{'='*70}")
    print(f"END OF CUBIC TRANSITION EXPLORATION")
    print(f"{'='*70}")


def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()

if __name__ == "__main__":
    main()
