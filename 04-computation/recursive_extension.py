#!/usr/bin/env python3
"""
opus-2026-05-27-S7: Recursive structure exploitation for fast sequence extension.

Key theorems used:
  THM-332 (Dominant Vertex): H(T) = H(T-Q) for absolute source Q
  THM-331 (King Increment): H(T) - H(T-Q) >= 2|N^-(Q)|
  Claim A / OCF: H(T) - H(T-v) = 2 * #{directed odd cycles through v}

Core idea: inserting vertex v_n into T_{n-1} changes H by exactly
  DeltaH = 2 * #{odd simple cycles through v_n in T_n}
         = 2 * #{odd-length simple paths from N+(v_n) to N-(v_n) in T_{n-1}}

So: a(n) >= H(T*_{n-1}) + 2 * max_{C subset V} OddPaths(T*_{n-1}, C, V\C)
where T*_{n-1} is the best known (n-1)-tournament.

This converts "find optimal n-tournament" to "find optimal 2-partition of (n-1)-tournament".
"""

import sys
import time
from itertools import combinations

# ============================================================
# PART 0: Tournament primitives
# ============================================================

def count_hp(n, adj):
    """Count Hamiltonian paths via bitmask DP. adj[v] = bitmask of vertices v beats."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    total = 0
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            d = dp[mask][v]
            if not d:
                continue
            if mask == full:
                total += d
                continue
            tmp = adj[v] & ~mask
            while tmp:
                lsb = tmp & (-tmp)
                w = lsb.bit_length() - 1
                dp[mask | lsb][w] += d
                tmp ^= lsb
    return total

def build_Tk(k):
    """
    Build all-0 staircase T_k on n=2k vertices.
    Base path:  n-1 -> n-2 -> ... -> 0  (each v beats v-1)
    All tiles (x,y) with x-y >= 2: arc y -> x (UPWARD, lower beats higher)
    Returns (n, adj) where adj[v] = bitmask of vertices v beats.
    """
    n = 2 * k
    adj = [0] * n
    for v in range(1, n):
        adj[v] |= (1 << (v - 1))           # base path: v -> v-1
    for x in range(2, n):
        for y in range(x - 1):             # x - y >= 2
            adj[y] |= (1 << x)             # y -> x (upward tile)
    return n, adj

def adj_from_str(n, s):
    """Parse adjacency from space-separated decimal bitmasks."""
    return list(map(int, s.split()))

def is_valid(n, adj):
    """Check tournament validity: exactly one arc per pair, no self-loops."""
    for i in range(n):
        if adj[i] >> i & 1:
            return False
        for j in range(i + 1, n):
            a = (adj[i] >> j) & 1
            b = (adj[j] >> i) & 1
            if a + b != 1:
                return False
    return True

# ============================================================
# PART 1: Staircase H values + Berlekamp-Massey
# ============================================================

def berlekamp_massey(seq):
    """Find minimal LFSR for sequence. Returns coefficients c s.t. s[n] = sum c[i]*s[n-1-i]."""
    n = len(seq)
    C, B = [1], [1]
    L, m, b = 0, 1, 1
    for i in range(n):
        d = seq[i]
        for j in range(1, L + 1):
            d += C[j] * seq[i - j]
        if d == 0:
            m += 1
            continue
        T = C[:]
        coef = -d * pow(b, -1, None)  # rational arithmetic
        # For integer sequences, do this with fractions
        from fractions import Fraction
        d_f = Fraction(d)
        b_f = Fraction(b)
        coef_f = -d_f / b_f
        if len(C) < len(B) + m:
            C += [Fraction(0)] * (len(B) + m - len(C))
        for j in range(len(B)):
            if j + m < len(C):
                C[j + m] += coef_f * B[j]
        if 2 * L <= i:
            L = i + 1 - L
            B = T
            b = d
            m = 1
        else:
            m += 1
    return C, L

def find_recurrence(seq, max_order=8):
    """Find minimal linear recurrence for integer sequence using fraction arithmetic."""
    from fractions import Fraction
    n = len(seq)
    S = [Fraction(x) for x in seq]

    for order in range(1, min(max_order + 1, n // 2 + 1)):
        # Try to find c[0],...,c[order-1] s.t. s[k] = sum c[i]*s[k-1-i]
        # Linear system: use first 2*order equations
        if n < 2 * order:
            break
        A = []
        b_vec = []
        for k in range(order, n):
            A.append([S[k - 1 - i] for i in range(order)])
            b_vec.append(S[k])

        # Solve overdetermined system via Gaussian elimination on first `order` rows
        mat = [A[i] + [b_vec[i]] for i in range(order)]
        # Gaussian elimination
        sol = None
        try:
            for col in range(order):
                pivot = None
                for row in range(col, order):
                    if mat[row][col] != 0:
                        pivot = row
                        break
                if pivot is None:
                    break
                mat[col], mat[pivot] = mat[pivot], mat[col]
                scale = mat[col][col]
                mat[col] = [x / scale for x in mat[col]]
                for row in range(order):
                    if row != col and mat[row][col] != 0:
                        factor = mat[row][col]
                        mat[row] = [mat[row][j] - factor * mat[col][j] for j in range(order + 1)]
            sol = [mat[i][order] for i in range(order)]
        except Exception:
            continue

        # Verify against ALL terms
        valid = True
        for k in range(order, n):
            pred = sum(sol[i] * S[k - 1 - i] for i in range(order))
            if pred != S[k]:
                valid = False
                break
        if valid:
            return order, [int(c) if c == int(c) else c for c in sol]
    return None, None

def staircase_extension():
    """Compute H(T_k) for k=1..8 and find recurrence."""
    print("=" * 60)
    print("PART 1: Staircase H values + recurrence detection")
    print("=" * 60)

    H_vals = []
    for k in range(1, 9):
        n, adj = build_Tk(k)
        t0 = time.time()
        h = count_hp(n, adj)
        elapsed = time.time() - t0
        H_vals.append(h)
        # Compute odd cycle count via 3-cycle formula (THM-321)
        # alpha_1 = #{odd cycles} in T_k is just sum of independent set size 1
        print(f"  H(T_{k}) = {h:>10}  (n={n}, {elapsed:.2f}s)")

    print(f"\nSequence: {H_vals}")
    print(f"Ratios:   {[round(H_vals[i+1]/H_vals[i], 3) for i in range(len(H_vals)-1)]}")

    # Find recurrence
    order, coeffs = find_recurrence(H_vals)
    if order:
        print(f"\nLinear recurrence found! Order {order}")
        print(f"  H(T_k) = " + " + ".join(f"{c}*H(T_{{k-{i+1}}})" for i, c in enumerate(coeffs)))
        # Verify and extend
        extended = H_vals[:]
        for k in range(len(H_vals), 30):
            next_val = sum(int(coeffs[i]) * extended[k - 1 - i] for i in range(order))
            extended.append(next_val)
        print(f"\nExtended sequence (k=1..30):")
        for i, v in enumerate(extended):
            print(f"  H(T_{i+1}) = {v}")
    else:
        print(f"\nNo recurrence found up to order 8 with {len(H_vals)} terms.")
        print("Need more terms (k=9,10) to detect higher-order recurrence.")

    return H_vals

# ============================================================
# PART 2: Dominant vertex analysis
# ============================================================

def dominant_vertices(n, adj):
    """Find absolute sources (d+=n-1) and sinks (d+=0) in T."""
    sources = [v for v in range(n) if bin(adj[v]).count('1') == n - 1]
    sinks = [v for v in range(n) if adj[v] == 0]
    return sources, sinks

def reduce_dominant(n, adj):
    """Repeatedly strip dominant vertices. Returns (reduced_n, reduced_adj, stripped_count).
    By THM-332: H is unchanged throughout."""
    adj = list(adj)
    stripped = 0
    while True:
        sources = [v for v in range(n) if bin(adj[v]).count('1') == n - 1]
        sinks = [v for v in range(n) if adj[v] == 0]
        if not sources and not sinks:
            break
        # Remove all sources and sinks
        to_remove = set(sources) | set(sinks)
        mapping = [v for v in range(n) if v not in to_remove]
        new_adj = []
        for v in mapping:
            new_mask = 0
            for w_idx, w in enumerate(mapping):
                if adj[v] >> w & 1:
                    new_mask |= (1 << w_idx)
            new_adj.append(new_mask)
        stripped += len(to_remove)
        n = len(mapping)
        adj = new_adj
        if n == 0:
            break
    return n, adj, stripped

def dominant_reduction_analysis():
    """Analyze dominant vertex reduction on staircase and known optimal tournaments."""
    print("\n" + "=" * 60)
    print("PART 2: Dominant vertex reduction analysis")
    print("=" * 60)

    # Test on staircase
    for k in range(1, 7):
        n, adj = build_Tk(k)
        rn, radj, stripped = reduce_dominant(n, adj)
        h_orig = count_hp(n, adj)
        h_red = count_hp(rn, radj) if rn > 0 else 1
        print(f"  T_{k} (n={n}): H={h_orig}, stripped={stripped} vertices, "
              f"core size={rn}, H_core={h_red}, match={h_orig==h_red}")

    # Test on best a(12) tournament from THM-329
    n12 = 12
    adj12_str = "3472 1585 2187 371 484 2629 1671 2858 3174 1309 2236 602"
    adj12 = adj_from_str(n12, adj12_str)
    if is_valid(n12, adj12):
        h12 = count_hp(n12, adj12)
        rn, radj, stripped = reduce_dominant(n12, adj12)
        h_red = count_hp(rn, radj) if rn > 0 else 1
        print(f"\n  Best a(12) tournament: H={h12}, stripped={stripped}, core={rn}, H_core={h_red}")

    # Test on best a(13) tournament
    n13 = 13
    adj13_str = "3666 5352 6411 3009 686 3213 6708 5701 243 3366 4444 5522 825"
    adj13 = adj_from_str(n13, adj13_str)
    if is_valid(n13, adj13):
        h13 = count_hp(n13, adj13)
        rn, radj, stripped = reduce_dominant(n13, adj13)
        h_red = count_hp(rn, radj) if rn > 0 else 1
        print(f"  Best a(13) tournament: H={h13}, stripped={stripped}, core={rn}, H_core={h_red}")

# ============================================================
# PART 3: Optimal vertex insertion via F_odd precomputation
# ============================================================

def compute_F_odd(n, adj):
    """
    For each pair (s, e) in {0,...,n-1}^2, compute:
      F_odd[s][e] = number of simple directed paths of ODD length from s to e in T.

    Uses bitmask DP:
      dp[mask][v] = number of simple directed paths starting at some
                    fixed source s, visiting vertices in mask, ending at v,
                    with odd path length (1, 3, 5, ...).

    This is computed separately for each starting vertex s.
    Total cost: O(n * 2^n * n) — same order as HP counting.
    """
    F_odd = [[0] * n for _ in range(n)]

    for s in range(n):
        # dp_even[mask][v] = # paths from s to v through exactly 'mask', even length (0,2,4,...)
        # dp_odd[mask][v]  = # paths from s to v through exactly 'mask', odd length (1,3,5,...)
        # Initialize: path of length 0 from s to s
        dp_even = [[0] * n for _ in range(1 << n)]
        dp_odd  = [[0] * n for _ in range(1 << n)]
        dp_even[1 << s][s] = 1

        for mask in range(1, 1 << n):
            if not (mask >> s & 1):
                continue
            for v in range(n):
                if not (mask >> v & 1):
                    continue
                de = dp_even[mask][v]
                do = dp_odd[mask][v]
                if not de and not do:
                    continue
                # Extend to unvisited neighbors
                tmp = adj[v] & ~mask
                while tmp:
                    lsb = tmp & (-tmp)
                    w = lsb.bit_length() - 1
                    new_mask = mask | lsb
                    dp_odd[new_mask][w]  += de   # even + 1 step = odd
                    dp_even[new_mask][w] += do   # odd + 1 step = even
                    tmp ^= lsb

        # Collect: F_odd[s][e] = sum over all masks of dp_odd[mask][e]
        for e in range(n):
            total = 0
            for mask in range(1, 1 << n):
                if mask >> e & 1 and mask >> s & 1:
                    total += dp_odd[mask][e]
            F_odd[s][e] = total

    return F_odd

def optimal_insertion(n_base, adj_base):
    """
    Given T_{n_base}, find the optimal orientation for inserting vertex n_base
    to maximize H(T_{n_base+1}) = H(T_{n_base}) + 2 * #{odd cycles through new vertex}.

    #{odd cycles through v_new | Court=C} = sum_{c in C, r not in C} F_odd[c][r]

    For n_base <= 14, bruteforce all 2^{n_base} partitions.
    Returns: (best_score, best_court, adj_new)
    """
    print(f"\n  Precomputing F_odd for n={n_base} tournament...", flush=True)
    t0 = time.time()
    F = compute_F_odd(n_base, adj_base)
    print(f"  F_odd computed in {time.time()-t0:.1f}s")

    # For each partition C (vertices v_new beats) vs R (vertices v_new loses to):
    best_score = -1
    best_court = None

    for court_mask in range(1 << n_base):
        score = 0
        for c in range(n_base):
            if not (court_mask >> c & 1):
                continue
            for r in range(n_base):
                if court_mask >> r & 1:
                    continue
                score += F[c][r]
        if score > best_score:
            best_score = score
            best_court = court_mask

    # Build new n_base+1 tournament
    n_new = n_base + 1
    adj_new = list(adj_base) + [0]
    v_new = n_base
    for v in range(n_base):
        if best_court >> v & 1:
            adj_new[v_new] |= (1 << v)   # v_new beats v (v in court)
        else:
            adj_new[v] |= (1 << v_new)   # v beats v_new (v in rivals)

    h_base = count_hp(n_base, adj_base)
    h_new = h_base + 2 * best_score

    return best_score, best_court, adj_new, h_base, h_new

def insertion_analysis():
    """Apply optimal insertion to best known tournaments to get new lower bounds."""
    print("\n" + "=" * 60)
    print("PART 3: Optimal vertex insertion analysis")
    print("=" * 60)

    # Known optimal tournaments from a038375 (best found)
    best_tournaments = {
        # n: (H, adj_str)
        4: (5, "5 9 2 4"),      # score [1,1,2,2]
        5: (15, "14 11 7 13 9"), # QR_5? Let me verify
    }

    # Build QR_5 (Paley on 5 vertices - but 5 ≡ 1 mod 4, NOT valid Paley!)
    # For n=5: QR_5 invalid. The optimal is the regular tournament.
    # Use the exhaustive n=5 result: all regular (score 2222) have H=15.
    # Build the unique (up to iso) regular 5-tournament.

    # Circulant(5, {1,2}): vertex i beats (i+1)%5, (i+2)%5
    n5 = 5
    adj5 = [0] * 5
    for i in range(5):
        adj5[i] |= (1 << ((i+1)%5))
        adj5[i] |= (1 << ((i+2)%5))
    h5 = count_hp(n5, adj5)
    print(f"\n  Circulant(5,{{1,2}}): H={h5}, valid={is_valid(n5,adj5)}")

    # Apply optimal insertion: n=5 -> n=6
    if h5 == 15 and is_valid(n5, adj5):
        t0 = time.time()
        score, court, adj6, h5_base, h6 = optimal_insertion(n5, adj5)
        print(f"  Optimal insertion n=5->6: odd_cycles={score}, H(T6)={h6}, "
              f"a(6)={45} {'✓' if h6==45 else '?'}, ({time.time()-t0:.1f}s)")

    # For n=7: use QR_7 (Paley, 7 ≡ 3 mod 4)
    n7 = 7
    adj7 = [0] * 7
    qr7 = {1, 2, 4}  # QR mod 7
    for i in range(7):
        for j in range(7):
            if i != j and (j - i) % 7 in qr7:
                adj7[i] |= (1 << j)
    h7 = count_hp(n7, adj7)
    print(f"\n  QR_7 (Paley): H={h7}, valid={is_valid(n7,adj7)}, expected=189")

    if h7 == 189 and is_valid(n7, adj7):
        t0 = time.time()
        score, court, adj8, h7_base, h8 = optimal_insertion(n7, adj7)
        print(f"  Optimal insertion n=7->8: odd_cycles={score}, H(T8)={h8}, "
              f"a(8)=661, ratio={h8/661:.3f}, ({time.time()-t0:.1f}s)")

        # Continue: n=8 -> n=9
        if h8 >= 600:
            t0 = time.time()
            score2, court2, adj9, h8_base, h9 = optimal_insertion(8, adj8)
            print(f"  Optimal insertion n=8->9: odd_cycles={score2}, H={h9}, "
                  f"a(9)=3357, ratio={h9/3357:.3f}, ({time.time()-t0:.1f}s)")

    # Analyze the BEST a(11) tournament: QR_11 (Paley)
    n11 = 11
    adj11 = [0] * 11
    # QR mod 11: 1,3,4,5,9
    qr11 = {1, 3, 4, 5, 9}
    for i in range(11):
        for j in range(11):
            if i != j and (j - i) % 11 in qr11:
                adj11[i] |= (1 << j)
    h11 = count_hp(n11, adj11)
    print(f"\n  QR_11 (Paley): H={h11}, valid={is_valid(n11,adj11)}, expected=95095")

    if h11 == 95095 and is_valid(n11, adj11):
        print(f"  Applying optimal insertion n=11->12...")
        t0 = time.time()
        score, court, adj12, h11_base, h12 = optimal_insertion(n11, adj11)
        elapsed = time.time() - t0
        print(f"  Optimal insertion n=11->12: odd_cycles={score}, H(T12)={h12}")
        print(f"  Current best a(12)=531205, our bound={h12}, ratio={h12/531205:.4f}")
        print(f"  Time: {elapsed:.1f}s")

        return adj12, h12

# ============================================================
# PART 4: Recursive lower bound table
# ============================================================

def compute_recursive_bounds():
    """
    Build a table of lower bounds using the optimal insertion chain.
    Start from QR_7 (exact) and insert upward.
    Compare with known a(n) values.
    """
    print("\n" + "=" * 60)
    print("PART 4: Recursive insertion lower bound chain")
    print("=" * 60)

    known = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}

    # Build QR_7
    n = 7
    adj = [0] * n
    qr7 = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%7 in qr7:
                adj[i] |= (1 << j)
    h = count_hp(n, adj)
    assert h == 189

    print(f"\n  n={n}: H={h} (QR_7, exact)")

    results = [(n, h, h, 1.0)]

    for target_n in range(8, 15):
        t0 = time.time()
        score, court, adj_new, h_prev, h_new = optimal_insertion(n, adj)
        elapsed = time.time() - t0

        known_h = known.get(target_n, None)
        ratio = h_new / known_h if known_h else None
        ratio_str = f"{ratio:.4f}" if ratio else "?"

        print(f"  n={target_n}: H_insert={h_new}, known={known_h or '?'}, "
              f"ratio={ratio_str}, score={score}, ({elapsed:.1f}s)")

        results.append((target_n, h_new, known_h, ratio))
        n = target_n
        adj = adj_new
        h = h_new

    return results

# ============================================================
# PART 5: Score sequence analysis — how dominant vertices affect search
# ============================================================

def analyze_score_vs_H():
    """
    For each n, compute correlation between score gap and H.
    Key insight: dominant vertex (gap = n-1) means H = H(T-v) [unchanged].
    This lets us PRUNE search space.
    """
    print("\n" + "=" * 60)
    print("PART 5: Score sequence and H reduction")
    print("=" * 60)

    for n in range(3, 8):
        # Enumerate all tournaments via bitmask
        num_pairs = n*(n-1)//2
        pairs = [(i, j) for i in range(n) for j in range(i+1, n)]

        H_by_gap = {}
        H_dominant = []
        H_nondominant = []

        for bits in range(1 << num_pairs):
            adj = [0] * n
            for idx, (i, j) in enumerate(pairs):
                if bits >> idx & 1:
                    adj[i] |= (1 << j)
                else:
                    adj[j] |= (1 << i)

            degrees = [bin(adj[v]).count('1') for v in range(n)]
            gap = max(degrees) - min(degrees)
            h = count_hp(n, adj)

            H_by_gap.setdefault(gap, []).append(h)
            if max(degrees) == n-1:
                H_dominant.append(h)
            else:
                H_nondominant.append(h)

        max_H = max(max(v) for v in H_by_gap.values())
        dominant_max = max(H_dominant) if H_dominant else 0
        nondominant_max = max(H_nondominant) if H_nondominant else 0

        print(f"\n  n={n}: max_H={max_H}")
        print(f"    With dominant vertex: max_H={dominant_max} "
              f"({'= global max' if dominant_max == max_H else '< global max'})")
        print(f"    Without dominant vertex: max_H={nondominant_max} "
              f"({'= global max' if nondominant_max == max_H else '< global max'})")
        print(f"    Gap distribution: " +
              ", ".join(f"gap={g}: max={max(v)}" for g, v in sorted(H_by_gap.items())))

# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Recursive Structure Analysis for Fast Sequence Extension")
    print("Session: opus-2026-05-27-S7")
    print("=" * 60)

    # Part 1: Staircase recurrence
    H_vals = staircase_extension()

    # Part 2: Dominant vertex reduction
    dominant_reduction_analysis()

    # Part 3: Optimal insertion
    result = insertion_analysis()

    # Part 4: Full recursive bound chain
    compute_recursive_bounds()

    # Part 5: Score analysis
    analyze_score_vs_H()

    print("\n\nDone.")
