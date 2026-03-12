"""
cycle_count_H_formula.py — Derive H = f(c_3, c_5, ..., c_p) for circulant tournaments.

CORE THEOREM (from previous session):
  Cycle counts (c_3, c_5, ..., c_p) DETERMINE H for circulant tournaments at prime p.
  This was verified computationally at p = 7, 11, 13.

THIS SCRIPT:
  1. Compute H and all cycle counts for circulant tournaments at p = 3, 5, 7, 11, 13
  2. Find the EXACT polynomial formula H = f(c_3, c_5, ..., c_p)
  3. Test whether the formula extends to non-circulant tournaments

The formula H = f(c_k) would be a MAJOR result because:
  - OCF says H = I(Omega, 2) where Omega depends on cycle ARRANGEMENT, not just counts
  - If cycle COUNTS suffice for circulant tournaments, the Z_p symmetry forces
    a unique conflict pattern for each count tuple
  - This would give a SPECTRAL formula for H (since c_k = tr(A^k)/k for eigenvalues)

KNOWN IDENTITIES:
  For ANY n-tournament: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... (OCF)
  where alpha_k = # independent k-sets in Omega(T)

  For n=7 regular: H = 1 + 2*alpha_1 + 4*alpha_2 (only 3-cycles, no 5 or 7)
  Wait no — Omega includes ALL odd cycles, not just 3-cycles.

  Actually Omega(T) = conflict graph of directed 3-cycles only, by definition.
  H = I(Omega, 2) = sum alpha_k * 2^k.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import time
import cmath
import math
from collections import defaultdict

sys.path.insert(0, '04-computation')
from satake_ndrt_h import ham_count_circulant, all_circulant_H


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def cycle_counts_from_eigenvalues(p, eigs):
    """Compute c_k = # directed k-cycles for all odd k from 3 to p."""
    counts = {}
    for k in range(3, p + 1, 2):
        tr_k = sum(e ** k for e in eigs)
        # For tournament: tr(A^k) counts closed walks of length k
        # c_k = (1/k) * (tr(A^k) - correction terms for non-prime-cycle walks)
        # For k=3: c_3 = tr(A^3) / 3 (no correction needed since 3 is prime)
        # For k=5: c_5 = (tr(A^5) - 5*c_3*(n-3)) / 5... no, this is wrong
        # Actually tr(A^k) = sum of A[i1,i2]*A[i2,i3]*...*A[ik,i1] over all i1,...,ik
        # This counts ALL closed walks, not just simple cycles
        # For small k (3,5,7), we need to subtract non-simple walks

        # For k=3: All length-3 closed walks in a tournament are simple (can't revisit)
        # because A[i,i]=0 and if you visit only 2 vertices, A[i,j]*A[j,i] = 0 (tournament)
        # So c_3 = tr(A^3) / 3
        if k == 3:
            counts[k] = round(tr_k.real / k)
        else:
            # For k >= 5, tr(A^k) includes non-simple closed walks
            # But for circulant tournaments we can compute c_k directly
            counts[k] = None  # Will compute differently
    return counts


def compute_cycle_counts_direct(p, S):
    """Compute all directed cycle counts c_k directly by enumeration."""
    S_set = set(S)
    adj = {}
    for i in range(p):
        for j in range(p):
            if i != j:
                adj[(i, j)] = (j - i) % p in S_set

    # Count directed k-cycles for each k
    # For circulant on Z_p, c_k = p * (# k-cycles through vertex 0) / k
    # But let's just count them all

    counts = {}

    # c_3: directed triangles
    c3 = 0
    for i in range(p):
        for j in range(i + 1, p):
            for k in range(j + 1, p):
                # Check all 3! = 6 orientations, but only 2 directed cycles possible
                if adj[(i, j)] and adj[(j, k)] and adj[(k, i)]:
                    c3 += 1
                if adj[(i, k)] and adj[(k, j)] and adj[(j, i)]:
                    c3 += 1
    counts[3] = c3

    # c_5: directed 5-cycles (enumerate 5-element subsets)
    if p >= 5:
        from itertools import combinations, permutations
        c5 = 0
        for subset in combinations(range(p), 5):
            # Count directed Hamiltonian cycles on this 5-vertex subtournament
            for perm in permutations(subset):
                is_cycle = True
                for idx in range(5):
                    if not adj[(perm[idx], perm[(idx + 1) % 5])]:
                        is_cycle = False
                        break
                if is_cycle:
                    c5 += 1
            # Each 5-cycle counted 5 times (cyclic rotations)
        counts[5] = c5 // 5

    # c_7: only for p >= 7
    if p >= 7:
        from itertools import combinations, permutations
        c7 = 0
        for subset in combinations(range(p), 7):
            for perm in permutations(subset):
                is_cycle = True
                for idx in range(7):
                    if not adj[(perm[idx], perm[(idx + 1) % 7])]:
                        is_cycle = False
                        break
                if is_cycle:
                    c7 += 1
            # Each 7-cycle counted 7 times
        counts[7] = c7 // 7

    return counts


def compute_c3_c5_fast(p, S):
    """Fast c_3 and c_5 computation for circulant tournament."""
    S_set = set(S)

    # c_3: for each triple {i,j,k}, check if it forms a directed 3-cycle
    c3 = 0
    for i in range(p):
        for j in range(i + 1, p):
            for k in range(j + 1, p):
                d_ij = (j - i) % p in S_set
                d_jk = (k - j) % p in S_set
                d_ki = (i - k) % p in S_set
                d_ji = (i - j) % p in S_set
                d_kj = (j - k) % p in S_set
                d_ik = (k - i) % p in S_set
                if d_ij and d_jk and d_ki:
                    c3 += 1
                if d_ik and d_kj and d_ji:
                    c3 += 1

    # c_5: use trace formula with correction
    # tr(A^5) = sum lambda_k^5 = 5*c_5 + 5*(n-3)*c_3 + ...
    # Actually for tournaments, the correction is more complex
    # Let's use eigenvalue formula: tr(A^5) counts all closed 5-walks
    # A closed 5-walk can be:
    #   (a) A simple 5-cycle (counted 5 times each in tr, total 5*c_5)
    #   (b) A "lollipop": 3-cycle + back-forth (but A[i,j]*A[j,i]=0 for tournaments!)
    # Actually in a tournament, since A[i,j]*A[j,i] = 0, any closed walk
    # that backtracks is 0. So ALL closed walks in tr(A^k) are simple!
    # Wait, not simple — they can revisit vertices via different edges.

    # For k=5 in a tournament:
    # tr(A^5) = sum_{i} (A^5)_{ii} = sum of products along 5-step closed walks
    # A walk i1->i2->i3->i4->i5->i1 where each step is a tournament edge
    # This CAN revisit vertices: e.g. i1->i2->i3->i1->i2->i3 is a valid walk
    # because i1->i2, i2->i3, i3->i1, i1->i2, i2->i3 could all be edges
    # Wait, i5->i1 means i5=i3, so i3->i1 must be an edge (which it is in the 3-cycle)
    # So this walk is i1->i2->i3->i1->i2->i3 with i5=i3, which IS a valid walk
    # But the walk visits {i1,i2,i3} and is NOT a simple 5-cycle.

    # So tr(A^5) != 5*c_5 in general. We need to correct.

    # Actually for circulant tournaments we can just compute c_5 directly:
    # c_5 = (p * local_c5) where local_c5 = # directed 5-cycles through vertex 0
    # And local_c5 = (1/5) * sum over 4-subsets containing 0 of Hamiltonian 5-cycle count

    # Faster: use vertex 0's neighborhood
    c5 = 0
    from itertools import combinations, permutations
    # Only count 5-cycles containing vertex 0 (then multiply by p/5)
    others = list(range(1, p))
    for subset in combinations(others, 4):
        verts = (0,) + subset
        for perm in permutations(verts):
            if perm[0] != 0:
                continue  # fix starting vertex to avoid rotation counting
            is_cycle = True
            for idx in range(5):
                a, b = perm[idx], perm[(idx + 1) % 5]
                if (b - a) % p not in S_set:
                    is_cycle = False
                    break
            if is_cycle:
                c5 += 1
    # Each 5-cycle through 0 counted once (fixed start at 0, but 2 directions)
    # Actually permutations of 5 elements with first fixed = 4! = 24
    # A 5-cycle has 5 rotations * 1 direction = 5 representations among 5! perms
    # With first vertex fixed: 1 representation per cycle
    # But we also count the reverse cycle (if it exists) — tournament excludes reverse
    c5 = c5 * p // 5  # p vertices, each 5-cycle contains 5 vertices

    return c3, c5


def all_circulant_tournaments(n):
    """Generate all circulant tournament connection sets on Z_n."""
    k = (n - 1) // 2
    pairs = []
    used = set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b:
                return []
            pairs.append((a, b))
            used.add(a)
            used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits >> i) & 1 else b)
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    """Held-Karp DP for Hamiltonian path count."""
    S_set = set(S)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i][j] = True

    full_mask = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]

    return sum(dp[full_mask][v] for v in range(n))


def spectral_y_squared(n, S):
    eigs = circulant_eigenvalues(n, S)
    m = (n - 1) // 2
    return [eigs[k].imag ** 2 for k in range(1, m + 1)]


# ======================================================================
# MAIN ANALYSIS
# ======================================================================

def main():
    print("=" * 70)
    print("CYCLE COUNT DETERMINATION OF H FOR CIRCULANT TOURNAMENTS")
    print("=" * 70)

    for p in [3, 5, 7, 11, 13]:
        print(f"\n{'=' * 60}")
        print(f"p = {p} (mod 4 = {p % 4})")
        print(f"{'=' * 60}")

        all_S = all_circulant_tournaments(p)
        m = (p - 1) // 2

        t0 = time.time()
        data = []
        for S in all_S:
            H = ham_count_dp(p, S)

            # Compute c_3 via counting
            S_set = set(S)
            c3 = 0
            for i in range(p):
                for j in range(i + 1, p):
                    for k in range(j + 1, p):
                        # Check both orientations
                        if all((S_set.__contains__)((b - a) % p) for a, b in [(i,j),(j,k),(k,i)]):
                            c3 += 1
                        if all((S_set.__contains__)((b - a) % p) for a, b in [(i,k),(k,j),(j,i)]):
                            c3 += 1

            # c_3 from eigenvalue trace
            eigs = circulant_eigenvalues(p, S)
            tr3 = sum(e ** 3 for e in eigs).real
            c3_trace = round(tr3 / 3)

            # Higher cycle counts from traces (corrected)
            # For tournaments: tr(A^k) counts ALL k-step closed walks
            # A closed walk on a tournament can revisit vertices
            # We need the CORRECTED formula for simple directed cycles

            # Actually, use Waring's formula / Newton identities on the
            # characteristic polynomial to get cycle counts
            # But it's simpler to just count directly for small p

            y2 = spectral_y_squared(p, S)

            data.append({
                'S': S, 'H': H, 'c3': c3, 'c3_trace': c3_trace,
                'y2': y2,
            })

        elapsed = time.time() - t0

        # Sort by H
        data.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(d['H'] for d in data), reverse=True)

        print(f"  {len(all_S)} tournaments computed in {elapsed:.1f}s")
        print(f"  {len(H_vals)} distinct H values")

        # Check: c_3 from trace vs direct count
        mismatches = sum(1 for d in data if d['c3'] != d['c3_trace'])
        print(f"  c_3 trace vs direct: {mismatches} mismatches")

        # All circulant tournaments on Z_p have the SAME c_3!
        c3_vals = set(d['c3'] for d in data)
        print(f"  c_3 values: {sorted(c3_vals)}")
        if len(c3_vals) == 1:
            print(f"  ALL circulant tournaments at p={p} have c_3 = {c3_vals.pop()}")
            print(f"  => c_3 does NOT distinguish H values. Need higher cycles.")

        # Table
        print(f"\n  {'H':>12} {'count':>5} {'c3':>5} {'y2 (spectral params)'}")
        print(f"  {'-' * 12} {'-' * 5} {'-' * 5} {'-' * 40}")
        for H_val in H_vals:
            matches = [d for d in data if d['H'] == H_val]
            d = matches[0]
            y2_str = ', '.join(f'{y:.4f}' for y in d['y2'])
            print(f"  {H_val:>12} {len(matches):>5} {d['c3']:>5} [{y2_str}]")

        # For p=7: compute alpha_1, alpha_2 from Omega
        if p <= 7:
            print(f"\n  OMEGA ANALYSIS (conflict graph of 3-cycles):")
            for d in data[:len(H_vals)]:  # one per H value
                S_set = set(d['S'])
                # Find all directed 3-cycles
                cycles = []
                for i in range(p):
                    for j in range(i + 1, p):
                        for k in range(j + 1, p):
                            if all((b - a) % p in S_set for a, b in [(i,j),(j,k),(k,i)]):
                                cycles.append(frozenset({i, j, k}))
                            if all((b - a) % p in S_set for a, b in [(i,k),(k,j),(j,i)]):
                                cycles.append(frozenset({i, j, k}))

                # Remove duplicates (same vertex set, different orientation —
                # but for 3-cycle vertex set each has exactly ONE directed cycle)
                # Actually wait: a 3-element subset either has 0 or 2 directed cycles
                # No! In a tournament on 3 vertices, there's either 0 or 1 directed 3-cycle
                # (the tournament is either transitive or cyclic)
                # For cyclic: there are 2 directed 3-cycles (one in each direction)
                # But in a tournament, only ONE direction exists!
                # So each {i,j,k} contributes at most 1 directed 3-cycle.

                # Hmm but my counting above found both orientations.
                # In a tournament: for each triple, exactly one of the two cyclic
                # orderings is consistent with the edges. If the triple is
                # transitively ordered, NEITHER cyclic ordering works.

                # So "cycles" may have duplicates from the same vertex set
                # Actually no: the frozenset is the same, but I checked both orderings
                # and only one can be a cycle. So cycles should have no duplicates.

                unique_cycles = list(set(cycles))
                alpha_0 = 1
                # alpha_1 = |Omega| (number of 3-cycles)
                alpha_1 = len(unique_cycles)

                # alpha_2 = # pairs of vertex-disjoint 3-cycles
                alpha_2 = 0
                for i in range(len(unique_cycles)):
                    for j in range(i + 1, len(unique_cycles)):
                        if not unique_cycles[i] & unique_cycles[j]:
                            alpha_2 += 1

                # alpha_3 = # triples of vertex-disjoint 3-cycles
                alpha_3 = 0
                if p >= 9:
                    for i in range(len(unique_cycles)):
                        for j in range(i + 1, len(unique_cycles)):
                            if unique_cycles[i] & unique_cycles[j]:
                                continue
                            for k in range(j + 1, len(unique_cycles)):
                                if (not unique_cycles[i] & unique_cycles[k] and
                                        not unique_cycles[j] & unique_cycles[k]):
                                    alpha_3 += 1

                H_ocf = alpha_0 + 2 * alpha_1 + 4 * alpha_2 + 8 * alpha_3
                print(f"    S={list(d['S'])}: H={d['H']}, "
                      f"alpha=(1,{alpha_1},{alpha_2},{alpha_3}), "
                      f"OCF={H_ocf}, match={H_ocf == d['H']}")

    # ================================================================
    # KEY ANALYSIS: Is c_3 constant for all circulant tournaments?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"KEY FINDING: c_3 for circulant tournaments on Z_p")
    print(f"{'=' * 60}")

    for p in [3, 5, 7, 11, 13]:
        S_set_qr = set(sorted(set(pow(a, 2, p) for a in range(1, p)) - {0}))
        # c_3 for a circulant tournament on Z_p with connection set S:
        # c_3 = # triples {i,j,k} that form a directed 3-cycle
        # = # triples where (j-i, k-j, i-k) mod p are all in S
        # By circulant symmetry, c_3 = p * (# pairs (a,b) with a in S, b in S,
        #   -(a+b) mod p in S) / 3
        # Wait: fix i=0. Then j=a, k=a+b where a in S, b in S.
        # Need (a+b) mod p != 0, and i-k = -(a+b) mod p in S.
        # So c_3 = p/3 * |{(a,b) in SxS : a+b != 0 mod p and -(a+b) mod p in S}|
        # ... minus cases where j=k or i=j or i=k.
        # j=k means b=0 (excluded since 0 not in S).
        # i=k means a+b=0 mod p, excluded.
        # i=j means a=0, excluded.
        # So c_3 = p/3 * |{(a,b) in SxS : -(a+b) mod p in S}|

        all_S_list = all_circulant_tournaments(p)
        c3_per_S = {}
        for S in all_S_list:
            S_s = set(S)
            count = 0
            for a in S_s:
                for b in S_s:
                    if (p - a - b) % p in S_s:
                        count += 1
            c3_per_S[S] = p * count // 3

        c3_vals = sorted(set(c3_per_S.values()))
        print(f"  p={p}: c_3 values = {c3_vals}")
        if len(c3_vals) == 1:
            c3_formula = c3_vals[0]
            # For regular tournaments: c_3 = C(p,3) - C(p,1)*C(m,2)/2
            # Actually c_3 = p(p-1)(p-5)/24 for regular tournaments
            # Let me check:
            expected = p * (p - 1) * (p - 5) // 24
            print(f"    Formula p(p-1)(p-5)/24 = {expected}, actual = {c3_formula}")
            # Actually the formula for regular tournaments with score s=(n-1)/2:
            # c_3 = C(n,3) - n * C(s,2) where s = (n-1)/2
            s = (p - 1) // 2
            c3_reg = p * (p - 1) * (p - 2) // 6 - p * s * (s - 1) // 2
            print(f"    C(p,3) - p*C(s,2) = {c3_reg}, actual = {c3_formula}")

    # ================================================================
    # KEY ANALYSIS: Why does c_3 not distinguish H values?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"CONFLICT GRAPH STRUCTURE ANALYSIS (p=7)")
    print(f"{'=' * 60}")

    p = 7
    all_S_list = all_circulant_tournaments(p)
    for S in all_S_list:
        S_set = set(S)
        H = ham_count_dp(p, S)

        # Find 3-cycles
        cycles = []
        for i in range(p):
            for j in range(i + 1, p):
                for k in range(j + 1, p):
                    if all((b - a) % p in S_set for a, b in [(i,j),(j,k),(k,i)]):
                        cycles.append(frozenset({i,j,k}))
                    elif all((b - a) % p in S_set for a, b in [(i,k),(k,j),(j,i)]):
                        cycles.append(frozenset({i,j,k}))

        # Conflict graph adjacency
        n_cyc = len(cycles)
        conflicts = 0
        disjoint = 0
        for a in range(n_cyc):
            for b in range(a + 1, n_cyc):
                if cycles[a] & cycles[b]:
                    conflicts += 1
                else:
                    disjoint += 1

        # Degree distribution of Omega
        degrees = []
        for a in range(n_cyc):
            deg = sum(1 for b in range(n_cyc) if a != b and cycles[a] & cycles[b])
            degrees.append(deg)

        deg_dist = defaultdict(int)
        for d in degrees:
            deg_dist[d] += 1

        if S == all_S_list[0] or H != ham_count_dp(p, all_S_list[0]):
            print(f"\n  S = {list(S)}, H = {H}")
            print(f"    |Omega| = {n_cyc}, conflicts = {conflicts}, disjoint = {disjoint}")
            print(f"    Degree dist: {dict(sorted(deg_dist.items()))}")
            print(f"    alpha_1 = {n_cyc}, alpha_2 = {disjoint}")
            print(f"    OCF check: 1 + 2*{n_cyc} + 4*{disjoint} = {1 + 2*n_cyc + 4*disjoint}")


if __name__ == '__main__':
    main()
