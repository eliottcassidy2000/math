import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
trace_ocf_bridge.py
kind-pasteur-2026-03-07-S39b

Investigate connections between trace formulas and OCF.

Key observation: tr(A^k)[v][v] = t_k(v) = # directed k-cycles through v
(for k = 3, 4, 5). And OCF involves the independence polynomial of
Omega(T), the conflict graph of odd directed cycles.

Question: Can we express Omega(T) or I(Omega(T), 2) using trace/matrix data?

The conflict graph adjacency: two odd cycles conflict iff they share a vertex.
This is equivalent to asking whether their vertex sets intersect.

For k=3: the number of pairs of 3-cycles sharing a vertex is related to
sum_v C(t_3(v), 2). And the independence number of Omega_3 is alpha_2.

Goal: Express OCF terms (alpha_1, alpha_2, ...) via trace computations.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from tournament_fast import c3_from_score, c5_fast, c4_fast
from itertools import combinations
from math import comb
import random


def find_all_directed_3_cycles(T):
    """Find all directed 3-cycles as vertex triples. O(n^3)."""
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles


def compute_alpha_2_3cycles(T):
    """Compute alpha_2 = max independent set size in 3-cycle conflict graph.
    More practically: count the number of DISJOINT pairs of 3-cycles."""
    cycles = find_all_directed_3_cycles(T)
    m = len(cycles)
    if m <= 1:
        return 0

    count = 0
    for i in range(m):
        for j in range(i+1, m):
            if not (cycles[i] & cycles[j]):
                count += 1
    return count


def trace_per_vertex(T, k):
    """Compute (A^k)[v][v] for each v. O(n^3)."""
    n = len(T)
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * T[l][j]
        Ak = new
    return [Ak[v][v] for v in range(n)]


def hamiltonian_paths_ocf(T):
    """Compute H(T) via OCF with full odd cycle enumeration."""
    n = len(T)
    if n <= 1:
        return 1

    # Find all odd cycles
    cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            v = list(verts)
            dp = [[0]*k for _ in range(1 << k)]
            dp[1][0] = 1
            for mask in range(1, 1 << k):
                for last in range(k):
                    if dp[mask][last] == 0 or not (mask & (1 << last)):
                        continue
                    for nxt in range(1, k):
                        if mask & (1 << nxt):
                            continue
                        if T[v[last]][v[nxt]]:
                            dp[mask | (1 << nxt)][nxt] += dp[mask][last]
            full = (1 << k) - 1
            for last in range(1, k):
                if T[v[last]][v[0]]:
                    cnt = dp[full][last]
                    for _ in range(cnt):
                        cycles.append(frozenset(verts))

    if not cycles:
        return 1

    m = len(cycles)
    # Build conflict graph and compute I(Omega, 2)
    adj = [0]*m
    for a in range(m):
        for b in range(a+1, m):
            if cycles[a] & cycles[b]:
                adj[a] |= 1 << b
                adj[b] |= 1 << a

    # Count independent sets of each size
    total = 1 + 2*m  # empty + singletons
    # Size 2
    pairs = []
    for a in range(m):
        for b in range(a+1, m):
            if not (adj[a] & (1 << b)):
                pairs.append((a,b))
    total += 4 * len(pairs)
    # Size 3
    for a, b in pairs:
        for c in range(b+1, m):
            if not (adj[a] & (1 << c)) and not (adj[b] & (1 << c)):
                total += 8
    return total


print("=" * 70)
print("TRACE-OCF BRIDGE: per-vertex trace data and OCF")
print("=" * 70)

# Key identity from THM-118:
# (A^3)[v][v] = t_3(v) = # directed 3-cycles through v
# (A^5)[v][v] = t_5(v) = # directed 5-cycles through v
# Sum: sum_v t_k(v) = k * c_k (each k-cycle counted k times)

# OCF: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
# where alpha_j = # independent sets of size j in Omega(T)
# = # collections of j pairwise vertex-disjoint odd directed cycles

# For small n, alpha_1 = c_3 + c_5 + ... (total # odd directed cycles)
# alpha_2 = # disjoint pairs of odd cycles

# Can we express alpha_2 via per-vertex trace data?

# Idea: Two 3-cycles {a,b,c} and {d,e,f} are disjoint iff they share
# no vertex. This means:
# alpha_2(3-cycles only) = C(c_3, 2) - E(Omega_3)
# where E(Omega_3) = # pairs sharing a vertex = (1/2) * sum_v C(t_3(v), 2)

# Wait: sum_v C(t_3(v), 2) counts ordered pairs sharing vertex v.
# But a pair might share 2 vertices, counted at both. In a tournament,
# two 3-cycles can share 0, 1, or 2 vertices (never 3 since unique
# directed cycle per triple).
# Sharing 2 vertices: cycles {a,b,c} and {a,b,d}. In a tournament,
# both directed cycles go through a,b,c and a,b,d with specific directions.
# This is possible: if T[a][b]=T[b][c]=T[c][a]=1 and T[a][d]=T[d][b]=T[b][a]=1...
# Wait, T[b][a]=1 and T[a][b]=1 is impossible in a tournament.
# So if T[a][b]=1 (in the first cycle a->b->c->a), then in the second
# cycle through {a,b,d}, we need a->b or b->a. Since T[a][b]=1, we have
# a->b in the second cycle too. So the second cycle has form a->b->d->a
# (T[b][d]=1 and T[d][a]=1). This is possible!
# These two cycles share vertices a and b, but not the third vertex.
# They conflict (share a and b).

# So pairs can share 1 or 2 vertices. The conflict graph edge count is:
# E = sum_{v} C(t_3(v), 2) - (double-counted pairs sharing 2 vertices)
# Hmm, this overcounts pairs sharing 2 vertices because they're counted
# at both shared vertices.

# Actually: E = # conflicting pairs = # pairs sharing >= 1 vertex.
# If I count sum_v C(t_3(v), 2), each pair sharing exactly 1 vertex
# is counted once, and each pair sharing 2 vertices is counted twice.
# Let s1 = # pairs sharing exactly 1 vertex, s2 = # sharing exactly 2.
# sum_v C(t_3(v), 2) = s1 + 2*s2
# E = s1 + s2
# So E = sum_v C(t_3(v), 2) - s2

# And alpha_2(3-only) = C(c_3, 2) - E = C(c_3, 2) - sum_v C(t_3(v), 2) + s2

# Can we compute s2 from trace data?
# s2 = # pairs of 3-cycles sharing exactly 2 vertices.
# Two 3-cycles sharing 2 vertices: {a,b,c} and {a,b,d}.
# The edge a-b is in both cycles. Since tournaments have unique
# direction on each edge, the number of 3-cycles containing edge (a,b)
# (where T[a][b]=1) is the number of c such that T[b][c]=1 and T[c][a]=1.
# This equals (A^2)[b][a] - T[b][a] ... hmm wait.
# Actually (A^3)[a][a] = sum paths a->?->?->a. This includes all 3-cycles
# through a. For edge a->b, the 3-cycles containing directed edge a->b
# are: a->b->c->a where T[b][c]=1 and T[c][a]=1.
# Let e(a,b) = # c with T[b][c]*T[c][a] = 1 = (A^2)[b][a] restricted to
# vertices != a, != b. Actually (A^2)[b][a] = sum_c T[b][c]*T[c][a], but
# c can equal a (T[b][a]*T[a][a]=0 since no self-loops) or c=b
# (T[b][b]*T[b][a]=0). So (A^2)[b][a] = sum_{c!=a,c!=b} T[b][c]*T[c][a].
# This IS exactly e(a,b)!

# So: for edge a->b (T[a][b]=1), the number of 3-cycles containing
# the directed edge a->b is (A^2)[b][a].

# The number of 2-vertex-sharing pairs for edge a->b is C(e(a,b), 2).
# Total s2 = sum_{(a,b): T[a][b]=1} C(e(a,b), 2).
# This is computable in O(n^3) since (A^2) takes O(n^3) and the sum
# is over O(n^2) edges.

# Let's verify this chain of formulas!

for n in [5, 6, 7]:
    m_bits = n * (n - 1) // 2
    total_tours = 1 << m_bits

    if n >= 7:
        random.seed(42)
        sample = [random.randint(0, total_tours - 1) for _ in range(500)]
        sample_type = "sampled"
    else:
        sample = range(total_tours)
        sample_type = "exhaustive"

    print(f"\nn={n} ({sample_type}, {len(list(sample)) if n < 7 else 500} tournaments):")

    formula_match = 0
    tested = 0

    for bits in (sample if n < 7 else [random.randint(0, total_tours-1) for _ in range(500)]):
        T = tournament_from_bits(n, bits)
        nn = len(T)

        # Direct computation of alpha_2 (3-cycles only: disjoint pairs)
        cycles3 = find_all_directed_3_cycles(T)
        c3 = len(cycles3)
        alpha2_direct = 0
        for i in range(c3):
            for j in range(i+1, c3):
                if not (cycles3[i] & cycles3[j]):
                    alpha2_direct += 1

        # Trace-based computation
        # t_3(v) = (A^3)[v][v]
        t3v = trace_per_vertex(T, 3)

        # A^2 matrix
        A2 = [[sum(T[i][k]*T[k][j] for k in range(nn)) for j in range(nn)] for i in range(nn)]

        # e(a,b) = (A^2)[b][a] for edge a->b
        s2 = 0
        for a in range(nn):
            for b in range(nn):
                if a != b and T[a][b]:
                    eab = A2[b][a]
                    s2 += eab * (eab - 1) // 2

        # Formula: alpha_2 = C(c3, 2) - sum_v C(t3v, 2) + s2
        sum_ct3 = sum(t3v[v] * (t3v[v] - 1) // 2 for v in range(nn))
        alpha2_formula = comb(c3, 2) - sum_ct3 + s2

        if alpha2_direct == alpha2_formula:
            formula_match += 1
        elif tested < 5:
            print(f"  MISMATCH bits={bits}: direct={alpha2_direct}, formula={alpha2_formula}, "
                  f"c3={c3}, sum_ct3={sum_ct3}, s2={s2}")
        tested += 1

    print(f"  alpha_2 formula (3-cycles only): {formula_match}/{tested} "
          f"({'PASS' if formula_match == tested else 'FAIL'})")


# ============================================================
# Now check: full alpha_2 including 5-cycles
# ============================================================
print("\n" + "=" * 70)
print("FULL alpha_2 including 5-cycles")
print("=" * 70)

for n in [5, 6]:
    m_bits = n * (n - 1) // 2
    total_tours = 1 << m_bits

    print(f"\nn={n} (exhaustive, {total_tours} tournaments):")

    # At n=5: max disjoint set of odd cycles has size 1 (5 vertices,
    # each 3-cycle uses 3, two disjoint need 6). So alpha_2=0 always.
    # Except: a 3-cycle (3 verts) and... nothing disjoint.
    # At n=5, alpha_2 = 0 always for Omega_3. For full Omega (including
    # 5-cycles), alpha_2 still 0 since a 5-cycle uses all 5 vertices.

    # At n=6: two disjoint 3-cycles are possible.
    # A 3-cycle and a disjoint 3-cycle: 6 vertices total, fits.
    # No 5-cycle disjoint from a 3-cycle (5+3=8 > 6).

    # At n=7: two disjoint 3-cycles (6 verts), or a 3-cycle and
    # a disjoint 5-cycle (but 3+5=8 > 7, impossible!). Wait, actually
    # the 5-cycle and 3-cycle share no vertex means 5+3=8 vertices needed.
    # At n=7 we only have 7 vertices. So 3+5 disjoint is impossible at n=7!
    # But 3+3 disjoint is possible (6 verts, 1 unused).

    # So at n<=7: alpha_2 comes ONLY from disjoint 3-cycle pairs.
    # At n=8: 3+5 needs 8 vertices, so disjoint 3+5 IS possible at n=8.
    # Also 3+3+3 needs 9 > 8, so alpha_3 from 3-cycles needs n >= 9.
    # But 5+3 = 8, so alpha_2 can include 5-cycle+3-cycle at n=8.

    # For n<=7: alpha_2 = alpha_2(3-cycles only).
    # This means: at n <= 7, OCF depends on alpha_2 which is fully
    # determined by 3-cycle data (computable from trace/A^2).

    all_match = True
    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        nn = len(T)

        # Direct: find all odd cycles
        all_odd_cycles = []
        for k in range(3, nn+1, 2):
            for verts in combinations(range(nn), k):
                v = list(verts)
                dp = [[0]*k for _ in range(1 << k)]
                dp[1][0] = 1
                for mask in range(1, 1 << k):
                    for last in range(k):
                        if dp[mask][last] == 0 or not (mask & (1 << last)):
                            continue
                        for nxt in range(1, k):
                            if mask & (1 << nxt):
                                continue
                            if T[v[last]][v[nxt]]:
                                dp[mask | (1 << nxt)][nxt] += dp[mask][last]
                full = (1 << k) - 1
                for last in range(1, k):
                    if T[v[last]][v[0]]:
                        cnt = dp[full][last]
                        for _ in range(cnt):
                            all_odd_cycles.append(frozenset(verts))

        # alpha_2 from full Omega
        alpha2_full = 0
        m_cyc = len(all_odd_cycles)
        for i in range(m_cyc):
            for j in range(i+1, m_cyc):
                if not (all_odd_cycles[i] & all_odd_cycles[j]):
                    alpha2_full += 1

        # alpha_2 from 3-cycles only
        cycles3 = find_all_directed_3_cycles(T)
        alpha2_3only = 0
        c3 = len(cycles3)
        for i in range(c3):
            for j in range(i+1, c3):
                if not (cycles3[i] & cycles3[j]):
                    alpha2_3only += 1

        if alpha2_full != alpha2_3only:
            all_match = False
            if bits < 10:
                print(f"  bits={bits}: alpha2_full={alpha2_full}, alpha2_3only={alpha2_3only}")

    if all_match:
        print(f"  alpha_2(full Omega) = alpha_2(3-cycles only) for ALL {total_tours} tournaments")
    else:
        print(f"  Some mismatches (5-cycles contribute to alpha_2 at n={n})")


# ============================================================
# Putting it together: H(T) from trace data at n<=7
# ============================================================
print("\n" + "=" * 70)
print("H(T) FROM TRACE DATA: all O(n^3)")
print("=" * 70)

for n in [5, 6, 7]:
    m_bits = n * (n - 1) // 2
    total_tours = 1 << m_bits

    if n >= 7:
        random.seed(42)
        sample = [random.randint(0, total_tours - 1) for _ in range(300)]
        sample_type = "sampled"
    else:
        sample = range(total_tours)
        sample_type = "exhaustive"

    print(f"\nn={n} ({sample_type}):")

    matches = 0
    tested = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)
        nn = len(T)

        # H(T) via DP (ground truth)
        from tournament_lib import hamiltonian_path_count
        h_true = hamiltonian_path_count(T)

        # alpha_1 = total odd directed cycles
        c3 = c3_from_score(T)
        c5 = c5_fast(T)
        alpha_1 = c3 + c5
        # At n<=7, no 7-cycle contributes to alpha_1 since 7 uses all vertices
        # and alpha_1 counts ALL directed cycles, but for n=7 the 7-cycles DO
        # contribute. Wait: at n=7, c_7 is the number of directed 7-cycles.
        # These are directed Hamiltonian cycles. They contribute to alpha_1.
        if n >= 7:
            from tournament_fast import count_3_cycles
            # Need c7 - can compute but it's expensive
            # Actually for alpha_1 at n=7 we need c3 + c5 + c7
            # c7 is the number of directed Hamiltonian cycles
            # This is H(T_on_7_vertices_as_cycle) which is hard to compute fast.
            # Skip the trace formula for alpha_1+c7 for now.
            pass

        # alpha_2 = disjoint pairs of odd cycles (only 3-cycles at n<=7)
        t3v = trace_per_vertex(T, 3)
        A2 = [[sum(T[i][k]*T[k][j] for k in range(nn)) for j in range(nn)] for i in range(nn)]
        s2 = 0
        for a in range(nn):
            for b in range(nn):
                if a != b and T[a][b]:
                    eab = A2[b][a]
                    s2 += eab * (eab - 1) // 2
        sum_ct3 = sum(t3v[v] * (t3v[v] - 1) // 2 for v in range(nn))
        alpha_2 = comb(c3, 2) - sum_ct3 + s2

        if n <= 6:
            # alpha_3 = 0 at n<=6 (need 9 vertices for 3 disjoint 3-cycles)
            # Actually: n=6 can have 2 disjoint 3-cycles. alpha_3 needs 3 disjoint,
            # using 9 vertices. At n=6,7,8: alpha_3=0 for 3-cycles only.
            # But at n=9+, alpha_3 can be positive.
            h_formula = 1 + 2 * alpha_1 + 4 * alpha_2
        elif n == 7:
            # Need c7 for alpha_1
            c7_dp = 0
            for verts in combinations(range(nn), 7):
                v = list(verts)
                dp = [[0]*7 for _ in range(1 << 7)]
                dp[1][0] = 1
                for mask in range(1, 1 << 7):
                    for last in range(7):
                        if dp[mask][last] == 0 or not (mask & (1 << last)):
                            continue
                        for nxt in range(1, 7):
                            if mask & (1 << nxt):
                                continue
                            if T[v[last]][v[nxt]]:
                                dp[mask | (1 << nxt)][nxt] += dp[mask][last]
                full = (1 << 7) - 1
                for last in range(1, 7):
                    if T[v[last]][v[0]]:
                        c7_dp += dp[full][last]
            alpha_1_full = c3 + c5 + c7_dp
            h_formula = 1 + 2 * alpha_1_full + 4 * alpha_2

        if h_true == h_formula:
            matches += 1
        elif tested < 5:
            print(f"  bits={bits}: H_true={h_true}, H_formula={h_formula}, "
                  f"c3={c3}, c5={c5}, alpha_2={alpha_2}")
        tested += 1

    print(f"  H(T) = 1 + 2*alpha_1 + 4*alpha_2: {matches}/{tested} "
          f"({'PASS' if matches == tested else 'FAIL'})")


print("\nDone.")
