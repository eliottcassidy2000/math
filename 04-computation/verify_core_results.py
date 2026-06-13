"""
Independent verification of core computational claims in the math research project.

Tests:
  1. H(T) = I(Omega(T), 2) for ALL tournaments on n=5 vertices
  2. H(T_11) = 95095 for the Paley tournament on 11 vertices
  3. Claim B: I(Omega(T),2) - I(Omega(T-v),2) = 2 * sum_{C through v} mu(C) for n=6

All code is self-contained -- no imports from project libraries.

Author: verification script (independent re-derivation from definitions)
"""

from itertools import permutations, combinations
import random
import time


# ============================================================
# Core primitives
# ============================================================

def all_tournaments(n):
    """Generate all tournaments on {0,...,n-1} as adjacency matrices."""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for bits in range(2**m):
        T = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T


def hamiltonian_path_count_dp(T, n):
    """Count Hamiltonian paths in tournament T on n vertices using bitmask DP.

    dp[mask][v] = number of Hamiltonian paths visiting exactly the vertices in mask,
    ending at vertex v.
    """
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    # Base: single vertex
    for v in range(n):
        dp[1 << v][v] = 1
    # Fill
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def hamiltonian_path_count_brute(T, vertices):
    """Brute force count of directed Hamiltonian paths on given vertex set."""
    count = 0
    for perm in permutations(vertices):
        valid = True
        for i in range(len(perm) - 1):
            if T[perm[i]][perm[i + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def find_all_odd_cycles(T, n):
    """Find all directed odd cycles in T. Return list of frozensets of vertex sets,
    along with a dict mapping frozenset -> list of cycle tuples."""
    cycles_by_vset = {}
    for length in range(3, n + 1, 2):
        for subset in combinations(range(n), length):
            # Try all directed cycles on this subset
            # Fix first vertex to avoid counting rotations multiple times
            first = subset[0]
            rest = list(subset[1:])
            for perm in permutations(rest):
                cycle = (first,) + perm
                # Check directed cycle: cycle[0]->cycle[1]->...->cycle[-1]->cycle[0]
                valid = True
                for i in range(length - 1):
                    if T[cycle[i]][cycle[i + 1]] != 1:
                        valid = False
                        break
                if valid and T[cycle[-1]][cycle[0]] == 1:
                    vset = frozenset(subset)
                    if vset not in cycles_by_vset:
                        cycles_by_vset[vset] = []
                    cycles_by_vset[vset].append(cycle)
    # Each directed cycle on k vertices is found (k-1)! / k times via rotation?
    # Actually we fixed first vertex, so each cycle (first, p1, ..., p_{k-1})
    # is found exactly once if it starts at first.
    # But a directed cycle has k rotations, and we only see those starting at first.
    # So we get exactly 1 representative per directed cycle per vertex set.
    # But we need to deduplicate: a cycle v0->v1->...->v_{k-1}->v0 is the same
    # directed cycle as v1->v2->...->v0->v1. Since we fix v0=first, we get
    # exactly one per directed cycle.
    # Actually wait -- we need DISTINCT directed cycles. Two directed cycles are
    # the same iff one is a rotation of the other. By fixing the first vertex,
    # each directed cycle that includes first is counted exactly once.
    # But cycles that don't include first... we chose subsets, and first is always
    # subset[0], so all subsets include first? No -- combinations iterates over
    # all subsets of size `length`.
    #
    # Correction: we iterate over all subsets of size `length` from range(n).
    # For each subset, first = subset[0] (smallest element).
    # Every directed cycle on that subset that starts at the smallest element
    # is counted exactly once. This is correct: fixing the minimum element
    # breaks rotational symmetry completely.

    # Return as list of (vertex_set_frozenset, representative_cycle_tuple)
    result = []
    for vset, cyc_list in cycles_by_vset.items():
        for cyc in cyc_list:
            result.append((vset, cyc))
    return result


def build_conflict_graph(cycles):
    """Build conflict graph: vertices = odd cycles, edge iff they share a vertex.

    cycles: list of (frozenset, cycle_tuple) pairs.
    Returns adjacency list as list of sets.
    """
    k = len(cycles)
    adj = [set() for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if cycles[i][0] & cycles[j][0]:  # shared vertex
                adj[i].add(j)
                adj[j].add(i)
    return adj


def independence_poly_at_x(adj, x):
    """Compute I(G, x) = sum over independent sets S of x^|S|.

    adj: adjacency list. Uses bitmask enumeration for small graphs.
    """
    k = len(adj)
    if k == 0:
        return 1  # Only the empty set

    # For large k, use inclusion via bitmask
    if k > 20:
        raise ValueError(f"Too many cycles ({k}) for bitmask enumeration")

    # Precompute neighbor masks
    nbr = [0] * k
    for i in range(k):
        for j in adj[i]:
            nbr[i] |= (1 << j)

    total = 0
    # Enumerate all subsets, check independence
    for mask in range(1 << k):
        # Check if mask is independent
        independent = True
        size = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            temp &= temp - 1
            size += 1
            if mask & nbr[v]:
                independent = False
                break
        if independent:
            total += x ** size
    return total


def subtournament(T, vertices):
    """Extract subtournament on given vertices. Returns new matrix and vertex list."""
    n = len(vertices)
    vlist = sorted(vertices)
    idx = {v: i for i, v in enumerate(vlist)}
    S = [[0] * n for _ in range(n)]
    for i, u in enumerate(vlist):
        for j, w in enumerate(vlist):
            if i != j:
                S[i][j] = T[u][w]
    return S, vlist


# ============================================================
# TEST 1: H(T) = I(Omega(T), 2) for all n=5 tournaments
# ============================================================

def test_ocf_identity_n5():
    print("=" * 60)
    print("TEST 1: H(T) = I(Omega(T), 2) for ALL n=5 tournaments")
    print("=" * 60)

    n = 5
    total = 0
    failures = 0
    t0 = time.time()

    for T in all_tournaments(n):
        total += 1

        # Compute H(T) via DP
        ht = hamiltonian_path_count_dp(T, n)

        # Find all odd cycles
        cycles = find_all_odd_cycles(T, n)

        if not cycles:
            io2 = 1  # Only empty set, contributes 2^0 = 1
        else:
            adj = build_conflict_graph(cycles)
            io2 = independence_poly_at_x(adj, 2)

        if ht != io2:
            failures += 1
            if failures <= 5:
                print(f"  FAILURE #{failures}: H(T)={ht}, I(Omega,2)={io2}, #cycles={len(cycles)}")

    elapsed = time.time() - t0
    print(f"  Checked {total} tournaments in {elapsed:.1f}s")
    if failures == 0:
        print(f"  PASSED: H(T) = I(Omega(T), 2) holds for all {total} tournaments on n={n}")
    else:
        print(f"  FAILED: {failures}/{total} failures")

    # Also verify small n
    for nn in [3, 4]:
        fails_nn = 0
        cnt = 0
        for T in all_tournaments(nn):
            cnt += 1
            ht = hamiltonian_path_count_dp(T, nn)
            cycles = find_all_odd_cycles(T, nn)
            if not cycles:
                io2 = 1
            else:
                adj = build_conflict_graph(cycles)
                io2 = independence_poly_at_x(adj, 2)
            if ht != io2:
                fails_nn += 1
        status = "PASSED" if fails_nn == 0 else f"FAILED ({fails_nn} failures)"
        print(f"  n={nn}: {cnt} tournaments -- {status}")

    return failures == 0


# ============================================================
# TEST 2: H(T_11) = 95095 for the Paley tournament
# ============================================================

def paley_tournament(p):
    """Build Paley tournament on Z/pZ.
    Arc i->j iff (j-i) is a quadratic residue mod p.
    p must be prime, p = 3 mod 4.
    """
    # Compute quadratic residues mod p
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)

    T = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j:
                if ((j - i) % p) in qr:
                    T[i][j] = 1
    return T


def test_paley_h_t11():
    print()
    print("=" * 60)
    print("TEST 2: H(T_11) = 95095 for the Paley tournament on 11 vertices")
    print("=" * 60)

    p = 11
    T = paley_tournament(p)

    # Verify it's a valid tournament
    for i in range(p):
        for j in range(p):
            if i != j:
                assert T[i][j] + T[j][i] == 1, f"Not a tournament at ({i},{j})"
    print(f"  T_11 is a valid tournament on {p} vertices")

    # Verify QR: residues mod 11 are {1,3,4,5,9}
    qr11 = set()
    for x in range(1, 11):
        qr11.add((x * x) % 11)
    print(f"  QR mod 11: {sorted(qr11)}")

    # Compute H(T_11) via DP
    t0 = time.time()
    ht = hamiltonian_path_count_dp(T, p)
    elapsed = time.time() - t0

    print(f"  H(T_11) = {ht}  (computed in {elapsed:.1f}s)")
    if ht == 95095:
        print(f"  PASSED: H(T_11) = 95095 confirmed")
    else:
        print(f"  FAILED: expected 95095, got {ht}")

    # Also check T_7 for sanity
    T7 = paley_tournament(7)
    ht7 = hamiltonian_path_count_dp(T7, 7)
    print(f"  Sanity check: H(T_7) = {ht7}")

    return ht == 95095


# ============================================================
# TEST 3: Claim B verification for random n=6 tournaments
# ============================================================

def find_odd_cycles_through_v(T, n, v):
    """Find all directed odd cycles through vertex v in T on {0,...,n-1}.
    Returns list of (vertex_frozenset, cycle_tuple) where cycle starts at v.
    """
    others = [u for u in range(n) if u != v]
    results = []
    for length in range(3, n + 1, 2):
        # Choose length-1 vertices from others
        for subset in combinations(others, length - 1):
            for perm in permutations(subset):
                # Cycle: v -> perm[0] -> ... -> perm[-1] -> v
                valid = True
                if T[v][perm[0]] != 1:
                    valid = False
                if valid:
                    for i in range(len(perm) - 1):
                        if T[perm[i]][perm[i + 1]] != 1:
                            valid = False
                            break
                if valid and T[perm[-1]][v] != 1:
                    valid = False
                if valid:
                    cycle = (v,) + perm
                    vset = frozenset(cycle)
                    results.append((vset, cycle))
    return results


def compute_mu(T, n, v, cycle_vset, cycle_tuple):
    """Compute mu(C) for cycle C through v.

    mu(C) = I(Omega(T-v)|_{avoid C\\{v}}, 2)

    That is: in T-v, find all odd cycles, restrict to those vertex-disjoint
    from C\\{v}, build their conflict graph, evaluate I at 2.
    """
    c_minus_v = cycle_vset - {v}

    # Build T-v
    remaining = [u for u in range(n) if u != v]
    n2 = len(remaining)

    # Find all odd cycles in T-v
    Tv, vlist = subtournament(T, remaining)
    # Map original vertices to new indices
    orig_to_new = {vlist[i]: i for i in range(n2)}

    all_cycles_tv = find_all_odd_cycles(Tv, n2)

    # Map back to original vertex labels
    all_cycles_orig = []
    for (vs, cyc) in all_cycles_tv:
        orig_vs = frozenset(vlist[i] for i in vs)
        orig_cyc = tuple(vlist[i] for i in cyc)
        all_cycles_orig.append((orig_vs, orig_cyc))

    # Filter: keep only cycles vertex-disjoint from C\{v}
    filtered = [(vs, cyc) for (vs, cyc) in all_cycles_orig if not (vs & c_minus_v)]

    if not filtered:
        return 1  # Only empty independent set

    adj = build_conflict_graph(filtered)
    return independence_poly_at_x(adj, 2)


def test_claim_b_n6():
    print()
    print("=" * 60)
    print("TEST 3: Claim B for random n=6 tournaments")
    print("  I(Omega(T),2) - I(Omega(T-v),2) = 2 * sum_{C through v} mu(C)")
    print("=" * 60)

    n = 6
    num_tests = 10
    random.seed(42)
    failures = 0
    t0 = time.time()

    for trial in range(num_tests):
        # Generate random tournament
        T = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    T[i][j] = 1
                else:
                    T[j][i] = 1

        # Pick a random vertex
        v = random.randint(0, n - 1)

        # Compute I(Omega(T), 2)
        cycles_T = find_all_odd_cycles(T, n)
        if not cycles_T:
            I_T = 1
        else:
            adj_T = build_conflict_graph(cycles_T)
            I_T = independence_poly_at_x(adj_T, 2)

        # Compute I(Omega(T-v), 2)
        remaining = [u for u in range(n) if u != v]
        Tv_sub, vlist = subtournament(T, remaining)
        n2 = len(remaining)
        cycles_Tv = find_all_odd_cycles(Tv_sub, n2)
        if not cycles_Tv:
            I_Tv = 1
        else:
            adj_Tv = build_conflict_graph(cycles_Tv)
            I_Tv = independence_poly_at_x(adj_Tv, 2)

        lhs = I_T - I_Tv

        # Compute RHS: 2 * sum_{C through v} mu(C)
        cycles_through_v = find_odd_cycles_through_v(T, n, v)

        # Deduplicate cycles: a directed cycle v->a->b->...->v and its rotations
        # are the same cycle. Since we fix start at v, but different directed cycles
        # on the same vertex set are distinct (different orderings).
        # Actually we need to be careful: each directed cycle on vertex set S
        # with |S|=k has k rotations. find_odd_cycles_through_v finds all cycles
        # starting at v, so each directed cycle through v is counted exactly once.
        # But wait -- we also need to count each UNDIRECTED cycle structure?
        # No -- the definitions say "directed odd cycles". So v->a->b->v and
        # v->b->a->v are two different directed cycles (if both exist).

        mu_sum = 0
        for (vset, cyc) in cycles_through_v:
            mu_c = compute_mu(T, n, v, vset, cyc)
            mu_sum += mu_c

        rhs = 2 * mu_sum

        status = "OK" if lhs == rhs else "FAIL"
        if lhs != rhs:
            failures += 1
        print(f"  Trial {trial+1}: v={v}, I(T)={I_T}, I(T-v)={I_Tv}, "
              f"LHS={lhs}, RHS={rhs}, #cycles_through_v={len(cycles_through_v)} -- {status}")

    elapsed = time.time() - t0
    print(f"  Completed {num_tests} trials in {elapsed:.1f}s")
    if failures == 0:
        print(f"  PASSED: Claim B verified for all {num_tests} random n=6 tournaments")
    else:
        print(f"  FAILED: {failures}/{num_tests} failures")

    return failures == 0


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("Core Results Verification Script")
    print("================================\n")

    results = {}

    results["OCF identity (n=3,4,5)"] = test_ocf_identity_n5()
    results["H(T_11) = 95095"] = test_paley_h_t11()
    results["Claim B (n=6, random)"] = test_claim_b_n6()

    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, passed in results.items():
        print(f"  {'PASS' if passed else 'FAIL'}: {name}")

    all_passed = all(results.values())
    print()
    if all_passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
