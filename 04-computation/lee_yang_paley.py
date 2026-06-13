"""
lee_yang_paley.py - Compute Lee-Yang zeros (roots of I(Omega(T_p), x)) for Paley tournaments.

The connection:
  H(T) = I(Omega(T), 2)   [OCF, proved]
  I(G, x) = grand partition function of hard-core lattice gas on G at fugacity x
  Lee-Yang zeros = roots of I(G, x) in complex plane

For claw-free G: all roots are real negative (Chudnovsky-Seymour 2007)
  - Proved for T_p with p<=8 (since Omega(T_p) is claw-free at n<=8)
  - T_7 has n=7<=8: should have all real negative roots
  - T_11 has n=11>8: Omega(T_11) may not be claw-free, roots might be complex

Key question: Does T_p have all-real-roots for ALL p? Or does it fail for large p?
This is the "Lee-Yang conjecture for Paley tournaments."

Method: enumerate all odd directed cycles of T_p, build Omega(T_p), compute I(G,x).
For small p (3,7,11) this is feasible.

Author: kind-pasteur-2026-03-10-S52
"""
import sys
import time
import numpy as np
from itertools import combinations
import sympy
from sympy import symbols, expand, factor, roots, Poly, re, im

sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def qr_set(p):
    return set((a * a) % p for a in range(1, p))


def tournament_adj(p):
    """Build adjacency matrix of T_p."""
    QR = qr_set(p)
    A = np.zeros((p, p), dtype=int)
    for u in range(p):
        for v in range(p):
            if u != v and (v - u) % p in QR:
                A[u, v] = 1
    return A


def find_odd_cycles(A, p):
    """Find all directed odd cycles in tournament T_p.
    Returns list of frozensets of vertices (undirected vertex sets).
    """
    n = p
    # DFS-based cycle finding
    # For efficiency: a directed cycle visits each vertex at most once
    # (simple cycles). Use Johnson's algorithm conceptually.
    # For small p: enumerate all subsets and check for cycles.

    odd_cycles = []
    # For each odd subset size k = 3, 5, 7, ..., p-2 (and p if complete)
    for k in range(3, p + 1, 2):
        for vset in combinations(range(p), k):
            # Check if vset supports at least one directed k-cycle
            # A k-cycle: permutation sigma of vset with A[sigma[i],sigma[(i+1)%k]]=1 for all i
            # For small k, enumerate permutations of vset
            if k <= 7:
                from itertools import permutations
                found = False
                vset_list = list(vset)
                # Fix first vertex to avoid counting rotations
                v0 = vset_list[0]
                rest = vset_list[1:]
                for perm in permutations(rest):
                    cycle = [v0] + list(perm)
                    if all(A[cycle[i], cycle[(i + 1) % k]] == 1 for i in range(k)):
                        found = True
                        break
                if found:
                    odd_cycles.append(frozenset(vset))
            else:
                # For larger k: use DP
                # Is there a Hamiltonian cycle in the induced subgraph?
                vset_list = list(vset)
                idx = {v: i for i, v in enumerate(vset_list)}
                # DP: dp[mask][i] = can we reach vertex i from vertex 0 using exactly the vertices in mask?
                sub_A = A[np.ix_(vset_list, vset_list)]
                m = k
                dp = np.zeros((1 << m, m), dtype=bool)
                dp[1][0] = True  # start at vertex 0 (= vset_list[0])
                for mask in range(1, 1 << m):
                    for last in range(m):
                        if not (mask >> last & 1):
                            continue
                        if not dp[mask][last]:
                            continue
                        for nxt in range(m):
                            if mask >> nxt & 1:
                                continue
                            if sub_A[last, nxt]:
                                dp[mask | (1 << nxt)][nxt] = True
                full_mask = (1 << m) - 1
                # Check if any last vertex can close back to vertex 0
                found = any(dp[full_mask][last] and sub_A[last, 0] for last in range(1, m))
                if found:
                    odd_cycles.append(frozenset(vset))

    return odd_cycles


def build_conflict_graph(cycles):
    """Build conflict graph Omega: cycles share a vertex iff adjacent."""
    n = len(cycles)
    G = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if cycles[i] & cycles[j]:  # share a vertex
                G[i, j] = G[j, i] = 1
    return G


def independence_polynomial_sympy(G, n_cycles):
    """Compute I(G, x) as a sympy polynomial."""
    x = symbols('x')
    n = n_cycles
    poly = sympy.Integer(0)

    # Enumerate all independent sets
    for size in range(n + 1):
        count = 0
        for indep_set in combinations(range(n), size):
            # Check if independent
            is_indep = True
            for i, j in combinations(indep_set, 2):
                if G[i, j]:
                    is_indep = False
                    break
            if is_indep:
                count += 1
        poly += count * x**size

    return poly


def compute_independence_polynomial_coeffs(G, n_cycles):
    """Compute coefficients of I(G, x) via dynamic programming."""
    # Use DP with bitmask for small n_cycles
    if n_cycles > 25:
        print(f"  Warning: {n_cycles} cycles, DP may be slow")

    # Direct enumeration
    coeffs = [0] * (n_cycles + 1)
    for size in range(n_cycles + 1):
        for indep_set in combinations(range(n_cycles), size):
            is_indep = all(not G[i, j] for i, j in combinations(indep_set, 2))
            if is_indep:
                coeffs[size] += 1

    return coeffs


def analyze_roots(poly_str, coeffs, p):
    """Analyze roots of the independence polynomial."""
    x = symbols('x')
    poly = sum(c * x**k for k, c in enumerate(coeffs))

    print(f"  I(Omega(T_{p}), x) = {poly}")
    print(f"  Coefficients: {coeffs}")
    print(f"  I(Omega(T_{p}), 2) = {sum(c * 2**k for k, c in enumerate(coeffs))}")

    # Find roots numerically
    try:
        p_obj = Poly(list(reversed(coeffs)), x)
        rts = p_obj.nroots(n=20)
        print(f"  Roots (numerical):")
        all_real = True
        all_negative = True
        for r in rts:
            r_re = float(re(r))
            r_im = float(im(r))
            is_real = abs(r_im) < 1e-10
            if not is_real:
                all_real = False
            if r_re >= 0:
                all_negative = False
            print(f"    {r_re:.6f} + {r_im:.6f}i  {'(real)' if is_real else '(COMPLEX!)'}")
        print(f"  All real? {all_real}")
        print(f"  All negative real? {all_real and all_negative}")
        if all_real and all_negative:
            print(f"  LEE-YANG CONFIRMED: All roots are real and negative.")
        elif all_real:
            print(f"  Roots real but some non-negative!")
        else:
            print(f"  LEE-YANG FAILS: Complex roots exist!")
        return all_real, all_negative, rts
    except Exception as e:
        print(f"  Root computation error: {e}")
        return None, None, []


def check_claw_free(G, n_cycles):
    """Check if G is claw-free (no induced K_{1,3} subgraph)."""
    for center in range(n_cycles):
        # Find neighbors of center
        nbrs = [j for j in range(n_cycles) if G[center, j] and j != center]
        # Check if any 3 neighbors are mutually non-adjacent
        for leaves in combinations(nbrs, 3):
            # Leaves are nbrs of center. Check if they're all non-adjacent to each other
            if all(not G[leaves[i], leaves[j]] for i, j in combinations(range(3), 2)):
                return False, (center, leaves)
    return True, None


def main():
    print("=" * 70)
    print("LEE-YANG ZEROS FOR PALEY TOURNAMENTS")
    print("=" * 70)
    print()
    print("Connection: H(T) = I(Omega(T), 2) [OCF]")
    print("Lee-Yang zeros = roots of I(Omega(T), x)")
    print("Theorem (Chudnovsky-Seymour): claw-free G => all roots real negative")
    print("Our result: T_p claw-free for p<=8, fails p=9 (first case n=11 is T_11)")
    print()

    for p in [3, 7, 11]:
        print(f"\n{'='*50}")
        print(f"T_{p} (n={p} vertices)")
        print(f"{'='*50}")
        t0 = time.time()

        A = tournament_adj(p)
        print(f"Adjacency matrix computed.")

        # Find odd cycles
        print(f"Finding odd directed cycles...")
        cycles = find_odd_cycles(A, p)
        print(f"  Found {len(cycles)} distinct odd cycle vertex-sets")

        # Sort cycles by size
        by_size = {}
        for c in cycles:
            k = len(c)
            by_size[k] = by_size.get(k, 0) + 1
        for k in sorted(by_size.keys()):
            print(f"  {k}-cycles: {by_size[k]}")

        # Build conflict graph
        G = build_conflict_graph(cycles)
        n_edges = G.sum() // 2
        print(f"  Conflict graph: {len(cycles)} vertices, {n_edges} edges")

        # Check claw-free
        is_cf, claw_witness = check_claw_free(G, len(cycles))
        print(f"  Omega(T_{p}) claw-free? {is_cf}")
        if not is_cf:
            center, leaves = claw_witness
            print(f"  Claw witness: center={center}, leaves={leaves}")

        # Compute independence polynomial
        print(f"  Computing I(Omega(T_{p}), x)...")
        if len(cycles) <= 30:
            coeffs = compute_independence_polynomial_coeffs(G, len(cycles))
            print()
            analyze_roots(f"T_{p}", coeffs, p)
        else:
            print(f"  Too many cycles ({len(cycles)}) for exact polynomial. Skipping.")

        print(f"  Time: {time.time()-t0:.1f}s")

    print("\n" + "=" * 70)
    print("SUMMARY & INTERPRETATION")
    print("=" * 70)
    print("""
The independence polynomial I(Omega(T_p), x) has:
  - All real negative roots for p=3,7 (claw-free Omega => Chudnovsky-Seymour)
  - For p=11: Omega(T_11) may or may not be claw-free (n=11 > 8)

The Lee-Yang theorem for hard-core models says:
  If all roots are real and negative, then the "pressure" log Z(lambda) is
  analytic in lambda for Re(lambda) > 0 (high-fugacity phase).

For the Hamiltonian path problem:
  H(T) = Z(Omega(T), 2) = partition function at fugacity 2
  Lee-Yang zeros at x = r_i (all real negative)
  => Z(Omega(T), lambda) = prod_i (1 + lambda/|r_i|) has NO complex zeros in Re(lambda) > 0
  => The "phase" of the lattice gas is well-defined at lambda=2

If T_p is a LEE-YANG SAFE tournament (all-real-roots for all p):
  This would explain WHY H(T_p) = I(Omega(T_p), 2) behaves so regularly.
  The H(T_p) values would lie on a "normal" branch of the partition function.

PALEY LEE-YANG CONJECTURE:
  For all Paley primes p == 3 mod 4, all roots of I(Omega(T_p), x) are real and negative.
""")


if __name__ == '__main__':
    main()
    print("DONE.")
