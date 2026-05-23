"""
opus-2026-05-23-S5: Full investigation of HYP-1732 and real-rootedness.

Goals:
1. Build all-0 staircase T_k; find ALL odd cycles (not just 3-cycles).
2. Compute I(Omega(T_k), x) and I_3(T_k, x), compare.
3. Check real-rootedness of both.
4. For alpha(Omega)=2 cases: HYP-1732 structure analysis.
   - Decompose alpha_2 = (C*-B pairs) + alpha_{A-B} + alpha_{A-A}
   - Check alpha_{A-A} = 0 (bipartite structure conjecture)
   - Verify alpha_{A-A} + alpha_{A-B} <= p*|A|

Arc structure of all-0 staircase T_k on n=2k vertices
(pairs: pair_i = {2i, 2i+1} for i=0,...,k-1):
  Type 1: 2i+1 -> 2i           (within pair: odd beats even)
  Type 2: 2a -> 2b   (a<b)     (even-a beats even-b)
  Type 3: 2b -> 2a+1 (a<b)     (even-b beats odd-a)
  Type 4: 2a+1 -> 2b+1 (a<b)   (odd-a beats odd-b)
  Type 5: 2a -> 2b+1 (a<b)     (even-a beats odd-b)
"""
import sys
import time
import numpy as np
from math import comb, factorial
from itertools import combinations

# =====================================================
# PART 1: Tournament construction
# =====================================================

def build_Tk(k):
    """Build adjacency matrix for all-0 staircase T_k."""
    n = 2 * k
    A = [[0]*n for _ in range(n)]
    for i in range(k):
        A[2*i+1][2*i] = 1  # type 1
    for a in range(k):
        for b in range(a+1, k):
            A[2*a][2*b] = 1       # type 2
            A[2*b][2*a+1] = 1     # type 3
            A[2*a+1][2*b+1] = 1   # type 4
            A[2*a][2*b+1] = 1     # type 5
    return A

def verify_tournament(A, n):
    """Check A is a tournament: exactly one of A[i][j], A[j][i] is 1."""
    for i in range(n):
        for j in range(i+1, n):
            assert A[i][j] + A[j][i] == 1, f"Not a tournament at ({i},{j})"
        assert A[i][i] == 0

# =====================================================
# PART 2: Cycle enumeration
# =====================================================

def find_all_directed_cycles(A, n, max_length=None):
    """Find all simple directed cycles in tournament A.
    Returns list of frozensets (vertex sets).
    Uses canonical representative (rotate to min vertex) to avoid duplicates.
    """
    if max_length is None:
        max_length = n

    found = set()  # canonical tuples
    result = []    # frozensets

    def dfs(start, current, path, in_path):
        for nxt in range(n):
            if not A[current][nxt]:
                continue
            if nxt == start:
                # Closed cycle
                L = len(path)
                if L >= 2:
                    # Canonical: rotate so path[0] = min vertex
                    mi = path.index(min(path))
                    canon = tuple(path[mi:] + path[:mi])
                    if canon not in found:
                        found.add(canon)
                        result.append(frozenset(path))
            elif nxt > start and nxt not in in_path and len(path) < max_length:
                # Only extend from vertices > start to avoid double-counting
                path.append(nxt)
                in_path.add(nxt)
                dfs(start, nxt, path, in_path)
                path.pop()
                in_path.remove(nxt)

    for start in range(n):
        dfs(start, start, [start], {start})

    return result

def find_odd_cycles(A, n, max_length=None):
    """Find all simple directed ODD cycles."""
    all_cycles = find_all_directed_cycles(A, n, max_length)
    return [c for c in all_cycles if len(c) % 2 == 1]

def find_3cycles(A, n):
    """Fast enumeration of 3-cycles only."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

# =====================================================
# PART 3: Independence polynomial
# =====================================================

def build_conflict_adj(cycles_list):
    """Build bitmask adjacency: conflict iff share a vertex."""
    m = len(cycles_list)
    adj = [0] * m
    for i in range(m):
        for j in range(i+1, m):
            if cycles_list[i] & cycles_list[j]:
                adj[i] |= (1 << j)
                adj[j] |= (1 << i)
    return adj

def compute_ip(adj, m):
    """Compute independence polynomial via bitmask backtracking."""
    ip = [0] * (m + 1)
    ip[0] = 1

    def backtrack(v, size, forbidden):
        for i in range(v, m):
            if (forbidden >> i) & 1:
                continue
            ip[size + 1] += 1
            backtrack(i + 1, size + 1, forbidden | adj[i])

    backtrack(0, 0, 0)
    while len(ip) > 1 and ip[-1] == 0:
        ip.pop()
    return ip

def compute_ip_from_cycles(cycles_list):
    """Convenience: build adj then compute IP."""
    m = len(cycles_list)
    adj = build_conflict_adj(cycles_list)
    return compute_ip(adj, m)

# =====================================================
# PART 4: Real-rootedness
# =====================================================

def check_real_roots(poly):
    """Check if polynomial (list of coefficients) has all real roots.
    poly = [a0, a1, ..., ad] means a0 + a1*x + ... + ad*x^d.
    Returns (is_real, roots_array).
    """
    d = len(poly) - 1
    if d <= 0:
        return True, []
    if d == 1:
        return True, [-poly[0] / poly[1]]

    # numpy roots: pass coefficients in descending order
    coeffs = list(reversed(poly))
    roots = np.roots(coeffs)

    imaginary_parts = [abs(r.imag) for r in roots]
    is_real = all(im < 1e-6 for im in imaginary_parts)
    real_parts = [r.real for r in roots]

    return is_real, sorted(real_parts)

def log_concave_check(poly):
    """Check log-concavity: alpha_i^2 >= alpha_{i-1} * alpha_{i+1}."""
    violations = []
    for i in range(1, len(poly)-1):
        if poly[i]**2 < poly[i-1] * poly[i+1]:
            violations.append(i)
    return violations

# =====================================================
# PART 5: HYP-1732 structure analysis
# =====================================================

def hyp1732_analysis(cycles_list, adj, ip):
    """
    For each valid pair-partner construction (C*, C**):
    - C* from a max IS NOT in the first max IS S
    - Compute A, B, p, alpha_{A-A}, alpha_{A-B}
    - Verify alpha_2 <= p*(m-p)
    """
    m = len(cycles_list)
    alpha = len(ip) - 1

    if alpha != 2:
        return {'alpha': alpha, 'applicable': False}

    alpha2 = ip[2]

    if alpha2 < 2:
        return {'alpha': 2, 'alpha2': alpha2, 'applicable': False,
                'note': 'alpha2 < 2, handled by Turan-ULC'}

    # Find all max ISs (size 2 = alpha)
    max_is_list = []
    for i in range(m):
        for j in range(i+1, m):
            if not ((adj[i] >> j) & 1):
                max_is_list.append((i, j))

    assert len(max_is_list) == alpha2, f"Expected {alpha2} max ISs, got {len(max_is_list)}"

    results = []
    violations = []

    # Apply pair-partner construction
    for s_idx, (c1, c2) in enumerate(max_is_list):
        S = frozenset([c1, c2])
        for t_idx, (c3, c4) in enumerate(max_is_list):
            if t_idx == s_idx:
                continue
            # C* = element of {c3,c4} not in S
            C_star_candidates = [c for c in [c3, c4] if c not in S]
            if not C_star_candidates:
                continue  # Both c3,c4 in S (only possible if c3=c1 and c4=c2, but t_idx != s_idx)
            C_star = C_star_candidates[0]
            C_partner = c3 if C_star == c4 else c4

            # Compute A, B sets
            B = [i for i in range(m) if i != C_star and not ((adj[C_star] >> i) & 1)]
            A_set = [i for i in range(m) if i != C_star and ((adj[C_star] >> i) & 1)]
            p = len(B)
            A_size = len(A_set)

            # alpha_{A-A}: disjoint pairs within A
            aa = 0
            for i in range(len(A_set)):
                for j in range(i+1, len(A_set)):
                    if not ((adj[A_set[i]] >> A_set[j]) & 1):
                        aa += 1

            # alpha_{A-B}: disjoint pairs between A and B
            ab = 0
            for a in A_set:
                for b in B:
                    if not ((adj[a] >> b) & 1):
                        ab += 1

            # Verify: alpha2 = p (C*-B pairs) + alpha_{A-B} + alpha_{A-A}
            check_sum = p + ab + aa
            assert check_sum == alpha2, f"Sum mismatch: p={p} ab={ab} aa={aa} sum={check_sum} != alpha2={alpha2}"

            # HYP-1732: alpha2 <= p*(m-p)
            bound = p * (m - p)
            holds = (alpha2 <= bound)

            # Bipartite check: A-A pairs should be 0
            bipartite = (aa == 0)

            rec = {
                'S': (c1, c2), 'C_star': C_star, 'C_partner': C_partner,
                'p': p, 'A_size': A_size,
                'alpha_AA': aa, 'alpha_AB': ab,
                'alpha2': alpha2, 'bound': bound, 'holds': holds,
                'bipartite': bipartite,
                'm': m
            }
            results.append(rec)

            if not holds:
                violations.append(rec)

    aa_vals = [r['alpha_AA'] for r in results]
    bipartite_all = all(r['bipartite'] for r in results)

    return {
        'alpha': 2, 'alpha2': alpha2, 'applicable': True,
        'm': m, 'n_tests': len(results),
        'violations': len(violations),
        'bipartite_always': bipartite_all,
        'max_alpha_AA': max(aa_vals) if aa_vals else 0,
        'results': results
    }

# =====================================================
# PART 6: 3-cycle IP formula verification
# =====================================================

def alpha_m_formula(k, m):
    """Compute alpha_m for 3-cycle IP of T_k using THM-321 formula."""
    def A_n(n):
        return factorial(2*n) // factorial(n)

    def c_coeff(m, j):
        if j > m // 2:
            return 0
        prod = 1
        for i in range(j):
            prod *= comb(2*m - j - 3*i, 3)
        return A_n(m - 2*j) * prod // factorial(j)

    return sum(c_coeff(m, j) * comb(k, 2*m - j) for j in range(m // 2 + 1))

def I3_formula(k):
    """Full 3-cycle IP for T_k using formula."""
    d = 2 * k // 3
    return [1] + [alpha_m_formula(k, m) for m in range(1, d+1)]

# =====================================================
# MAIN
# =====================================================

def main():
    print("=" * 70)
    print("FULL ODD-CYCLE INVESTIGATION: ALL-0 STAIRCASE T_k")
    print("=" * 70)

    print("\n=== PART A: Cycle counts and independence polynomials ===\n")
    print(f"{'k':>3} {'n':>3} {'#3cyc':>7} {'#5cyc':>7} {'#7cyc':>7} {'#total':>8} "
          f"{'I_3':>20} {'I_full':>20} {'real_I3':>8} {'real_full':>10}")
    print("-" * 100)

    all_data = {}

    for k in range(2, 7):
        n = 2 * k
        t0 = time.time()

        A = build_Tk(k)
        verify_tournament(A, n)

        # Find 3-cycles
        three_cycles = find_3cycles(A, n)

        # Find all odd cycles
        all_odd = find_odd_cycles(A, n)

        # Categorize by length
        len_count = {}
        for c in all_odd:
            L = len(c)
            len_count[L] = len_count.get(L, 0) + 1

        n3 = len_count.get(3, 0)
        n5 = len_count.get(5, 0)
        n7 = len_count.get(7, 0)

        assert n3 == len(three_cycles), f"k={k}: 3-cycle count mismatch {n3} vs {len(three_cycles)}"

        # Compute IPs
        adj3 = build_conflict_adj(three_cycles)
        ip3 = compute_ip(adj3, len(three_cycles))

        if len(all_odd) <= 28:
            adj_full = build_conflict_adj(all_odd)
            ip_full = compute_ip(adj_full, len(all_odd))
        else:
            ip_full = None

        # Check real roots
        real3, roots3 = check_real_roots(ip3)
        if ip_full is not None:
            real_full, roots_full = check_real_roots(ip_full)
        else:
            real_full = None

        # Log-concavity
        lc3 = log_concave_check(ip3)

        elapsed = time.time() - t0

        ip3_str = str(ip3) if len(str(ip3)) <= 20 else str(ip3)[:17]+"..."
        ipf_str = str(ip_full)[:20] if ip_full else "N/A"

        print(f"k={k}: n={n:2d}, #3={n3:4d}, #5={n5:4d}, #7={n7:4d}, tot={len(all_odd):4d}  "
              f"I3={ip3}  full={ip_full}  "
              f"realI3={'Y' if real3 else 'N'}  realFull={'Y' if real_full else ('N' if real_full is False else '?')}")

        # Formula check
        ip3_formula = I3_formula(k)
        if ip3 != ip3_formula:
            print(f"  WARNING: formula gives {ip3_formula} but computed {ip3}")
        else:
            print(f"  [Formula check: OK, I_3 = {ip3}]")

        if real3 and ip3 != [1]:
            print(f"  I_3 roots: {[f'{r:.4f}' for r in roots3]}")

        if ip_full and ip_full != ip3:
            diff = [ip_full[i] - (ip3[i] if i < len(ip3) else 0) for i in range(len(ip_full))]
            print(f"  FULL - I_3 correction: {diff}")
            if real_full:
                print(f"  FULL roots: {[f'{r:.4f}' for r in roots_full]}")
            else:
                print(f"  FULL NOT real-rooted (complex roots found)")

        print()

        all_data[k] = {
            'A': A, 'n': n, 'three_cycles': three_cycles, 'all_odd': all_odd,
            'len_count': len_count, 'ip3': ip3, 'ip_full': ip_full,
            'real3': real3, 'real_full': real_full
        }

    print("\n=== PART B: HYP-1732 Structure Analysis (alpha=2 cases) ===\n")

    def fast_hyp1732(all_odd, n_vertices):
        """Fast HYP-1732 analysis: count pairs directly, no full IP needed."""
        m = len(all_odd)
        # Build conflict graph
        adj = build_conflict_adj(all_odd)

        # Find alpha (independence number) -- for alpha=2 check
        max_is = [(i, j) for i in range(m) for j in range(i+1, m)
                  if not ((adj[i] >> j) & 1)]

        # Check no triple IS
        alpha2 = len(max_is)

        # Verify alpha <= 2: check no triple of mutually disjoint cycles
        alpha_ge3 = False
        for a, b in max_is:
            for c in range(b+1, m):
                if not ((adj[a]>>c)&1) and not ((adj[b]>>c)&1):
                    alpha_ge3 = True
                    break
            if alpha_ge3:
                break

        alpha = 1 if not max_is else (3 if alpha_ge3 else 2)
        if not max_is:
            # Check if any cycle exists
            alpha = 1 if m > 0 else 0

        if alpha != 2:
            return {'alpha': alpha, 'alpha2': alpha2, 'applicable': False}

        if alpha2 < 2:
            return {'alpha': 2, 'alpha2': alpha2, 'applicable': False,
                    'note': 'alpha2<2, Turan handles', 'm': m}

        results = []
        violations = []

        for s_idx, (c1, c2) in enumerate(max_is):
            S = frozenset([c1, c2])
            for t_idx, (c3, c4) in enumerate(max_is):
                if t_idx == s_idx:
                    continue
                C_star_candidates = [c for c in [c3, c4] if c not in S]
                if not C_star_candidates:
                    continue
                C_star = C_star_candidates[0]

                B = [i for i in range(m) if i != C_star and not ((adj[C_star]>>i)&1)]
                A_s = [i for i in range(m) if i != C_star and ((adj[C_star]>>i)&1)]
                p = len(B)

                aa = sum(1 for i in range(len(A_s)) for j in range(i+1, len(A_s))
                        if not ((adj[A_s[i]]>>A_s[j])&1))
                ab = sum(1 for a in A_s for b in B if not ((adj[a]>>b)&1))

                bound = p * (m - p)
                holds = (alpha2 <= bound)
                bipartite = (aa == 0)

                check = p + ab + aa
                assert check == alpha2, f"Decomp error: {p}+{ab}+{aa}={check} != {alpha2}"

                rec = {'S': (c1,c2), 'C_star': C_star, 'p': p,
                       'A_size': len(A_s), 'alpha_AA': aa, 'alpha_AB': ab,
                       'alpha2': alpha2, 'bound': bound, 'holds': holds,
                       'bipartite': bipartite, 'm': m}
                results.append(rec)
                if not holds:
                    violations.append(rec)

        aa_vals = [r['alpha_AA'] for r in results]
        return {
            'alpha': 2, 'alpha2': alpha2, 'm': m, 'applicable': True,
            'n_tests': len(results), 'violations': len(violations),
            'bipartite_always': all(r['bipartite'] for r in results),
            'max_alpha_AA': max(aa_vals) if aa_vals else 0,
            'results': results
        }

    for k in range(2, 7):
        data = all_data[k]
        all_odd = data['all_odd']
        n_v = data['n']

        if len(all_odd) > 200:
            print(f"k={k}: {len(all_odd)} cycles, too large for full HYP-1732 analysis.")
            continue

        result = fast_hyp1732(all_odd, n_v)

        if not result.get('applicable', False):
            if 'alpha' in result:
                print(f"k={k}: alpha(Omega)={result['alpha']}, alpha2={result.get('alpha2','?')}, "
                      f"not applicable. {result.get('note','')}")
            continue

        print(f"k={k} (n={n_v}):")
        print(f"  m={result['m']} total odd cycles, alpha_2={result['alpha2']}")
        print(f"  Pair-partner tests: {result['n_tests']}, violations: {result['violations']}")
        print(f"  Always bipartite G[A,B] (alpha_AA=0 always): {result['bipartite_always']}")
        print(f"  Max alpha_AA observed: {result['max_alpha_AA']}")

        if result['violations'] == 0:
            print(f"  HYP-1732 HOLDS for all {result['n_tests']} pair-partner choices. ✓")
        else:
            print(f"  HYP-1732 VIOLATED {result['violations']} times! ✗")

        aa_examples = [r for r in result['results'] if r['alpha_AA'] > 0]
        if aa_examples:
            r = aa_examples[0]
            print(f"  Example alpha_AA>0: C*={r['C_star']}, p={r['p']}, "
                  f"|A|={r['A_size']}, aa={r['alpha_AA']}, ab={r['alpha_AB']}, "
                  f"bound={r['bound']}, holds={r['holds']}")
        print()

    print("\n=== PART C: Real-Rootedness Analysis ===\n")

    print("I_3(T_k, x) roots:")
    for k, data in all_data.items():
        ip3 = data['ip3']
        real3, roots3 = check_real_roots(ip3)
        d = len(ip3) - 1
        if real3:
            # Check all roots are negative
            all_neg = all(r < 0 for r in roots3) if roots3 else True
            print(f"  k={k}: I_3={ip3}, d={d}, REAL-ROOTED, roots={[f'{r:.3f}' for r in roots3]}, all_neg={all_neg}")
        else:
            print(f"  k={k}: I_3={ip3}, NOT real-rooted!")

    print()
    print("I(Omega(T_k), x) roots (full):")
    for k, data in all_data.items():
        ip_full = data['ip_full']
        if ip_full is None:
            print(f"  k={k}: Too many cycles for full IP.")
            continue
        real_full, roots_full = check_real_roots(ip_full)
        d = len(ip_full) - 1
        if real_full:
            all_neg = all(r < 0 for r in roots_full) if roots_full else True
            print(f"  k={k}: I_full={ip_full}, d={d}, REAL-ROOTED, roots={[f'{r:.3f}' for r in roots_full]}, all_neg={all_neg}")
        else:
            print(f"  k={k}: I_full={ip_full}, NOT real-rooted!")

    print("\n=== PART D: Bipartite Structure Theorem Investigation ===\n")

    print("Testing: For ALL choices of C* (not just pair-partner), is G[A,B] bipartite?")
    print("(i.e., is alpha_AA = 0 for every C*?)")
    print()

    for k in range(2, 7):
        data = all_data[k]
        all_odd = data['all_odd']
        m = len(all_odd)

        if m > 200:
            print(f"k={k}: {m} cycles, too large.")
            continue

        adj_full = build_conflict_adj(all_odd)

        # Find alpha and alpha2
        max_is_list = [(i,j) for i in range(m) for j in range(i+1,m)
                       if not ((adj_full[i]>>j)&1)]
        alpha2 = len(max_is_list)

        alpha_ge3 = any(
            not ((adj_full[a]>>c)&1) and not ((adj_full[b]>>c)&1)
            for a,b in max_is_list for c in range(b+1,m)
        )
        alpha = 1 if not max_is_list else (3 if alpha_ge3 else 2)
        if not max_is_list:
            alpha = 1 if m > 0 else 0

        if alpha != 2:
            print(f"k={k}: alpha={alpha}, skipping.")
            continue

        all_bipartite_any_Cstar = True
        aa_any_nonzero = []

        for C_star in range(m):
            B = [i for i in range(m) if i != C_star and not ((adj_full[C_star]>>i)&1)]
            A_s = [i for i in range(m) if i != C_star and ((adj_full[C_star]>>i)&1)]

            aa = sum(1 for i in range(len(A_s)) for j in range(i+1, len(A_s))
                    if not ((adj_full[A_s[i]]>>A_s[j])&1))

            if aa > 0:
                all_bipartite_any_Cstar = False
                aa_any_nonzero.append((C_star, len(A_s), len(B), aa))

        print(f"k={k} (m={m} cycles, alpha2={alpha2}): alpha_AA=0 for ALL C*? {all_bipartite_any_Cstar}")
        if aa_any_nonzero:
            print(f"  Counterexamples (C*, |A|, p, alpha_AA): {aa_any_nonzero[:3]}")
            for C_star, A_size, p, aa_val in aa_any_nonzero[:3]:
                A_s = [i for i in range(m) if i != C_star and ((adj_full[C_star]>>i)&1)]
                B = [i for i in range(m) if i != C_star and not ((adj_full[C_star]>>i)&1)]
                ab = sum(1 for a in A_s for b in B if not ((adj_full[a]>>b)&1))
                bound = p * (m - p)
                holds = (alpha2 <= bound)
                print(f"    C*={C_star}: p={p}, |A|={A_size}, aa={aa_val}, ab={ab}, "
                      f"alpha2={alpha2}, bound={bound}, HYP-1732={holds}")
        else:
            print(f"  G[A union B] is bipartite for EVERY C*. ✓")
        print()

    print("\n=== PART E: Key ratios and patterns ===\n")

    for k, data in all_data.items():
        ip3 = data['ip3']
        ip_full = data['ip_full']
        if ip_full is None:
            continue

        real3, roots3 = check_real_roots(ip3)
        real_full, roots_full = check_real_roots(ip_full)

        print(f"k={k}:")
        if roots3:
            print(f"  I_3 roots:  {[f'{r:.4f}' for r in roots3]}")
        if roots_full and real_full:
            print(f"  full roots: {[f'{r:.4f}' for r in roots_full]}")

        # Are I_3 roots a subset of I_full roots?
        if real3 and real_full and roots3 and roots_full:
            # Check interlacing
            sorted3 = sorted(roots3)
            sorted_full = sorted(roots_full)
            print(f"  Roots of I_full vs I_3: {sorted_full} vs {sorted3}")
        print()

if __name__ == '__main__':
    main()
