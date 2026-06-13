"""
Explore the all-0 interleaved staircase tournaments and their H values.

The all-0 staircase at n=2k:
  - Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
  - Global ranking: 0 > 2 > 4 > ... > 2k-2 > 1 > 3 > ... > 2k-1
  - Crossing: follows ranking (lower rank beats higher rank)
  - Within-pair: all-0 means RECESSIVE beats DOMINANT (bit=0 -> recessive wins)

Known values: H(k=2)=5, H(k=3)=29, H(k=4)=233
Question: is this a Markov-related sequence?
"""

from itertools import combinations, product
from collections import defaultdict
import numpy as np

# -----------------------------------------------------------------
# 1. Build all-0 staircase tournament
# -----------------------------------------------------------------

def build_staircase_allzero(k):
    """
    Return adjacency dict: arc_set[(i,j)] = True means i -> j.
    n = 2k vertices: 0, 1, ..., 2k-1
    Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
    Global ranking: 0 > 2 > 4 > ... > 2k-2 > 1 > 3 > ... > 2k-1
    Within-pair: recessive (odd of pair) beats dominant (even of pair).
    """
    n = 2 * k
    # ranking[v] = position in global order (0 = best)
    ranking = {}
    for p in range(k):
        ranking[2*p] = p          # dominants: positions 0..k-1
        ranking[2*p+1] = k + p   # recessives: positions k..2k-1

    arcs = set()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # Check if i and j are in the same pair
            same_pair = (i // 2 == j // 2)
            if same_pair:
                # All-0: recessive beats dominant
                pair_p = i // 2
                dom = 2 * pair_p
                rec = 2 * pair_p + 1
                # rec -> dom  (recessive wins)
                arcs.add((rec, dom))
            else:
                # Crossing: lower ranking wins
                if ranking[i] < ranking[j]:
                    arcs.add((i, j))
    return arcs, n

def get_out_neighbors(arcs, n):
    out = defaultdict(set)
    for (i, j) in arcs:
        out[i].add(j)
    return out

# -----------------------------------------------------------------
# 2. Find all directed odd cycles
# -----------------------------------------------------------------

def find_all_directed_cycles(arcs, n):
    """Find all directed cycles in tournament (not just odd ones)."""
    out = get_out_neighbors(arcs, n)

    all_cycles = []

    # DFS from each vertex, track path
    def dfs(start, current, path, visited):
        for nxt in out[current]:
            if nxt == start and len(path) >= 2:
                # Found a cycle
                all_cycles.append(tuple(path))
            elif nxt > start and nxt not in visited:
                # Only visit vertices with index > start to avoid duplicates
                visited.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, visited)
                path.pop()
                visited.remove(nxt)

    for v in range(n):
        dfs(v, v, [v], {v})

    return all_cycles

def find_odd_directed_cycles(arcs, n):
    """Find all directed odd cycles (length 1 mod 2, i.e., 3,5,7,...)."""
    all_cycles = find_all_directed_cycles(arcs, n)
    odd_cycles = [c for c in all_cycles if len(c) % 2 == 1]
    return odd_cycles

# -----------------------------------------------------------------
# 3. Build conflict graph Omega(T)
# -----------------------------------------------------------------

def build_omega(odd_cycles):
    """
    Conflict graph: vertices = odd cycles, edges = sharing a vertex.
    Returns: list of cycle vertex-sets, adjacency matrix.
    """
    m = len(odd_cycles)
    cycle_sets = [frozenset(c) for c in odd_cycles]

    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:  # share a vertex
                adj[i][j] = adj[j][i] = 1

    return cycle_sets, adj

# -----------------------------------------------------------------
# 4. Compute independence polynomial I(Omega, x)
# -----------------------------------------------------------------

def independence_polynomial(cycle_sets, adj):
    """
    Compute I(Omega, x) = sum_k alpha_k * x^k as a list of coefficients.
    Uses exhaustive enumeration of independent sets.
    Returns list [alpha_0, alpha_1, alpha_2, ...] where alpha_0 = 1.
    """
    m = len(cycle_sets)
    coeffs = defaultdict(int)

    # Enumerate all independent sets via DFS
    def dfs(candidates, current_set, size):
        coeffs[size] += 1
        for idx, v in enumerate(candidates):
            # v can be added if it doesn't conflict with current_set
            new_candidates = [u for u in candidates[idx+1:]
                             if adj[v][u] == 0]
            dfs(new_candidates, current_set + [v], size + 1)

    dfs(list(range(m)), [], 0)

    max_k = max(coeffs.keys()) if coeffs else 0
    return [coeffs[k] for k in range(max_k + 1)]

def eval_poly(coeffs, x):
    """Evaluate polynomial with given coefficients at x."""
    return sum(c * x**k for k, c in enumerate(coeffs))

# -----------------------------------------------------------------
# 5. Find real roots of polynomial
# -----------------------------------------------------------------

def real_roots(coeffs):
    """Find roots of polynomial using numpy."""
    if len(coeffs) <= 1:
        return []
    # numpy poly1d uses descending powers, so reverse
    poly = np.poly1d(list(reversed(coeffs)))
    roots = poly.roots
    return roots

# -----------------------------------------------------------------
# Main computation
# -----------------------------------------------------------------

if __name__ == '__main__':
    print("=" * 70)
    print("All-0 Interleaved Staircase: H values and independence polynomial")
    print("=" * 70)

    known_markov = {1, 2, 5, 13, 29, 34, 89, 169, 194, 233, 433, 610, 985,
                    1325, 1597, 2897, 4181, 5741, 6466, 7561, 9077, 10946,
                    14701, 17711, 28657, 37666, 43261, 51641, 62210, 75025,
                    96557, 135721, 196418, 294685, 426389, 499393, 514229,
                    646018, 925765}

    results = []

    for k in range(2, 7):  # n = 4, 6, 8, 10, 12
        n = 2 * k
        print(f"\n--- k={k}, n={n} ---")

        arcs, _ = build_staircase_allzero(k)

        # Find out-degrees
        out = get_out_neighbors(arcs, n)
        scores = sorted([len(out[v]) for v in range(n)])
        print(f"Score sequence: {scores}")

        odd_cycles = find_odd_directed_cycles(arcs, n)
        print(f"Number of directed odd cycles: {len(odd_cycles)}")

        cycle_lens = sorted(set(len(c) for c in odd_cycles))
        for L in cycle_lens:
            cnt = sum(1 for c in odd_cycles if len(c) == L)
            print(f"  {L}-cycles: {cnt}")

        cycle_sets, adj_omega = build_omega(odd_cycles)
        print(f"Omega(T): {len(cycle_sets)} vertices")

        coeffs = independence_polynomial(cycle_sets, adj_omega)
        print(f"I(Omega, x) coefficients [alpha_0..alpha_k]: {coeffs}")

        H = eval_poly(coeffs, 2)
        print(f"H = I(Omega, 2) = {H}")

        is_markov = H in known_markov
        print(f"Is H a known Markov number? {is_markov}")

        if len(coeffs) > 1:
            roots = real_roots(coeffs)
            all_real = all(abs(r.imag) < 1e-8 for r in roots)
            real_part = sorted([r.real for r in roots])
            print(f"Roots of I(Omega, x): {[f'{r:.4f}' for r in real_part]}")
            print(f"All roots real? {all_real}")
            if all_real:
                print(f"All roots negative? {all(r < 0 for r in real_part)}")

        results.append({
            'k': k, 'n': n, 'H': H,
            'alpha': coeffs,
            'num_odd_cycles': len(odd_cycles),
            'is_markov': is_markov
        })

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'k':>4} {'n':>4} {'H':>12} {'Markov?':>10} {'alpha_k coefficients'}")
    print("-" * 70)
    for r in results:
        print(f"{r['k']:>4} {r['n']:>4} {r['H']:>12} {str(r['is_markov']):>10} {r['alpha']}")

    # Check Markov equation for consecutive H values
    print("\n--- Markov equation check (x, y, z): x²+y²+z² vs 3xyz ---")
    H_vals = [r['H'] for r in results]
    for i in range(len(H_vals) - 2):
        x, y, z = H_vals[i], H_vals[i+1], H_vals[i+2]
        lhs = x**2 + y**2 + z**2
        rhs = 3*x*y*z
        print(f"H[{i}]={x}, H[{i+1}]={y}, H[{i+2}]={z}: LHS={lhs}, RHS={rhs}, diff={lhs-rhs}")

    # Check if consecutive pairs satisfy any Markov triple relation
    print("\n--- Pairwise Markov triple search ---")
    for i in range(len(H_vals) - 1):
        a, b = H_vals[i], H_vals[i+1]
        # If (a, b, c) is a Markov triple: c² - 3ab*c + (a²+b²) = 0
        # c = [3ab ± sqrt(9a²b² - 4(a²+b²))]/2
        discriminant = 9*a**2*b**2 - 4*(a**2 + b**2)
        if discriminant >= 0:
            sqrt_d = int(discriminant**0.5)
            if sqrt_d**2 == discriminant:
                c1 = (3*a*b + sqrt_d) // 2
                c2 = (3*a*b - sqrt_d) // 2
                print(f"H[{i}]={a}, H[{i+1}]={b}: Markov partners c = {c1}, {c2}")
                print(f"  Check ({a},{b},{c1}): {a**2+b**2+c1**2} vs {3*a*b*c1}")
                print(f"  Check ({a},{b},{c2}): {a**2+b**2+c2**2} vs {3*a*b*c2}")
