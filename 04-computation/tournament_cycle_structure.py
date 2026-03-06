#!/usr/bin/env python3
"""
Tournament cycle structure and hook positivity mechanism.
Instance: opus-2026-03-06-S10

KEY DISCOVERY: For tournaments (with canonical path labeling), the ONLY
non-zero p-coefficients in U_T come from ODD cycle types.

Specifically, p_{mu} = 0 unless all parts of mu are odd.

This is because:
- A k-cycle sigma must be a directed k-cycle in T or T^op
- For even k: in a tournament, a directed k-cycle and its reverse
  can't both fail to exist (one must work), but the SIGN contribution
  from T^op is (-1)^{k-1}. For even k, (-1)^{k-1} = -1.

Wait, that doesn't make them zero. Let me re-examine.

Actually for TRANSPOSITIONS (2-cycles): sigma = (i,j) means sigma(i)=j, sigma(j)=i.
For this to be a directed cycle in T, we need BOTH i->j AND j->i.
But in a tournament, exactly one direction exists. So transpositions
NEVER form directed cycles in T or T^op. Hence p_{(...,2,...)} = 0.

For 4-cycles: sigma = (i,j,k,l) needs i->j->k->l->i.
The reverse (i,l,k,j) needs i->l->k->j->i.
In a tournament, one of each pair {i->j, j->i} exists.
It's possible for NEITHER the 4-cycle nor its reverse to be directed in T.

So p_{(4,...)} can be zero or nonzero depending on T's structure.

Let's test this computationally and see which cycle types actually appear.
"""

from itertools import permutations
from math import factorial
from fractions import Fraction
from collections import defaultdict, Counter

CHARACTER_TABLES = {
    3: {
        'partitions': [(3,), (2,1), (1,1,1)],
        'classes': [(1,1,1), (2,1), (3,)],
        'table': [
            [1, 1, 1],
            [2, 0, -1],
            [1, -1, 1],
        ]
    },
    4: {
        'partitions': [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)],
        'classes': [(1,1,1,1), (2,1,1), (2,2), (3,1), (4,)],
        'table': [
            [1, 1, 1, 1, 1],
            [3, 1, -1, 0, -1],
            [2, 0, 2, -1, 0],
            [3, -1, -1, 0, 1],
            [1, -1, 1, 1, -1],
        ]
    },
    5: {
        'partitions': [(5,), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)],
        'classes': [(1,1,1,1,1), (2,1,1,1), (2,2,1), (3,1,1), (3,2), (4,1), (5,)],
        'table': [
            [1, 1, 1, 1, 1, 1, 1],
            [4, 2, 0, 1, -1, 0, -1],
            [5, 1, 1, -1, 1, -1, 0],
            [6, 0, -2, 0, 0, 0, 1],
            [5, -1, 1, -1, -1, 1, 0],
            [4, -2, 0, 1, 1, 0, -1],
            [1, -1, 1, 1, -1, -1, 1],
        ]
    },
}

def z_lambda(partition):
    cnt = Counter(partition)
    result = 1
    for k, mk in cnt.items():
        result *= (k ** mk) * factorial(mk)
    return result

def is_hook(partition):
    return sum(1 for p in partition if p > 1) <= 1

def tournament_from_bits(bits, n):
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j: adj[(i,j)] = 0
    for i in range(n-1):
        adj[(i,i+1)] = 1; adj[(i+1,i)] = 0
    idx = 0
    for gap in range(2, n):
        for i in range(n - gap):
            j = i + gap
            if (bits >> idx) & 1: adj[(i,j)] = 1; adj[(j,i)] = 0
            else: adj[(j,i)] = 1; adj[(i,j)] = 0
            idx += 1
    return adj

def get_cycles(perm):
    n = len(perm); visited = [False]*n; cycles = []
    for start in range(n):
        if visited[start]: continue
        cycle = []; v = start
        while not visited[v]:
            visited[v] = True; cycle.append(v); v = perm[v]
        cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle_in(cycle, adj):
    k = len(cycle)
    if k == 1: return True
    return all(adj.get((cycle[i], cycle[(i+1)%k]), 0) == 1 for i in range(k))

def compute_p_expansion(n, adj):
    adj_op = {(j,i): adj.get((i,j),0) for i in range(n) for j in range(n) if i!=j}
    coeffs = defaultdict(int)
    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple); cycles = get_cycles(perm)
        phi = 0; valid = True
        for cycle in cycles:
            if len(cycle) == 1: continue
            in_T = is_directed_cycle_in(cycle, adj)
            in_Top = is_directed_cycle_in(cycle, adj_op)
            if not in_T and not in_Top: valid = False; break
            if in_Top and not in_T: phi += len(cycle) - 1
        if not valid: continue
        partition = tuple(sorted([len(c) for c in cycles], reverse=True))
        coeffs[partition] += (-1)**phi
    return dict(coeffs)

def count_ham_paths(n, adj):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj.get((v, u), 0):
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_directed_kcycles(n, adj, k):
    """Count directed k-cycles in tournament T."""
    count = 0
    from itertools import combinations
    for vertices in combinations(range(n), k):
        # Check all cyclic orderings
        from itertools import permutations as perms
        for perm in perms(vertices):
            # Check if canonical (smallest element first)
            if perm[0] != min(perm): continue
            # Check if directed cycle
            if all(adj.get((perm[i], perm[(i+1)%k]), 0) == 1 for i in range(k)):
                count += 1
    return count

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def main():
    print("=== TOURNAMENT CYCLE STRUCTURE & HOOK POSITIVITY MECHANISM ===\n")

    # Part 1: Which cycle types have non-zero p-coefficients?
    print("PART 1: Non-zero p-coefficient cycle types for tournaments")
    print("="*70)

    for n in range(3, 7):
        print(f"\nn = {n}")
        m = n*(n-1)//2 - (n-1)
        all_nonzero_types = set()
        seen = {}

        for bits in range(2**m):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))

            key = (scores, H)
            if key in seen: continue
            seen[key] = bits

            p_coeffs = compute_p_expansion(n, adj)
            for mu, c in p_coeffs.items():
                if c != 0:
                    all_nonzero_types.add(mu)

        print(f"  Non-zero cycle types: {sorted(all_nonzero_types)}")
        all_odd = [mu for mu in all_nonzero_types if all(p % 2 == 1 for p in mu)]
        has_even = [mu for mu in all_nonzero_types if any(p % 2 == 0 for p in mu)]
        print(f"  All-odd-parts types: {sorted(all_odd)}")
        print(f"  Has-even-part types: {sorted(has_even)}")

    # Part 2: Hook character signs at relevant cycle types
    print("\n\nPART 2: Hook character signs at tournament-relevant cycle types")
    print("="*70)

    for n in [3, 4, 5]:
        ct = CHARACTER_TABLES[n]
        parts = ct['partitions']; classes = ct['classes']; table = ct['table']
        hooks = [p for p in parts if is_hook(p)]
        non_hooks = [p for p in parts if not is_hook(p)]

        # Identify which classes have all odd parts
        odd_classes = [c for c in classes if all(p % 2 == 1 for p in c)]

        print(f"\nn = {n}")
        print(f"  Tournament-relevant classes (all odd parts): {[partition_str(c) for c in odd_classes]}")

        print(f"\n  Hook characters at these classes:")
        for lam_idx, lam in enumerate(parts):
            if not is_hook(lam): continue
            vals = []
            for c in odd_classes:
                c_idx = classes.index(c)
                vals.append(table[lam_idx][c_idx])
            all_pos = all(v >= 0 for v in vals)
            print(f"    chi^{partition_str(lam)} at {[partition_str(c) for c in odd_classes]}: "
                  f"{vals}  {'ALL >= 0' if all_pos else 'HAS NEGATIVE'}")

        if non_hooks:
            print(f"\n  Non-hook characters at these classes:")
            for lam_idx, lam in enumerate(parts):
                if is_hook(lam): continue
                vals = []
                for c in odd_classes:
                    c_idx = classes.index(c)
                    vals.append(table[lam_idx][c_idx])
                all_pos = all(v >= 0 for v in vals)
                print(f"    chi^{partition_str(lam)} at {[partition_str(c) for c in odd_classes]}: "
                      f"{vals}  {'ALL >= 0' if all_pos else 'HAS NEGATIVE'}")

    # Part 3: The hook positivity theorem at n=4
    print("\n\nPART 3: Hook positivity theorem at n=4")
    print("="*70)
    print("""
At n=4, tournament-relevant p-types are: (1,1,1,1) and (3,1).
(p_{(2,1,1)} = p_{(2,2)} = p_{(4)} = 0 for ALL tournaments.)

Hook characters at these types:
  chi^{(4)}((1^4))=1, chi^{(4)}((3,1))=1
  chi^{(3,1)}((1^4))=3, chi^{(3,1)}((3,1))=0
  chi^{(2,1,1)}((1^4))=3, chi^{(2,1,1)}((3,1))=0
  chi^{(1^4)}((1^4))=1, chi^{(1^4)}((3,1))=1

ALL entries are >= 0!

Since p_{(1^4)} = 1 >= 0 and p_{(3,1)} = 2*c3 >= 0 (OCF),
every hook Schur coefficient is a sum of non-negative terms.

THEOREM (n=4): [s_{hook}] U_T >= 0 for all tournaments on 4 vertices.

Proof: p_{mu} = 0 for mu with even parts (no even directed cycles).
       chi^{hook}(mu) >= 0 for mu with all odd parts.
       OCF gives p_{mu} >= 0.
       [s_{hook}] = sum_mu chi^{hook}(mu) * p_{mu} / z_{mu} >= 0.  QED

Non-hook coefficient:
  chi^{(2,2)}((1^4))=2, chi^{(2,2)}((3,1))=-1
  [s_{(2,2)}] = 2/24 + (-1)/3 * 2*c3 = 1/12 - 2c3/3
  Negative whenever c3 >= 1 (all non-transitive tournaments).
""")

    # Part 4: Does this extend to n=5?
    print("PART 4: Testing the mechanism at n=5")
    print("="*70)

    n = 5
    ct = CHARACTER_TABLES[n]
    parts = ct['partitions']; classes = ct['classes']; table = ct['table']
    hooks = [p for p in parts if is_hook(p)]
    odd_classes = [c for c in classes if all(p % 2 == 1 for p in c)]

    print(f"Tournament-relevant classes: {[partition_str(c) for c in odd_classes]}")
    print(f"  = (1^5), (3,1,1), (5)")

    print("\nHook characters at tournament-relevant types:")
    for lam_idx, lam in enumerate(parts):
        if not is_hook(lam): continue
        vals = []
        for c in odd_classes:
            c_idx = classes.index(c)
            vals.append(table[lam_idx][c_idx])
        print(f"  chi^{partition_str(lam)}: {vals}")

    print("""
chi^{(4,1)} at (5) = -1 < 0!
chi^{(2,1,1,1)} at (5) = -1 < 0!

So the n=4 mechanism (all hook chars non-negative at tournament types)
FAILS at n=5. But hook positivity still holds empirically!

This means there must be a QUANTITATIVE constraint:
  p_{(3,1,1)}/3 - p_{(5)}/5 >= -1/30  (for s_{(4,1)} >= 0)
  i.e., 2*c3/3 >= 2*(#5-cycles)/5 - 1/30

Is this guaranteed by tournament structure?
""")

    # Check the constraint for all n=5 tournaments
    print("Verifying quantitative constraint for all n=5 tournaments:")
    m = n*(n-1)//2 - (n-1)
    seen = {}

    for bits in range(2**m):
        adj = tournament_from_bits(bits, n)
        H = count_ham_paths(n, adj)
        scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))

        key = (scores, H)
        if key in seen: continue
        seen[key] = bits

        p_coeffs = compute_p_expansion(n, adj)
        p311 = p_coeffs.get((3,1,1), 0)
        p5 = p_coeffs.get((5,), 0)

        # s_{(4,1)} = 4/120 + 1/3 * p311 + (-1)/5 * p5
        s41 = Fraction(4, 120) + Fraction(p311, 3) + Fraction(-1, 5) * Fraction(p5)

        c3 = p311 // 2  # p_{(3,1,1)} = 2 * c3
        c5 = p5 // 2    # p_{(5)} = 2 * #5-cycles

        margin = float(s41)
        print(f"  c3={c3}, #5cyc={c5}, p(3,1,1)={p311}, p(5)={p5}: "
              f"[s(4,1)]={float(s41):.6f} (margin={margin:.6f})")

    # Part 5: Relationship to grid geometry
    print("\n\nPART 5: Geometric interpretation")
    print("="*70)
    print("""
The perpendicular grid symmetry (M[a,b] = M[b,a]) corresponds to U_T = U_{T^op}.
This means U_T is fixed by the omega involution at the level of Schur functions:
  [s_lambda] U_T = [s_{lambda'}] U_T

For hooks: (k, 1^{n-k})' = (n-k+1, 1^{k-1}), so hook coefficients are palindromic.

The triangular pin grid has a PERPENDICULAR AXIS OF SYMMETRY.
Tilings symmetric under this axis correspond to SC (self-converse) tournaments.

For non-SC tournaments, the grid has NO perpendicular symmetry,
but the GENERATING FUNCTION still has it (because omega(U_T) = U_T).

This means: even though individual tilings break perpendicular symmetry,
the WEIGHTED COUNT respects it.

Hook Schur functions correspond to EXTERIOR POWERS of the standard rep.
In the grid picture, these correspond to counting tilings where
the "level sets" of the tiling form a specific pattern related to
descents of permutations.

The fact that hook coefficients are non-negative (at n>=4) means:
for each exterior power, the weighted tiling count is non-negative.
This is a form of "signed tiling positivity" restricted to hook shapes.

Question: does the perpendicular symmetry of the grid EXPLAIN why
non-hook coefficients carry all the negativity?
""")

if __name__ == '__main__':
    main()
