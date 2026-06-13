#!/usr/bin/env python3
"""
Hook positivity and geometric interpretation.
Instance: opus-2026-03-06-S10

DISCOVERY: At n>=4, all hook Schur coefficients of U_T are non-negative.
Only non-hook coefficients can be negative.

This script investigates WHY by connecting to:
1. The exterior power / transfer matrix eigenvalue interpretation
2. The perpendicular grid symmetry (M[a,b] = M[b,a])
3. The relationship between hooks and the pin grid tiling

Hook characters have the formula:
  chi^{(n-j, 1^j)}(sigma) = (-1)^j * e_j(x_1-1, ..., x_r-1)
where x_i are cycle INDICATORS (not lengths).

Actually, the key identity: for the standard representation V of S_n,
  s_{(n-j, 1^j)} = Lambda^j(V) where V = s_{(n-1,1)}

So [s_{(n-j,1^j)}] U_T = trace of Lambda^j action on U_T's representation.

For U_T = sum_sigma f(sigma) sigma, the hook coefficients are:
  [s_{(n-j,1^j)}] U_T = (1/n!) sum_sigma f(sigma) chi^{(n-j,1^j)}(sigma)

Key question: can we express this as a POSITIVE quantity using the
transfer matrix / tiling interpretation?

Also investigates: the n=3 failure (s(2,1) < 0) and why n>=4 succeeds.
"""

from itertools import permutations
from math import factorial, comb
from fractions import Fraction
from collections import defaultdict, Counter

# Character tables
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

def p_to_schur_full(n, p_coeffs):
    ct = CHARACTER_TABLES.get(n)
    if ct is None: return None
    parts = ct['partitions']; classes = ct['classes']; table = ct['table']
    class_idx = {c: i for i, c in enumerate(classes)}
    s_coeffs = {}
    for i, lam in enumerate(parts):
        coeff = Fraction(0)
        for mu, c_mu in p_coeffs.items():
            j = class_idx.get(mu)
            if j is None: continue
            coeff += Fraction(table[i][j]) * Fraction(c_mu) / Fraction(z_lambda(mu))
        s_coeffs[lam] = coeff
    return s_coeffs

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

def count_directed_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                    if i < j and i < k: count += 1
    return count

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def main():
    print("=== HOOK POSITIVITY: GEOMETRIC & ALGEBRAIC ANALYSIS ===\n")

    # Part 1: Hook coefficient structure
    print("PART 1: Hook coefficient structure as function of H and c3")
    print("="*70)

    for n in [3, 4, 5]:
        ct = CHARACTER_TABLES[n]
        parts = ct['partitions']
        hooks = [p for p in parts if is_hook(p)]
        non_hooks = [p for p in parts if not is_hook(p)]

        print(f"\nn = {n}")
        print(f"Hook partitions: {[partition_str(h) for h in hooks]}")
        print(f"Non-hook partitions: {[partition_str(h) for h in non_hooks]}")

        m = n*(n-1)//2 - (n-1)
        seen = {}

        for bits in range(2**m):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
            c3 = count_directed_3cycles(n, adj)
            key = (scores, H, c3)
            if key in seen: continue
            seen[key] = bits

            p_coeffs = compute_p_expansion(n, adj)
            s_coeffs = p_to_schur_full(n, p_coeffs)

            # Express hook coefficients as function of p-coefficients
            # This reveals what COMBINATIONS of cycle counts determine hooks
            print(f"\n  H={H}, c3={c3}, scores={scores}")
            print(f"  p-expansion: ", end='')
            for mu, c in sorted(p_coeffs.items()):
                print(f"  {c}*p{partition_str(mu)}", end='')
            print()

            for h in hooks:
                c = s_coeffs[h]
                print(f"    [s{partition_str(h)}] = {c} = {float(c):.6f}")

            for h in non_hooks:
                c = s_coeffs[h]
                print(f"    [s{partition_str(h)}] = {c} = {float(c):.6f}  (NON-HOOK)")

    # Part 2: Analyze WHY hooks are positive at n>=4
    print("\n\n" + "="*70)
    print("PART 2: What makes hook coefficients positive?")
    print("="*70)

    print("""
Key algebraic fact:
  [s_{(n-j, 1^j)}] f = (1/n!) sum_sigma f(sigma) * chi^{(n-j,1^j)}(sigma)

For hooks, chi^{(n-j,1^j)} at cycle type mu = (mu_1,...,mu_r) is:
  = sum over subsets S of {mu_1,...,mu_r} with sum = j of (-1)^{j - |S|}

For U_T, f(sigma) counts signed directed decompositions.

The geometric interpretation via transfer matrix M:
  U_T(x_1,...,x_n) = sum over (a,b) M[a,b] * x^a * x^b
where M[a,b] is the transfer matrix entry counting tilings.

Since M[a,b] = M[b,a] (the symmetry we want to prove = OCF),
U_T is symmetric under the perpendicular grid reflection.

For hooks specifically:
  [s_{(n-j, 1^j)}] = [s_{(n-j, 1^j)}] of a symmetric polynomial
The hook characters are related to EXTERIOR POWERS of the standard rep.
""")

    # Part 3: Ratio analysis - hook coefficients vs H
    print("PART 3: Hook coefficients normalized by H (Hamiltonian path count)")
    print("="*70)

    for n in [4, 5]:
        ct = CHARACTER_TABLES[n]
        parts = ct['partitions']
        hooks = [p for p in parts if is_hook(p)]

        print(f"\nn = {n}")
        m = n*(n-1)//2 - (n-1)
        seen = {}

        for bits in range(2**m):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
            c3 = count_directed_3cycles(n, adj)
            key = (scores, H, c3)
            if key in seen: continue
            seen[key] = bits

            p_coeffs = compute_p_expansion(n, adj)
            s_coeffs = p_to_schur_full(n, p_coeffs)

            ratios = []
            for h in hooks:
                c = s_coeffs[h]
                ratio = float(c) / H if H > 0 else 0
                ratios.append(ratio)

            ratio_strs = [f"{r:.6f}" for r in ratios]
            print(f"  H={H:>3}, c3={c3}: ratios = [{', '.join(ratio_strs)}]")

    # Part 4: Connection to p-coefficients
    print("\n\nPART 4: Hook Schur from p-coefficients (explicit formula)")
    print("="*70)
    print("""
At n=4, the hooks are (4), (3,1), (2,1,1), (1,1,1,1).
The character table restricted to hooks gives:

     p(1^4)  p(2,1,1)  p(2,2)  p(3,1)  p(4)
(4)     1       1        1       1       1
(3,1)   3       1       -1       0      -1
(2,1,1) 3      -1       -1       0       1
(1^4)   1      -1        1       1      -1

So: [s_{(4)}] = p1/24 + p21/4 + p22/8 + p31/3 + p4/4
    [s_{(3,1)}] = 3*p1/24 + p21/4 - p22/8 + 0 - p4/4
    = p1/8 + p21/4 - p22/8 - p4/4

Let's verify with actual data.
""")

    # Verify the formula explicitly for n=4
    n = 4
    ct = CHARACTER_TABLES[n]
    parts = ct['partitions']; classes = ct['classes']; table = ct['table']
    hooks = [p for p in parts if is_hook(p)]

    for i, lam in enumerate(parts):
        if not is_hook(lam): continue
        print(f"  [s{partition_str(lam)}] = ", end='')
        terms = []
        for j, mu in enumerate(classes):
            chi_val = table[i][j]
            z = z_lambda(mu)
            if chi_val != 0:
                terms.append(f"({chi_val}/{z})*p{partition_str(mu)}")
        print(" + ".join(terms))

    print("\nFor any tournament T, OCF gives p_{mu} >= 0 for all mu.")
    print("So positivity of hooks reduces to: is the weighted sum positive?")
    print()

    # Check: for each hook, what is the minimum possible value given p >= 0?
    # The issue is that some chi values are negative, e.g., chi^{(3,1)}(2,2) = -1
    # So we need enough positive terms to outweigh negative ones.

    # The key observation: at n=4, for non-transitive tournaments,
    # p_{(2,2)} and p_{(4)} tend to be small relative to p_{(1,1,1,1)} and p_{(2,1,1)}

    print("OCF p-coefficient ranges across all n=4 tournaments:")
    m = n*(n-1)//2 - (n-1)
    all_p = defaultdict(list)
    seen = {}
    for bits in range(2**m):
        adj = tournament_from_bits(bits, n)
        H = count_ham_paths(n, adj)
        scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
        c3 = count_directed_3cycles(n, adj)
        key = (scores, H, c3)
        if key in seen: continue
        seen[key] = bits
        p_coeffs = compute_p_expansion(n, adj)
        for mu in classes:
            all_p[mu].append(p_coeffs.get(mu, 0))

    for mu in classes:
        vals = all_p[mu]
        print(f"  p{partition_str(mu)}: min={min(vals)}, max={max(vals)}, "
              f"ratio to z={min(vals)}/{z_lambda(mu)}..{max(vals)}/{z_lambda(mu)}")

    # Part 5: Why n=3 fails but n>=4 succeeds
    print("\n\nPART 5: Why n=3 fails but n>=4 succeeds")
    print("="*70)
    print("""
At n=3, ALL partitions are hooks: (3), (2,1), (1,1,1).
The standard representation (2,1) has character [2, 0, -1] at classes [(1^3), (2,1), (3)].

For the 3-cycle tournament (the unique non-transitive tournament at n=3):
  p-expansion: p_{(1,1,1)} = 0, p_{(3)} = 2 (two 3-cycles: T and T^op)
  (Actually p_{(1,1,1)} counts identity = always 1 for any tournament)

  [s_{(2,1)}] = chi^{(2,1)}((1^3))/z_{(1^3)} * p_{(1^3)}
              + chi^{(2,1)}((3))/z_{(3)} * p_{(3)}
              = (2/6)*? + (-1/3)*?

The issue: at n=3, the 3-cycle contributes NEGATIVELY to s_{(2,1)} and there's
not enough positive contribution from smaller cycle types to compensate.

At n>=4, the additional structure (more cycle types, larger z_lambda denominators)
allows positive contributions to dominate for hooks.

Key structural reason: hook characters have the property that
chi^{hook}(identity) = dim(hook rep) which is ALWAYS positive and large
(it equals C(n-1, j) for hook (n-j, 1^j)).
The identity contribution is dim/n! * p_{(1^n)} = C(n-1,j)/n! * n! = C(n-1,j).
Wait, p_{(1^n)} = n! for the identity (since all perms have (1^n) as their p-contribution).
Actually no: p_{(1^n)} is the coefficient of the power sum p_{1}^n in U_T.

Let me just compute the identity contribution carefully.
""")

    for n in [3, 4, 5]:
        ct = CHARACTER_TABLES[n]
        parts = ct['partitions']; hooks = [p for p in parts if is_hook(p)]

        print(f"\nn={n}:")
        # p_{(1^n)} = n! for the identity permutation (always contributes +1)
        # Actually p_{(1^n)} coefficient in U_T is the number of permutations
        # that contribute to U_T with all fixed points... that's just 1 (the identity)
        # times n! from the expansion... wait.
        #
        # compute_p_expansion sums over ALL perms. Identity has cycle type (1^n).
        # For identity: all cycles are fixed points, so phi=0, valid=True.
        # Contribution: +1 to p_{(1^n)}.
        # So [p_{(1^n)}] U_T = 1 always (just the identity).

        id_part = tuple([1]*n)
        id_z = z_lambda(id_part)

        print(f"  z(1^{n}) = {id_z}")
        print(f"  [p(1^{n})]U_T = 1 (always, from identity perm)")
        print(f"  Identity contribution to [s_hook]:")

        for i, lam in enumerate(parts):
            if not is_hook(lam): continue
            j_idx = [ii for ii, c in enumerate(ct['classes']) if c == id_part][0]
            chi_id = ct['table'][i][j_idx]
            contrib = Fraction(chi_id, id_z)
            print(f"    [s{partition_str(lam)}] += {chi_id}/{id_z} = {float(contrib):.6f}")

if __name__ == '__main__':
    main()
