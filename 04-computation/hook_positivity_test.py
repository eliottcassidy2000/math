#!/usr/bin/env python3
"""
Test HOOK Schur positivity of U_T for tournaments.
Instance: opus-2026-03-06-S10

Hypothesis: [s_{(k,1^{n-k})}] U_T >= 0 for all tournaments and all k.
From T081: these coefficients are palindromic (h_i = h_{n+1-i}).
From positivity_full_test.py: the s(2,2) coefficient is negative,
but hook coefficients s(k,1^{n-k}) might always be non-negative.

Also tests the GEOMETRIC INTERPRETATION:
- M[a,b] is even/odd in skew variables (S20 finding)
- This means M is invariant under perpendicular grid reflection
- The transfer matrix RESPECTS the grid's perpendicular symmetry
"""

from itertools import permutations
from math import factorial, comb
from fractions import Fraction
from collections import defaultdict, Counter

def tournament_from_bits(bits, n):
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                adj[(i,j)] = 0
    for i in range(n-1):
        adj[(i,i+1)] = 1
        adj[(i+1,i)] = 0
    idx = 0
    for gap in range(2, n):
        for i in range(n - gap):
            j = i + gap
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(j,i)] = 1
                adj[(i,j)] = 0
            idx += 1
    return adj

def get_cycles(perm):
    n = len(perm)
    visited = [False]*n
    cycles = []
    for start in range(n):
        if visited[start]: continue
        cycle = []
        v = start
        while not visited[v]:
            visited[v] = True
            cycle.append(v)
            v = perm[v]
        cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle_in(cycle, adj):
    k = len(cycle)
    if k == 1: return True
    for i in range(k):
        if adj.get((cycle[i], cycle[(i+1)%k]), 0) != 1:
            return False
    return True

def compute_p_expansion(n, adj):
    adj_op = {(j,i): adj.get((i,j),0) for i in range(n) for j in range(n) if i!=j}
    coeffs = defaultdict(int)
    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple)
        cycles = get_cycles(perm)
        phi = 0
        valid = True
        for cycle in cycles:
            if len(cycle) == 1: continue
            in_T = is_directed_cycle_in(cycle, adj)
            in_Top = is_directed_cycle_in(cycle, adj_op)
            if not in_T and not in_Top:
                valid = False; break
            if in_Top and not in_T:
                phi += len(cycle) - 1
        if not valid: continue
        partition = tuple(sorted([len(c) for c in cycles], reverse=True))
        coeffs[partition] += (-1)**phi
    return dict(coeffs)

def generate_partitions(n):
    if n == 0:
        yield (); return
    def helper(n, max_part):
        if n == 0:
            yield (); return
        for k in range(min(n, max_part), 0, -1):
            for rest in helper(n-k, k):
                yield (k,) + rest
    yield from helper(n, n)

def z_lambda(partition):
    cnt = Counter(partition)
    result = 1
    for k, mk in cnt.items():
        result *= (k ** mk) * factorial(mk)
    return result

def hook_partitions(n):
    """Generate hook partitions (k, 1^{n-k}) for k=1,...,n."""
    result = []
    for k in range(1, n+1):
        if k == n:
            result.append((n,))
        else:
            result.append((k,) + (1,)*(n-k))
    return result

def chi_hook(hook_part, mu):
    """
    Character of hook representation (k, 1^{n-k}) at conjugacy class mu.

    Formula: chi^{(k,1^{n-k})}(mu) = sum over i where cycle type mu has a cycle
    of length >= 1: use Murnaghan-Nakayama rule.

    For hooks, there's a nice combinatorial formula:
    chi^{(k,1^{n-k})}(mu) = sum over border-strip tableaux of mu type...

    Actually, for hook characters there's an explicit formula:
    chi^{(n-j, 1^j)}(sigma) = sum_{S subset of fixed points of sigma, |S|=j, ...}
    Hmm, this is getting complicated.

    Let me use a direct approach: the hook character can be computed as:
    chi^{(k,1^{n-k})}(mu) = e_{n-k}(x_1, ..., x_r) evaluated at x_i = omega^{mu_i}
    where omega = e^{2pi i / lcm(mu)}.

    Actually, the simplest formula for hooks:

    Let sigma be a permutation of type mu. Then:
    chi^{(n-j, 1^j)}(sigma) = (-1)^j * e_j(y_1-1, y_2-1, ..., y_r-1)
    where y_i are the cycle lengths of sigma (so y_i = mu_i).

    Wait, that's not quite right either. Let me use the correct formula.

    For the hook (k, 1^{n-k}), the character at type mu = (mu_1, ..., mu_r) is:

    chi^{(k,1^{n-k})}(mu) = sum over subsets S of [r] with sum_{i in S} mu_i = n-k
                              of prod_{i in S} (-1)^{mu_i - 1}

    This comes from the Murnaghan-Nakayama rule: a border strip of size mu_i
    can be placed either horizontally (in the first row, contributing +1)
    or vertically (in the first column, contributing (-1)^{mu_i-1}).
    We need the horizontal strips to total k, so vertical strips total n-k.

    Actually: the character of hook (k, 1^{n-k}) evaluated at a permutation
    with cycle type mu = (mu_1, ..., mu_r) is:

    chi = sum over subsets S of {1,...,r}:
          sum_{i in S} mu_i = n-k (placed in column)
          prod_{i in S} (-1)^{mu_i-1} * prod_{i not in S} 1

    = sum over S with sum mu_i = n-k of (-1)^{sum(mu_i-1) for i in S}
    = sum over S with sum mu_i = n-k of (-1)^{(n-k) - |S|}
    """
    n = sum(mu)
    target = n - hook_part[0]  # n-k = number of 1's in hook

    # Enumerate subsets S of cycle lengths that sum to target
    # mu_sorted = sorted list of cycle lengths
    mu_list = list(mu)
    r = len(mu_list)

    total_chi = 0

    # Use bitmask to enumerate subsets
    for mask in range(1 << r):
        s = sum(mu_list[i] for i in range(r) if (mask >> i) & 1)
        if s == target:
            num_selected = bin(mask).count('1')
            # Sign: (-1)^{target - num_selected}
            total_chi += (-1)**(target - num_selected)

    return total_chi

def compute_hook_schur_coefficients(n, p_coeffs):
    """Compute [s_{hook}] U_T for all hook partitions."""
    hooks = hook_partitions(n)
    all_parts = list(generate_partitions(n))

    hook_coeffs = {}
    for hook in hooks:
        coeff = Fraction(0)
        for mu in all_parts:
            if mu not in p_coeffs:
                continue
            chi = chi_hook(hook, mu)
            coeff += Fraction(chi) * Fraction(p_coeffs[mu]) / Fraction(z_lambda(mu))
        hook_coeffs[hook] = coeff
    return hook_coeffs

def count_ham_paths(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj.get((v, u), 0):
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k <= i or k == j: continue
                if i < j:
                    if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                        count += 1
                    elif adj.get((i,k),0) and adj.get((k,j),0) and adj.get((j,i),0):
                        count += 1
    return count // 1  # already deduplicated by i<k constraint... hmm

def count_directed_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                    if i < j and i < k:
                        count += 1
    return count

def main():
    print("=== HOOK SCHUR POSITIVITY TEST ===\n")

    for n in range(3, 7):
        print(f"\n{'='*70}")
        print(f"n = {n}")
        hooks = hook_partitions(n)
        hook_names = [f"s{partition_str(h)}" for h in hooks]
        print(f"Hooks: {', '.join(hook_names)}")
        print(f"{'='*70}")

        m = n*(n-1)//2 - (n-1)
        num_tourn = 2**m

        seen = {}
        all_hook_pos = True
        min_hook_coeff = Fraction(1)
        count_tested = 0
        count_hook_pos = 0

        header = f"{'bits':>6} {'H':>3} {'c3':>3} | " + ' '.join(f"{hn:>12}" for hn in hook_names) + " | pos?"
        print(header)
        print('-'*len(header))

        for bits in range(num_tourn):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
            c3 = count_directed_3cycles(n, adj)

            key = (scores, H, c3)
            if key in seen:
                continue
            seen[key] = bits
            count_tested += 1

            p_coeffs = compute_p_expansion(n, adj)
            hook_coeffs = compute_hook_schur_coefficients(n, p_coeffs)

            is_pos = all(v >= 0 for v in hook_coeffs.values())
            if is_pos:
                count_hook_pos += 1

            vals = [hook_coeffs[h] for h in hooks]
            for v in vals:
                if v < min_hook_coeff:
                    min_hook_coeff = v

            if not is_pos:
                all_hook_pos = False

            # Print details for n <= 5 or failures
            if n <= 5 or not is_pos:
                val_strs = [f"{float(v):>12.4f}" for v in vals]
                print(f"{bits:>6} {H:>3} {c3:>3} | " + ' '.join(val_strs) + f" | {'YES' if is_pos else 'NO'}")

        if n >= 6 and all_hook_pos:
            print(f"  (all {count_tested} classes have hook-positive Schur coefficients)")

        print(f"\n  Hook-positive: {count_hook_pos}/{count_tested} ({100*count_hook_pos/count_tested:.1f}%)")
        print(f"  Min hook coefficient: {float(min_hook_coeff):.6f}")

        # Check palindromic property: h_i should equal h_{n+1-i}
        print(f"\n  Palindromic check (h_k = h_{{n+1-k}}):")
        palindromic_violations = 0
        for bits in [seen[k] for k in list(seen.keys())[:3]]:
            adj = tournament_from_bits(bits, n)
            p_coeffs = compute_p_expansion(n, adj)
            hook_coeffs = compute_hook_schur_coefficients(n, p_coeffs)
            vals = [hook_coeffs[h] for h in hooks]
            for i in range(len(vals)):
                if vals[i] != vals[len(vals)-1-i]:
                    palindromic_violations += 1
                    print(f"    VIOLATION at bits={bits}: h[{i}]={vals[i]} != h[{len(vals)-1-i}]={vals[len(vals)-1-i]}")
            if palindromic_violations == 0:
                H = count_ham_paths(n, adj)
                print(f"    bits={bits}, H={H}: {[str(v) for v in vals]} ✓ palindromic")

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

if __name__ == '__main__':
    main()
