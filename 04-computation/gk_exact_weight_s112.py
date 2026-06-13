#!/usr/bin/env python3
"""
gk_exact_weight_s112.py — Exact weight per cluster configuration
kind-pasteur-2026-03-15-S112

From gk_cluster_corr_s112.py data:
- Separated pairs: 2^r/(n)_{2r}
- Contiguous block L: 2/(n)_L
- Mixed cluster (L1,...,Ls): weight = ?

Hypothesis: weight = prod_i w(Li) * 1/(n)_{sum Li}
where w(Li) = 2/2^{Li/2} = 2^{1-Li/2} for a cluster of size Li

Actually from the data:
- (2,): ratio to 2/(n)_2 = 1. Weight factor = 2.
- (4,): ratio to 2^2/(n)_4 = 1/2. So actual = 2/(n)_4. Weight factor = 2.
- (6,): ratio to 2^3/(n)_6 = 1/4. So actual = 2/(n)_6. Weight factor = 2.
- (2,2) gap=any: ratio to 2^2/(n)_4 = 1. So actual = 4/(n)_4. Weight factor = 4.
- (2,4) gap=1: ratio to 2^3/(n)_6 = 1/2. So actual = 4/(n)_6. Weight factor = 4.
- (4,2) gap=1: same = 4/(n)_6. Weight factor = 4.

Pattern: weight factor = 2^(number of clusters + number of contiguous-pair-overlaps?)
No: single cluster of any size -> factor 2. Two separated clusters -> factor 4 = 2^2.

CLEAN: weight = 2^(number of connected components) / (n)_{total size}

So E[prod Z_j for positions in config] = 2^c / (n)_L
where c = number of connected components, L = total number of positions.

Let me verify this.
"""

from itertools import permutations, combinations
from fractions import Fraction
from math import factorial
from collections import defaultdict

def falling_factorial(n, k):
    r = 1
    for i in range(k):
        r *= (n - i)
    return r

def compute_expectations_n8():
    """Return all nonzero E[prod Z_j] for domino subsets at n=8."""
    n = 8
    m = n - 1
    nfact = factorial(n)
    adj_pairs = [(j, j+1) for j in range(m - 1)]

    z_all = []
    for perm in permutations(range(n)):
        z = []
        for j in range(m):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z.append(x - y)
        z_all.append(z)

    results = {}
    max_r = 3

    for r in range(1, max_r + 1):
        for combo in combinations(range(len(adj_pairs)), r):
            chosen = [adj_pairs[i] for i in combo]
            positions = set()
            for p in chosen:
                positions.update(p)
            S = frozenset(positions)
            if S in results:
                continue
            total = Fraction(0)
            for z in z_all:
                prod = 1
                for j in S:
                    prod *= z[j]
                total += prod
            results[S] = total / nfact

    return results

def get_components(S):
    """Get connected components of position set S."""
    positions = sorted(S)
    if not positions:
        return []
    components = []
    current = [positions[0]]
    for p in positions[1:]:
        if p == current[-1] + 1:
            current.append(p)
        else:
            components.append(current)
            current = [p]
    components.append(current)
    return components

print("TESTING: E[prod Z_j] = 2^c / (n)_L")
print("where c = #components, L = #positions")
print("="*60)

for n in [6, 7, 8]:
    print(f"\nn={n}:")

    m = n - 1
    nfact = factorial(n)
    adj_pairs = [(j, j+1) for j in range(m - 1)]

    z_all = []
    for perm in permutations(range(n)):
        z = []
        for j in range(m):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z.append(x - y)
        z_all.append(z)

    max_r = min(3, len(adj_pairs))

    for r in range(1, max_r + 1):
        for combo in combinations(range(len(adj_pairs)), r):
            chosen = [adj_pairs[i] for i in combo]
            positions = set()
            for p in chosen:
                positions.update(p)
            S = frozenset(positions)

            total = Fraction(0)
            for z in z_all:
                prod = 1
                for j in S:
                    prod *= z[j]
                total += prod
            E = total / nfact
            if E == 0:
                continue

            comps = get_components(S)
            c = len(comps)
            L = len(S)

            pred = Fraction(2**c, falling_factorial(n, L))
            match = "OK" if E == pred else f"FAIL (E={E}, pred={pred}, ratio={E/pred})"

            comp_sizes = tuple(len(comp) for comp in comps)
            print(f"  S={sorted(S)}, comps={comp_sizes}, c={c}, L={L}: {match}")

# Now derive g_k from this formula
print("\n" + "="*60)
print("DERIVE g_k FROM WEIGHT FORMULA")
print("="*60)
print("E_{2k}/E_0 = sum over configs of r pairs from m=n-2 slots")
print("         = sum over domino subsets S with |S|=2k of 2^c/(n)_L")
print("where c=#components, L=|S|=2k")
print("So factor out 1/(n)_{2k}: g_k(m)/(n)_{2k} * 2 = sum 2^c/(n)_{2k}")
print("=> contribution = (sum 2^{c-1}) * 2/(n)_{2k}")
print("=> g_k(m) = sum over configs of k pairs in m+2k-1 positions of 2^{c-1}")
print()
print("Wait, that's not right. Let me reconsider.")
print()
print("Actually: the sum over all k-pair domino configs in positions 0..n-2")
print("of E[prod Z] = sum of 2^c/(n)_{|S|}.")
print("But |S| is NOT always 2k! Overlapping pairs share positions.")
print()

# Let me just directly compute g_k by counting
# g_k(m) is defined by: sum over configs = 2*g_k(m)/(n)_{2k}
# where m = n - 2k.
# But this assumes |S| = 2k for all configs (pairs don't overlap).
# When pairs DO overlap (contiguous clusters), |S| < 2k.

# So the actual formula is:
# Total = sum over k-pair configs of 2^c / (n)_{|S|}
# NOT simply 2*g_k(m)/(n)_{2k}

# This means THM-216's g_k is NOT just a count — it absorbs the
# different denominators from different |S| values.

# Let me compute more carefully.

print("CAREFUL DECOMPOSITION:")
print()

for n in range(5, 9):
    m = n - 1  # positions 0..m-1
    adj_pairs = [(j, j+1) for j in range(m - 1)]

    # For each k (number of pairs chosen)
    for k in range(1, min(4, len(adj_pairs) + 1)):
        total = Fraction(0)

        by_type = defaultdict(lambda: Fraction(0))

        for combo in combinations(range(len(adj_pairs)), k):
            chosen = [adj_pairs[i] for i in combo]
            positions = set()
            for p in chosen:
                positions.update(p)
            S = frozenset(positions)
            L = len(S)
            comps = get_components(S)
            c = len(comps)
            comp_sizes = tuple(len(comp) for comp in comps)

            weight = Fraction(2**c, falling_factorial(n, L))
            total += weight
            by_type[comp_sizes] += weight

        # Compare to 2*g_k(m)/(n)_{2k}
        m_val = n - 2*k
        if k == 1:
            gk = m_val
        elif k == 2:
            gk = m_val**2
        elif k == 3:
            gk = Fraction(2*m_val**3 + m_val, 3)
        else:
            gk = None

        if gk is not None:
            pred = 2 * gk / falling_factorial(n, 2*k)
            match = "OK" if total == pred else f"FAIL"
        else:
            pred = "?"
            match = "?"

        print(f"n={n}, k={k}: total={float(total):.8f}, "
              f"pred={float(pred) if pred != '?' else '?':.8f}, {match}")

        for cs, w in sorted(by_type.items()):
            print(f"  type {cs}: {float(w):.10f}")

print("\nDone!")
