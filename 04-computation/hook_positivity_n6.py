#!/usr/bin/env python3
"""
Hook Schur positivity test at n=6.
Instance: opus-2026-03-06-S10

Uses the Murnaghan-Nakayama rule for hooks (which has a clean formula)
to avoid needing the full S_6 character table.

For hook (n-j, 1^j), the MN rule gives:
  chi^{(n-j,1^j)}(mu) = sum over subsets S of cycles of mu
                          with sum of cycle lengths = j
                          of (-1)^{j - |S|}

This IS the correct formula (verified against character tables at n=3,4,5).
The bug in hook_positivity_test.py was a different formula.

But wait - let me re-check. The formula above was in hook_positivity_test.py
and produced wrong results. Let me verify against known tables first.

At n=4, hook (3,1): j=1 (one 1 in the tail).
chi^{(3,1)}((3,1)) = subsets S of {3,1} with sum = 1:
  S = {1}: (-1)^{1-1} = 1
  Total = 1.
But actual chi^{(3,1)}((3,1)) = 0.

So the formula is WRONG. The correct MN formula for hooks is:

For border strip tableaux of shape (n-j, 1^j):
- Remove strips of sizes mu_r, mu_{r-1}, ..., mu_1 from the hook
- At each step, the strip can be placed in the arm or leg
- Sign = (-1)^{height of strip - 1}

For a hook, each strip must be a connected border strip.
If a strip of size k is placed in the arm (row 1), height = 1, sign = +1.
If placed in the leg (column 1), height = k, sign = (-1)^{k-1}.

The key constraint: the strips must cover the entire hook exactly.

Let me implement MN for hooks correctly by recursive enumeration.
"""

from itertools import permutations
from math import factorial
from fractions import Fraction
from collections import defaultdict, Counter

def chi_hook_mn(hook_part, mu):
    """
    Compute chi^{(n-j,1^j)}(mu) using Murnaghan-Nakayama for hooks.

    Shape: hook (n-j, 1^j) has arm of length n-j and leg of length j.
    We remove border strips of sizes mu_1, mu_2, ..., mu_r in order.

    For hooks, at each step the remaining shape is still a hook (or empty).
    A border strip of size k on a hook (a, 1^b) can be:
    - From the arm: removes k cells from arm. New shape: (a-k, 1^b).
      Height = 1. Sign = +1. Only valid if k <= a.
    - From the leg: removes k cells from leg. New shape: (a, 1^{b-k}).
      Height = k. Sign = (-1)^{k-1}. Only valid if k <= b.
    - Corner: removes cells from both arm and leg.
      If we take p from arm and q from leg with p+q=k and p>=1, q>=1:
      New shape: (a-p, 1^{b-q}). But this is only a valid border strip
      if we take ALL remaining arm cells (p=a) or ALL remaining leg cells...
      Actually for a border strip to be connected on a hook:
      - Pure arm: take k cells from position (1,a-k+1) to (1,a). Height=1.
      - Pure leg: take k cells from position (b-k+2,1) to (b+1,1). Height=k.
      - Corner: take the corner cell (1,1) plus some from arm and leg.
        This forms an L-shape. Take p cells from arm (including corner)
        and q cells from leg (not including corner), p+q=k, p>=1.
        Height = q+1. Valid if p <= a and q <= b.

    Wait, I need to be more careful. A border strip is a connected skew shape
    with no 2x2 square. On a hook (a, 1^b):

    The hook has cells: (1,1), (1,2), ..., (1,a) [arm]
                        (2,1), (3,1), ..., (b+1,1) [leg]
    Total = a + b cells.

    Border strips that can be removed from the BORDER of this hook:
    1. Pure arm strip: cells (1, a-k+1) through (1, a). Size k, k <= a-1
       (can't take the corner unless we also take leg cells).
       Actually if k <= a and we take (1,a-k+1)...(1,a), this is valid
       even if it includes the corner (1,1) only when k=a.
       If k < a: remaining shape is (a-k, 1^b). Height = 1. Sign = +1.
       If k = a: remaining shape is (0, 1^b) = (1^b). Height = 1. Sign = +1.
       Wait, but if k=a we remove the entire arm including corner.
       Then remaining is column (2,1)...(b+1,1) which is (1^b). Valid if b >= 0.

    2. Pure leg strip: cells (b+1-k+1, 1) through (b+1, 1). Size k, k <= b.
       Remaining shape is (a, 1^{b-k}). Height = k. Sign = (-1)^{k-1}.

    3. L-shaped strip: corner (1,1) plus some arm cells AND some leg cells.
       Take arm cells (1,1)...(1,p) and leg cells (2,1)...(q+1,1).
       Size = p + q, p >= 1, q >= 1.
       Height = q + 1. Sign = (-1)^q.
       Remaining shape: must be a valid Young diagram.
       If p < a: remaining has arm cells (1,p+1)...(1,a) — but row 1 starts
       at column p+1, while rows 2...(b-q+1) have cells at column 1.
       This is NOT a valid Young diagram unless p = a (all arm removed)
       or q = b (all leg removed).

       Wait: if we remove cells (1,1)...(1,p) and (2,1)...(q+1,1),
       the remaining cells are (1,p+1)...(1,a) and (q+2,1)...(b+1,1).
       For this to be a valid Young diagram (or connected shape that's
       still a valid shape to continue MN), we need...

       Actually this remaining shape is NOT connected (arm fragment above,
       leg fragment below, with the L-corner removed). So it's not a valid
       Young diagram.

    So for hooks, only TWO types of border strips are valid:
    1. Pure arm strip (size k, 1 <= k <= a): height 1, sign +1
       Remaining: (a-k, 1^b) if k < a, or (1^b) if k = a
    2. Pure leg strip (size k, 1 <= k <= b): height k, sign (-1)^{k-1}
       Remaining: (a, 1^{b-k})

    Exception: k = a + b (entire hook). This is one strip of size n.
    Valid only when mu = (n). Height = b+1. Wait, the entire hook as a border
    strip: it goes along the arm then down the leg. Height = b+1.
    Sign = (-1)^b.

    Actually, I realize I should also consider: when k = a, "pure arm strip"
    removes the entire arm, leaving (1^b). Then when k > a... hmm, there's
    also the case where the strip wraps around the corner.

    Let me reconsider: a border strip of size k from the border of hook (a, 1^b):
    - It must be connected and contain no 2x2 square
    - It must include the OUTER border cells

    Actually in MN rule, we remove from the BOTTOM-RIGHT border.
    For a hook, the outer border cells are:
    - Bottom of leg: (b+1, 1)
    - Right of arm: (1, a)
    - And the corner (1,1)

    A border strip of the hook shape is a connected set of cells on the
    border that, when removed, leaves a valid Young diagram.

    Options:
    - Take k cells from right end of arm: remaining (a-k, 1^b), height 1
      Valid for 1 <= k <= a-1 (not the corner)
    - Take k cells from bottom of leg: remaining (a, 1^{b-k}), height k
      Valid for 1 <= k <= b
    - Take all of arm (k1 cells) + some of leg (k2 cells from bottom):
      k1 + k2 = k. Take the rightmost k1 of arm and bottommost k2 of leg.
      Wait, these aren't connected unless we also take the corner.
      Actually: arm cells (1,1)...(1,a) and leg cells (2,1)...(b+1,1).
      If we take (1, a-k1+1)...(1,a) from arm and (b+1-k2+1,1)...(b+1,1) from leg:
      these are only connected if the arm strip reaches the corner (k1=a) OR the
      leg strip reaches the corner... no, (1,a-k1+1)...(1,a) and
      (b+1-k2+1,1)...(b+1,1) are never adjacent unless both include corner.

      Hmm, actually (1,1) connects to (2,1). So arm cell (1,1) is adjacent to
      leg cell (2,1). So if arm strip includes (1,1) (k1=a) and leg strip
      starts at (2,1), they ARE connected via the (1,1)-(2,1) adjacency.

      k1 = a (take entire arm including corner) + k2 cells from top of leg.
      Take cells (1,1)...(1,a) and (2,1)...(k2+1,1).
      Total = a + k2. Remaining: (0, 1^{b-k2}) = (1^{b-k2}).
      Height of this strip: the strip goes from (1,a) left to (1,1) then
      down to (k2+1,1). Height = k2 + 1. Sign = (-1)^{k2}.
      Valid for k2 = 0 (just arm, height 1) through k2 = b (entire hook).

    So the complete list:
    A. k cells from right end of arm, 1 <= k <= a-1: remaining (a-k, 1^b), height 1, sign +1
    B. k cells from bottom of leg, 1 <= k <= b: remaining (a, 1^{b-k}), height k, sign (-1)^{k-1}
    C. entire arm + k2 cells from top of leg, k = a + k2, 0 <= k2 <= b:
       remaining (1^{b-k2}) if k2 < b, or empty if k2 = b
       height = k2 + 1, sign = (-1)^{k2}

    Note: A with k=a is the same as C with k2=0 (take entire arm, remaining (1^b), height 1, sign +1).
    And C with k2=b is the entire hook.

    So: A covers k=1..a (taking from arm, height 1, sign +1, remaining (a-k, 1^b))
        B covers k=1..b (taking from leg, height k, sign (-1)^{k-1}, remaining (a, 1^{b-k}))
        C covers k=a+1..a+b (taking arm + some leg, height k-a+1, sign (-1)^{k-a}, remaining (1^{a+b-k}))

    After removing a strip, the remaining shape is again a hook, so we recurse.
    """
    n = sum(mu)
    arm = hook_part[0]  # a
    leg = n - arm       # b = j

    # mu is a tuple of cycle lengths. We process them one at a time.
    # Recursion: process mu[0], then recurse on remaining mu and remaining shape.

    def mn_recurse(arm, leg, parts_remaining):
        if not parts_remaining:
            if arm == 0 and leg == 0:
                return 1
            else:
                return 0

        k = parts_remaining[0]
        rest = parts_remaining[1:]
        total = 0

        # Option A: take k cells from RIGHT end of arm.
        # When leg > 0: k must be < arm (can't take corner, would orphan leg).
        # When leg = 0: k can be <= arm.
        max_arm_k = arm - 1 if leg > 0 else arm
        if k <= max_arm_k and k >= 1:
            new_arm = arm - k
            new_leg = leg
            sign = 1  # height = 1
            total += sign * mn_recurse(new_arm, new_leg, rest)

        # Option B: take k cells from BOTTOM of leg (k <= leg)
        if k <= leg and k >= 1:
            new_arm = arm
            new_leg = leg - k
            sign = (-1) ** (k - 1)  # height = k
            total += sign * mn_recurse(new_arm, new_leg, rest)

        # Option C: take entire hook (all arm + leg cells).
        # Only valid when leg > 0 (otherwise same as option A with k=arm).
        # Height = leg + 1. Sign = (-1)^leg.
        if k == arm + leg and leg > 0:
            sign = (-1) ** leg  # height = leg + 1
            total += sign * mn_recurse(0, 0, rest)

        return total

    return mn_recurse(arm, leg, list(mu))

# Verify against known character tables
def verify_chi_hook():
    """Verify the MN formula against hardcoded character tables."""
    tables = {
        3: {
            'parts': [(3,), (2,1), (1,1,1)],
            'classes': [(1,1,1), (2,1), (3,)],
            'table': [[1,1,1],[2,0,-1],[1,-1,1]]
        },
        4: {
            'parts': [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)],
            'classes': [(1,1,1,1), (2,1,1), (2,2), (3,1), (4,)],
            'table': [[1,1,1,1,1],[3,1,-1,0,-1],[2,0,2,-1,0],[3,-1,-1,0,1],[1,-1,1,1,-1]]
        },
        5: {
            'parts': [(5,), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)],
            'classes': [(1,1,1,1,1), (2,1,1,1), (2,2,1), (3,1,1), (3,2), (4,1), (5,)],
            'table': [
                [1,1,1,1,1,1,1],[4,2,0,1,-1,0,-1],[5,1,1,-1,1,-1,0],
                [6,0,-2,0,0,0,1],[5,-1,1,-1,-1,1,0],[4,-2,0,1,1,0,-1],
                [1,-1,1,1,-1,-1,1]
            ]
        },
    }

    def is_hook(p):
        return sum(1 for x in p if x > 1) <= 1

    print("Verifying MN hook formula against character tables:")
    all_ok = True
    for n, ct in tables.items():
        for i, lam in enumerate(ct['parts']):
            if not is_hook(lam): continue
            for j, mu in enumerate(ct['classes']):
                expected = ct['table'][i][j]
                computed = chi_hook_mn(lam, mu)
                if computed != expected:
                    print(f"  MISMATCH n={n}: chi^{lam}({mu}) = {computed}, expected {expected}")
                    all_ok = False
    if all_ok:
        print("  All values match!")
    return all_ok

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
    result = []
    for k in range(1, n+1):
        if k == n:
            result.append((n,))
        else:
            result.append((k,) + (1,)*(n-k))
    return result

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

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def main():
    # First verify the formula
    if not verify_chi_hook():
        print("Formula verification FAILED. Aborting.")
        return

    print("\n=== HOOK POSITIVITY TEST AT n=6 ===\n")

    n = 6
    hooks = hook_partitions(n)
    hook_names = [f"s{partition_str(h)}" for h in hooks]
    print(f"Hooks: {', '.join(hook_names)}")

    # First, check which cycle types appear
    print("\nStep 1: Identify tournament-relevant cycle types")
    all_parts = list(generate_partitions(n))
    odd_parts = [mu for mu in all_parts if all(p % 2 == 1 for p in mu)]
    print(f"  Odd-part partitions of 6: {[partition_str(mu) for mu in odd_parts]}")

    # Check hook characters at these types
    print("\nStep 2: Hook characters at odd-part types")
    for h in hooks:
        vals = {}
        for mu in odd_parts:
            chi = chi_hook_mn(h, mu)
            vals[mu] = chi
        neg = [mu for mu, v in vals.items() if v < 0]
        print(f"  chi^{partition_str(h):12s}: {[vals[mu] for mu in odd_parts]}  "
              f"{'ALL >= 0' if not neg else 'NEG at ' + str([partition_str(mu) for mu in neg])}")

    # Now test all n=6 tournaments
    print(f"\nStep 3: Testing all n=6 tournaments")
    m = n*(n-1)//2 - (n-1)
    print(f"  Free bits: {m}, total: {2**m}")

    seen = {}
    count_tested = 0
    count_hook_pos = 0
    min_hook_val = Fraction(1)
    failures = []

    for bits in range(2**m):
        adj = tournament_from_bits(bits, n)
        H = count_ham_paths(n, adj)
        scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))

        key = (scores, H)
        if key in seen: continue
        seen[key] = bits
        count_tested += 1

        p_coeffs = compute_p_expansion(n, adj)

        # Compute hook Schur coefficients
        hook_coeffs = {}
        for h in hooks:
            coeff = Fraction(0)
            for mu in all_parts:
                if mu not in p_coeffs: continue
                chi = chi_hook_mn(h, mu)
                coeff += Fraction(chi) * Fraction(p_coeffs[mu]) / Fraction(z_lambda(mu))
            hook_coeffs[h] = coeff

        is_pos = all(v >= 0 for v in hook_coeffs.values())
        if is_pos:
            count_hook_pos += 1
        else:
            failures.append((bits, H, scores, hook_coeffs))

        for v in hook_coeffs.values():
            if v < min_hook_val:
                min_hook_val = v

        if count_tested % 100 == 0:
            print(f"  ... tested {count_tested}, hook-positive so far: {count_hook_pos}")

    print(f"\n  RESULT: Hook-positive {count_hook_pos}/{count_tested}")
    print(f"  Min hook coefficient: {float(min_hook_val):.6f}")

    if failures:
        print(f"\n  FAILURES:")
        for bits, H, scores, hc in failures[:5]:
            print(f"    bits={bits}, H={H}, scores={scores}")
            for h in hooks:
                if hc[h] < 0:
                    print(f"      [{partition_str(h)}] = {hc[h]} = {float(hc[h]):.6f}")
    else:
        print(f"\n  ALL TOURNAMENTS HOOK-POSITIVE AT n=6!")

    # Check palindromic property
    print(f"\n  Palindromic check (first 3 tournaments):")
    for bits in [seen[k] for k in list(seen.keys())[:3]]:
        adj = tournament_from_bits(bits, n)
        p_coeffs = compute_p_expansion(n, adj)
        vals = []
        for h in hooks:
            coeff = Fraction(0)
            for mu in all_parts:
                if mu not in p_coeffs: continue
                chi = chi_hook_mn(h, mu)
                coeff += Fraction(chi) * Fraction(p_coeffs[mu]) / Fraction(z_lambda(mu))
            vals.append(coeff)

        palindromic = all(vals[i] == vals[len(vals)-1-i] for i in range(len(vals)))
        print(f"    bits={bits}: {[str(v) for v in vals]} {'PALINDROMIC' if palindromic else 'NOT PALINDROMIC'}")

if __name__ == '__main__':
    main()
