#!/usr/bin/env python3
"""
Test the CORRECTED involution conjecture:
    Every SC tournament has at least one involution anti-automorphism.

The original claim (kind-pasteur S18e) that ALL anti-auts are involutions is
FALSE at n=6 (two counterexample classes with |Aut|>1 have order-6 anti-auts).

The CORRECTED claim: at least one anti-aut is an involution (sigma^2=id).
This is what matters for the SC maximizer mechanism, since we only need ONE
involution to create the orbit pairing structure.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import permutations
from collections import defaultdict

def canonical_form(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def opposite_tournament(T):
    n = len(T)
    return [[1 - T[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

def is_self_converse(T):
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

def find_anti_auts(T):
    n = len(T)
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def compose(s, t):
    return tuple(s[t[i]] for i in range(len(s)))

for n in [4, 5, 6, 7]:
    m = n * (n - 1) // 2
    print(f"\n{'='*70}")
    print(f"n={n}: Exhaustive involution-existence test")
    print("=" * 70)

    seen = {}
    has_involution = 0
    no_involution = 0
    total_sc = 0
    all_are_involutions = 0
    some_are_not = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cf = canonical_form(T)
        if cf in seen:
            continue
        seen[cf] = True

        if not is_self_converse(T):
            continue
        total_sc += 1

        anti_auts = find_anti_auts(T)
        identity = tuple(range(n))

        involutions = [s for s in anti_auts if compose(s, s) == identity]
        non_involutions = [s for s in anti_auts if compose(s, s) != identity]

        if involutions:
            has_involution += 1
            if not non_involutions:
                all_are_involutions += 1
            else:
                some_are_not += 1
                # Analyze the structure
                aut_orders = set()
                for s in non_involutions:
                    k = 1
                    x = s
                    while compose(x, x if k == 1 else x) != identity and k < 100:
                        # compute order properly
                        break
                    s2 = compose(s, s)
                    # order of s2 as automorphism
                    y = s2
                    for order in range(1, n+1):
                        if y == identity:
                            aut_orders.add(order)
                            break
                        y = compose(y, s2)
        else:
            no_involution += 1
            print(f"  NO INVOLUTION! bits={bits}")
            print(f"    |AntiAut| = {len(anti_auts)}")
            for s in anti_auts[:3]:
                s2 = compose(s, s)
                k = 1
                x = s
                while True:
                    k += 1
                    x = compose(x, s)
                    if x == identity:
                        break
                    if k > 100:
                        break
                print(f"    sigma={s}, sigma^2={s2}, order={k}")

    print(f"\nn={n} summary:")
    print(f"  SC iso classes: {total_sc}")
    print(f"  Has ≥1 involution anti-aut: {has_involution}")
    print(f"    All anti-auts are involutions: {all_are_involutions}")
    print(f"    Some anti-auts not involutions: {some_are_not}")
    print(f"  No involution anti-aut: {no_involution}")
    if no_involution == 0:
        print(f"  CONJECTURE HOLDS at n={n}: every SC tournament has an involution anti-aut")

print("\nDone.")
