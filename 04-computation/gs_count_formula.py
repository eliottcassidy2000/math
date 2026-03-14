"""
gs_count_formula.py -- kind-pasteur-2026-03-14-S78
Find a formula for #GS tilings per isomorphism class.

KEY OBSERVATIONS from S77 data:
  n=5: GS counts per SC class: [1, 0, 0, 1, 1, 0, 3, 0, 3, 1, 3, 3]
  Classes with GS > 0 are the "GS-reachable" classes.
  Not all SC classes have GS tilings!

HYPOTHESIS: #GS tilings in class C =
  |{sigma in S_n : sigma.T is self-converse under tau AND sigma.T has base path P_0}|
  where tau: v -> n+1-v.

Actually: #GS tilings = #tilings in C that are GS
  = #{bit strings b : tournament(b) is in class C AND b is GS}

A GS tiling corresponds to a tournament that:
1. Contains the base path P_0: n -> n-1 -> ... -> 1
2. Is self-converse under v -> n+1-v

The GS count per class C is:
  #GS(C) = #{sigma : sigma.T in C, sigma.T contains P_0, sigma.T = tau(sigma.T)}

This is a counting problem involving:
- The automorphism group Aut(T) of the class representative T
- The normalizer of tau in S_n
- The stabilizer of the base path P_0

ALSO: Why do NSC classes have 0 GS tilings?
Because a GS tiling is self-converse, so it can only be in an SC class.

Let me verify this and compute the formula.
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def build_tiling_data(n):
    tiles = []
    for b in range(1, n-1):
        for a in range(b+2, n+1):
            tiles.append((a, b))
    m = len(tiles)
    tile_idx = {t: i for i, t in enumerate(tiles)}

    trans_map = []
    for a, b in tiles:
        a2, b2 = n+1-b, n+1-a
        if a2 < b2: a2, b2 = b2, a2
        trans_map.append(tile_idx.get((a2, b2), -1))

    def is_gs(mask):
        for i in range(m):
            j = trans_map[i]
            if j != i and ((mask >> i) & 1) != ((mask >> j) & 1):
                return False
        return True

    def bits_to_adj(mask):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for idx, (a, b) in enumerate(tiles):
            a0, b0 = a-1, b-1
            if (mask >> idx) & 1: A[b0][a0] = 1
            else: A[a0][b0] = 1
        return A

    perms = list(permutations(range(n)))
    def canonicalize(A):
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    mask_data = {}
    for mask in range(1 << m):
        A = bits_to_adj(mask)
        canon = canonicalize(A)
        gs = is_gs(mask)
        mask_data[mask] = {'adj': A, 'canon': canon, 'gs': gs}

    groups = defaultdict(list)
    for mask, d in mask_data.items():
        groups[d['canon']].append(mask)

    class_list = sorted(groups.keys())
    class_idx = {c: i for i, c in enumerate(class_list)}

    return {'n': n, 'm': m, 'tiles': tiles, 'trans_map': trans_map,
            'mask_data': mask_data, 'groups': groups,
            'class_list': class_list, 'class_idx': class_idx}

def main():
    print("=" * 70)
    print("GS COUNT PER CLASS — FORMULA DERIVATION")
    print("kind-pasteur-2026-03-14-S78")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        data = build_tiling_data(n)
        m = data['m']
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        # For each class: compute GS count, automorphism group size, SC status
        tau = list(range(n-1, -1, -1))  # tau: v -> n-1-v (0-indexed)

        perms = list(permutations(range(n)))
        results = []

        for canon in data['class_list']:
            masks = data['groups'][canon]
            size = len(masks)
            aut_size = math.factorial(n) // size

            # GS count
            gs_count = sum(1 for mask in masks if data['mask_data'][mask]['gs'])

            # Is SC? Check if tau(T) is isomorphic to T
            rep_A = data['mask_data'][masks[0]]['adj']
            tau_A = [[rep_A[tau[i]][tau[j]] for j in range(n)] for i in range(n)]
            tau_canon = None
            for p in perms:
                s = ''.join(str(tau_A[p[i]][p[j]]) for i in range(n) for j in range(n))
                if tau_canon is None or s < tau_canon: tau_canon = s

            is_sc = (tau_canon == canon)

            # Is T self-converse (= T^op isomorphic to T)?
            op_A = [[rep_A[j][i] for j in range(n)] for i in range(n)]
            op_canon = None
            for p in perms:
                s = ''.join(str(op_A[p[i]][p[j]]) for i in range(n) for j in range(n))
                if op_canon is None or s < op_canon: op_canon = s

            is_self_converse = (op_canon == canon)

            # Count anti-automorphisms (sigma with sigma.T = T^op under vertex perm)
            anti_aut_count = 0
            for p in perms:
                # Check: p applied to T gives T^op?
                # i.e., T[p[i]][p[j]] = T[j][i] for all i,j
                is_anti = all(rep_A[p[i]][p[j]] == rep_A[j][i] for i in range(n) for j in range(n))
                if is_anti:
                    anti_aut_count += 1

            ci = data['class_idx'][canon]
            results.append({
                'ci': ci, 'size': size, 'aut': aut_size,
                'gs': gs_count, 'sc': is_sc, 'self_conv': is_self_converse,
                'anti_aut': anti_aut_count,
            })

        print(f"  {'CI':>3} {'Size':>5} {'|Aut|':>5} {'GS':>3} {'SC':>3} {'SConv':>5} {'AntiAut':>7} {'GS*|Aut|':>8}")
        for r in results:
            sc_str = 'Y' if r['sc'] else ''
            conv_str = 'Y' if r['self_conv'] else ''
            print(f"  {r['ci']:3d} {r['size']:5d} {r['aut']:5d} {r['gs']:3d} {sc_str:>3s} {conv_str:>5s} "
                  f"{r['anti_aut']:7d} {r['gs']*r['aut']:8d}")

        # KEY: Is GS > 0 iff SC (self-converse)?
        gs_iff_sc = all((r['gs'] > 0) == r['self_conv'] for r in results)
        print(f"\n  GS > 0 iff self-converse? {gs_iff_sc}")

        # Actually check: GS > 0 iff SC (under our tau symmetry)?
        gs_iff_sc_tau = all((r['gs'] > 0) == r['sc'] for r in results)
        print(f"  GS > 0 iff SC (under tau)? {gs_iff_sc_tau}")

        # Check: GS * |Aut| formula?
        # Does GS * |Aut| = anti_aut for SC classes?
        for r in results:
            if r['gs'] > 0:
                ratio = r['anti_aut'] / r['gs'] if r['gs'] > 0 else 'N/A'
                print(f"    class {r['ci']}: anti_aut / GS = {ratio}")

        # FORMULA: GS count = anti_aut_count / n! * size ?
        # Actually by orbit-stabilizer: #GS tilings in class C =
        # #{sigma : sigma(T) contains P_0 AND sigma(T) = (sigma(T))^{tau}} / ...
        # This is complicated. Let me just look for a simple pattern.

        # Does GS = anti_aut * size / n! ?
        print(f"\n  Checking: GS = anti_aut * size / n!?")
        for r in results:
            if r['gs'] > 0:
                predicted = r['anti_aut'] * r['size'] // math.factorial(n)
                match = (predicted == r['gs'])
                print(f"    class {r['ci']}: predicted={predicted}, actual={r['gs']}, match={match}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
