#!/usr/bin/env python3
"""
sl_mine_direct_s20ej.py — DIRECT computation of SL_mine
kind-pasteur-2026-03-23-S20ej

SL_mine = sum over iso classes C of #{arc orbits O of C where flipping preserves C}

Also computes:
  D = sum over C of neutral_labeled(C) = #{labeled arcs e in C such that [C flip e] = [C]}
  T = sum over C of #{arc orbits of C} (total orbit-level transitions)
  non_SL = T - SL_mine
  E = #{unordered pairs {C,C'} connected by some flip}
  multi_edge_surplus = non_SL - 2*E

KEY QUESTION: Is SL_mine = T - 2*E? Or is T - 2*E an overestimate?
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DIRECT SL_mine COMPUTATION")
print("  kind-pasteur-2026-03-23-S20ej")
print("=" * 80)

for n in range(3, 8):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    def canon(bits):
        best = None
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_ARCS):
                # Arc (i,j): bit k tells direction.
                # Under permutation p: i->p[i], j->p[j]
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                    if bits & (1 << k):
                        nb |= (1 << nk)
                else:
                    nk = ALL_ARCS.index((pj, pi))
                    if not (bits & (1 << k)):
                        nb |= (1 << nk)
            if best is None or nb < best:
                best = nb
        return best

    # Build canonical map
    canon_map = {}
    classes = defaultdict(list)  # canon -> [bits]
    for bits in range(1 << m):
        cn = canon(bits)
        canon_map[bits] = cn
        classes[cn].append(bits)

    num_classes = len(classes)

    # For each class: find automorphisms, arc orbits, self-loop orbits
    total_SL_mine = 0  # orbit-level self-loops
    total_D = 0  # labeled-level self-loops (= neutral_labeled sum)
    total_T = 0  # total orbit-level transitions
    total_E_set = set()  # edges (unordered class pairs)
    total_multi_edges = 0  # count of orbit-level non-SL transitions

    class_details = []

    for cn in sorted(classes.keys()):
        members = classes[cn]
        rep = cn  # canonical form is the representative

        # Find automorphisms (perms that fix the representative)
        auts = []
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_ARCS):
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                    if rep & (1 << k): nb |= (1 << nk)
                else:
                    nk = ALL_ARCS.index((pj, pi))
                    if not (rep & (1 << k)): nb |= (1 << nk)
            if nb == rep:
                auts.append(p)

        aut_size = len(auts)

        # Find arc orbits under automorphism group
        arc_visited = [False] * m
        arc_orbits = []
        for start_arc in range(m):
            if arc_visited[start_arc]:
                continue
            orbit = set()
            orbit.add(start_arc)
            for p in auts:
                i, j = ALL_ARCS[start_arc]
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                else:
                    nk = ALL_ARCS.index((pj, pi))
                orbit.add(nk)
            # Closure
            changed = True
            while changed:
                changed = False
                new_arcs = set()
                for a in orbit:
                    for p in auts:
                        i, j = ALL_ARCS[a]
                        pi, pj = p[i], p[j]
                        if pi < pj:
                            nk = ALL_ARCS.index((pi, pj))
                        else:
                            nk = ALL_ARCS.index((pj, pi))
                        new_arcs.add(nk)
                if not new_arcs.issubset(orbit):
                    orbit.update(new_arcs)
                    changed = True

            for a in orbit:
                arc_visited[a] = True
            arc_orbits.append(sorted(orbit))

        num_orbits = len(arc_orbits)

        # For each orbit: check if flipping preserves the class
        sl_orbits = 0
        neutral_labeled = 0
        orbit_targets = {}

        for oi, orbit in enumerate(arc_orbits):
            # Pick any arc in the orbit, flip it, check class
            a = orbit[0]
            flipped = rep ^ (1 << a)
            target_cn = canon_map[flipped]
            orbit_targets[oi] = target_cn

            if target_cn == cn:
                sl_orbits += 1
                neutral_labeled += len(orbit)
            else:
                total_E_set.add((min(cn, target_cn), max(cn, target_cn)))

        total_SL_mine += sl_orbits
        total_D += neutral_labeled
        total_T += num_orbits

        class_details.append({
            'cn': cn, 'size': len(members), 'aut': aut_size,
            'orbits': num_orbits, 'sl_orbits': sl_orbits,
            'neutral_labeled': neutral_labeled
        })

    total_E = len(total_E_set)
    non_SL = total_T - total_SL_mine

    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")
    print(f"  Classes: {num_classes}")
    print(f"  T (total arc orbits): {total_T}")
    print(f"  SL_mine (orbit-level self-loops): {total_SL_mine}")
    print(f"  D (labeled self-loops / class): {total_D}")
    print(f"  E (metagraph edges): {total_E}")
    print(f"  non_SL = T - SL_mine = {non_SL}")
    print(f"  2*E = {2*total_E}")
    print(f"  T - 2E = {total_T - 2*total_E}")
    print(f"  multi_edge_surplus = non_SL - 2E = {non_SL - 2*total_E}")
    print(f"  D == SL_mine? {total_D == total_SL_mine}")
    print(f"  T-2E == SL_mine? {total_T - 2*total_E == total_SL_mine}")

    # Check known values
    T_known = {3:4, 4:16, 5:88, 6:704, 7:8912}
    E_known = {3:1, 4:5, 5:30, 6:290, 7:4086}
    if n in T_known:
        print(f"\n  T_known = {T_known[n]}, match: {total_T == T_known[n]}")
    if n in E_known:
        print(f"  E_known = {E_known[n]}, match: {total_E == E_known[n]}")

    # Show classes with |Aut| > 1 that have self-loops
    print(f"\n  Classes with SL orbits:")
    for d in class_details:
        if d['sl_orbits'] > 0:
            tag = " <--- D != SL" if d['neutral_labeled'] != d['sl_orbits'] else ""
            print(f"    |Aut|={d['aut']:3d}, orbits={d['orbits']:2d}, SL={d['sl_orbits']:2d}, neutral_labeled={d['neutral_labeled']:3d}{tag}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

    if elapsed > 120:
        print("  Skipping larger n")
        break

print("\nDONE.")
print("=" * 80)
