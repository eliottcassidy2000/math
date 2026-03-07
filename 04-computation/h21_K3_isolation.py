#!/usr/bin/env python3
"""
Can 3 pairwise-intersecting 3-cycles be isolated in Omega?

For I(C,2)=7: need K_3 in Omega (3 cycles, all pairwise adjacent, no other
cycle adjacent to any of them).

Two cases:
A) Common vertex: all 3 cycles share a vertex v. THM-029 handles this.
B) No common vertex: {a,b,c}, {a,d,e}, {b,d,f} (sunflower-free).

For case B: check whether the 6 vertices always create additional 3-cycles
or 5-cycles that connect to the component.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations, permutations
from collections import defaultdict


def check_isolated_K3(n_max=8):
    """For each n, find tournaments where some triple of pairwise-intersecting
    3-cycles is isolated (no other 3-cycle shares a vertex with any of them)."""

    for n in range(5, n_max+1):
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(edges)
        total = 1 << m

        isolated_K3_count = 0
        examined = 0

        for bits in range(total):
            adj_bits = [0]*n
            A = [[0]*n for _ in range(n)]
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    A[j][i] = 1
                    adj_bits[j] |= (1 << i)
                else:
                    A[i][j] = 1
                    adj_bits[i] |= (1 << j)

            # Find all 3-cycles
            cycles = []
            for a in range(n):
                for b in range(a+1, n):
                    for c in range(b+1, n):
                        if (A[a][b] and A[b][c] and A[c][a]) or \
                           (A[a][c] and A[c][b] and A[b][a]):
                            cycles.append(frozenset({a, b, c}))

            if len(cycles) < 3:
                continue

            examined += 1

            # Check all triples of pairwise-intersecting 3-cycles
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if not (cycles[i] & cycles[j]):
                        continue
                    for k in range(j+1, len(cycles)):
                        if not (cycles[i] & cycles[k]):
                            continue
                        if not (cycles[j] & cycles[k]):
                            continue

                        # Triple (i,j,k) is pairwise intersecting
                        triple_verts = cycles[i] | cycles[j] | cycles[k]

                        # Check isolation: no OTHER 3-cycle shares a vertex
                        isolated = True
                        for other in range(len(cycles)):
                            if other in (i, j, k):
                                continue
                            if cycles[other] & triple_verts:
                                isolated = False
                                break

                        if isolated:
                            # Also check: no 5-cycle through these vertices
                            has_5cycle = False
                            for verts5 in combinations(range(n), 5):
                                vset5 = frozenset(verts5)
                                if not (vset5 & triple_verts):
                                    continue
                                # Check if this 5-set has a directed cycle
                                v0 = verts5[0]
                                rest = verts5[1:]
                                for p in permutations(rest):
                                    seq = (v0,) + p
                                    ok = all(A[seq[q]][seq[(q+1)%5]] for q in range(5))
                                    if ok:
                                        has_5cycle = True
                                        break
                                if has_5cycle:
                                    break

                            if not has_5cycle:
                                isolated_K3_count += 1
                                common = cycles[i] & cycles[j] & cycles[k]
                                print(f"  n={n}: ISOLATED K3! bits={bits}")
                                print(f"    cycles: {set(cycles[i])}, {set(cycles[j])}, {set(cycles[k])}")
                                print(f"    common vertex: {set(common) if common else 'NONE'}")
                                print(f"    all verts used: {set(triple_verts)}")

        print(f"n={n}: {examined} tournaments with >=3 cycles, {isolated_K3_count} isolated K3 found")


if __name__ == '__main__':
    check_isolated_K3(7)
