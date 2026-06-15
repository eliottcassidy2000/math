#!/usr/bin/env python3
"""
alpha_lattice_mean_field_gap_codex.py
=====================================
Exact n<=6 audit of the low-level OCF alpha-lattice in the "mean-field board vs
interacting Omega" reframe.

For n<=6:
  - alpha_3 = 0,
  - alpha_1 = total odd cycles = c3 + c5,
  - alpha_2 = number of vertex-disjoint odd-cycle pairs = disjoint 3-3 pairs only.

So H = 1 + 2*alpha_1 + 4*alpha_2 exactly, and the first interaction correction is alpha_2.

The script computes:
  1. exact achievable fibers alpha_2(alpha_1) for n=5,6;
  2. the missing lattice points on the H=7 and H=21 lines;
  3. the finer spectral->interaction fibers alpha_2(c3,c5), showing where the
     spectral mean-field data still leaves interaction holes.

This is the "baby Hodge" version of the permanent gaps:
  - moments / traces = the spectral mean-field shadow,
  - alpha_2 = the first interaction correction,
  - forbidden H values = missing alpha-lattice points.

Author: codex-2026-06-15.
"""
from collections import defaultdict
from itertools import combinations, permutations


def tournament_from_bits(n, bits):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def is_directed_3cycle(A, triple):
    i, j, k = triple
    s = A[i][j] + A[j][k] + A[k][i]
    return s == 3 or s == 0


def c3_and_triangles(A):
    n = len(A)
    triangles = []
    for triple in combinations(range(n), 3):
        if is_directed_3cycle(A, triple):
            triangles.append(triple)
    return len(triangles), triangles


def c5(A):
    n = len(A)
    total = 0
    for subset in combinations(range(n), 5):
        start = subset[0]
        others = subset[1:]
        count = 0
        for perm in permutations(others):
            cyc = (start,) + perm
            ok = True
            for i in range(5):
                if A[cyc[i]][cyc[(i + 1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                count += 1
        total += count
    return total


def alpha2_n_le_6(triangles):
    tri_set = set(triangles)
    verts = sorted({v for tri in triangles for v in tri})
    if not verts:
        return 0
    n = max(verts) + 1
    total = 0
    for triple in combinations(range(n), 3):
        comp = tuple(sorted(set(range(n)) - set(triple)))
        if len(comp) != 3:
            continue
        triple = tuple(sorted(triple))
        if triple < comp and triple in tri_set and comp in tri_set:
            total += 1
    return total


def audit_n(n):
    fibers_alpha = defaultdict(set)
    fibers_spec = defaultdict(set)
    h_values = set()
    total = 1 << (n * (n - 1) // 2)
    for bits in range(total):
        A = tournament_from_bits(n, bits)
        cc3, triangles = c3_and_triangles(A)
        cc5 = c5(A)
        a1 = cc3 + cc5
        a2 = alpha2_n_le_6(triangles)
        h = 1 + 2 * a1 + 4 * a2
        fibers_alpha[a1].add(a2)
        fibers_spec[(cc3, cc5)].add(a2)
        h_values.add(h)
    return fibers_alpha, fibers_spec, sorted(h_values)


def format_h_line(target_r):
    out = []
    for a1 in range(target_r, -1, -2):
        a2 = (target_r - a1) // 2
        out.append((a1, a2))
    return out


def main():
    print("=" * 72)
    print(" ALPHA-LATTICE / MEAN-FIELD GAP AUDIT")
    print("=" * 72)
    print()
    print("For n<=6: H = 1 + 2*alpha_1 + 4*alpha_2 exactly.")
    print("The free mean-field envelope at fixed alpha_1=m is 0 <= alpha_2 <= C(m,2).")
    print("The realized tournament lattice is much thinner.")
    print()

    for n in (5, 6):
        alpha_fibers, spec_fibers, h_values = audit_n(n)
        print(f"n={n}")
        print(f"  Distinct H values: {h_values}")
        print("  Achievable alpha_2 fibers at fixed alpha_1:")
        for a1 in sorted(alpha_fibers):
            vals = sorted(alpha_fibers[a1])
            envelope_max = a1 * (a1 - 1) // 2
            holes = [v for v in range(vals[0], vals[-1] + 1) if v not in alpha_fibers[a1]]
            print(
                f"    alpha_1={a1:2d}: alpha_2 in {vals}"
                f"  | mean-field envelope [0,{envelope_max}]"
                + (f"  | interior holes {holes}" if holes else "")
            )

        print("  Spectral -> interaction fibers (c3,c5) -> alpha_2:")
        for (cc3, cc5) in sorted(spec_fibers):
            vals = sorted(spec_fibers[(cc3, cc5)])
            print(f"    (c3,c5)=({cc3:2d},{cc5:2d}) -> alpha_2 in {vals}")
        print()

    print("Target H-lines in alpha-coordinates:")
    for h in (7, 21):
        r = (h - 1) // 2
        print(f"  H={h}: alpha_1 + 2*alpha_2 = {r} -> candidates {format_h_line(r)}")
    print()

    n6_alpha, _, h6 = audit_n(6)
    print("n=6 exact verdict on the first permanent gaps:")
    for h in (7, 21):
        r = (h - 1) // 2
        candidates = format_h_line(r)
        missing = [(a1, a2) for (a1, a2) in candidates if a2 not in n6_alpha.get(a1, set())]
        print(f"  H={h}: all candidate alpha-points are missing -> {missing}")
    print()

    print("Interpretation:")
    print("  - alpha_1 is the mean-field odd-cycle census (first spectral shadow).")
    print("  - alpha_2 is the first interaction correction: disjoint-cycle compatibility.")
    print("  - H=7 and H=21 fail because their required low-level alpha vectors are absent,")
    print("    not because of a scalar congruence accident.")
    print("  - This is the alpha-lattice version of the baby-Hodge hole picture.")


if __name__ == "__main__":
    main()
