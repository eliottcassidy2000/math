#!/usr/bin/env python3
"""
beta2_common_inneighbor.py — Test "common in-neighbor" condition for Z₂ filling

KEY INSIGHT FROM CONE ANALYSIS:
- Type 1 defects (v,a,c) with c→a ALWAYS vanish (Ω₂ constraint on z). PROVED.
- Type 2 defects (v,b,c) with b→v vanish iff v→b.

So if we cone path (a,b,c) from vertex v with BOTH v→a AND v→b,
there's NO defect. The resulting (v,a,b,c) is automatically in Ω₃!

Question: For each path (a,b,c) in a Z₂ cycle, does there always exist
a vertex v ∉ {a,b,c} with v→a AND v→b?

Such v ∈ in(a) ∩ in(b) ∩ (V\{a,b,c}) = "common in-neighbor of a and b, different from c"

If YES for all Z₂ support paths → path-dependent coning proves β₂=0.

ALSO: Even if some paths lack such v, the filling might still work if
the defect terms cancel across the Z₂ cycle.

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


print("=" * 70)
print("COMMON IN-NEIGHBOR CONDITION FOR Z₂ FILLING")
print("=" * 70)

for n in [4, 5, 6]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    if n == 6:
        import random
        random.seed(42)
        samples = random.sample(range(1 << m), 200)
    else:
        samples = range(1 << m)

    total_cycles = 0
    total_paths_in_cycles = 0
    paths_with_cin = 0  # paths (a,b,c) with ∃v: v→a, v→b, v∉{a,b,c}
    paths_without_cin = 0
    cycles_all_cin = 0
    cycles_some_no_cin = 0

    # Also: for paths with v→a,v→b,v→c (DT condition: all three coned)
    paths_with_full_cone = 0

    for bits in samples:
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        d2 = dim_om(om2)
        if d2 == 0:
            continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk2
        if z2_dim == 0:
            continue

        z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
        ap2_list = [tuple(x) for x in ap2]

        for z_idx in range(z2_dim):
            z_om2 = z2_basis[:, z_idx]
            z_a2 = om2 @ z_om2

            total_cycles += 1
            all_cin = True

            for i, (a,b,c) in enumerate(ap2_list):
                if abs(z_a2[i]) < 1e-10:
                    continue

                total_paths_in_cycles += 1

                # Check: ∃v ∈ V\{a,b,c} with v→a AND v→b?
                cin = [v for v in range(n) if v not in (a,b,c) and A[v][a] and A[v][b]]

                if cin:
                    paths_with_cin += 1
                else:
                    paths_without_cin += 1
                    all_cin = False

                # Check full cone: v→a, v→b, v→c
                full_cin = [v for v in cin if A[v][c]]
                if full_cin:
                    paths_with_full_cone += 1

            if all_cin:
                cycles_all_cin += 1
            else:
                cycles_some_no_cin += 1

    print(f"\n{'='*50}")
    print(f"n = {n} ({len(samples)} tournaments)")
    print(f"{'='*50}")
    print(f"Total Z₂ basis elements: {total_cycles}")
    print(f"Total paths in Z₂ support: {total_paths_in_cycles}")
    print(f"Paths WITH common in-neighbor: {paths_with_cin}")
    print(f"Paths WITHOUT common in-neighbor: {paths_without_cin}")
    print(f"Paths with FULL cone (v→a,v→b,v→c): {paths_with_full_cone}")
    print(f"Z₂ cycles where ALL paths have CIN: {cycles_all_cin}/{total_cycles}")
    print(f"Z₂ cycles with SOME paths missing CIN: {cycles_some_no_cin}/{total_cycles}")

    if paths_without_cin == 0:
        print("✓ EVERY Z₂ support path has a common in-neighbor!")
        print("⟹ Path-dependent coning proves β₂=0 at this n!")
    else:
        pct = 100 * paths_without_cin / total_paths_in_cycles if total_paths_in_cycles > 0 else 0
        print(f"✗ {pct:.1f}% of Z₂ support paths lack common in-neighbor")

# Now: for paths without CIN, analyze what's happening
print(f"\n{'='*70}")
print("DETAILED ANALYSIS: PATHS WITHOUT COMMON IN-NEIGHBOR")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

no_cin_examples = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
    ap2_list = [tuple(x) for x in ap2]

    for z_idx in range(z2_dim):
        z_om2 = z2_basis[:, z_idx]
        z_a2 = om2 @ z_om2

        for i, (a,b,c) in enumerate(ap2_list):
            if abs(z_a2[i]) < 1e-10: continue
            cin = [v for v in range(n) if v not in (a,b,c) and A[v][a] and A[v][b]]
            if not cin and len(no_cin_examples) < 10:
                # What vertices are available?
                d_out_a = sum(A[a])
                d_out_b = sum(A[b])
                in_a = [v for v in range(n) if v != a and A[v][a]]
                in_b = [v for v in range(n) if v != b and A[v][b]]
                no_cin_examples.append({
                    'bits': bits, 'path': (a,b,c), 'coeff': z_a2[i],
                    'in_a': in_a, 'in_b': in_b, 'd_out_a': d_out_a, 'd_out_b': d_out_b,
                    'scores': tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
                })

print(f"\nExamples of paths in Z₂ without common in-neighbor:")
for ex in no_cin_examples:
    a,b,c = ex['path']
    print(f"\n  T#{ex['bits']} scores={ex['scores']}")
    print(f"  Path ({a},{b},{c}): coeff={ex['coeff']:.4f}")
    print(f"  d_out(a={a})={ex['d_out_a']}, d_out(b={b})={ex['d_out_b']}")
    print(f"  in(a)={ex['in_a']}, in(b)={ex['in_b']}")
    print(f"  in(a) ∩ in(b) \\ {{a,b,c}} = {set(ex['in_a']) & set(ex['in_b']) - {a,b,c}}")
    # Why is this empty?
    # in(a) = vertices beating a. in(b) = vertices beating b.
    # a→b (given), so a ∈ in(b) but a ∉ eligible (a is the path start).
    # Common in-neighbors must beat BOTH a and b.
    # If a has high out-degree, few beat a. If b also has moderate out-degree...

# NEW TEST: "modified cone" — cone from RIGHT instead of LEFT
# h'_v(a,b,c) = (a,b,c,v) with c→v
# Then face (a,b,v): a→b ✓, b→v needed.
# Face (b,c,v): b→c ✓, c→v ✓.
# Face (a,c,v): a→c needed (type 1 issue).
#
# So for right cone: need c→v AND b→v (common out-neighbor of b and c)
# to avoid type 2 defect from face (a,b,v).
#
# Or more flexible: use LEFT cone for some paths, RIGHT cone for others.

print(f"\n{'='*70}")
print("LEFT vs RIGHT CONE COVERAGE")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

total_cycles = 0
cycles_covered_left = 0  # all paths have left-cone vertex
cycles_covered_right = 0  # all paths have right-cone vertex
cycles_covered_either = 0  # all paths have at least one (left or right)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T
    ap2_list = [tuple(x) for x in ap2]

    for z_idx in range(z2_dim):
        z_om2 = z2_basis[:, z_idx]
        z_a2 = om2 @ z_om2

        total_cycles += 1
        all_left = True
        all_right = True
        all_either = True

        for i, (a,b,c) in enumerate(ap2_list):
            if abs(z_a2[i]) < 1e-10: continue

            # Left cone: v→a AND v→b, v∉{a,b,c}
            left = any(v not in (a,b,c) and A[v][a] and A[v][b] for v in range(n))
            # Right cone: c→v AND b→v, v∉{a,b,c}
            right = any(v not in (a,b,c) and A[c][v] and A[b][v] for v in range(n))

            if not left: all_left = False
            if not right: all_right = False
            if not (left or right): all_either = False

        if all_left: cycles_covered_left += 1
        if all_right: cycles_covered_right += 1
        if all_either: cycles_covered_either += 1

print(f"Total Z₂ cycles: {total_cycles}")
print(f"All paths have LEFT cone: {cycles_covered_left}/{total_cycles}")
print(f"All paths have RIGHT cone: {cycles_covered_right}/{total_cycles}")
print(f"All paths have LEFT or RIGHT cone: {cycles_covered_either}/{total_cycles}")

if cycles_covered_either == total_cycles:
    print("\n✓ Every Z₂ path can be coned from LEFT or RIGHT!")
else:
    print(f"\n✗ {total_cycles - cycles_covered_either} cycles have paths with NEITHER cone")

    # What about "middle insertion"? (a,v,b,c) or (a,b,v,c)?
    # (a,v,b,c): a→v, v→b, b→c. Face issues: (a,b,c)-(a,v,c)+(a,v,b)-(v,b,c)... no wait
    # ∂₃(a,v,b,c) = (v,b,c)-(a,b,c)+(a,v,c)-(a,v,b)
    # Face (a,b,c) with sign -1. That's what we want.
    # Face (v,b,c): needs v→b ✓ (a→v→b). In A₂? v→b and b→c. ✓
    # Face (a,v,c): needs a→v ✓ and v→c? If c→v, then NOT in A₂.
    # Face (a,v,b): needs a→v ✓ and v→b ✓. In A₂.
    #
    # So middle-insert (a,v,b,c) with a→v, v→b, b→c:
    # Type 1 defect: (a,v,c) when c→v — vanishes by same Ω₂ argument? NO...
    # The Ω₂ constraint applies to the ORIGINAL z's non-arc terms.
    # For middle insertion, the non-A₂ face is (a,v,c) where v→c might not hold.
    # This is a DIFFERENT non-arc, not related to z's Ω₂ constraint.

    # So middle insertion has its own issues. Let's not go there.
    pass

print("\nDone.")
