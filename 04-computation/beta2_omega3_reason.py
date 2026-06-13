#!/usr/bin/env python3
"""beta2_omega3_reason.py - WHY is cone filling automatically in Omega_3?

Omega_3 = ker(d_2 composed with projection A_3 -> A_2/Omega_2)
equivalently: w in Omega_3 iff for every non-allowed pair (u,w),
sum of coefficients of 3-paths having (u,w) as a face = 0.

For the cone filling w = sum alpha_{abc} * [(v,a,b,c) + (a,b,c,v)]:
- Faces of (v,a,b,c): (a,b,c), (v,b,c), (v,a,c), (v,a,b)
- Faces of (a,b,c,v): (b,c,v), (a,c,v), (a,b,v), (a,b,c)

The T'-internal face (a,b,c) appears with coeff +1 from front and -1 from back,
so it cancels.

For Omega_3, we need: for each non-allowed 2-path (x,y,z),
sum_{3-paths p having (x,y,z) as face} sign * w[p] = 0.

Test: which non-allowed 2-paths get contributions from the cone,
and do the contributions automatically cancel?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


print("=" * 70)
print("WHY IS CONE FILLING IN OMEGA_3?")
print("=" * 70)
print()
print("For each non-allowed 2-path (x,y,z), check what contributions")
print("the cone 3-paths make. A 3-path (p0,p1,p2,p3) has 2-faces:")
print("  +face_0 = (p1,p2,p3)")
print("  -face_1 = (p0,p2,p3)")
print("  +face_2 = (p0,p1,p3)")
print("  -face_3 = (p0,p1,p2)")

# Boundary of a 3-path
def boundary_faces_3(path):
    """Return list of (face, sign) for a 3-path."""
    p0, p1, p2, p3 = path
    return [
        ((p1, p2, p3), +1),
        ((p0, p2, p3), -1),
        ((p0, p1, p3), +1),
        ((p0, p1, p2), -1),
    ]


n = 5
example_count = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_set = set(tuple(p) for p in paths2)
    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    # Find non-allowed 2-paths
    all_2paths = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] == 1 and A[j][k] == 1:
                    all_2paths.append((i, j, k))

    non_allowed = [p for p in all_2paths if p not in path2_set]
    if not non_allowed:
        continue

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]
        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        if not rows:
            continue
        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C
        if ker_dim == 0:
            continue

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        if not valid_Tp:
            continue

        # For EACH non-allowed 2-path, track contributions from
        # each T' path's cone pair
        if example_count < 3:
            example_count += 1
            print(f"\nT#{tidx}, v={v}, P={P}, Q={Q}")
            print(f"  Non-allowed 2-paths: {len(non_allowed)}")
            print(f"  Valid T' paths: {len(valid_Tp)}")

            # For each valid T' path, list which non-allowed faces it creates
            for abc in valid_Tp[:3]:
                a, b, c = abc
                front = (v, a, b, c)
                back = (a, b, c, v)

                front_faces = boundary_faces_3(front)
                back_faces = boundary_faces_3(back)

                na_front = [(f, s) for f, s in front_faces if f not in path2_set]
                na_back = [(f, s) for f, s in back_faces if f not in path2_set]

                print(f"\n  T' path ({a},{b},{c}):")
                print(f"    Front ({v},{a},{b},{c}) non-allowed faces: "
                      f"{[(f, '+' if s > 0 else '-') for f, s in na_front]}")
                print(f"    Back ({a},{b},{c},{v}) non-allowed faces: "
                      f"{[(f, '+' if s > 0 else '-') for f, s in na_back]}")

                # Check if same non-allowed faces appear in both
                na_f_set = {f for f, _ in na_front}
                na_b_set = {f for f, _ in na_back}
                common = na_f_set & na_b_set
                if common:
                    print(f"    SHARED non-allowed faces: {common}")
                    for cf in common:
                        sf = [s for f, s in front_faces if f == cf][0]
                        sb = [s for f, s in back_faces if f == cf][0]
                        print(f"      {cf}: front sign={'+' if sf > 0 else '-'}, "
                              f"back sign={'+' if sb > 0 else '-'}, "
                              f"{'CANCEL' if sf + sb == 0 else 'ADD'}")


# ============================================================
# PART 2: Systematic analysis - do non-allowed contributions always cancel?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: SYSTEMATIC - non-allowed face contributions")
print("=" * 70)

n = 5
na_contribution_types = set()

for A in all_tournaments_gen(n):
    paths2 = enumerate_allowed_paths(A, n, 2)
    path2_set = set(tuple(p) for p in paths2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path3_idx = {tuple(p): i for i, p in enumerate(paths3)}

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        for a, b, c in valid_Tp:
            front = (v, a, b, c)
            back = (a, b, c, v)

            # All faces of front and back
            f_faces = boundary_faces_3(front)
            b_faces = boundary_faces_3(back)

            # Combine: each face gets coeff alpha_{abc} * sign
            # For Omega_3, we need: for each non-allowed (x,y,z),
            # sum over all (a,b,c) of alpha_{abc} * total_sign = 0
            for face, sign in f_faces + b_faces:
                if face not in path2_set:
                    # This is a non-allowed face
                    # Classify: which positions of (v,a,b,c) / (a,b,c,v) is it from?
                    if face in [f for f, _ in f_faces]:
                        src = 'front'
                        face_idx = [i for i, (f, _) in enumerate(f_faces) if f == face][0]
                    else:
                        src = 'back'
                        face_idx = [i for i, (f, _) in enumerate(b_faces) if f == face][0]

                    # Does the OTHER cone also have this face?
                    other_faces = b_faces if src == 'front' else f_faces
                    partner = [(f, s) for f, s in other_faces if f == face]

                    if partner:
                        p_sign = partner[0][1]
                        na_contribution_types.add(
                            (src, face_idx, sign, 'partner', p_sign,
                             'cancel' if sign + p_sign == 0 else 'add'))
                    else:
                        na_contribution_types.add(
                            (src, face_idx, sign, 'no_partner', 0, 'exposed'))

print("\nNon-allowed face contribution patterns:")
for pattern in sorted(na_contribution_types):
    print(f"  {pattern}")


# ============================================================
# PART 3: Classify non-allowed faces by type
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: CLASSIFY NON-ALLOWED FACES")
print("=" * 70)
print("Faces of (v,a,b,c): (a,b,c) [T'], (v,b,c), (v,a,c), (v,a,b)")
print("Faces of (a,b,c,v): (b,c,v), (a,c,v), (a,b,v), (a,b,c) [T']")
print()
print("Non-allowed = not in A_2 = some consecutive pair missing arc")
print("For tournament: (x,y,z) non-allowed iff A[x][y]=0 or A[y][z]=0")
print()
print("v-faces of front: (v,b,c), (v,a,c), (v,a,b)")
print("  (v,x,y) non-allowed iff A[v][x]=0 or A[x][y]=0")
print("  Since v -> a (always), v -> b? v -> c?")
print()
print("v-faces of back: (b,c,v), (a,c,v), (a,b,v)")
print("  (x,y,v) non-allowed iff A[x][y]=0 or A[y][v]=0")
print("  Since c -> v (always), b -> v? a -> v?")

# At n=5, analyze which v-faces are non-allowed
n = 5
face_analysis = {}

for A in all_tournaments_gen(n):
    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_set = set(tuple(p) for p in paths2)

        for a, b, c in valid_Tp:
            # Classify each face
            front_v_faces = [(v, b, c), (v, a, c), (v, a, b)]
            back_v_faces = [(b, c, v), (a, c, v), (a, b, v)]

            for i, face in enumerate(front_v_faces):
                key = f"front_face_{i}"
                is_allowed = face in path2_set
                if key not in face_analysis:
                    face_analysis[key] = {'allowed': 0, 'not': 0}
                if is_allowed:
                    face_analysis[key]['allowed'] += 1
                else:
                    face_analysis[key]['not'] += 1

            for i, face in enumerate(back_v_faces):
                key = f"back_face_{i}"
                is_allowed = face in path2_set
                if key not in face_analysis:
                    face_analysis[key] = {'allowed': 0, 'not': 0}
                if is_allowed:
                    face_analysis[key]['allowed'] += 1
                else:
                    face_analysis[key]['not'] += 1

print(f"\nn=5 face allowedness:")
labels = {
    'front_face_0': '(v,b,c) sign=-1',
    'front_face_1': '(v,a,c) sign=+1',
    'front_face_2': '(v,a,b) sign=-1',
    'back_face_0': '(b,c,v) sign=+1',
    'back_face_1': '(a,c,v) sign=-1',
    'back_face_2': '(a,b,v) sign=+1',
}
for key in sorted(face_analysis.keys()):
    d = face_analysis[key]
    total = d['allowed'] + d['not']
    pct = 100 * d['not'] / total if total > 0 else 0
    print(f"  {labels.get(key, key)}: allowed={d['allowed']}, "
          f"non-allowed={d['not']} ({pct:.1f}%)")


# ============================================================
# PART 4: When a face IS non-allowed, does the partner cancel it?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: NON-ALLOWED FACE CANCELLATION ANALYSIS")
print("=" * 70)

n = 5
cancel_count = 0
expose_count = 0
na_face_details = []

for A in all_tournaments_gen(n):
    paths2 = enumerate_allowed_paths(A, n, 2)
    path2_set = set(tuple(p) for p in paths2)

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]

        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_all = [tuple(vlist[x] for x in p) for p in paths2_Tp]
        valid_Tp = [(a, b, c) for (a, b, c) in Tp_all
                    if A[v][a] == 1 and A[c][v] == 1]

        for a, b, c in valid_Tp:
            front = (v, a, b, c)
            back = (a, b, c, v)
            f_faces = boundary_faces_3(front)
            b_faces = boundary_faces_3(back)

            # Check each non-allowed face from front
            for face, sign in f_faces:
                if face in path2_set or face == (a, b, c):
                    continue  # allowed or T' internal
                # Does back also have this face?
                partner = [(f, s) for f, s in b_faces if f == face]
                if partner and sign + partner[0][1] == 0:
                    cancel_count += 1
                else:
                    expose_count += 1
                    na_face_details.append(('front', face, sign, v, a, b, c))

            # Check each non-allowed face from back
            for face, sign in b_faces:
                if face in path2_set or face == (a, b, c):
                    continue
                partner = [(f, s) for f, s in f_faces if f == face]
                if partner and sign + partner[0][1] == 0:
                    pass  # already counted
                elif not partner:
                    expose_count += 1
                    na_face_details.append(('back', face, sign, v, a, b, c))

print(f"n=5: Non-allowed face pairs that cancel: {cancel_count}")
print(f"n=5: Exposed non-allowed faces: {expose_count}")

if na_face_details:
    print("\nExposed non-allowed face examples:")
    for src, face, sign, v, a, b, c in na_face_details[:10]:
        print(f"  {src} ({v},{a},{b},{c}): face={face}, sign={'+' if sign > 0 else '-'}")
else:
    print("\nNO exposed non-allowed faces! All cancel in pairs.")
    print("This means the front+back cone is INTRINSICALLY in Omega_3.")
    print("Proof sketch: every non-allowed face of (v,a,b,c) appears")
    print("in (a,b,c,v) with opposite sign, and vice versa.")


print("\n\nDone.")
