#!/usr/bin/env python3
"""
Debug: Check whether flip (bit complement) is well-defined on isomorphism classes.
Instance: opus-2026-03-06-S11
"""
from itertools import permutations
from collections import defaultdict

def analyze_flip_consistency(n):
    """Check if all members of a class flip to the same class."""
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    flip_mask = (1 << m) - 1

    # Build adjacency from bits
    def bits_to_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1):
            A[k][k+1] = 1
        for i, (xL, yL) in enumerate(tiles):
            xi = verts.index(xL)
            yi = verts.index(yL)
            if (bits >> i) & 1 == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return tuple(tuple(row) for row in A)

    def canonicalize(A):
        best = None
        for p in permutations(range(n)):
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    # Build classes
    classes = {}
    class_members = defaultdict(list)
    tiling_class = {}
    cid_counter = 0

    for bits in range(1 << m):
        A = bits_to_adj(bits)
        canon = canonicalize(A)
        if canon not in classes:
            classes[canon] = cid_counter
            cid_counter += 1
        cid = classes[canon]
        class_members[cid].append(bits)
        tiling_class[bits] = cid

    num_classes = cid_counter
    print(f"\nn={n}: {num_classes} classes, {1 << m} tilings, {m} tiles")

    # For each class, check where all members flip to
    inconsistent = 0
    class_flip_targets = {}
    for cid in range(num_classes):
        targets = set()
        for bits in class_members[cid]:
            fbits = bits ^ flip_mask
            targets.add(tiling_class[fbits])
        class_flip_targets[cid] = targets
        if len(targets) > 1:
            inconsistent += 1
            print(f"  Class {cid} ({len(class_members[cid])} members) flips to MULTIPLE classes: {targets}")

    if inconsistent == 0:
        print(f"  Flip is well-defined on all classes!")
    else:
        print(f"  {inconsistent} classes have inconsistent flip targets!")

    # Now do the proper categorization
    class_flip = {}
    for cid in range(num_classes):
        if len(class_flip_targets[cid]) == 1:
            class_flip[cid] = next(iter(class_flip_targets[cid]))

    # Check involution
    involution_ok = True
    for cid in class_flip:
        fcid = class_flip[cid]
        if fcid in class_flip and class_flip[fcid] != cid:
            print(f"  NOT involution: flip({cid})={fcid} but flip({fcid})={class_flip[fcid]}")
            involution_ok = False

    if involution_ok and inconsistent == 0:
        # Count self-flip and pairs
        self_flip = [cid for cid in range(num_classes) if class_flip[cid] == cid]
        seen = set()
        pairs = []
        for cid in range(num_classes):
            if cid in seen:
                continue
            fcid = class_flip[cid]
            if fcid == cid:
                seen.add(cid)
            else:
                pairs.append((cid, fcid))
                seen.add(cid)
                seen.add(fcid)
        print(f"  Self-flip: {len(self_flip)}, Pairs: {len(pairs)}")
        print(f"  Check: {len(self_flip)} + 2*{len(pairs)} = {len(self_flip) + 2*len(pairs)} vs {num_classes}")

for n in range(3, 7):
    analyze_flip_consistency(n)
