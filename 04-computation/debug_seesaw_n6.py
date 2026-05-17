"""Debug the n=6 mismatch in the formula H3(T) = H2tilde(I(H_v))."""
import sys, itertools
from collections import defaultdict
import numpy as np

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from beta1_beta3_seesaw_experiments import *

def rank_snf(M):
    """Better rank computation via SVD with higher tolerance."""
    if M.size == 0:
        return 0
    s = np.linalg.svd(M.astype(float), compute_uv=False)
    tol = max(M.shape) * np.finfo(float).eps * (s[0] if s.size > 0 else 1)
    return int(np.sum(s > tol))

n = 6
mismatches = []

for arc_set in all_tournaments(n):
    betti = compute_tournament_betti(arc_set, n, max_p=4)
    b3 = betti.get(3, 0)

    for v in range(n):
        h2 = compute_H2tilde_I_Hv(arc_set, n, v)
        if h2 != b3:
            mismatches.append((arc_set, v, b3, h2))
            if len(mismatches) <= 3:
                print(f"Mismatch: v={v}, beta3={b3}, H2tilde={h2}")
                left, right, edges = compute_H_v(arc_set, n, v)
                print(f"  H_v: left={left}, right={right}, edges={sorted(edges)}")
                faces = build_independence_complex(edges, left+right)
                hom = compute_homology_of_complex(faces)
                print(f"  I(H_v) homology: {hom}")
                print(f"  Number of faces: {len(faces)}")
                # Also compute betti by dimension
                by_dim = defaultdict(int)
                for f in faces:
                    by_dim[len(f)-1] += 1
                print(f"  f-vector: {dict(sorted(by_dim.items()))}")
                print()

print(f"Total mismatches: {len(mismatches)}")

# Sample a mismatch to examine
if mismatches:
    arc_set, v, b3, h2 = mismatches[0]
    print(f"\nDetailed examination of first mismatch:")
    print(f"b3={b3}, H2tilde={h2}")

    # Check if b3 is computed correctly
    betti2 = compute_tournament_betti(arc_set, n, max_p=5)
    print(f"Full betti: {betti2}")

    # Check H_v more carefully
    left, right, edges = compute_H_v(arc_set, n, v)
    print(f"H_v: left={left}, right={right}")
    print(f"Edges: {sorted([sorted(e) for e in edges])}")

    # Build faces manually
    faces = build_independence_complex(edges, left + right)
    print(f"Independent sets: {sorted([sorted(f) for f in faces])}")

    # Compute homology step by step
    hom = compute_homology_of_complex(faces)
    print(f"Homology: {hom}")
