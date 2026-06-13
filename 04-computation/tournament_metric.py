"""
tournament_metric.py — opus-2026-04-04-S6

PRACTICAL APPLICATION: Use the quadratic form Q to define a
"tournament distance metric" and "tournament information content".

The quadratic form Q encodes pairwise tile interactions.
Two tournaments that differ in their "frustrated" (same-end) tiles
will have large quadratic distance, while those differing in
"synergistic" (disjoint) tiles will be closer.

Applications:
1. Tournament fingerprinting — classify tournaments by quadratic projection
2. Tournament information content — how many bits does H encode?
3. Tournament similarity search — find similar tournaments quickly
"""
import numpy as np
from collections import defaultdict
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def compute_Q_metric(n):
    """Build Q and use it as a metric on tiling space."""
    print(f"\n{'='*60}")
    print(f"TOURNAMENT Q-METRIC n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m
    coeffs = mobius_inversion(H, m)

    # Build Q
    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        Q[bits[0]][bits[1]] = c
        Q[bits[1]][bits[0]] = c

    # Also add linear term on diagonal (not part of Q, but useful for metric)
    # c_k = 2^{skip-1}
    linear = np.zeros(m)
    for k in range(m):
        skip = tiles[k][0] - tiles[k][1]
        linear[k] = 2**(skip-1)

    # Q-distance between tilings t1, t2:
    # d_Q(t1,t2) = (x1-x2)^T Q (x1-x2) where xi = tiling vector
    # This is a signed quadratic form — not a metric in the usual sense.
    # But |d_Q| measures "quadratic interaction difference".

    # Compute Q-distance between all pairs and correlate with |H(t1)-H(t2)|
    # (sampling for n=7)
    if N <= 1024:
        # All pairs
        dQ_vals = []
        dH_vals = []
        for t1 in range(0, N, max(1, N//200)):
            for t2 in range(t1+1, min(t1+50, N)):
                x1 = np.array([(t1>>k)&1 for k in range(m)], dtype=float)
                x2 = np.array([(t2>>k)&1 for k in range(m)], dtype=float)
                diff = x1 - x2
                dQ = diff @ Q @ diff
                dH = abs(int(H[t1]) - int(H[t2]))
                dQ_vals.append(dQ)
                dH_vals.append(dH)

        corr = np.corrcoef(dQ_vals, dH_vals)[0,1]
        print(f"  Corr(Q-distance, |ΔH|): {corr:.4f}")
    else:
        np.random.seed(42)
        dQ_vals = []
        dH_vals = []
        for _ in range(5000):
            t1 = np.random.randint(N)
            t2 = np.random.randint(N)
            x1 = np.array([(t1>>k)&1 for k in range(m)], dtype=float)
            x2 = np.array([(t2>>k)&1 for k in range(m)], dtype=float)
            diff = x1 - x2
            dQ = diff @ Q @ diff
            dH = abs(int(H[t1]) - int(H[t2]))
            dQ_vals.append(dQ)
            dH_vals.append(dH)
        corr = np.corrcoef(dQ_vals, dH_vals)[0,1]
        print(f"  Corr(Q-distance, |ΔH|) [sampled]: {corr:.4f}")

    # INFORMATION CONTENT: how much of H's variance is explained by
    # the linear and quadratic terms alone?
    h_pred_lin = np.zeros(N)
    h_pred_quad = np.zeros(N)
    for t in range(N):
        x = np.array([(t>>k)&1 for k in range(m)], dtype=float)
        h_pred_lin[t] = 1 + linear @ x
        h_pred_quad[t] = 1 + linear @ x + x @ Q @ x / 2  # include cross-terms

    var_H = np.var(H.astype(float))
    var_lin_resid = np.var(H.astype(float) - h_pred_lin)
    var_quad_resid = np.var(H.astype(float) - h_pred_quad)

    r2_lin = 1 - var_lin_resid / var_H
    r2_quad = 1 - var_quad_resid / var_H

    print(f"\n  H variance explained:")
    print(f"    Linear only (degree 1): R² = {r2_lin:.4f}")
    print(f"    Linear + Quadratic (degree 1+2): R² = {r2_quad:.4f}")
    print(f"    Remaining variance fraction: {var_quad_resid/var_H:.4f}")

    # How many "effective bits" does the quadratic form encode?
    # Information = log2(1/R²_residual) ≈ bits of H explained by quadratic
    if r2_quad > 0 and r2_quad < 1:
        bits_explained = -np.log2(1 - r2_quad)
        print(f"    Bits of H explained by quadratic: {bits_explained:.2f}")

    # TOURNAMENT FINGERPRINT: project into top-3 Q eigenvectors
    eigvals, eigvecs = np.linalg.eigh(Q)
    order = sorted(range(m), key=lambda i: abs(eigvals[i]), reverse=True)

    print(f"\n  Top-3 Q eigenvalues: {[f'{eigvals[order[k]]:.2f}' for k in range(min(3,m))]}")

    # Cluster H values in fingerprint space
    if N <= 32768:
        from collections import Counter
        h_vals = [int(H[t]) for t in range(N)]
        h_counter = Counter(h_vals)
        h_unique = sorted(h_counter.keys())

        # For each H value, compute mean fingerprint
        top_k = min(3, m)
        print(f"\n  Mean fingerprint by H value (top-{top_k} projections):")
        for h in h_unique[:10]:
            tilings = [t for t in range(N) if int(H[t]) == h]
            fps = []
            for t in tilings:
                x = np.array([(t>>k)&1 for k in range(m)], dtype=float) - 0.5
                fp = [eigvecs[:, order[k]].dot(x) for k in range(top_k)]
                fps.append(fp)
            mean_fp = np.mean(fps, axis=0)
            std_fp = np.std(fps, axis=0)
            print(f"    H={h:3d} ({len(tilings):5d} tilings): "
                  f"mean=[{', '.join(f'{m:.3f}' for m in mean_fp)}], "
                  f"std=[{', '.join(f'{s:.3f}' for s in std_fp)}]")

if __name__ == "__main__":
    for n in [5, 6, 7]:
        compute_Q_metric(n)
