"""
hessian_eigenvalues.py — opus-2026-04-04-S6

THE KEY CLAIM: The quadratic form Q(t) = Σ_{i<j} c_{ij} t_i t_j
has EXACTLY n-2 negative eigenvalues.

The Hessian matrix H_ij = c_{ij} (the quadratic multilinear coefficients)
determines the curvature of H(t) at the transitive point t=0.

The n-2 negative eigenvalues come from the TWO RECURSION MODES:
  Mode A: add hypotenuse strip (n → n+1)
  Mode B: add both legs (n → n+2)

And: Mode A ∘ Mode A ≅ Mode B in some isomorphic sense.

Let's verify the eigenvalue count and understand the structure.
"""
import numpy as np
from collections import defaultdict

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def compute_quadratic_form(n):
    """Compute the quadratic form matrix Q_{ij} = c_{ij} (multilinear coefficients)."""
    tiles = make_tiles(n)
    m = len(tiles)

    # Compute H for all tilings
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)

    # Extract quadratic coefficients via Möbius
    Q = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(i, m):
            S = (1 << i) | (1 << j)
            if i == j:
                # Diagonal: this is the coefficient of t_i^2 = 0 in multilinear
                # But for the Hessian we use ∂²H/∂t_i∂t_j
                # For i=j in multilinear: there's no t_i² term.
                # The "Hessian" at t=1/2 would include the linear terms.
                # Let me use the proper quadratic form: the off-diagonal entries only.
                continue
            else:
                # c_{ij} via Möbius: c_{ij} = H(ij) - H(i) - H(j) + H(0)
                h_ij = int(H[(1 << i) | (1 << j)])
                h_i = int(H[1 << i])
                h_j = int(H[1 << j])
                h_0 = int(H[0])
                c_ij = h_ij - h_i - h_j + h_0
                Q[i][j] = c_ij
                Q[j][i] = c_ij

    return Q, tiles, m, H

def analyze_eigenvalues():
    print("HESSIAN EIGENVALUE ANALYSIS")
    print("=" * 60)

    for n in range(3, 9):
        tiles = make_tiles(n)
        m = len(tiles)

        if m > 21:
            print(f"\nn={n}: m={m} too large, skipping")
            continue

        Q, tiles_list, m, H = compute_quadratic_form(n)

        # Eigenvalues
        eigvals = np.linalg.eigvalsh(Q)
        eigvals_sorted = sorted(eigvals)

        n_neg = sum(1 for e in eigvals if e < -1e-10)
        n_zero = sum(1 for e in eigvals if abs(e) < 1e-10)
        n_pos = sum(1 for e in eigvals if e > 1e-10)

        print(f"\nn={n}: m={m} tiles")
        print(f"  Eigenvalues of Q (sorted): {[f'{e:.2f}' for e in eigvals_sorted]}")
        print(f"  Signature: ({n_pos}+, {n_zero}0, {n_neg}-)")
        print(f"  Expected n-2 negative: {n-2}")
        print(f"  Actual negative: {n_neg}")
        print(f"  {'✓ MATCH' if n_neg == n-2 else '✗ MISMATCH'}")

        # Trace and determinant
        print(f"  Trace(Q) = {np.trace(Q):.2f} (sum of eigenvalues)")
        print(f"  Sum of c_{'{ij}'} = {np.sum(Q)/2:.0f}")

        # The positive eigenvalues: how many?
        # m - (n-2) - n_zero positive eigenvalues
        print(f"  Positive eigenvalues: {n_pos}")
        print(f"  m - (n-2) = {m - (n-2)}")

        # MODE A vs MODE B:
        # The staircase δ_{n-2} has tiles organized by "strip" (anti-diagonal r+c=k).
        # Mode A (hypotenuse) adds the outermost strip.
        # Mode B (legs) adds the bottom column and right row.

        # Let me identify which eigenvectors are negative.
        eigvals_full, eigvecs = np.linalg.eigh(Q)
        neg_idx = [i for i in range(m) if eigvals_full[i] < -1e-10]

        if n <= 7 and neg_idx:
            print(f"\n  Negative eigenvectors (showing dominant tiles):")
            for idx in neg_idx:
                vec = eigvecs[:, idx]
                # Find the tiles with largest |coefficient|
                top = sorted(range(m), key=lambda k: abs(vec[k]), reverse=True)[:4]
                print(f"    λ={eigvals_full[idx]:.2f}: "
                      f"{[(tiles_list[k], f'{vec[k]:.3f}') for k in top]}")

        # CLASSIFY tiles by position on staircase
        if n <= 7:
            # Strip (anti-diagonal): r+c where r=x-y-1, c=y (1-indexed)
            # So strip = (x-y-1)+y = x-1
            # Hmm, that gives strip = x-1 which is just the upper vertex index.
            # Let me use the standard strip: r+c = a-b-1+b = a-1 where a=x, b=y (1-indexed).
            # So strip(tile (x,y)) = x-1.

            # The hypotenuse tiles have strip = n-1 (the outermost anti-diagonal).
            # These are tiles (n, y) for y=1,...,n-2. Count: n-2.

            # The leg tiles:
            # Bottom leg: tiles (x, 1) for x=3,...,n. Strip values: 2,...,n-1.
            # Top leg: tiles (n, y) for y=1,...,n-2. Strip values: n-1,...,n-1.
            # Wait, tiles (n, y) have strip = n-1 for all y. So ALL top-leg tiles
            # are on the HYPOTENUSE strip!

            # This means:
            # Hypotenuse strip = {(n, y) : y=1,...,n-2} = top-leg tiles
            # This is the same as Mode A's addition!

            hyp_tiles = [i for i, (x,y) in enumerate(tiles_list) if x == n]
            bottom_tiles = [i for i, (x,y) in enumerate(tiles_list) if y == 1]
            inner_tiles = [i for i, (x,y) in enumerate(tiles_list)
                          if x < n and y > 1]

            print(f"\n  Tile classification:")
            print(f"    Hypotenuse (x=n={n}): {len(hyp_tiles)} tiles "
                  f"{[tiles_list[i] for i in hyp_tiles]}")
            print(f"    Bottom leg (y=1): {len(bottom_tiles)} tiles "
                  f"{[tiles_list[i] for i in bottom_tiles]}")
            print(f"    Inner (x<n, y>1): {len(inner_tiles)} tiles")

            # The OVERLAP between hypotenuse and bottom: tile (n, 1) = apex
            apex = [i for i in hyp_tiles if i in bottom_tiles]
            print(f"    Apex (overlap): {len(apex)} tiles")

            # PROJECT negative eigenvectors onto tile groups
            for idx in neg_idx:
                vec = eigvecs[:, idx]
                hyp_weight = sum(vec[i]**2 for i in hyp_tiles)
                bot_weight = sum(vec[i]**2 for i in bottom_tiles)
                inn_weight = sum(vec[i]**2 for i in inner_tiles)
                print(f"    λ={eigvals_full[idx]:.2f}: "
                      f"hyp={hyp_weight:.3f} bot={bot_weight:.3f} "
                      f"inner={inn_weight:.3f}")

def mode_a_vs_mode_b():
    """
    Mode A: n → n+1 by adding hypotenuse strip.
    This means: add tiles (n+1, y) for y=1,...,n-1. These are n-1 new tiles.
    The new staircase goes from δ_{n-2} (m = C(n-1,2)) to δ_{n-1} (m' = C(n,2) = m + n-1).

    Mode B: n → n+2 by adding both legs.
    Add bottom (x, 1) for x=3,...,n+2 and top (n+2, y) for y=2,...,n plus apex.
    Total: 2n-1 new tiles. Goes from δ_{n-2} to δ_n (m'' = C(n+1,2) = m + 2n-1).

    Mode A twice: n → n+1 → n+2.
    First add n-1 tiles (hypotenuse of n+1), then n tiles (hypotenuse of n+2).
    Total: 2n-1 new tiles. Same count as Mode B!

    So Mode A ∘ Mode A adds the SAME NUMBER of tiles as Mode B.
    But the tiles are DIFFERENT:
    - Mode A²: tiles (n+1, y) for y=1..n-1 UNION tiles (n+2, y) for y=1..n
    - Mode B: tiles (x, 1) for x=3..n+2 UNION tiles (n+2, y) for y=2..n UNION apex (n+2,1)

    Wait — Mode A² adds tiles involving vertices n+1 and n+2 from the TOP.
    Mode B adds tiles involving vertices 1 and n+2 from the LEGS.

    These are DIFFERENT sets of tiles! They're not the same staircase decomposition.

    But the resulting staircase δ_n is the same in both cases. The difference is in
    which INTERMEDIATE STRUCTURE is created:
    - Mode A²: creates δ_{n-1} as intermediate
    - Mode B: goes directly to δ_n

    The "isomorphism" might be about the TRANSITION MATRICES, not the tile identities.
    """
    print("\n" + "=" * 60)
    print("MODE A vs MODE B: TWO PATHS TO δ_n")
    print("=" * 60)

    for n in [5, 6, 7]:
        print(f"\nn={n}: δ_{n-2} → δ_n")

        tiles_n = make_tiles(n)
        tiles_np1 = make_tiles(n+1) if n+1 <= 8 else None
        tiles_np2 = make_tiles(n+2)

        m_n = len(tiles_n)
        m_np2 = len(tiles_np2)

        # Mode A step 1: n → n+1 adds hypotenuse strip
        # Tiles of (n+1)-tournament not in n-tournament:
        # (n+1)-tournament tiles (x,y) with x-y≥2, x,y ∈ {1,...,n+1}
        # n-tournament tiles (x,y) with x-y≥2, x,y ∈ {1,...,n}
        # New tiles: those involving vertex n+1.
        # = (n+1, y) for y=1,...,n-1 → n-1 tiles

        hyp1 = [(n+1, y) for y in range(1, n)]
        print(f"  Mode A step 1 (hypotenuse): {len(hyp1)} tiles: {hyp1[:5]}...")

        # Mode A step 2: n+1 → n+2 adds another hypotenuse strip
        # New tiles: (n+2, y) for y=1,...,n → n tiles
        hyp2 = [(n+2, y) for y in range(1, n+1)]
        print(f"  Mode A step 2 (hypotenuse): {len(hyp2)} tiles: {hyp2[:5]}...")

        mode_a_total = set(hyp1) | set(hyp2)
        print(f"  Mode A² total new tiles: {len(mode_a_total)}")

        # Mode B: n → n+2 directly (add legs)
        # New tiles: those in (n+2)-tournament but not in n-tournament
        tiles_n_set = set(tiles_n)
        # But wait: n-tournament tiles have vertices 1..n
        # (n+2)-tournament tiles have vertices 1..n+2
        # "Same" inner tiles: (x+1, y+1) for (x,y) in n-tournament
        # are tiles of n+2-tournament connecting {2,...,n+1}

        # The n+2 tournament has tiles connecting {1,...,n+2}.
        # The "overlap" with n-tournament (shifted to {2,...,n+1}) is
        # tiles (x+1, y+1) for (x,y) ∈ tiles_n.

        overlap = set((x+1, y+1) for x, y in tiles_n)
        mode_b_new = set(tiles_np2) - overlap
        print(f"  Mode B new tiles (legs): {len(mode_b_new)}")

        # KEY: what's the DIFFERENCE between Mode A² tiles and Mode B tiles?
        # Mode A² adds tiles involving vertices n+1 and n+2 (from the top of the staircase)
        # Mode B adds tiles involving vertices 1 and n+2 (from the legs)

        # But both get to the SAME staircase δ_n.
        # The isomorphism must be a RELABELING of vertices!

        # Mode A² keeps vertices {1,...,n} fixed, adds n+1 then n+2 at the top.
        # Mode B adds vertex 1 (bottom) and vertex n+2 (top).

        # The difference: Mode A² adds TWO top vertices.
        # Mode B adds one bottom + one top vertex.

        # Relabeling: shift all old vertices up by 1 (so {1,...,n} → {2,...,n+1})
        # then add vertex 1 (bottom) and vertex n+2 (top).
        # This IS Mode B!

        # While Mode A² adds vertex n+1 and vertex n+2 to the TOP.

        # The "isomorphism" between A² and B is the vertex relabeling:
        # A² uses {1,...,n,n+1,n+2} with old={1,...,n} and new={n+1,n+2}
        # B uses {1,2,...,n+1,n+2} with old={2,...,n+1} and new={1,n+2}

        # Under the shift i → i+1:
        # A²'s old vertices {1,...,n} map to B's inner vertices {2,...,n+1} ✓
        # A²'s new vertex n+1 maps to... not 1 or n+2.

        # So it's NOT a simple shift. The isomorphism is more subtle.
        print(f"\n  Mode A² tiles (involving n+1={n+1} or n+2={n+2}):")
        for t in sorted(mode_a_total):
            print(f"    {t}")
        print(f"  Mode B tiles (legs + apex):")
        for t in sorted(mode_b_new):
            print(f"    {t}")

if __name__ == "__main__":
    analyze_eigenvalues()
    mode_a_vs_mode_b()
