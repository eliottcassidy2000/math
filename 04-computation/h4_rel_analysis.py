"""
h4_rel_analysis.py — H_4(T,T\v) and the connecting homomorphism δ

From the LES: H_4(T) → H_4(T,T\v) →δ H_3(T\v) →i_* H_3(T)
- If β_4(T)=0: ker(δ)=0, so im(δ)=H_4^rel, rank(i_*)=1-dim(H_4^rel)
- rank(i_*)=0 iff H_4^rel≥1

Key question: why is H_4^rel = 0 universally at n=7 (b3=1)?

This script computes:
1. β_4(T) for b3=1 tournaments
2. dim(H_4^rel) directly from relative complex
3. How ψ(ker) compares to im(d_4^{T\v}) dimension-by-dimension

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, deletion_adj,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def analyze_h4_rel(A, n, v):
    """Compute H_4(T,T\\v) dimension and related invariants."""
    max_p = min(n - 1, 6)
    cc = full_chain_complex_modp(A, n, max_p)
    b4 = cc['bettis'].get(4, 0)
    b3 = cc['bettis'].get(3, 0)

    if b3 != 1:
        return None

    # Get T\v data
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}

    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))
    b3_Tv = cc_Tv['bettis'].get(3, 0)
    b4_Tv = cc_Tv['bettis'].get(4, 0)

    if b3_Tv != 1:
        return None

    # Now compute chain-level data for H_4^rel
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    def classify(ps, v):
        tv = [i for i, p in enumerate(ps) if v in p]
        old = [i for i, p in enumerate(ps) if v not in p]
        return tv, old

    tv4, old4 = classify(paths.get(4, []), v)
    tv3, old3 = classify(paths.get(3, []), v)
    tv5, old5 = classify(paths.get(5, []), v) if paths.get(5) else ([], [])

    if not tv4:
        return None

    # Build boundary maps
    idx3_T = {p: i for i, p in enumerate(paths.get(3, []))}
    idx4_T = {p: i for i, p in enumerate(paths.get(4, []))}

    bd4_T = np.zeros((len(paths.get(3, [])), len(paths.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    bd5_T = np.zeros((len(paths.get(4, [])), len(paths.get(5, []))), dtype=np.int64)
    for j, path in enumerate(paths.get(5, [])):
        for sign, face in boundary_faces(path):
            if face in idx4_T:
                bd5_T[idx4_T[face], j] = (bd5_T[idx4_T[face], j] + sign) % PRIME

    # Get Omega bases
    def get_omega(ap, deg):
        ps = ap.get(deg, [])
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    ob4_T = get_omega(ap_T, 4)
    ob5_T = get_omega(ap_T, 5)

    # d_4^{tv→tv}: the relative boundary
    d4_tvtv = bd4_T[np.ix_(tv3, tv4)] % PRIME if tv3 and tv4 else np.zeros((0, 0), dtype=np.int64)

    # Map d_4^{tv→tv} through Omega_4 basis
    d4_tvtv_omega = (d4_tvtv @ ob4_T[:, tv4].T % PRIME) if ob4_T.shape[0] > 0 and tv4 else np.zeros((len(tv3), 0), dtype=np.int64)

    # ker(d_4^{tv→tv} restricted to Omega_4)
    if d4_tvtv_omega.size > 0 and d4_tvtv_omega.shape[1] > 0:
        _, kv = _gauss_nullbasis_modp(
            d4_tvtv_omega.astype(np.int32), d4_tvtv_omega.shape[0], d4_tvtv_omega.shape[1], PRIME)
        dim_ker_rel4 = len(kv) if kv else 0
    else:
        dim_ker_rel4 = ob4_T.shape[0]

    # im(d_5^{tv→tv} restricted to Omega_5) — this is im(d_5^rel)
    d5_tvtv = bd5_T[np.ix_(tv4, tv5)] % PRIME if tv4 and tv5 else np.zeros((len(tv4), 0), dtype=np.int64)
    d5_tvtv_omega = (d5_tvtv @ ob5_T[:, tv5].T % PRIME) if ob5_T.shape[0] > 0 and tv5 else np.zeros((len(tv4), 0), dtype=np.int64)

    # But we need im(d_5^rel) inside ker(d_4^rel). The image of d_5^rel
    # lands in Omega_4 space. Its tv4-components are what enters d_4^{tv→tv}.
    # Actually d_5^rel maps from C_5^rel to C_4^rel, where C_p^rel ≅ Ω_p^{tv-supported}.
    # Hmm, this is the quotient complex. Let me think...

    # Actually, for the relative complex approach:
    # C_p^rel = C_p(T) / C_p(T\v)
    # The relative d_4 is induced by d_4 on the quotient.
    # If we represent C_4^rel by the tv4 components (since old4 ≅ C_4(T\v)):
    #   d_4^rel([w]) = [d_4(w)] where w has tv4 components only
    #   [d_4(w)] = tv3-components of d_4(w) (mod old3 ≅ C_3(T\v))
    #   = d_4^{tv→tv}(w_tv4)

    # So d_4^rel = d_4^{tv→tv} on C_4^{tv} / C_4^{tv} ∩ C_4(T\v) = C_4^{tv}.
    # And d_5^rel = d_5^{tv→tv} on C_5^{tv}.

    # But we need to work in Omega coords. In Omega_4(T), the tv4-part of
    # an Omega basis vector may be partial. The correct relative Omega_4 is:
    # Ω_4^rel = Ω_4(T) / i(Ω_4(T\v))

    # Let me just compute H_4^rel numerically.
    # H_4^rel = ker(d_4^rel) / im(d_5^rel) in the relative chain complex.

    # The relative chain complex has:
    # C_p^rel = A_p(T) / i(A_p(T\v)) at the allowed-path level.
    # But we should work at the Omega level for the actual chain complex.

    # Alternative: compute directly from ranks.
    # From LES: ... → H_4(T) → H_4(T,T\v) → H_3(T\v) → H_3(T) → ...
    # So dim(H_4^rel) = dim(ker i_*:H_3(T\v)→H_3(T)) + dim(coker H_4(T)→H_4^rel)
    # With b4=0: dim(H_4^rel) = dim(ker i_*)
    # And dim(ker i_*) = b3_Tv - rank(i_*) = 1 - rank(i_*)

    # So H_4^rel = 0 iff rank(i_*) = 1, H_4^rel = 1 iff rank(i_*) = 0.
    # This is circular... I need to compute i_* directly.

    # Compute rank(i_*) directly
    # Embed H_3(T\v) generator into Ω_3(T) and check if it's in im(d_4^T)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)
    ob3_T = get_omega(ap_T, 3)

    # d_3^{T\v} in Omega coords
    idx3_Tv = {p: i for i, p in enumerate(paths_Tv.get(3, []))}
    idx2_Tv = {p: i for i, p in enumerate(paths_Tv.get(2, []))}

    bd3_Tv = np.zeros((len(paths_Tv.get(2, [])), len(paths_Tv.get(3, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(3, [])):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3_Tv = np.array(kv_Tv, dtype=np.int64) if kv_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3_Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3_Tv_A = ker_d3_Tv @ ob3_Tv % PRIME  # in A_3(T\v) coords

    # im(d_4^{T\v}) in A_3(T\v) coords
    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    rk_im_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

    # Full boundary image of d_4^T in Omega coords (restricted to old3 rows)
    # This is ψ(ker(d_4^{tv→tv})) in the correct sense
    if dim_ker_rel4 > 0:
        ker_rel4 = np.array(kv if kv else [], dtype=np.int64).reshape(-1, ob4_T.shape[0]) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
        ker_A4 = ker_rel4 @ ob4_T % PRIME  # in A_4(T) coords
        # Full old3-projection of boundary
        psi_full = (bd4_T[old3, :] @ ker_A4.T % PRIME).T
    else:
        psi_full = np.zeros((0, len(old3)), dtype=np.int64)

    # Map psi_full to T\v coords
    old3_to_Tv = {}
    for local_j, global_i in enumerate(old3):
        path_T = paths[3][global_i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_to_Tv[local_j] = idx3_Tv[path_Tv]

    psi_Tv = np.zeros((psi_full.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
    for local_j in range(len(old3)):
        if local_j in old3_to_Tv:
            j_Tv = old3_to_Tv[local_j]
            psi_Tv[:, j_Tv] = (psi_Tv[:, j_Tv] + psi_full[:, local_j]) % PRIME

    rk_psi = int(_gauss_rank_np(psi_Tv.copy(), PRIME)) if psi_Tv.shape[0] > 0 else 0

    # Codimension of ψ(ker) in ker(d_3^{T\v}) and in im(d_4^{T\v})
    if psi_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        psi_in_im = (rk_combined == rk_im_d4_Tv)
        extra = rk_combined - rk_im_d4_Tv
    elif psi_Tv.shape[0] > 0:
        psi_in_im = (rk_psi == 0)
        extra = rk_psi
    else:
        psi_in_im = True
        extra = 0

    # Also compute: dim(Ω_4^{tv only}) and dim(Ω_5^{tv only})
    dim_omega4 = ob4_T.shape[0]
    dim_omega5 = ob5_T.shape[0]

    return {
        'b3': b3, 'b4': b4, 'b3_Tv': b3_Tv, 'b4_Tv': b4_Tv,
        'dim_omega4': dim_omega4, 'dim_omega5': dim_omega5,
        'dim_tv4': len(tv4), 'dim_tv5': len(tv5),
        'dim_old4': len(old4), 'dim_old5': len(old5),
        'dim_ker_rel4': dim_ker_rel4,
        'rk_psi': rk_psi,
        'rk_im_d4_Tv': rk_im_d4_Tv,
        'rk_ker_d3_Tv': ker_d3_Tv_A.shape[0],
        'psi_in_im': psi_in_im,
        'extra': extra,
        'codim_psi_in_imd4': rk_im_d4_Tv - rk_psi if psi_in_im else None,
        'codim_psi_in_kerd3': ker_d3_Tv_A.shape[0] - rk_psi,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"H_4^rel AND δ ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 400

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)

            for v_cand in range(n):
                r = analyze_h4_rel(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs with b3=b3_Tv=1, {elapsed:.1f}s")

        # β_4 distribution
        b4_dist = Counter(r['b4'] for r in results)
        print(f"\n  β_4(T) distribution: {dict(sorted(b4_dist.items()))}")

        b4_Tv_dist = Counter(r['b4_Tv'] for r in results)
        print(f"  β_4(T\\v) distribution: {dict(sorted(b4_Tv_dist.items()))}")

        # ψ containment
        in_im = sum(1 for r in results if r['psi_in_im'])
        print(f"\n  ψ(ker) ⊂ im(d_4^Tv): {in_im}/{len(results)}")

        # Codimension of ψ(ker) in im(d_4^{T\v})
        codims = [r['codim_psi_in_imd4'] for r in results if r['codim_psi_in_imd4'] is not None]
        if codims:
            codim_dist = Counter(codims)
            print(f"  codim(ψ(ker), im d_4^Tv): {dict(sorted(codim_dist.items()))}")

        # Codimension of ψ(ker) in ker(d_3^{T\v})
        codims2 = [r['codim_psi_in_kerd3'] for r in results]
        codim2_dist = Counter(codims2)
        print(f"  codim(ψ(ker), ker d_3^Tv): {dict(sorted(codim2_dist.items()))}")

        # Dimensional stats
        for key in ['dim_omega4', 'dim_omega5', 'dim_tv4', 'dim_tv5',
                     'dim_ker_rel4', 'rk_psi', 'rk_im_d4_Tv', 'rk_ker_d3_Tv']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # The key relationship: rk_psi vs rk_im_d4_Tv
        same = sum(1 for r in results if r['rk_psi'] == r['rk_im_d4_Tv'])
        print(f"\n  rk(ψ(ker)) = rk(im d_4^Tv): {same}/{len(results)}")

        # Cases where ψ escapes
        escapes = [r for r in results if not r['psi_in_im']]
        if escapes:
            print(f"\n  --- {len(escapes)} cases where ψ escapes im(d_4^Tv) ---")
            for r in escapes[:10]:
                print(f"    b4={r['b4']}, b4_Tv={r['b4_Tv']}, rk_psi={r['rk_psi']}, "
                      f"rk_d4Tv={r['rk_im_d4_Tv']}, dim_ker_rel={r['dim_ker_rel4']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
