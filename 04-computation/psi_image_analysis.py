"""
psi_image_analysis.py — Why does ψ(ker(d_4^{tv→tv})) ⊂ im(d_4^{T\v})?

The chain criterion says rank(i_*) = 0 iff ψ(ker(d_4^{tv→tv})) escapes
im(d_4^{T\v}) into H_3(T\v). At n=7 this never happens.

Key question: Is there a chain-complex reason ψ(ker) ⊂ im(d_4^{T\v})?

Note: ψ(ker) ⊂ ker(d_3^{T\v}) always (proven). And ker/im = H_3(T\v)
has dimension beta_3(T\v) = 1. So ψ(ker) avoids H_3 iff
ψ(ker) ⊂ im(d_4^{T\v}).

Dimensional check: if dim(ψ(ker)) is large compared to dim(im(d_4^{T\v})),
it SHOULD escape. If dim(ψ(ker)) << dim(ker(d_3^{T\v})), it might not.

Also: the relative complex C_*(T,T\v) has a connecting homomorphism
  δ: H_4(T,T\v) → H_3(T\v)
from the LES. This δ maps relative 4-cycles to their "boundary in T\v".
For a relative cycle [w_tv] with d_4^{tv→tv}(w_tv) = 0,
δ([w_tv]) = [ψ(w_tv)] in H_3(T\v).

So: δ(H_4^rel) = im(ψ on ker(d_4^{tv→tv})) projected to H_3(T\v).
And rank(i_*) = 0 iff H_4^rel ≠ 0 iff δ ≠ 0.

The LES says: H_4(T) → H_4(T,T\v) →[δ] H_3(T\v) →[i_*] H_3(T)
with exact sequence. When H_4(T) = 0:
  ker(δ) = 0, so δ is injective.
  im(δ) = ker(i_*).
So: rank(i_*) = dim(H_3(T\v)) - dim(im(δ)) = 1 - dim(H_4^rel).

rank(i_*) = 0 iff H_4^rel = 1 iff δ is nonzero.

But H_4^rel is the homology of the RELATIVE chain complex at degree 4.
The relative complex is C_*(T)/i(C_*(T\v)).

Let's compute H_4^rel directly and verify it matches.

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


def relative_h4_and_delta(A, n, v):
    """Compute H_4(T,T\\v) and the connecting homomorphism δ."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1 - 1) + 1)}

    # Classify paths as tv or old for degrees 3,4,5
    def classify(ps, v):
        tv = [i for i, p in enumerate(ps) if v in p]
        old = [i for i, p in enumerate(ps) if v not in p]
        return tv, old

    tv3, old3 = classify(paths.get(3, []), v)
    tv4, old4 = classify(paths.get(4, []), v)
    tv5, old5 = classify(paths.get(5, []), v) if 5 in paths and paths[5] else ([], [])

    if not tv4:
        return None

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
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)

    # Build boundary maps
    idx3_T = {p: i for i, p in enumerate(paths.get(3, []))}
    idx4_T = {p: i for i, p in enumerate(paths.get(4, []))}
    idx3_Tv = {p: i for i, p in enumerate(paths_Tv.get(3, []))}
    idx2_Tv = {p: i for i, p in enumerate(paths_Tv.get(2, []))}

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

    # d_4^{tv→tv}: bd4_T[tv3, tv4]
    d4_tvtv = bd4_T[np.ix_(tv3, tv4)] % PRIME if tv3 and tv4 else np.zeros((0, 0), dtype=np.int64)
    # ψ = d_4^{tv→old}: bd4_T[old3, tv4]
    psi = bd4_T[np.ix_(old3, tv4)] % PRIME if old3 and tv4 else np.zeros((0, 0), dtype=np.int64)

    # d_5^{tv→tv}: bd5_T[tv4, tv5]
    d5_tvtv = bd5_T[np.ix_(tv4, tv5)] % PRIME if tv4 and tv5 else np.zeros((len(tv4), 0), dtype=np.int64)

    # Relative d_4^rel = d_4^{tv→tv} on Omega_4^T
    # But need to work in Omega coords
    d4_tvtv_omega = d4_tvtv @ ob4_T[:, tv4].T % PRIME if ob4_T.shape[0] > 0 and tv4 else np.zeros((len(tv3), 0), dtype=np.int64)

    # ker(d_4^rel) = ker(d_4^{tv→tv} on Omega_4)
    if d4_tvtv_omega.size > 0:
        _, kv = _gauss_nullbasis_modp(
            d4_tvtv_omega.astype(np.int32), d4_tvtv_omega.shape[0], d4_tvtv_omega.shape[1], PRIME)
        ker_d4_rel = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
    else:
        ker_d4_rel = np.eye(ob4_T.shape[0], dtype=np.int64)

    dim_ker_d4_rel = ker_d4_rel.shape[0]

    # im(d_5^rel) = image of d_5^{tv→tv} on Omega_5, in Omega_4 coords
    # d_5^rel: Omega_5 → C_4^tv via d_5^{tv→tv}
    d5_tvtv_omega = d5_tvtv @ ob5_T[:, tv5].T % PRIME if ob5_T.shape[0] > 0 and tv5 else np.zeros((len(tv4), 0), dtype=np.int64)

    # Actually, d_5^rel maps Omega_5 vectors to A_4^tv vectors.
    # We need to see this in Omega_4 coords. But im(d_5) is in A_4 coords.
    # We need to express im(d_5^rel) as vectors in Omega_4 coords to
    # intersect with ker(d_4^rel).
    #
    # This is getting complicated. Let me just compute H_4^rel directly:
    # H_4^rel = ker(d_4^rel in Omega) / im(d_5^rel in Omega)
    #
    # d_4^rel sends Omega_4 vectors to C_3^tv = A_3^tv (tv part of boundary)
    # d_5^rel sends Omega_5 vectors to C_4^tv = A_4^tv (tv part of boundary)
    # But these live in different coordinate spaces...

    # SIMPLER: use the full relative complex approach.
    # The relative complex C_p^rel = Ω_p(T) / i(Ω_p(T\v))
    # with d^rel([x]) = [d(x)].
    # Since old→tv = 0, the tv-part of d(x) depends only on x's tv part.
    # And i(Ω_p(T\v)) = old part of Ω_p(T) (approximately).

    # For H_4^rel: we need the relative complex at degree 4.
    # C_4^rel = Ω_4(T) / i(Ω_4(T\v)).
    # d_4^rel: C_4^rel → C_3^rel = Ω_3(T) / i(Ω_3(T\v)).
    # d_5^rel: C_5^rel → C_4^rel.

    # Since i(Ω_p(T\v)) embeds as old-supported Omega vectors:
    # C_4^rel ≅ Ω_4(T) / old-Ω_4.
    # d_4^rel on this quotient: [w] ↦ [d_4(w)] in C_3^rel.
    # For [w] where w is an old Omega vector: d_4(w) is old, so [d_4(w)] = 0.
    # For [w] where w is tv: d_4(w) = tv_part + old_part, and [d_4(w)] = [tv_part].
    # So d_4^rel ≅ d_4^{tv→tv} on quotient representatives.

    # Let's just compute directly via matrix algebra.
    # H_4^rel = dim(ker of d_4^rel quotient map) - dim(im of d_5^rel quotient map)

    # For the quotient: pick a basis for Ω_4 / i(Ω_4(T\v)).
    # i(Ω_4(T\v)) in Ω_4(T) coords: embed ob4_Tv into Ω_4(T).

    # Map Ω_4(T\v) to A_4(T) coords then project to Ω_4(T) coords
    # Ω_4(T\v) → A_4(T\v) → A_4(T) via embedding → Ω_4(T) coords
    embed_4 = np.zeros((len(paths_Tv.get(4, [])), len(paths.get(4, []))), dtype=np.int64)
    for jv, pv in enumerate(paths_Tv.get(4, [])):
        mapped = tuple(remaining[x] for x in pv)
        if mapped in idx4_T:
            embed_4[jv, idx4_T[mapped]] = 1

    if ob4_Tv.shape[0] > 0:
        embedded_ob4_Tv_A = ob4_Tv @ embed_4 % PRIME  # in A_4(T) coords
    else:
        embedded_ob4_Tv_A = np.zeros((0, len(paths.get(4, []))), dtype=np.int64)

    # Now express in Ω_4(T) coords... this requires solving ob4_T @ x = embedded_A
    # This is expensive. Let me just compute im(d_4) and im(d_4^{T\v}) and compare.

    # Actually, let me just check: rank of ψ-image in H_3(T\v).
    # ψ applied to ker(d_4^{tv→tv}) vectors
    ker_A4 = ker_d4_rel @ ob4_T % PRIME
    # BUG FIX: use FULL boundary map restricted to old3 rows, not just tv→old block.
    # Kernel vectors have both tv4 and old4 components; the old4 components
    # contribute via d_4^{old→old} which the previous code missed.
    psi_img = (bd4_T[old3, :] @ ker_A4.T % PRIME).T if ker_d4_rel.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)

    # Map to T\v
    old3_map = {}
    for i_idx in old3:
        path_T = paths[3][i_idx]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_map[i_idx] = idx3_Tv[path_Tv]

    psi_Tv = np.zeros((psi_img.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
    for local_j, global_i in enumerate(old3):
        if global_i in old3_map:
            j_Tv = old3_map[global_i]
            psi_Tv[:, j_Tv] = (psi_Tv[:, j_Tv] + psi_img[:, local_j]) % PRIME

    rk_psi = int(_gauss_rank_np(psi_Tv.copy(), PRIME)) if psi_Tv.shape[0] > 0 else 0

    # im(d_4^{T\v})
    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

    # ker(d_3^{T\v})
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

    ker_d3_Tv_A = ker_d3_Tv @ ob3_Tv % PRIME
    beta3_Tv = ker_d3_Tv_A.shape[0] - rk_d4_Tv

    if beta3_Tv != 1:
        return None

    # Is ψ image in im(d_4^{T\v})?
    if psi_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        psi_in_im = (rk_combined == rk_d4_Tv)
        extra = rk_combined - rk_d4_Tv
    elif psi_Tv.shape[0] > 0:
        rk_psi_only = int(_gauss_rank_np(psi_Tv.copy(), PRIME))
        psi_in_im = (rk_psi_only == 0)
        extra = rk_psi_only
    else:
        psi_in_im = True
        extra = 0

    return {
        'dim_ker_d4_rel': dim_ker_d4_rel,
        'rk_psi': rk_psi,
        'rk_d4_Tv': rk_d4_Tv,
        'rk_ker_d3_Tv': ker_d3_Tv_A.shape[0],
        'beta3_Tv': beta3_Tv,
        'psi_in_im_d4Tv': psi_in_im,
        'extra': extra,
        'ratio_psi_to_im': rk_psi / max(rk_d4_Tv, 1),
        'ratio_psi_to_ker': rk_psi / max(ker_d3_Tv_A.shape[0], 1),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"ψ IMAGE vs im(d_4^Tv) AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(12345)
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 800

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = relative_h4_and_delta(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD (T,v) pairs, {elapsed:.1f}s")

        in_im = sum(1 for r in results if r['psi_in_im_d4Tv'])
        print(f"  ψ(ker) ⊂ im(d_4^Tv): {in_im}/{len(results)}")

        for key in ['dim_ker_d4_rel', 'rk_psi', 'rk_d4_Tv', 'rk_ker_d3_Tv']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # Ratio: how full is ψ(ker) relative to im(d_4^{T\v})?
        ratios_im = [r['ratio_psi_to_im'] for r in results]
        ratios_ker = [r['ratio_psi_to_ker'] for r in results]
        print(f"\n  rk(ψ(ker))/rk(im d_4^Tv): min={min(ratios_im):.3f}, max={max(ratios_im):.3f}, mean={np.mean(ratios_im):.3f}")
        print(f"  rk(ψ(ker))/rk(ker d_3^Tv): min={min(ratios_ker):.3f}, max={max(ratios_ker):.3f}, mean={np.mean(ratios_ker):.3f}")

        # Extra directions (= rank in H_3(T\v))
        extras = Counter(r['extra'] for r in results)
        print(f"\n  ψ(ker) directions in H_3(T\\v): {dict(sorted(extras.items()))}")

        if any(r['extra'] > 0 for r in results):
            print(f"\n  --- Cases where ψ reaches H_3 ---")
            for r in results:
                if r['extra'] > 0:
                    print(f"    rk_psi={r['rk_psi']}, rk_d4Tv={r['rk_d4_Tv']}, dim_ker_rel={r['dim_ker_d4_rel']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
