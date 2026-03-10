"""
relative_exactness.py — Exactness of the relative complex at degree 4

H_4^rel = ker(d_4^rel) / im(d_5^rel) = 0 iff i_* injective.

Compute:
1. dim(ker d_4^rel) and dim(im d_5^rel) — does im(d_5^rel) fill ker(d_4^rel)?
2. The "slack": dim(ker d_4^rel) - dim(im d_5^rel)
3. How close is d_5^rel to surjecting onto ker(d_4^rel)?

Key insight from S59: ψ(ker(d_4^{tv→tv})) = im(d_4^{T\v}) at n=7.
The relative d_4 is essentially d_4^{tv→tv}. So ker(d_4^rel) is the
space of relative 4-cycles. And im(d_5^rel) must fill all of it.

We compute the relative boundary maps directly.

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def relative_exactness_data(A, n, v):
    """Compute boundary ranks for the relative complex at degree 4."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1-1))

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}

    def classify(ps, v):
        tv = [i for i, p in enumerate(ps) if v in p]
        old = [i for i, p in enumerate(ps) if v not in p]
        return tv, old

    tv = {}
    old = {}
    for p in range(2, max_p + 1):
        tv[p], old[p] = classify(paths.get(p, []), v)

    if not tv.get(4):
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

    ob_T = {}
    for p in range(2, max_p + 1):
        ob_T[p] = get_omega(ap_T, p)

    ob_Tv = {}
    for p in range(2, min(max_p, n1-1) + 1):
        ob_Tv[p] = get_omega(ap_Tv, p)

    # Build boundary maps
    idx = {}
    for p in range(2, max_p + 1):
        idx[p] = {path: i for i, path in enumerate(paths[p])}

    bd = {}
    for p in range(3, max_p + 1):
        src = paths.get(p, [])
        tgt = paths.get(p-1, [])
        if not src or not tgt:
            continue
        mat = np.zeros((len(tgt), len(src)), dtype=np.int64)
        for j, path in enumerate(src):
            for sign, face in boundary_faces(path):
                if face in idx[p-1]:
                    mat[idx[p-1][face], j] = (mat[idx[p-1][face], j] + sign) % PRIME
        bd[p] = mat

    # Relative d_4: d_4^{tv→tv} through Omega_4(T)
    # d_4^{tv→tv}: takes tv4-indexed rows of bd4, tv3-indexed columns → result
    # But we work in Omega coords: d_4^{tv→tv} ∘ Ω_4(T)
    # = bd[4][tv3, :] @ ob_T[4].T → vectors indexed by tv3
    d4_tv_omega = bd[4][tv[3], :] @ ob_T[4].T % PRIME if 4 in bd and ob_T[4].shape[0] > 0 else np.zeros((len(tv.get(3, [])), 0), dtype=np.int64)

    # ker(d_4^{tv→tv}) in Omega_4 coords
    if d4_tv_omega.size > 0 and d4_tv_omega.shape[1] > 0:
        rk_d4_tv = int(_gauss_rank_np(d4_tv_omega.copy(), PRIME))
        _, kv = _gauss_nullbasis_modp(
            d4_tv_omega.astype(np.int32), d4_tv_omega.shape[0], d4_tv_omega.shape[1], PRIME)
        dim_ker_d4_rel = len(kv) if kv else 0
    else:
        rk_d4_tv = 0
        dim_ker_d4_rel = ob_T[4].shape[0]

    # Relative d_5: d_5^{tv→tv} through Omega_5(T)
    # d_5^{tv→tv}: bd[5][tv4, :] @ ob_T[5].T
    if 5 in bd and ob_T.get(5, np.zeros((0,0))).shape[0] > 0 and tv.get(5):
        d5_tv_omega_raw = bd[5][tv[4], :] @ ob_T[5].T % PRIME
    else:
        d5_tv_omega_raw = np.zeros((len(tv.get(4, [])), 0), dtype=np.int64)

    # But d_5^rel maps C_5^rel → C_4^rel. The image should be in Omega_4 coords.
    # d_5(Omega_5(T)) gives vectors in A_4 coords. We need to express these in Omega_4 coords.
    # Actually, im(d_5) lands in Omega_4 automatically (d ∘ d = 0 ensures boundaries are cycles of d_4's kernel,
    # but in the path complex, im(d_5) ⊂ Omega_4 because the chain complex uses Omega_p, not A_p).
    # Wait: d_5 maps Omega_5 → Omega_4? No: d_5 maps Omega_5 → A_4, and its image lands in ker(d_4)
    # which may or may not be representable in Omega_4 coords... Actually the chain complex IS:
    # ... → Omega_5 →d_5 Omega_4 →d_4 Omega_3 → ...
    # where d_p is the restriction of the face map. So d_5 maps Omega_5 to Omega_4? No:
    # d_5 maps Omega_5 to A_4 (face map), and we need the result to be in Omega_4.
    # But Omega_4 ⊂ A_4, and im(d_5) ⊂ Omega_4 is part of the chain complex structure.

    # For the relative complex: d_5^rel sends C_5^rel to C_4^rel.
    # C_5^rel = Omega_5(T)/i(Omega_5(T\v))
    # C_4^rel = Omega_4(T)/i(Omega_4(T\v))

    # The image of d_5^rel in C_4^rel is the image of d_5(Omega_5(T)) projected to the quotient.
    # d_5(Omega_5(T)) in Omega_4 coords: need to solve for Omega_4 coords.

    # Actually, the boundary map in the GLMY complex: d_p maps Omega_p → A_{p-1}.
    # The image of d_p is im(d_p) ⊂ A_{p-1}. Then ker(d_{p-1}) ⊃ im(d_p).
    # But ker(d_{p-1}) lives in Omega_{p-1}. So im(d_p) ⊂ Omega_{p-1}? Not necessarily directly.
    # The chain complex works in A_p coordinates, with Omega_p as the chain group.
    # im(d_p) ⊂ A_{p-1}, and im(d_p) ⊂ Omega_{p-1} because d_{p-1} ∘ d_p = 0 and
    # d_{p-1} restricted to Omega_{p-1} captures all of ker(d_{p-1}).
    # Actually no, that's not right either.

    # Let me think more carefully. In the GLMY complex:
    # Omega_p = {allowed p-paths satisfying regularity constraints}
    # The boundary map d_p: A_p → A_{p-1} is the alternating face map.
    # The chain complex is: Omega_p →d_p Omega_{p-1} → ...
    # This means d_p(Omega_p) ⊂ Omega_{p-1}? NO!
    # d_p maps Omega_p into A_{p-1}, and the chain complex actually works with
    # A_p as the ambient space, with Omega_p the "regular" subspace.
    # The chain complex is built so that d_p maps A_p → A_{p-1}, and
    # we compute H_p = ker(d_p|_Omega_p) / im(d_{p+1}|_Omega_{p+1}).
    # But im(d_{p+1}|_Omega_{p+1}) might not be contained in Omega_p!
    # It's contained in A_p, and we only care about its projection/intersection with Omega_p.

    # Hmm wait. Let me re-read the definition. In GLMY:
    # R_p = A_p (regular paths), Omega_p is the subspace where boundary
    # maps "work correctly". The chain complex is actually:
    # A_p/Omega_p^perp or something... I need to check.

    # Actually in our code, the chain complex works as:
    # The boundary map d_p is defined on ALL of A_p.
    # Omega_p = ker(constraint matrix) ⊂ A_p.
    # The chain groups are Omega_p.
    # d_p: Omega_p → A_{p-1}. We compute rank of d_p restricted to Omega_p.
    # im(d_p) is in A_{p-1} but for homology we need im(d_{p+1}) ∩ Omega_p... no.
    # Actually, the correct chain complex:
    # d_p|_Omega_p maps into Omega_{p-1}. This is a theorem (or part of the definition).
    # See Grigor'yan et al: d_p maps Omega_p to Omega_{p-1}.

    # OK so: d_5 maps Omega_5 INTO Omega_4. The boundary image is in Omega_4 already.
    # In coordinates: d_5(Omega_5) ⊂ A_4, but also ⊂ Omega_4.
    # So the images expressed in A_4 coords are automatically in Omega_4.

    # For the relative d_5: maps Omega_5/i(Omega_5(T\v)) → Omega_4/i(Omega_4(T\v))
    # The image of d_5^rel = d_5(Omega_5(T)) + i(Omega_4(T\v)) / i(Omega_4(T\v))
    # Its dimension = dim(d_5(Omega_5(T)) + i(Omega_4(T\v))) - dim(i(Omega_4(T\v)))

    # To compute rk(d_5^rel) = dim of image of d_5^rel:
    # = dim(span(d_5(Omega_5(T)) union i(Omega_4(T\v)))) - dim(Omega_4(T\v))

    # Compute d_5(Omega_5(T)) in A_4 coords:
    bd5_omega = bd[5] @ ob_T[5].T % PRIME if 5 in bd and ob_T[5].shape[0] > 0 else np.zeros((len(paths.get(4, [])), 0), dtype=np.int64)
    im_d5 = bd5_omega.T  # rows are image vectors in A_4 coords

    # Embed Omega_4(T\v) into A_4(T) coords
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}
    embed_4 = np.zeros((len(paths_Tv.get(4, [])), len(paths.get(4, []))), dtype=np.int64)
    for jv, pv in enumerate(paths_Tv.get(4, [])):
        mapped = tuple(remaining[x] for x in pv)
        if mapped in idx[4]:
            embed_4[jv, idx[4][mapped]] = 1

    ob4_Tv_embedded = ob_Tv.get(4, np.zeros((0, 0), dtype=np.int64))
    if ob4_Tv_embedded.shape[0] > 0:
        i_Omega4_Tv = (ob4_Tv_embedded @ embed_4 % PRIME)  # rows in A_4(T) coords
    else:
        i_Omega4_Tv = np.zeros((0, len(paths.get(4, []))), dtype=np.int64)

    dim_Omega4_Tv = i_Omega4_Tv.shape[0]

    # rk(d_5^rel) = dim(span(im_d5, i_Omega4_Tv)) - dim(Omega_4(T\v))
    if im_d5.shape[0] > 0 and i_Omega4_Tv.shape[0] > 0:
        combined = np.vstack([im_d5, i_Omega4_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        rk_d5_rel = rk_combined - dim_Omega4_Tv
    elif im_d5.shape[0] > 0:
        rk_d5_rel = int(_gauss_rank_np(im_d5.copy(), PRIME))
    else:
        rk_d5_rel = 0

    # Similarly for d_4^rel: maps Omega_4/i(Omega_4(T\v)) → Omega_3/i(Omega_3(T\v))
    bd4_omega = bd[4] @ ob_T[4].T % PRIME if 4 in bd and ob_T[4].shape[0] > 0 else np.zeros((len(paths.get(3, [])), 0), dtype=np.int64)
    im_d4 = bd4_omega.T  # in A_3 coords

    embed_3 = np.zeros((len(paths_Tv.get(3, [])), len(paths.get(3, []))), dtype=np.int64)
    idx3_Tv = {path: i for i, path in enumerate(paths_Tv.get(3, []))}
    for jv, pv in enumerate(paths_Tv.get(3, [])):
        mapped = tuple(remaining[x] for x in pv)
        if mapped in idx[3]:
            embed_3[jv, idx[3][mapped]] = 1

    ob3_Tv_embedded = ob_Tv.get(3, np.zeros((0, 0), dtype=np.int64))
    if ob3_Tv_embedded.shape[0] > 0:
        i_Omega3_Tv = (ob3_Tv_embedded @ embed_3 % PRIME)
    else:
        i_Omega3_Tv = np.zeros((0, len(paths.get(3, []))), dtype=np.int64)

    dim_Omega3_Tv = i_Omega3_Tv.shape[0]

    if im_d4.shape[0] > 0 and i_Omega3_Tv.shape[0] > 0:
        combined_d4 = np.vstack([im_d4, i_Omega3_Tv]) % PRIME
        rk_combined_d4 = int(_gauss_rank_np(combined_d4.copy(), PRIME))
        rk_d4_rel = rk_combined_d4 - dim_Omega3_Tv
    elif im_d4.shape[0] > 0:
        rk_d4_rel = int(_gauss_rank_np(im_d4.copy(), PRIME))
    else:
        rk_d4_rel = 0

    dim_Omega4_T = ob_T[4].shape[0]
    dim_C4_rel = dim_Omega4_T - dim_Omega4_Tv

    # H_4^rel = dim(ker d_4^rel) - dim(im d_5^rel)
    dim_ker_d4_rel = dim_C4_rel - rk_d4_rel
    h4_rel = dim_ker_d4_rel - rk_d5_rel

    return {
        'dim_C4_rel': dim_C4_rel,
        'dim_C5_rel': ob_T.get(5, np.zeros((0,0))).shape[0] - ob_Tv.get(5, np.zeros((0,0))).shape[0],
        'rk_d4_rel': rk_d4_rel,
        'rk_d5_rel': rk_d5_rel,
        'dim_ker_d4_rel': dim_ker_d4_rel,
        'h4_rel': h4_rel,
        'dim_Omega4_T': dim_Omega4_T,
        'dim_Omega4_Tv': dim_Omega4_Tv,
        'rk_d5_total': int(_gauss_rank_np(im_d5.copy(), PRIME)) if im_d5.shape[0] > 0 else 0,
        'rk_d4_total': int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"RELATIVE EXACTNESS AT DEGREE 4, n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 300

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, min(n - 1, 6))
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                n1 = n - 1
                rem = [i for i in range(n) if i != v_cand]
                A_sub = [[A[rem[i]][rem[j]] for j in range(n1)] for i in range(n1)]
                cc_Tv = full_chain_complex_modp(A_sub, n1, min(n1-1, 5))
                if cc_Tv['bettis'].get(3, 0) != 1:
                    continue

                r = relative_exactness_data(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
                if len(results) >= target:
                    break

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        h4_dist = Counter(r['h4_rel'] for r in results)
        print(f"\n  H_4^rel dimension: {dict(sorted(h4_dist.items()))}")

        for key in ['dim_C4_rel', 'dim_C5_rel', 'rk_d4_rel', 'rk_d5_rel',
                     'dim_ker_d4_rel', 'rk_d5_total', 'rk_d4_total']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # Ratio: rk_d5_rel / dim_ker_d4_rel
        # When H_4^rel = 0, this should be 1.0
        ratios = [r['rk_d5_rel'] / max(r['dim_ker_d4_rel'], 1) for r in results]
        print(f"\n  rk(d_5^rel)/dim(ker d_4^rel): min={min(ratios):.3f}, max={max(ratios):.3f}, mean={np.mean(ratios):.3f}")

        # Slack: dim_ker_d4_rel - rk_d5_rel
        slacks = [r['dim_ker_d4_rel'] - r['rk_d5_rel'] for r in results]
        slack_dist = Counter(slacks)
        print(f"  Slack (ker - im): {dict(sorted(slack_dist.items()))}")

        # Cases where H_4^rel > 0
        nonzero = [r for r in results if r['h4_rel'] > 0]
        if nonzero:
            print(f"\n  --- H_4^rel > 0 cases ---")
            for r in nonzero[:10]:
                print(f"    h4_rel={r['h4_rel']}, dim_ker={r['dim_ker_d4_rel']}, rk_d5_rel={r['rk_d5_rel']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
