"""
les_codim1_proof.py — Prove HYP-408 via the Long Exact Sequence

The LES of the pair (T, T\\v):
  ... → H_4(T) → H_4(T,T\\v) → H_3(T\\v) →[i_*] H_3(T) → H_3(T,T\\v) → H_2(T\\v) → ...

Key facts for b3(T)=1, b3(T\\v)=1:
  b4(T)=0, b4(T\\v)=0 (HYP-406, confirmed)
  b2(T)=0, b2(T\\v)=0 (HYP-322, confirmed)

So the LES simplifies to:
  0 → H_4(T,T\\v) → H_3(T\\v) →[i_*] H_3(T) → H_3(T,T\\v) → 0

This gives us:
  dim H_4(T,T\\v) = dim ker(i_*) = b3(T\\v) - rank(i_*)
  dim H_3(T,T\\v) = b3(T) - rank(i_*) = 1 - rank(i_*)

When rank(i_*)=1: H_4^rel=0, H_3^rel=0
When rank(i_*)=0: H_4^rel=1, H_3^rel=1

Now, HYP-408 says: codim(im(d_4^T)_old, ker(d_3^T)_old) = 1.

Let's see if this follows from the LES. The "old projection" of ker(d_3^T)
is related to the restriction map ker(d_3^T) → C_3(T\\v) (not exactly,
since the restriction is to subcomplex, not quotient).

Actually, let me think about this differently. The inclusion i: T\\v → T
gives a chain map i_#: C_*(T\\v) → C_*(T). The "old projection"
π: C_3(T) → C_3(T) restricted to old-indexed coords is NOT a chain map.

But there IS a natural projection: the restriction map
r: C_*(T) → C_*(T\\v) defined by r(σ) = σ if σ doesn't use v, = 0 otherwise.
This IS a chain map! Because d commutes with restriction to subcomplex.

So r ∘ i = id on C_*(T\\v), making (i, r) a retraction pair.

This means: H_*(T\\v) is a RETRACT of H_*(T).
In particular: i_*: H_3(T\\v) → H_3(T) is always INJECTIVE!

Wait — that can't be right, because i_* FAILS at n=8.
Let me check this more carefully...

The issue is whether r is actually a chain map. Let's verify:
  d_p ∘ r = r ∘ d_p?

For a p-chain σ = (v_0, ..., v_p):
  d_p(σ) = Σ (-1)^i (v_0, ..., hat{v_i}, ..., v_p)

If σ doesn't use v: r(σ) = σ, d(r(σ)) = d(σ), and r(d(σ)) = d(σ)
  since all faces also don't use v. ✓

If σ uses v at position k: r(σ) = 0.
  d(r(σ)) = 0.
  r(d(σ)) = r(Σ (-1)^i face_i)
           = Σ_{i ≠ k} (-1)^i r(face_i) + (-1)^k r(face_k)
  face_k = (v_0, ..., hat{v}, ..., v_p) doesn't use v, so r(face_k) = face_k.
  face_i (i ≠ k) still uses v, so r(face_i) = 0.
  Therefore r(d(σ)) = (-1)^k face_k ≠ 0 in general!

So d ∘ r ≠ r ∘ d. The restriction r is NOT a chain map for GLMY path homology!

This is the KEY difference from simplicial homology (where r IS a chain map).
In GLMY homology, the Omega constraints (the quotient by non-allowed paths)
break the chain map property.

Wait — actually, in GLMY path homology, we work with ALLOWED paths only.
The chain groups are C_p = span of allowed p-paths. A p-path σ is "allowed"
if it satisfies certain tournament-dependent conditions.

The restriction r sends an allowed path in T to the same path if it doesn't
use v. But the face of a through-v path, with v deleted, may or may not be
allowed in T\\v. Actually from S57 we know:
- All old 3-paths ARE allowed in T\\v (no orphans)
- But v-deletion faces of 4-paths are NOT always allowed in T\\v (22% miss)

So even on the chain level, r doesn't send C_*(T) to C_*(T\\v) cleanly.

Let me verify this computationally: is the restriction a chain map?

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
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


def check_retraction_chain_map(A, n, v):
    """Check if restriction r: C_*(T) -> C_*(T\\v) is a chain map."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    results = {}

    for p in range(2, min(6, n)):
        paths_p_T = ap_T.get(p, [])
        paths_pm1_T = ap_T.get(p - 1, [])
        paths_p_Tv = set(ap_Tv.get(p, []))
        paths_pm1_Tv = set(ap_Tv.get(p - 1, []))

        idx_pm1_T = {path: i for i, path in enumerate(paths_pm1_T)}
        idx_pm1_Tv_list = ap_Tv.get(p - 1, [])
        idx_pm1_Tv = {path: i for i, path in enumerate(idx_pm1_Tv_list)}

        # For through-v p-paths, check: r(d(sigma)) vs d(r(sigma)) = 0
        n_through_v = 0
        n_rd_nonzero = 0  # cases where r(d(sigma)) ≠ 0
        n_rd_in_Tv = 0  # cases where r(d(sigma)) is in C_{p-1}(T\\v)

        for path in paths_p_T:
            if v not in path:
                continue
            n_through_v += 1

            # d(sigma) = sum of signed faces
            # r(d(sigma)) = keep only faces not containing v
            old_faces = []
            for sign, face in boundary_faces(path):
                if v not in face:
                    old_faces.append((sign, face))

            if old_faces:
                n_rd_nonzero += 1
                # Check if all old faces are allowed in T\\v
                all_in_Tv = True
                for sign, face in old_faces:
                    tv_face = tuple(remaining_inv[x] for x in face)
                    if tv_face not in paths_pm1_Tv:
                        all_in_Tv = False
                        break
                if all_in_Tv:
                    n_rd_in_Tv += 1

        results[p] = {
            'n_through_v': n_through_v,
            'n_rd_nonzero': n_rd_nonzero,
            'n_rd_in_Tv': n_rd_in_Tv,
            'frac_nonzero': n_rd_nonzero / max(n_through_v, 1),
            'frac_in_Tv': n_rd_in_Tv / max(n_rd_nonzero, 1) if n_rd_nonzero > 0 else 1.0,
        }

    return results


def les_dimension_check(A, n, v):
    """Verify LES dimensions for the pair (T, T\\v)."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    bettis_T = {p: cc_T['bettis'].get(p, 0) for p in range(n)}
    bettis_Tv = {p: cc_Tv['bettis'].get(p, 0) for p in range(n - 1)}

    if bettis_T.get(3, 0) < 1 or bettis_Tv.get(3, 0) < 1:
        return None

    # Euler characteristics
    chi_T = sum((-1)**p * cc_T.get('omega_dims', {}).get(p, 0) for p in range(n))
    chi_Tv = sum((-1)**p * cc_Tv.get('omega_dims', {}).get(p, 0) for p in range(n - 1))
    chi_rel = chi_T - chi_Tv  # Euler char of relative complex

    # From LES: sum(-1)^p b_p^rel = chi_rel
    # b_p^rel from LES when consecutive groups vanish:
    # H_4^rel = ker(i_*: H_3(T\\v) -> H_3(T)) when H_4(T)=0
    # H_3^rel = coker(i_*: H_3(T\\v) -> H_3(T)) when H_2(T\\v)=0

    return {
        'bettis_T': bettis_T,
        'bettis_Tv': bettis_Tv,
        'chi_T': chi_T,
        'chi_Tv': chi_Tv,
        'chi_rel': chi_rel,
    }


def main():
    print("=" * 70)
    print("LES AND RETRACTION ANALYSIS")
    print("=" * 70)

    # Part 1: Check if restriction is a chain map
    print("\n--- RETRACTION CHAIN MAP CHECK ---")
    for n in [7, 8]:
        print(f"\nn={n}:")
        rng = np.random.RandomState(42)
        total_results = {}
        n_checked = 0

        for trial in range(500):
            if n_checked >= 30:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = check_retraction_chain_map(A, n, v_cand)
                for p, data in r.items():
                    if p not in total_results:
                        total_results[p] = []
                    total_results[p].append(data)
                n_checked += 1
                if n_checked >= 30:
                    break

        for p in sorted(total_results.keys()):
            data_list = total_results[p]
            frac_nz = np.mean([d['frac_nonzero'] for d in data_list])
            frac_tv = np.mean([d['frac_in_Tv'] for d in data_list])
            avg_thru = np.mean([d['n_through_v'] for d in data_list])
            print(f"  p={p}: {len(data_list)} checks, avg through-v={avg_thru:.0f}, "
                  f"r(d(σ))≠0: {frac_nz:.3f}, all faces in T\\v: {frac_tv:.3f}")

    # Part 2: LES dimension analysis
    print("\n\n--- LES DIMENSION ANALYSIS ---")
    for n in [7, 8]:
        print(f"\nn={n}:")
        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()

        for trial in range(5000):
            if len(results) >= 100:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = les_dimension_check(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD vertices, {elapsed:.1f}s")

        # Betti analysis
        for p in range(max(n, 8)):
            bT = [r['bettis_T'].get(p, 0) for r in results]
            bTv = [r['bettis_Tv'].get(p, 0) for r in results if p < n - 1]
            if any(b > 0 for b in bT):
                print(f"  b_{p}(T): {dict(sorted(set((b, bT.count(b)) for b in bT)))}")
            if bTv and any(b > 0 for b in bTv):
                print(f"  b_{p}(T\\v): {dict(sorted(set((b, bTv.count(b)) for b in bTv)))}")

        # Chi analysis
        chi_rels = [r['chi_rel'] for r in results]
        from collections import Counter
        print(f"  chi_rel distribution: {dict(sorted(Counter(chi_rels).items()))}")

        # KEY: relative Betti numbers from LES
        # When b2(T\\v)=0, b4(T)=0, b4(T\\v)=0:
        #   H_4^rel = ker(i_*)
        #   H_3^rel = coker(i_*)
        #   Alternating sum: -b_4^rel + b_3^rel = -ker(i_*) + coker(i_*)
        #   = -(b3Tv - rk_i*) + (b3T - rk_i*) = b3T - b3Tv
        # Since b3T=b3Tv=1: this gives 0 always.
        # So chi_rel contribution from p=3,4 is 0.
        # The full chi_rel = 0 (confirmed S57) constrains higher degrees too.

        print(f"\n  LES prediction: when b4=b2=0,")
        print(f"  b_4^rel = b3(Tv) - rk(i_*), b_3^rel = b3(T) - rk(i_*)")
        print(f"  Since b3(T)=b3(Tv)=1: b_4^rel = b_3^rel = 1-rk(i_*)")
        print(f"  So b_4^rel=b_3^rel ∈ {{0,1}} and they're equal.")


if __name__ == '__main__':
    main()
    print("\nDONE.")
