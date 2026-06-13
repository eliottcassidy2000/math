"""
omega_dim_comparison.py — Compare Omega dimensions T vs T\v at n=7 and n=8

For the i_*-injectivity mechanism, the key is:
  New 5-chains (through v) have boundary d_5 that creates im(d_5^new) in Ω_4
  New 4-chains (through v) have boundary d_4^new in Ω_3
  If im(d_4^new) doesn't reach the H_3(T\v) generator, i_* is injective.

The dimensional picture:
  dim(Ω_p^new) = dim(Ω_p(T)) - dim(Ω_p(T\v))
  These are the "new" chains introduced by vertex v.

At n=7: does the relative chain complex have dimensional constraints
that FORCE H_4^rel = 0?

Author: opus-2026-03-09-S57
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, deletion_adj,
    full_chain_complex_modp,
    RANK_PRIME
)


def main():
    print("=" * 70)
    print("OMEGA DIMENSION COMPARISON: T vs T\\v")
    print("=" * 70)

    for n in [7, 8]:
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        rng = np.random.RandomState(42)
        max_p = n - 1

        # Collect: for b3=1 tournaments, dim(Omega_p(T)) and dim(Omega_p(T\v))
        # for BAD vertices only
        data_bad = []
        found = 0
        target = 200 if n == 7 else 100
        t0 = time.time()

        for trial in range(20000):
            A = random_tournament(n, rng)
            cc_T = full_chain_complex_modp(A, n, max_p)
            if cc_T['bettis'].get(3, 0) != 1:
                continue

            omega_T = {p: cc_T.get('omega_dims', {}).get(p, 0) for p in range(n)}
            ranks_T = {p: cc_T.get('ranks', {}).get(p, 0) for p in range(n)}
            bettis_T = {p: cc_T['bettis'].get(p, 0) for p in range(n)}

            for v in range(n):
                Av, nv = deletion_adj(A, n, v)
                cc_Tv = full_chain_complex_modp(Av, nv, nv - 1)
                b3v = cc_Tv['bettis'].get(3, 0)

                if b3v == 0:
                    continue  # GOOD vertex

                omega_Tv = {p: cc_Tv.get('omega_dims', {}).get(p, 0) for p in range(n - 1)}
                ranks_Tv = {p: cc_Tv.get('ranks', {}).get(p, 0) for p in range(n - 1)}
                bettis_Tv = {p: cc_Tv['bettis'].get(p, 0) for p in range(n - 1)}

                # Relative dimensions
                omega_rel = {p: omega_T.get(p, 0) - omega_Tv.get(p, 0) for p in range(n)}
                rank_rel_d = {}
                for p in range(n):
                    # rank of d_p^rel = rank(d_p^T) - rank(d_p^{T\v})
                    # (this is an approximation — true only under certain conditions)
                    rank_rel_d[p] = ranks_T.get(p, 0) - ranks_Tv.get(p, 0)

                data_bad.append({
                    'omega_T': omega_T, 'omega_Tv': omega_Tv, 'omega_rel': omega_rel,
                    'ranks_T': ranks_T, 'ranks_Tv': ranks_Tv, 'rank_rel': rank_rel_d,
                    'bettis_T': bettis_T, 'bettis_Tv': bettis_Tv,
                    'trial': trial, 'v': v,
                })
                found += 1

            if found >= target:
                break

        elapsed = time.time() - t0
        print(f"  {found} BAD vertices collected, {elapsed:.1f}s")

        # Summary statistics
        print(f"\n  Omega_p(T) for b3=1 tournaments:")
        for p in range(n):
            vals = [d['omega_T'][p] for d in data_bad]
            if vals:
                print(f"    Ω_{p}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        print(f"\n  Omega_p(T\\v) for BAD vertices:")
        for p in range(n - 1):
            vals = [d['omega_Tv'][p] for d in data_bad]
            if vals:
                print(f"    Ω_{p}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        print(f"\n  Relative dimensions Ω_p^rel = Ω_p(T) - Ω_p(T\\v):")
        for p in range(n):
            vals = [d['omega_rel'][p] for d in data_bad]
            if vals:
                print(f"    Ω_{p}^rel: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        print(f"\n  Relative ranks (d_p^T - d_p^Tv):")
        for p in range(1, n):
            vals = [d['rank_rel'][p] for d in data_bad]
            if vals:
                print(f"    rk_rel(d_{p}): min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # KEY: relative Euler characteristic
        chi_rel_vals = []
        for d in data_bad:
            chi_rel = sum((-1)**p * d['omega_rel'].get(p, 0) for p in range(n))
            chi_rel_vals.append(chi_rel)

        chi_rel_dist = Counter(chi_rel_vals)
        print(f"\n  Relative Euler char chi(T,T\\v) = chi(T) - chi(T\\v):")
        print(f"    Distribution: {dict(sorted(chi_rel_dist.items()))}")

        # KEY: H_4 relative from LES
        # When b4(T)=0, b4(T\v)=0: H_4^rel = ker(i_*)
        # When b4(T)=0, b4(T\v)>0: more complex
        print(f"\n  Betti distribution:")
        b4_T_dist = Counter(d['bettis_T'].get(4, 0) for d in data_bad)
        b4_Tv_dist = Counter(d['bettis_Tv'].get(4, 0) for d in data_bad)
        b5_T_dist = Counter(d['bettis_T'].get(5, 0) for d in data_bad)
        print(f"    b4(T):   {dict(sorted(b4_T_dist.items()))}")
        print(f"    b4(T\\v): {dict(sorted(b4_Tv_dist.items()))}")
        print(f"    b5(T):   {dict(sorted(b5_T_dist.items()))}")

        # Ratio analysis: dim(Ω_4^rel) / dim(Ω_3^rel) ?
        print(f"\n  Dimension ratios:")
        for p in [3, 4, 5]:
            vals = [d['omega_rel'].get(p, 0) / d['omega_rel'].get(p - 1, 1)
                    if d['omega_rel'].get(p - 1, 0) > 0 else float('inf')
                    for d in data_bad]
            finite_vals = [v for v in vals if v != float('inf')]
            if finite_vals:
                print(f"    Ω_{p}^rel / Ω_{p-1}^rel: min={min(finite_vals):.2f}, "
                      f"max={max(finite_vals):.2f}, mean={np.mean(finite_vals):.2f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
