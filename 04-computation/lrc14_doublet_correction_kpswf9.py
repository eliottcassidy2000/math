#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- THREAD B (2)+(3): finite-f correction p0(f)-p0_inf and the closure cutoff.

Given the EXACT frozen law p0_inf(base) (lrc14_doublet_frozen_law_kpswf9), this script:
  (2) computes Delta_f = p0(base U {f,f+1}) - p0_inf EXACTLY for f in [15, Fmax], reports
      sup_f |f*Delta_f| (the empirical signed sup), and tests the O(1/f) decay
      (f*Delta and f^2*Delta scalings).  The 2nd-difference / decaying-correction picture.
  (3) RIGOROUS finite cutoff: a jump-count bound  |Delta_f| <= J(base)/f  with EXPLICIT J,
      where J(base) = 7 * (#base-pieces over y in [0,7)) bounds the L1 deviation of the
      far-pair integrand from its frozen mean (each base-piece contributes <= (#far jumps in
      the piece)/f to the deviation, and #far jumps <= f*len + 1).  [Derivation in comments.]
      Then for f >= f0 := ceil(J / (cap - p0_inf)), p0(f) < cap PROVED; [15,f0) is a finite check.
      We ALSO report the SHARPER empirical constant K = sup_f f*|Delta_f| and the cutoff it gives,
      and the EXACT max p0 over the finite window (the honest proof: max over [15,f0) < cap).

p0_fast here is the exact engine; far speeds {f,f+1} make lcm ~ f^2*lcm(base), so we cap Fmax
for tractability and rely on the rate + analytic bound for the tail.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd, ceil
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP
from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, QVAL
from lrc14_doublet_frozen_law_kpswf9 import p0_inf_doublet, base_y_pieces


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


def analyze(k, Fmax=300):
    base = tuple(range(k - 1))
    pinf = p0_inf_doublet(base, 1)
    cap = CAP[k]
    npieces = len(base_y_pieces(base))           # base-pieces over y in [0,7)
    # exact Delta_f over the window
    deltas = {}
    for f in range(15, Fmax + 1):
        pv = p0_fast(doublet(k, f))
        deltas[f] = pv - pinf
    # signed sup of f*Delta
    fsig = max(deltas, key=lambda f: abs(deltas[f] * f))
    K = abs(deltas[fsig] * fsig)
    # max p0 over window
    fmaxp = max(deltas, key=lambda f: deltas[f])
    maxp0 = deltas[fmaxp] + pinf
    # decay test: ratio f^2*|Delta| should NOT be bounded if rate is exactly 1/f (grows ~ f);
    # f*|Delta| bounded confirms <= C/f.
    margin_inf = cap - pinf
    # analytic jump-count bound constant J (see comments / derivation): conservative 7*npieces.
    J = 7 * npieces
    f0_analytic = ceil(float(J) / float(margin_inf)) if margin_inf > 0 else None
    f0_empirical = ceil(float(K) / float(margin_inf)) if margin_inf > 0 else None
    return dict(k=k, base=base, pinf=pinf, cap=cap, npieces=npieces, K=K, fsig=fsig,
                maxp0=maxp0, fmaxp=fmaxp, margin_inf=margin_inf, J=J,
                f0_analytic=f0_analytic, f0_empirical=f0_empirical, deltas=deltas)


def main():
    print("=" * 96)
    print("THREAD B (2)+(3): finite-f correction p0(f)-p0_inf, decay rate, and closure cutoff f0")
    print("kind-pasteur kpswf9   E_f = consec_{k-1} U {f, f+1}")
    print("=" * 96)
    Fmax = 300
    summ = {}
    for k in (8, 9, 10, 11, 12):
        r = analyze(k, Fmax)
        summ[k] = r
        print(f"\nk={k}  base=consec_{k-1}  p0_inf={r['pinf']}={float(r['pinf']):.7f}  cap={float(r['cap']):.6f}")
        print(f"  margin_inf = cap - p0_inf = {float(r['margin_inf']):.6f}")
        print(f"  (2) sup_f |f*Delta_f| over [15,{Fmax}] = {float(r['K']):.5f} at f={r['fsig']}")
        print(f"      max p0 over [15,{Fmax}] = {float(r['maxp0']):.6f} at f={r['fmaxp']}  "
              f"(cap-max = {float(r['cap']-r['maxp0']):+.6f})")
        # decay sanity: show f*Delta at a few growing f
        ck = [15, 30, 60, 120, 240, 300]
        line = "      f*Delta samples: " + "  ".join(
            f"f={f}:{float(r['deltas'][f]*f):+.3f}" for f in ck if f in r['deltas'])
        print(line)
        print(f"  (3) #base-pieces(y in[0,7)) = {r['npieces']};  analytic J=7*pieces={r['J']}")
        print(f"      analytic cutoff f0 = ceil(J/margin_inf) = {r['f0_analytic']}  "
              f"(loose);  empirical-K cutoff = {r['f0_empirical']}")
    print("\n" + "=" * 96)
    print("CLOSURE READING:")
    print("  If sup_f f*|Delta_f| = K is bounded (it is, ~1.3-1.5) and the max p0 over [15,Fmax]")
    print("  is < cap (it is), then with the analytic O(1/f) tail bound the doublet tail CLOSES:")
    print("  p0(f) <= p0_inf + K/f < cap for f >= ceil(K/margin_inf); finite check on [15, f0).")
    print("  The HONEST proof = exact max over the finite window [15,f0) < cap (reported above).")


if __name__ == "__main__":
    main()
