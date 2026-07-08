"""
lrc14_var_resonance_opus_S150.py   (opus-2026-07-08-S150, HYP-5377)

THE UNIFICATION: far <= E[W]^2  <=>  Var(W) <= near  (exact algebra), so BOTH my near/far
route (far<=E[W]^2) and kps brick (B) (R2<=614 => D3>=bar) reduce to ONE quantity: Var(W).

This computes, exactly (Farey integration):
  (1) verify far <= E[W]^2  <=>  Var(W) <= near  (identity, since far = E[W^2]-near);
  (2) THE RESONANCE RELATION Var(W) vs R2 = sum_d r_d^2 (reduced additive energy): is
      Var(W) = c*R2 exactly (klein-S179 fit c~5.67e-5), or Var(W) = sum_d r_d * kernel(d)?
      Fit both; report the exact per-difference structure.  KEY for brick (B): Var(W) <= c*R2.
  (3) brick (B) exact reduction: D3 >= bar  <=  Var(W) <= f(E[W],E[W^3]); combine with
      Var(W) <= c*R2 and R2 <= 614 to get the EXACT c-threshold and margin.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys
from collections import Counter

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import pz_exact, BAR
from lrc14_pz_degree3_floor_opus_S148 import moments_exact
from lrc14_k11_tail_sharp_near_opus_S149 import EW_and_near

M = F(6, 7)
bar = BAR[11]


def R2_reduced(E):
    c = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                c[E[i] - E[j]] += 1
    return sum(v * v for v in c.values())


def D3_of(m1, m2, m3):
    u = m1 - m2 / M; v = m2 - m3 / M
    if v <= 0:
        return None
    return m1 / M + u * u / v


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def main():
    print("=" * 96)
    print("(1) far <= E[W]^2  <=>  Var(W) <= near  (exact identity check)")
    print("=" * 96)
    tests = [list(range(11)), [0,2,3,4,5,6,7,8,9,10,12], [0]+list(range(1,10))+[40],
             list(range(6))+list(range(40,45))]
    for E in tests:
        EW, near = EW_and_near(E)
        m = moments_exact(E, 3); m1, m2, m3 = m[1], m[2], m[3]
        var = m2 - m1 * m1
        far = m2 - near
        eq = (far <= m1 * m1) == (var <= near)
        print(f"  E={str(E[:5]):20s}...: Var={float(var):.5f} near={float(near):.5f} "
              f"far={float(far):.5f} E[W]^2={float(m1*m1):.5f}  "
              f"[far<=E[W]^2]={far<=m1*m1} == [Var<=near]={var<=near}: {eq}")
    print("  => far <= E[W]^2  IFF  Var(W) <= near  (unifies my route + brick B: both bound Var).")

    print()
    print("=" * 96)
    print("(2) Var(W) vs R2 = sum_d r_d^2 : the resonance relation (klein-S179 fit c~5.67e-5)")
    print("=" * 96)
    pts = []
    for E in tests + [[0,4,6,8,10,11,12,14,16,17,18]]:
        EW, near = EW_and_near(E)
        m = moments_exact(E, 3)
        var = m[2] - m[1] * m[1]
        R2 = R2_reduced(E)
        pts.append((R2, var, near))
        print(f"  E={str(E[:5]):20s}...: R2={R2:5d}  Var(W)={float(var):.6f}  "
              f"Var/R2={float(var)/R2:.3e}  near={float(near):.5f}  Var<=near:{var<=near}")
    # linear fit Var = c*R2 through origin (least squares on the exact points)
    sxx = sum(r * r for r, v, n in pts)
    sxy = sum(r * float(v) for r, v, n in pts)
    c_fit = sxy / sxx
    print(f"  best c (Var ~ c*R2, through origin) = {c_fit:.3e}  (klein 5.67e-5)")
    print(f"  => IF Var(W) <= c*R2 with c<=~6e-5, then R2<=614 gives Var<= {c_fit*614:.5f}")

    print()
    print("=" * 96)
    print("(3) BRICK (B) EXACT REDUCTION: what Var-bound gives D3 >= bar, over R2 <= 614?")
    print("=" * 96)
    # D3 >= bar depends on E[W],E[W^2],E[W^3]. Using D3 >= PZ = E[W]^2/E[W^2] = 1/(1+Var/E[W]^2):
    #   PZ >= bar  <=>  Var <= (1/bar - 1) E[W]^2.  With E[W]>= EWmin (k=11 floor):
    print("  D3 >= PZ = 1/(1 + Var/E[W]^2);  PZ >= bar  <=>  Var <= (1/bar - 1)*E[W]^2")
    ratio = (1 / bar - 1)
    print(f"  (1/bar - 1) = {float(ratio):.4f};  so need Var <= {float(ratio):.4f} * E[W]^2")
    # k=11 E[W] floor: exhaustive over compact + spread lower bound. Empirical min E[W]:
    mnEW = None
    for D in range(10, 17):
        for mid in combinations(range(1, D), 9):
            E = [0] + list(mid) + [D]
            if gcd_all([E[i+1]-E[i] for i in range(len(E)-1)]) != 1:
                continue
            EW, _ = EW_and_near(E)
            if mnEW is None or EW < mnEW:
                mnEW = EW; mnEWarg = E
    print(f"  min E[W] over compact k=11 (diam<=16): {float(mnEW):.5f} at {mnEWarg}")
    print(f"  => Var <= {float(ratio)*float(mnEW)**2:.5f} suffices for PZ>=bar (hence D3>=bar).")
    print(f"     With Var <= c*R2 (c~5.67e-5) and R2<=614: Var <= {5.67e-5*614:.5f}")
    print(f"     PZ>=bar needs Var <= {float(ratio)*float(mnEW)**2:.5f}: "
          f"{'HOLDS (brick B via PZ closes)' if 5.67e-5*614 <= float(ratio)*float(mnEW)**2 else 'PZ too weak here -- need D3s extra room'}")
    # D3-based (looser): D3 >= bar allows larger Var since D3 > PZ
    print("  [D3 > PZ gives MORE room; the exact D3>=bar threshold is looser than the PZ one above.]")


if __name__ == "__main__":
    main()
