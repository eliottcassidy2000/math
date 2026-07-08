"""
lrc14_degree3_closed_form_floor_opus_S148.py   (opus-2026-07-08-S148, HYP-5327)

THE EXACT CLOSED-FORM DEGREE-3 COVERING FLOOR -- a sharpening of THM-660 (PZ = degree 2)
within mac-mini's THM-661 degree-d moment framework, applied to DE-RISK the binding k=11 leg.

mac-mini's THM-661: mu >= B_d(E) = max{ sum c_i E[W^i] : poly p = sum c_i w^i <= 1_{w>0}
on [0,6/7] }.  They use B_2 (=PZ) for k=11,12,13 -- leaving the k=11 margin THIN (+0.0159),
and B_4 (an LP) for k=8,9,10.

THIS SESSION identifies the OPTIMAL DEGREE-3 polynomial in CLOSED FORM (no LP):

    p(t) = 1 - (1 - t/r)^2 (1 - t/M),   M = 6/7,   r = r* := (m2 - m3/M)/(m1 - m2/M),

which satisfies p(0)=0, p(t) <= 1 on [0,M] for ANY r (both factors >= 0), so it is a valid
minorant of 1_{w>0}; optimizing r gives the EXACT rational floor

    D3(E) = E[p(W)] = m1/M + (m1 - m2/M)^2 / (m2 - m3/M)      (m_i = E[W^i], valid when
                                                               m2 - m3/M > 0),

strictly above B_2 = m1^2/m2 (Paley-Zygmund).  For the block:
    k=11: PZ +0.0159 -> D3 +0.0735   (4.6x thicker -- the binding leg de-risked)
    k=12: +0.109 -> +0.159 ;  k=13: +0.216 -> +0.257.
Exhaustive (diam <= 14): the BLOCK is the exact D3-minimizer, min D3 = 0.404751 >= bar,
margin +0.0735 uniformly.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_degree3_floor_opus_S148 import moments_exact
from lrc14_pz_general_integrator_opus_S148 import BAR

M = F(6, 7)


def D3_exact(E):
    """Exact degree-3 covering floor (closed form). Returns (D3, PZ, r*, valid)."""
    m = moments_exact(E, 3)
    m1, m2, m3 = m[1], m[2], m[3]
    u = m1 - m2 / M
    v = m2 - m3 / M
    PZ = m1 * m1 / m2
    if v <= 0 or u <= 0:
        return None, PZ, None, False
    rstar = v / u
    D3 = m1 / M + u * u / v
    return D3, PZ, rstar, (0 < rstar <= M)


def gcd_all(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g


def main():
    print("=" * 96)
    print("THE EXACT CLOSED-FORM DEGREE-3 COVERING FLOOR (sharpens THM-660 PZ, within THM-661)")
    print("  D3(E) = E[W]/M + (E[W]-E[W^2]/M)^2/(E[W^2]-E[W^3]/M),  M=6/7,  optimal deg-3 poly")
    print("  p(t)=1-(1-t/r*)^2(1-t/M),  r*=(E[W^2]-E[W^3]/M)/(E[W]-E[W^2]/M)  -- exact, no LP")
    print("=" * 96)
    print(f"  {'k':>3} {'bar':>10} {'PZ=B_2':>10} {'PZ margin':>10} {'D3':>10} {'D3 margin':>10} {'gain':>9}")
    for k in (8, 9, 10, 11, 12, 13):
        E = list(range(k))
        D3, PZ, rstar, valid = D3_exact(E)
        bar = BAR[k]
        print(f"  {k:>3} {float(bar):>10.6f} {float(PZ):>10.6f} {float(PZ-bar):>+10.6f} "
              f"{float(D3):>10.6f} {float(D3-bar):>+10.6f} {float(D3-PZ):>+9.6f}"
              f"  {'r* in (0,M]' if valid else 'r* OOR!'}")
    print()
    print("  NOTE: mac-mini's THM-661 clears k=8,9,10 with degree-4 (D3 for the block also clears")
    print("  k=10 at +0.101; k=8,9 need degree 4).  The DELTA here is the binding k=11 leg:")
    print("  degree-2 PZ margin +0.0159 (razor-thin) -> exact closed-form degree-3 +0.0735 (4.6x).")

    print()
    print("=" * 96)
    print("EXHAUSTIVE k=11: the block is the exact D3-minimizer over the compact regime")
    print("=" * 96)
    K = 11
    bar = BAR[K]
    omin = None
    oarg = None
    for D in range(K - 1, 15):
        mn = None
        arg = None
        cnt = 0
        for mid in combinations(range(1, D), K - 2):
            E = [0] + list(mid) + [D]
            if gcd_all([E[i + 1] - E[i] for i in range(len(E) - 1)]) != 1:
                continue
            cnt += 1
            D3, _, _, valid = D3_exact(E)
            if D3 is None:
                continue
            if mn is None or D3 < mn:
                mn = D3
                arg = E
        if mn is None:
            continue
        if omin is None or mn < omin:
            omin = mn
            oarg = arg
        print(f"  diam {D:3d}: {cnt:5d} shapes, min D3 = {float(mn):.6f} margin +{float(mn-bar):.5f}"
              f"  {'CLEARS' if mn >= bar else 'BELOW'}")
    print()
    print(f"  EXHAUSTIVE min D3 (k=11, diam {K-1}..14) = {omin}")
    print(f"    = {float(omin):.6f} at {oarg}, margin +{float(omin-bar):.5f} (vs PZ min 0.346788 +0.0156)")
    print(f"  => the binding k=11 leg clears the honest bar by +{float(omin-bar):.4f} at degree 3,")
    print(f"     block-minimized, exact rational -- the thin PZ margin is comfortably thickened.")


if __name__ == "__main__":
    main()
