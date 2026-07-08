"""
lrc14_pz_degree3_floor_opus_S148.py   (opus-2026-07-08-S148, HYP-5327)

THE DEGREE-3 COVERING FLOOR -- strengthening the thin k=11 Paley-Zygmund margin (+0.0159).

PZ = E[W]^2/E[W^2] is the OPTIMAL *2-moment* lower bound on P(W>0), W in [0,6/7] (THM-660).
A *3-moment* bound (adding E[W^3]) is strictly stronger.  The optimal n-moment lower bound
on P(X>0), X in [0,M], with moments m1..mn is the LP:
    P(X>0) >= max over deg-n poly p:  E[p(X)] = sum c_i m_i,  s.t. p(0) <= 0 and
              p(t) <= 1 for all t in (0, M].
(Any such p satisfies p(X) <= 1_{X>0} pointwise, so E[p] <= P(X>0).)  We solve it for
the covering W with EXACT E[W], E[W^2], E[W^3] (Farey-cell integration, W^3 -> quartic),
and compare the degree-3 floor to the degree-2 PZ at the block for k=11,12,13.

If degree-3 substantially thickens the k=11 margin, the binding leg is de-risked.
"""
from fractions import Fraction as F
from math import floor
import sys
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import BAR

TH = F(1, 7)
M = F(6, 7)   # W max


def moments_exact(E, upto=3):
    """Exact E[W^p], p=1..upto, for integer family E via Farey-cell integration."""
    E = sorted(E); k = len(E); e0 = E[0]
    Es = [e - e0 for e in E]
    dens = set(e for e in Es if e > 0)
    for i in range(k):
        for j in range(i + 1, k):
            dens.add(Es[j] - Es[i])
    bps = set([F(0), F(1)])
    for d in dens:
        for p in range(0, d + 1):
            bps.add(F(p, d))
    bps = sorted(bps)
    mom = [F(0)] * (upto + 1)
    for c in range(len(bps) - 1):
        x0, x1 = bps[c], bps[c + 1]; xm = (x0 + x1) / 2
        lin = [(F(-floor(e * xm)), F(e)) for e in Es]
        order = sorted(range(k), key=lambda j: lin[j][0] + lin[j][1] * xm)
        sp = [lin[j] for j in order]
        gaps = [(sp[i + 1][0] - sp[i][0], sp[i + 1][1] - sp[i][1]) for i in range(k - 1)]
        gaps.append((F(1) + sp[0][0] - sp[k - 1][0], sp[0][1] - sp[k - 1][1]))
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                xs = (TH - a) / b
                if x0 < xs < x1:
                    subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs) - 1):
            u0, u1 = subs[s], subs[s + 1]; um = (u0 + u1) / 2
            A = F(0); B = F(0)
            for (a, b) in gaps:
                if a + b * um > TH:
                    A += (a - TH); B += b
            # integrate (A+Bx)^p over [u0,u1] exactly, p=1..upto
            for p in range(1, upto + 1):
                # (A+Bx)^p ; expand and integrate
                acc = F(0)
                from math import comb
                for r in range(p + 1):
                    coef = comb(p, r) * A**(p - r) * B**r
                    # int x^r dx = (u1^{r+1}-u0^{r+1})/(r+1)
                    acc += coef * (u1**(r + 1) - u0**(r + 1)) / (r + 1)
                mom[p] += acc
    return mom  # mom[1], mom[2], mom[3]


def degree3_floor(m1, m2, m3, Mval, ngrid=4000):
    """LP: max c0 + c1 m1 + c2 m2 + c3 m3 s.t. p(0)=c0 <= 0 and p(t) <= 1 for t in (0,M].
       Returns the degree-3 lower bound on P(W>0). (Numeric; valid because constraints are
       enforced on a fine grid + endpoints; conservative refinement possible.)"""
    m1f, m2f, m3f, Mf = float(m1), float(m2), float(m3), float(Mval)
    # variables c0,c1,c2,c3; maximize c0 + c1 m1 + c2 m2 + c3 m3  => minimize negative
    cobj = [-1.0, -m1f, -m2f, -m3f]
    Aub = []; bub = []
    # p(0) = c0 <= 0
    Aub.append([1, 0, 0, 0]); bub.append(0.0)
    # p(t) <= 1 for t in (0,M] on a grid (include small t near 0 and endpoint M)
    ts = np.concatenate([np.linspace(1e-6, Mf, ngrid), [Mf]])
    for t in ts:
        Aub.append([1, t, t * t, t**3]); bub.append(1.0)
    # bound variables to keep LP bounded
    bounds = [(-1e4, 1e4)] * 4
    res = linprog(cobj, A_ub=np.array(Aub), b_ub=np.array(bub), bounds=bounds, method="highs")
    if not res.success:
        return None, None
    c = res.x
    val = c[0] + c[1] * m1f + c[2] * m2f + c[3] * m3f
    return val, c


def main():
    print("=" * 96)
    print("DEGREE-3 COVERING FLOOR vs degree-2 PZ (block k=11,12,13); W in [0,6/7]")
    print("=" * 96)
    for k in (11, 12, 13):
        E = list(range(k))
        m = moments_exact(E, 3)
        m1, m2, m3 = m[1], m[2], m[3]
        pz2 = m1 * m1 / m2   # exact degree-2
        f3, c = degree3_floor(m1, m2, m3, M)
        bar = float(BAR[k])
        print(f"  k={k} block:")
        print(f"    E[W]={float(m1):.6f}  E[W^2]={float(m2):.6f}  E[W^3]={float(m3):.6f}")
        print(f"    degree-2 PZ  = {float(pz2):.6f}   margin {float(pz2)-bar:+.6f}")
        print(f"    degree-3 LP  = {f3:.6f}   margin {f3-bar:+.6f}   "
              f"(improvement +{f3-float(pz2):.6f})   bar {bar:.6f}")
        print()
    print("  If degree-3 thickens k=11's margin substantially, the binding leg is de-risked")
    print("  (the exact degree-3 bound = a rational LP with the exact E[W^i]; the numeric LP")
    print("   here is the indicative value -- an exact 3-point Markov-Krein rep pins it).")


if __name__ == "__main__":
    main()
