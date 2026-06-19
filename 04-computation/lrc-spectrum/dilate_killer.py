"""
Dilation + killer constructions, exact M.

Idea that AVOIDS the dense-prefix collapse:
 - A "dilated AP" S0 = {d, 2d, ..., k d} has the SAME M as {1..k} = 1/(k+1) (scaling
   t by 1/d). So pure dilation keeps tight. To break tightness we must add gcd-1
   structure: replace ONE multiple by a nearby non-multiple (a "defect"), or add a
   killer at a coprime residue.

 - KILLER family: S = {1, 2, ..., k-1} is dense (collapses). Instead use a SPARSE AP
   with a defect: S = {d, 2d, ..., (k-1)d, X} where X breaks gcd to 1. Choose d, X so
   the optimum t stays near 1/(d(k)) -> M near 1/(?) ... explore.

 - Actually the cleanest near-floor non-tight known: take the AP {1,...,k} and SCALE
   the WHOLE problem: consider {1,...,k} but with one speed multiplied. Let's just
   exhaustively perturb the AP by adding +q to ONE element where q is a multiple of
   (k+1) (so at t=1/(k+1) that element's ||.|| is unchanged), making M still attainable
   near floor but the set non-tight at the OPTIMAL t (which shifts slightly).
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def ev(S, k):
    S = sorted(set(int(x) for x in S))
    if len(S) != k or min(S) <= 0 or setgcd(S) != 1:
        return None
    M, t = M_exact_fast(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, g=g, gk2=g * k * k, tight=(g == 0),
                a=(level_a(M, k) if g > 0 else None), maxS=max(S))


def fmt(r):
    if r is None:
        return "INVALID/tight"
    return (f"M={r['M']} (~{float(r['M']):.10f}) g*k^2={float(r['gk2']):.6f} "
            f"a={float(r['a']):.4f} maxS={r['maxS']}")


def perturb_AP(k, max_mult=8):
    """AP {1..k}, add j*(k+1) to one element (keeps that element's ||.|| at t=1/(k+1)).
    Several elements perturbed simultaneously also tested (single + pairs)."""
    base = list(range(1, k + 1))
    floor = Fraction(1, k + 1)
    best = None
    cand = []
    # single-element shift by j*(k+1)
    for idx in range(k):
        for j in range(1, max_mult + 1):
            S = base[:]
            S[idx] = base[idx] + j * (k + 1)
            r = ev(S, k)
            if r and not r['tight'] and r['g'] > 0:
                cand.append(r)
    cand.sort(key=lambda r: r['gk2'])
    return cand[:8]


if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1:]] or [12, 20, 30, 31, 40]
    for k in ks:
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} ====", flush=True)
        cl = perturb_AP(k)
        for r in cl[:5]:
            print("   ", r['S'])
            print("      ", fmt(r), flush=True)
        if not cl:
            print("   no non-tight perturbation", flush=True)
