"""
Independent verification of the "two-far curvature" decorrelation claim.

Setup (my own reading of the spec):
  - Sectors: [j/7, (j+1)/7) for j=0..6 on the circle [0,1).
  - For a finite set E of nonnegative integers and x in [0,1):
        residues(x) = { floor(7 * frac(e*x)) mod 7 : e in E }   (the sectors hit)
  - coverage indicator C_E(x) = 1 iff residues(x) == {0,1,...,6} (all 7 sectors hit).
  - p0(E) = Lebesgue measure { x in [0,1) : C_E(x) = 1 }.

  - I_B(u,v) = p0(B u {u,v}) - p0(B u {u}) - p0(B u {v}) + p0(B)
  - p1(B) = meas{ x : the set {frac(b*x): b in B} misses EXACTLY 1 of the 7 sectors }
  - p2(B) = meas{ x : ...                          misses EXACTLY 2 of the 7 sectors }
  - Phi_2(B) = (2*p2(B) - p1(B)) / 49

We compute EVERYTHING exactly with fractions.Fraction.

Breakpoints: the sector of frac(e*x) changes exactly when e*x crosses k/7 (mod 1),
i.e. at x = m/(7e) for integers m.  So the union over e in E of { m/(7e) : 0 <= m <= 7e }
gives a partition of [0,1) into intervals on which every e's sector is constant; we
evaluate the indicator at the midpoint of each interval and weight by its length.
"""

from fractions import Fraction as F
from math import floor, gcd


def sector_of(e, x):
    """Return which 7-sector frac(e*x) lands in, as integer 0..6. Exact."""
    t = e * x                      # exact Fraction
    t = t - floor(t)               # frac, exact
    s = floor(t * 7)               # 0..6 (t<1 so t*7<7)
    return s


def breakpoints(E):
    """All x in [0,1) where some e in E changes sector: x = m/(7e)."""
    pts = set()
    pts.add(F(0))
    pts.add(F(1))
    for e in E:
        if e == 0:
            continue  # frac(0*x)=0 always in sector 0, never changes
        # k/(7e) for k = 1 .. 7e-1, restricted to (0,1)
        denom = 7 * e
        for k in range(0, denom + 1):
            x = F(k, denom)
            if 0 <= x <= 1:
                pts.add(x)
    return sorted(pts)


def sectors_hit(E, x):
    """Set of sectors covered by E at point x."""
    return frozenset(sector_of(e, x) for e in E)


def measure_by_misscount(E):
    """
    Return dict: miss_count -> total measure of x in [0,1) where the set
    {frac(e*x): e in E} misses exactly that many of the 7 sectors.
    """
    pts = breakpoints(E)
    out = {}
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        if b == a:
            continue
        mid = (a + b) / 2
        hit = sectors_hit(E, mid)
        miss = 7 - len(hit)
        out[miss] = out.get(miss, F(0)) + (b - a)
    return out


def p0(E):
    """Measure of x where E hits ALL 7 sectors (miss count == 0)."""
    return measure_by_misscount(E)[0] if 0 in measure_by_misscount(E) else F(0)


# small optimization: compute miss-distribution once
def p0_fast(E):
    d = measure_by_misscount(E)
    return d.get(0, F(0))


def I_B(B, u, v):
    Bu = tuple(B) + (u,)
    Bv = tuple(B) + (v,)
    Buv = tuple(B) + (u, v)
    return p0_fast(Buv) - p0_fast(Bu) - p0_fast(Bv) + p0_fast(B)


def Phi2(B):
    d = measure_by_misscount(B)
    p1 = d.get(1, F(0))
    p2 = d.get(2, F(0))
    return (2 * p2 - p1) / 49, p1, p2


def show(label, val):
    print(f"  {label} = {val}   (~ {float(val):.10f})")


if __name__ == "__main__":
    print("=" * 70)
    print("PART (a): closed-form Phi_2(B) = (2 p2 - p1)/49")
    print("=" * 70)

    for B, expect in [
        ((0, 1, 2, 3, 4, 5, 6, 7), F(47, 24010)),
        ((0, 1, 2, 4, 8), F(1, 98)),
    ]:
        print(f"\nB = {B}")
        phi, p1, p2 = Phi2(B)
        show("p1(B)", p1)
        show("p2(B)", p2)
        show("Phi_2(B) = (2 p2 - p1)/49", phi)
        show("EXPECTED", expect)
        print(f"  MATCH: {phi == expect}")

    print()
    print("=" * 70)
    print("PART (b): I_B(u,v) -> Phi_2(B) for large coprime (u,v)")
    print("=" * 70)

    for B in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 2, 4, 8)]:
        phi, p1, p2 = Phi2(B)
        print(f"\nB = {B}")
        show("Phi_2(B)", phi)
        for (u, v) in [(101, 211), (211, 401)]:
            g = gcd(u, v)
            iv = I_B(B, u, v)
            diff = iv - phi
            print(f"  (u,v)=({u},{v}) gcd={g}:")
            show("    I_B(u,v)", iv)
            show("    I_B - Phi_2", diff)
