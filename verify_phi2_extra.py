"""
(1) Push generic pairs further to confirm I_B -> Phi_2.
(2) Independent cross-check of p0 via an exact rational fine grid that does NOT
    use the breakpoint method (it samples deterministic rationals and refines),
    to make sure the exact-breakpoint p0 is right.
(3) Average I_B over a band of generic pairs (Cesaro-style) -> should track Phi_2.
"""

from fractions import Fraction as F
from math import floor, gcd

SECTORS = 7


def measure_by_misscount(E):
    pts = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        denom = SECTORS * e
        for k in range(0, denom + 1):
            pts.add(F(k, denom))
    pts = sorted(pts)
    out = {}
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        if b == a:
            continue
        mid = (a + b) / 2
        hit = set()
        for e in E:
            t = e * mid
            t = t - floor(t)
            hit.add(floor(t * SECTORS))
        miss = SECTORS - len(hit)
        out[miss] = out.get(miss, F(0)) + (b - a)
    return out


def p0(E):
    return measure_by_misscount(E).get(0, F(0))


def I_B(B, u, v):
    Bt = tuple(B)
    return p0(Bt + (u, v)) - p0(Bt + (u,)) - p0(Bt + (v,)) + p0(Bt)


def Phi2(B):
    d = measure_by_misscount(B)
    return (2 * d.get(2, F(0)) - d.get(1, F(0))) / 49


# ---- Independent p0 via uniform deterministic grid of N points (no breakpoints)
def p0_grid_lowerbound(E, N):
    """Riemann lower estimate: sample midpoints of N equal cells.  Converges to p0
    but is NOT exact; used only as a sanity range-check vs the exact value."""
    cnt = 0
    for i in range(N):
        x = F(2 * i + 1, 2 * N)
        hit = set()
        for e in E:
            t = e * x
            t = t - floor(t)
            hit.add(floor(t * SECTORS))
        if len(hit) == SECTORS:
            cnt += 1
    return F(cnt, N)


if __name__ == "__main__":
    print("CROSS-CHECK: exact breakpoint p0 vs grid sampling")
    for E in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 2, 4, 8), (0, 1, 2, 4, 8, 144, 233)]:
        exact = p0(E)
        g = p0_grid_lowerbound(E, 20000)
        print(f"  E={E}: exact={float(exact):.8f}  grid20000={float(g):.8f}  "
              f"|diff|={abs(float(exact)-float(g)):.2e}")

    print()
    bigfib = [28657, 46368, 75025, 121393, 196418, 317811]
    pairs = [(bigfib[i], bigfib[i + 1]) for i in range(len(bigfib) - 1)]
    pairs = [(u, v) for (u, v) in pairs if gcd(u, v) == 1 and u % 7 and v % 7]
    for B in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 2, 4, 8)]:
        phi = Phi2(B)
        print(f"\nB={B}  Phi_2={float(phi):.12f}  (deeper Fibonacci pairs)")
        for (u, v) in pairs:
            diff = I_B(B, u, v) - phi
            print(f"   (u,v)=({u},{v})  I_B-Phi_2 = {float(diff):+.3e}")

    # Cesaro average over a band of generic coprime pairs near a fixed large scale
    print("\nCESARO AVERAGE of I_B over a band of coprime pairs (scale ~ 1000-1200):")
    for B in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 2, 4, 8)]:
        phi = Phi2(B)
        tot = F(0)
        ct = 0
        for u in range(1000, 1030):
            for v in range(1130, 1160):
                if gcd(u, v) == 1 and u % 7 and v % 7 and gcd(u, 7) == 1:
                    tot += I_B(B, u, v)
                    ct += 1
        avg = tot / ct
        print(f"   B={B}: avg I_B over {ct} pairs = {float(avg):.10f}  "
              f"Phi_2={float(phi):.10f}  diff={float(avg-phi):+.2e}")
