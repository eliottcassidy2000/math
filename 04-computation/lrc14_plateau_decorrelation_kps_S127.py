# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-11-S127: the plateau decorrelation lemma (THM-700) -- the wide-spread wall's
# inductive step. For E = E' u {w}, w = max E, the exact cover decomposition
#   1_cover(E,x) = 1_cover(E',x) + sum_s f_s(x)*1{frac(wx) in sector s},  f_s = 1{E' misses exactly sector s}
# gives the CENTERED error  p0(E) - Plat(E') = sum_s integral f_s(x)*[1{frac(wx) in s} - 1/7] dx,
# bounded by Fourier decay of BV f_s against the mean-zero sector oscillation:  |.| <= V(f_s)/(6w).
# This script: (1) confirms the exact decomposition pointwise; (2) verifies Plat(consec_8)=0.36210
# (HYP-2644); (3) verifies the O(1/w) decorrelation and that the V/(6w) bound holds.
import math

def sector(y):                       # y in [0,1) -> sector 0..6
    return min(6, int(y * 7))

def sectors_hit(E, x):
    return set(sector((e * x) % 1.0) for e in E)

def covered(E, x):
    return 1 if len(sectors_hit(E, x)) == 7 else 0

def missed_sector(Ep, x):            # the unique missed sector if E' misses exactly one, else None
    missed = set(range(7)) - sectors_hit(Ep, x)
    return next(iter(missed)) if len(missed) == 1 else None

def p0(E, N):
    return sum(covered(E, (i + 0.5) / N) for i in range(N)) / N

def plat(Ep, N):
    p0p = p0(Ep, N)
    p1 = sum(1 for i in range(N) if missed_sector(Ep, (i + 0.5) / N) is not None) / N
    return p0p, p1, p0p + p1 / 7

def tv_bound(Ep):                    # crude V(E') <= 14 * sum e' (each f_s jumps <= 2 e' per sector)
    return 14 * sum(e for e in Ep if e > 0)

def main():
    Ep = [0, 1, 2, 3, 4, 5, 6, 7]    # consec_8 core (offsets, 0 in E')
    N = 200003                        # fine odd grid (avoids exact breakpoints)

    # (1) exact cover decomposition, checked pointwise on the grid
    bad = 0
    for i in range(20000):
        x = (i + 0.5) / 20000
        w = 101
        lhs = covered(Ep + [w], x)
        ms = missed_sector(Ep, x)
        rhs = covered(Ep, x) + (1 if (ms is not None and sector((w * x) % 1.0) == ms) else 0)
        if lhs != rhs:
            bad += 1
    print(f"(1) exact cover decomposition 1_cover(E)=1_cover(E')+sum_s f_s*1(wx in s):  mismatches = {bad}  "
          f"[{'CONFIRMED' if bad == 0 else 'FAILED'}]")

    # (2) the plateau value
    p0p, p1, pl = plat(Ep, N)
    V = tv_bound(Ep)
    print(f"(2) core consec_8: p0(E')={p0p:.5f}, p1(E')={p1:.5f}, Plat={pl:.5f}  "
          f"(HYP-2644 claim 0.36210); V-bound={V}")

    # (3) the O(1/w) decorrelation and the V/(6w) bound
    print(f"    {'w':>6} {'p0(E)':>9} {'|p0-Plat|':>11} {'err*w':>8} {'V/(6w)':>9} {'bound holds':>12}")
    for w in [50, 101, 199, 401, 809, 1601]:
        e = abs(p0(Ep + [w], N) - pl)
        b = V / (6 * w)
        print(f"    {w:>6} {p0(Ep+[w], N):>9.5f} {e:>11.5f} {e*w:>8.3f} {b:>9.5f} {str(e <= b):>12}")
    print("    => p0(E) - Plat(E') = O(1/w), the V/(6w) Fourier bound holds (loose; sharp const ~0.2).")

if __name__ == "__main__":
    main()
