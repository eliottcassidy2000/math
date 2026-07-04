#!/usr/bin/env python3
"""
klein-2026-07-03-S122 (HYP-4071) - IS THE KERSHNER (2D HEXAGONAL) BRIDGE A REDUCTION OR AN ANALOGY?

opus's standing candidate for NEW leverage on the covering-min lower bound M(covering) >= n/Phi6(n):
the covering-min witness is the zeta_6 hexagonal rotation, Phi6(n)=N(n-omega) (Eisenstein norm),
so MAYBE Kershner's theorem (hexagonal lattice = thinnest covering of the plane) proves the bound.

opus's own honest barrier: M(covering) >= n/Phi6(n) IS LRC(n) (it exceeds 1/n). So the question is
whether Kershner REDUCES it (a geometric shortcut) or is only a structural ANALOGY.

DECISIVE TEST: Kershner bounds a 2D EUCLIDEAN covering radius (the Eisenstein NORM metric). LRC
optimizes a 1D DIOPHANTINE min-distance (the residue/phase metric ||r/Phi6||). If these two metrics
on Z/Phi6(n) = Z[omega]/(n-omega) induce DIFFERENT orderings, Kershner's 2D optimality cannot bound
the 1D LRC quantity -- the bridge is an analogy, not a reduction. This script tests that, exactly.

Also: the exact decomposition n/Phi6(n) = (1/(n-1))*(1 - 1/Phi6(n)), and a Kershner numeric check.
"""
from fractions import Fraction as F
from math import gcd, sqrt
import itertools

def phi6(n): return n*n - n + 1

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def eisen_norm(x, y):
    # N(x + y*omega), omega = primitive 6th root (omega^2 = omega - 1): x^2 - x y + y^2
    return x*x - x*y + y*y

def short_eisen(r, q, n, B=40):
    """minimal Eisenstein norm over lifts x+y*omega with x + n*y == r (mod q)."""
    best = None
    for y in range(-B, B+1):
        x0 = (r - n*y) % q
        for x in (x0, x0 - q):
            nm = eisen_norm(x, y)
            if best is None or nm < best: best = nm
    return best

print("="*78)
print("PART 1 (exact): the covering-min decomposition")
print("  n/Phi6(n) = 1/(n-1) * (1 - 1/Phi6(n))  -- the (n-1)-point value, arithmetic-discounted")
print("="*78)
for n in [4,7,14,20]:
    q = phi6(n)
    lhs = F(n, q)
    rhs = F(1, n-1) * (1 - F(1, q))
    print(f"  n={n:>3}: n/Phi6={lhs} ({float(lhs):.6f}); 1/(n-1)*(1-1/Phi6)={rhs} ({float(rhs):.6f}); "
          f"equal? {lhs==rhs}; 1/n={float(F(1,n)):.6f}, 1/(n-1)={float(F(1,n-1)):.6f}")
print("  => covering-min sits in (1/n, 1/(n-1)), a 1D circle scale, NOT a 2D-density scale.")

print()
print("="*78)
print("PART 2 (THE DECISIVE TEST): 1D phase-distance vs 2D Eisenstein-norm on Z/Phi6(14)=Z/183")
print("  Kershner lives on the Eisenstein NORM (2D). LRC lives on the phase DISTANCE (1D).")
print("  If the orderings disagree, Kershner's 2D optimum does NOT bound LRC's 1D min.")
print("="*78)
n = 14; q = phi6(n)
rows = []
for r in range(1, q):
    d = cdist_q(r, q)          # 1D phase distance (numerator; actual = d/q)
    nm = short_eisen(r, q, n)  # 2D min Eisenstein norm of the lift
    rows.append((r, d, nm))
# Are the two metrics monotonically related?  Check: do small-norm residues = small-distance?
# The Eisenstein UNITS (norm 1) are the 6 closest 2D points; where do they land in 1D distance?
units = [(r,d,nm) for (r,d,nm) in rows if nm == 1]
print(f"  Eisenstein UNITS (2D norm 1, the 6 closest 2D vectors) and their 1D phase distances d/q:")
for r,d,nm in units:
    print(f"    r={r:>4}  2D-norm={nm}  1D-dist={d}/{q} ({float(F(d,q)):.4f})")
print(f"  -> the 6 nearest 2D points have 1D distances {sorted(set(d for _,d,_ in units))}/{q} "
      f"(NOT a single value): the 2D and 1D metrics DISAGREE.")
# quantify: rank correlation between 1D-dist and 2D-norm
import statistics
ds = [d for _,d,_ in rows]; nms = [nm for _,_,nm in rows]
# Spearman-ish: fraction of pairs where the two metrics ORDER oppositely (inversions)
N = len(rows); samp = rows[::1]
inv = 0; tot = 0
for i in range(0, len(samp), 3):
    for j in range(i+1, len(samp), 7):
        a,b = samp[i], samp[j]
        tot += 1
        if (a[1]-b[1])*(a[2]-b[2]) < 0: inv += 1   # opposite order on the two metrics
print(f"  metric-disagreement rate (sampled pairs ordered OPPOSITELY by 1D-dist vs 2D-norm): "
      f"{inv}/{tot} = {inv/tot:.2%}")
print("  A genuine reduction needs the SAME ordering (0% inversions). High inversion => Kershner's")
print("  2D covering radius is a DIFFERENT functional from LRC's 1D min-distance.")

print()
print("="*78)
print("PART 3: does the deep-well phase set = the 2D-hexagonally-optimal spread, or a 1D-AP spread?")
print("  deep-well phases at t*=14/183 (from S119): compare to (a) 1D equal-spacing, (b) 2D-nearest")
print("="*78)
DW = list(range(1,13)) + [182]
phases = sorted(cdist_q(v*n, q) for v in DW)   # 1D distances-from-0 of the 13 runner phases
print(f"  deep-well 1D distances-from-0 (sorted): {phases} (/{q})")
print(f"  min = {min(phases)}/{q} = M*Phi6 = {min(phases)} (= n = {n}); these are an AP-like 1D comb")
# 1D equal spacing of 13 points: distances would be ~q/26.. ; the point is min-dist ~ q/(2*13)
print(f"  1D equal-spacing of 13 points on the circle: min-dist-from-0 ~ 1/(2*13)*q ={q/26:.1f}/{q}; ")
print(f"  the covering-min config instead pins min at n={n} (arithmetic), NOT the 2D hexagonal radius.")

print()
print("="*78)
print("PART 4: Kershner numeric -- the 2D hexagonal covering radius does NOT reproduce n/Phi6(n)")
print("="*78)
# Kershner: hexagonal lattice A2 covering density = 2*pi/sqrt(27) ~ 1.209 (dimensionless, a DENSITY).
# It is not a length that yields n/Phi6(n). The covering-min is a rational Farey value n/Phi6(n),
# an ARITHMETIC (Eisenstein-Farey) distance, not the transcendental Kershner density.
kersh = 2*3.141592653589793/sqrt(27)
print(f"  Kershner hexagonal covering density = 2pi/sqrt(27) = {kersh:.6f} (transcendental, a density)")
print(f"  covering-min = n/Phi6(n) = {float(F(n,q)):.6f} = the RATIONAL Eisenstein-Farey value 14/183")
print(f"  No dimensional/numeric identity links them: Kershner gives a density, LRC a Farey length.")

print()
print("VERDICT (to be read with the reflection): the zeta_6/hexagonal structure is the SYMMETRY GROUP")
print("of the extremal witness (real, provable: n order 6 mod Phi6, Phi6=N(n-omega)). But the covering-")
print("min is a 1D Diophantine max-min on the residue metric, NOT a 2D Euclidean covering radius on the")
print("Eisenstein-norm metric. The two metrics disagree (Part 2), the value is a rational Farey length")
print("not a Kershner density (Part 4). So Kershner's 2D optimality does NOT reduce the covering-min:")
print("the bridge is a STRUCTURAL ANALOGY (right group, wrong theorem), not a proof route. The lower")
print("bound stays LRC-equivalent (opus's honest barrier).")
print("DONE")
