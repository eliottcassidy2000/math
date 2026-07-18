#!/usr/bin/env python3
"""
IS THE SELF-SIMILAR RESONANCE A GENUINE RECURSION? (boxeph-2026-07-18-S97)

S96 found: the density route's forced peel w = d (far element) is maximally
resonant, and "w plays a Lonely Runner against the mode lattice tZ --
self-similarity one level up" (THM-886 V). Owner: test whether that
self-similarity is a GENUINE RECURSION (LRC(14) resonance reducible to a
LOWER-count LRC, inductable / closable by settled LRC(<=13)) or a WALL
(a 1-D scaling law terminating at a fixed base = the 7-section comb).

Two decisive tests:

 TEST A -- SCALING LAW vs TOWER.  Compute Error(t) = |S(t)|/t at the resonant
   peel for the family E_t = {1..6, t}, t growing.  A GENUINE recursion would
   descend into a SHRINKING sequence of LRC instances (no fixed limit / a tower).
   A WALL/scaling-law gives Error(t) -> a FIXED constant c* = |nu_hat(1)|, the
   amplitude of a t-INDEPENDENT 7-comb.  Also |S(at)|/t -> |nu_hat(a)| for a=1,2,3.

 TEST B -- IS THE BASE FIXED (one step) OR DOES IT RECURSE?  Compute the section
   measure nu_c = N_c / t (signed far-element endpoint density in section
   c = j mod 7).  If nu_c converges to a FIXED 7-vector independent of t, the
   "recursion" terminates in ONE step at a fixed finite base -- NOT a tower.
   Report nu_hat = DFT(nu); check the empty-section (LRC(6) gap) structure.

 TEST C -- GENERATOR COUNT of the mode comb.  Over n in [1,P), find the large
   |S(n)| set; is it (essentially) the single lattice tZ (=> 1 generator =>
   Dirichlet, NOT a multi-runner LRC), or a multi-generator lattice (=> a genuine
   LRC of rank = #generators)?  A genuine LRC recursion needs rank in [2,13].

Pure Python.  Reuses the S25 endpoint machinery (same R_s convention).
"""

from fractions import Fraction as Fr
from math import gcd, lcm, pi, sin, cos, sqrt
import cmath


def endpoints(E, s):
    """R_s endpoints (position Fr in [0,1), sign, owner). R_s = {x : the 7e-section
    occupancy has len==6 and s not in occ}. +1 entering, -1 leaving."""
    bps = sorted(set(Fr(k, 7 * e) for e in E for k in range(7 * e)) | {Fr(0), Fr(1)})
    pts = []
    prev_in = None
    first_in = None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if prev_in is None:
            first_in = cur
        else:
            if cur and not prev_in:
                pts.append([bps[i], +1])
            elif prev_in and not cur:
                pts.append([bps[i], -1])
        prev_in = cur
    if prev_in != first_in:
        pts.append([Fr(0), +1 if (first_in and not prev_in) else -1])
    out = []
    for p, sg in pts:
        own = [e for e in E if (p * 7 * e).denominator == 1]
        out.append((p, sg, min(own) if own else 0))
    return out


def S_at(pts, P, n):
    """S(n) = sum sign * e(n * pos / P), pos = int(p*P)."""
    return sum(sg * cmath.exp(2j * pi * n * int(p * P) / P) for p, sg, o in pts)


def worst_section(E):
    """return the section s with the largest |S(t)|^2 (the resonant worst case)."""
    P = 7 * lcm(*E)
    t = max(E)
    best = None
    for s in range(7):
        pts = endpoints(E, s)
        if not pts:
            continue
        v = abs(S_at(pts, P, t)) ** 2
        if best is None or v > best[1]:
            best = (s, v, pts, P)
    return best


def section_measure(E, s):
    """nu_c signed density of the FAR element's endpoints in class c = j mod 7,
    where a far-owned endpoint at p = j/(7t) has class j mod 7.  Returns N_c list
    (signed counts) so nu_c = N_c / t."""
    t = max(E)
    P = 7 * lcm(*E)
    pts = endpoints(E, s)
    N = [0] * 7
    for p, sg, o in pts:
        if o == t:                     # far-element-owned endpoint
            j = int(p * 7 * t)         # p = j/(7t)  ->  j
            N[j % 7] += sg
    return N


def dft7(N):
    return [sum(N[c] * cmath.exp(2j * pi * a * c / 7) for c in range(7)) for a in range(7)]


# ---------------------------------------------------------------------------
print("=" * 78)
print("TEST A + B -- scaling law vs tower; is the base a FIXED 7-comb?")
print("family E_t = {1..6, t};  Error(t) = |S(t)|/t  at the forced resonant peel")
print("=" * 78)
print(f"{'t':>6} {'s*':>3} {'M':>6} | {'|S(t)|/t':>9} {'|S(2t)|/t':>10} {'|S(3t)|/t':>10}"
      f" | {'nu_hat(1)':>9} {'nu_hat(2)':>9} {'nu_hat(3)':>9}")
prev = None
for t in [6, 12, 24, 60, 120, 240, 480, 960, 1920]:
    E = [1, 2, 3, 4, 5, 6, t]
    s, v, pts, P = worst_section(E)
    M = len(pts)
    e1 = abs(S_at(pts, P, t)) / t
    e2 = abs(S_at(pts, P, 2 * t)) / t
    e3 = abs(S_at(pts, P, 3 * t)) / t
    N = section_measure(E, s)
    nu = [x / t for x in N]
    nh = dft7(nu)
    print(f"{t:>6} {s:>3} {M:>6} | {e1:>9.4f} {e2:>10.4f} {e3:>10.4f} | "
          f"{abs(nh[1]):>9.4f} {abs(nh[2]):>9.4f} {abs(nh[3]):>9.4f}")

print()
print("Interpretation: if |S(at)|/t and |nu_hat(a)| CONVERGE to fixed constants")
print("(t-independent), the self-similarity is a SCALING LAW terminating at the")
print("fixed 7-comb nu_hat -- ONE step, NOT a tower.  Error -> c* = |nu_hat(1)| fixed.")
print("A genuine recursion would instead give a SHRINKING sequence with no fixed base.")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("TEST B' -- the fixed base nu: empty-section (LRC(6) gap) structure")
print("=" * 78)
t = 1920
E = [1, 2, 3, 4, 5, 6, t]
s, v, pts, P = worst_section(E)
N = section_measure(E, s)
print(f"t={t}, worst s*={s}:  signed section counts N_c =", N)
print(f"  nu_c = N_c/t =", [round(x / t, 4) for x in N])
print(f"  section s*={s} count = {N[s]}  (the LRC(6)-frame gap the far element sits in)")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("TEST C -- generator count of the mode comb: single lattice tZ (Dirichlet)")
print("          or multi-generator (genuine LRC of rank = #generators)?")
print("=" * 78)
for t in [24, 60, 120]:
    E = [1, 2, 3, 4, 5, 6, t]
    s, v, pts, P = worst_section(E)
    M = len(pts)
    # large-|S| positions over a full period (P manageable for these t)
    thr = 0.35 * sqrt(v)          # 35% of the peak amplitude |S(t)|
    comb = [n for n in range(1, P) if abs(S_at(pts, P, n)) > thr]
    on_t = [n for n in comb if n % t == 0]
    off_t = [n for n in comb if n % t != 0]
    g = 0
    for n in comb:
        g = gcd(g, n)
    print(f"t={t:>4} P={P:>6} M={M:>5} peak|S|={sqrt(v):.1f} thr={thr:.1f} | "
          f"comb size={len(comb):>4}  on tZ={len(on_t):>4}  off tZ={len(off_t):>4}  "
          f"gcd(comb)={g}")
    if off_t:
        # residues of off-lattice teeth mod t -- do they form a fixed sub-structure?
        res = sorted(set(n % t for n in off_t))
        print(f"        off-lattice residues mod t (first 12): {res[:12]}"
              f"{'...' if len(res) > 12 else ''}  (count {len(res)})")
print()
print("If off-tZ teeth are absent or a FIXED small residue set mod t (e.g. multiples")
print("of t/7), the comb = single lattice tZ modulated by the fixed 7 -> 1 nontrivial")
print("generator -> DIRICHLET, not a multi-runner LRC.  Rank >= 2 with growing/coprime")
print("generators would be needed for a genuine LRC recursion.")
