#!/usr/bin/env python3
"""
resonant_crossing_chord_spectrum_monad.py
monad-explorer-2026-06-13  (OPEN-Q-057 / HYP-2461 / THM-493 / THM-494 follow-up)

GOAL: pin WHICH Moser norm t produces the n=28 unit-distance 3N-crossing, and decide
whether t=3 (sqrt(-11)) is arithmetically UNIQUE or merely FIRST.

THESIS (to verify, exact integer arithmetic only):
  For a resonant product P = G  (+ omega_t)  H  (THM-493),
      U(P) = e(G)|H| + |G|e(H) + Delta_t,
      Delta_t(G,H) = (1/2) sum_{N(alpha)=t} m_alpha(G) m_alpha(H),
  the bonus is nonzero ONLY IF  t  is a shared CHORD NORM of BOTH factors:
      t in ChordSpec(G) and ChordSpec(H),   ChordSpec(X)={N(x-x'): x,x' in X, x!=x'}.
  Hence the crossing t is confined to the (small) factor's chord spectrum.
  At n=28 the only dense lattice factorization is R (+) W7 (4 x 7); ChordSpec(R)={1,3},
  so t in {1,3}; t=1 is degenerate => t=3 is FORCED, unique.

All arithmetic is exact in Eisenstein integers a+b*zeta6, zeta6^2 = zeta6 - 1,
norm N(a,b) = a^2 + a*b + b^2.
"""

from itertools import combinations
from collections import defaultdict

# ---------- exact Eisenstein integer arithmetic ----------
def esub(p, q):           # p - q
    return (p[0]-q[0], p[1]-q[1])
def enorm(p):             # N(a+b zeta6) = a^2+ab+b^2
    a, b = p
    return a*a + a*b + b*b

def chord_norms(pts):
    """multiset of nonzero pairwise difference norms (unordered)"""
    out = defaultdict(int)
    for p, q in combinations(pts, 2):
        out[enorm(esub(p, q))] += 1
    return dict(out)

def edges(pts):
    return sum(1 for p, q in combinations(pts, 2) if enorm(esub(p, q)) == 1)

def m_alpha(pts):
    """ordered displacement census: alpha -> #{(p,q): p-q = alpha}"""
    out = defaultdict(int)
    for p in pts:
        for q in pts:
            if p != q:
                out[esub(p, q)] += 1
    return out

def norm_t_alphas(t, R=6):
    """all Eisenstein alpha with N(alpha)=t inside a box (R big enough for small t)"""
    res = []
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            if a*a + a*b + b*b == t:
                res.append((a, b))
    return res

def delta_t(G, H, t):
    """Delta_t(G,H) = 1/2 sum_{N(alpha)=t} m_alpha(G) m_alpha(H)  (exact, returns int)"""
    mG, mH = m_alpha(G), m_alpha(H)
    s = 0
    for alpha in norm_t_alphas(t):
        s += mG.get(alpha, 0) * mH.get(alpha, 0)
    assert s % 2 == 0, "Delta_t must be an integer"
    return s // 2

def resonant_total(G, H, t):
    return edges(G)*len(H) + len(G)*edges(H) + delta_t(G, H, t)

# ---------- the canonical small factors (triangular-lattice realizations) ----------
W7 = [(0,0),(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]   # hub + 6 units; e=12
R  = [(0,0),(1,0),(0,1),(1,1)]                          # rhombus K4-e = {0,1,z,1+z}; e=5

# Loeschian numbers (norms representable as a^2+ab+b^2), used for "is t admissible"
def is_loeschian(t):
    return len(norm_t_alphas(t)) > 0

print("="*72)
print("PART A — verify the W7 (+)_3 R = u(28)>=85 resonant product, exactly")
print("="*72)
print(f"  e(W7)={edges(W7)}, e(R)={edges(R)}, |W7|={len(W7)}, |R|={len(R)}")
print(f"  ChordSpec(W7) = {chord_norms(W7)}")
print(f"  ChordSpec(R)  = {chord_norms(R)}")
csW, csR = set(chord_norms(W7)), set(chord_norms(R))
shared = sorted(csW & csR)
print(f"  shared chord norms ChordSpec(W7) & ChordSpec(R) = {shared}")
nontriv_shared = [t for t in shared if t >= 2]
print(f"  shared NON-unit (t>=2) norms = {nontriv_shared}   <-- the admissible resonant t")
print()
P28 = edges(W7)*len(R) + len(W7)*edges(R)
print(f"  generic product cap P(28) = e(W7)*|R| + |W7|*e(R) = {edges(W7)}*4 + 7*{edges(R)} = {P28}")
print(f"  3N = 84")
print()
for t in [3,4,7,9,12,13]:
    d = delta_t(R, W7, t)
    tot = P28 + d
    tag = ""
    if t == 3: tag = "  <-- Moser sqrt(-11)"
    cross = ">3N CROSS" if tot > 84 else ("=3N" if tot==84 else "<3N")
    print(f"   t={t:2d}: loeschian={str(is_loeschian(t)):5s}  Delta_t={d}  U={tot:3d}  {cross}{tag}")

print()
print("  ===> uniqueness check: Delta_t(R,W7) is NONZERO only for t in", end=" ")
nz = [t for t in range(2, 60) if delta_t(R, W7, t) > 0]
print(nz, " (scanned t=2..59)")
assert nz == [3], "EXPECTED t=3 unique!"
print("  CONFIRMED: t=3 is the UNIQUE resonant norm for the R x W7 factorization.")
print("  Reason: ChordSpec(R)={1,3}; R's only non-edge chord is sqrt3 (norm 3).")
print(f"  And Delta_3 = {delta_t(R,W7,3)} = the exact crossing 83 < 84 < 85.")
