#!/usr/bin/env python3
"""hn_rotation_field_law_s688.py — the rotation-field law for Eisenstein
unit-distance gadgets: a norm-N rhombus rotation lives in ℚ(√−(4N−1)); √−11
(N=3, the Moser spindle) is the SMALLEST field-escape from Eisenstein.

Follows up HYP-2276 (the Moser spindle lives in ℚ(√−3,√−11); χ≥4 needs a 2nd
imaginary-quadratic field). Here we make the "which field?" map precise.

SETUP. A spindle-type χ=4 gadget joins two Eisenstein lattice points that are each
at squared-distance N from a shared pivot, by rotating one rhombus about the pivot
through θ until the two far vertices are a UNIT distance apart. Two points at
distance √N from the pivot, separated by angle θ, are unit-distance iff
    2√N · sin(θ/2) = 1  ⟹  cos θ = 1 − 1/(2N) = (2N−1)/(2N).
The rotation e^{iθ} then satisfies z² − 2cos θ·z + 1 = 0, i.e.
    N·z² − (2N−1)·z + N = 0,
with discriminant  (2N−1)² − 4N² = −(4N−1).  So:

    THE ROTATION FIELD IS  ℚ(√−(4N−1)).                              (★)

It is a NEW field (≠ the Eisenstein field ℚ(√−3)) unless the squarefree part of
4N−1 is 3 — i.e. unless 4N−1 = 3·(odd square), N ∈ {1,7,19,37,61,...}, where the
rotation stays inside ℚ(√−3) and forces no escape.

N must be a representable EISENSTEIN norm: N = a²−ab+b² ∈ {1,3,4,7,9,12,13,...}.
The smallest representable N that ESCAPES (N ∉ {1,7,19,...}) is N=3 ⟹ √−11 ⟹ the
MOSER SPINDLE. (N=2 would give the smaller √−7, but 2 is not an Eisenstein norm.)

Session: claude-2026-06-06-S688 (hn-rotation-field-law)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
from math import gcd, isqrt

def squarefree_part(n):
    """squarefree part of positive integer n."""
    n = abs(n); sf = 1; d = 2
    while d*d <= n:
        cnt = 0
        while n % d == 0: n //= d; cnt += 1
        if cnt % 2 == 1: sf *= d
        d += 1
    if n > 1: sf *= n
    return sf

def eisenstein_norms(limit):
    """representable norms a²−ab+b² ≤ limit (a,b ∈ ℤ)."""
    S = set()
    R = isqrt(limit) + 2
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            v = a*a - a*b + b*b
            if 1 <= v <= limit: S.add(v)
    return sorted(S)

HEEGNER = {1, 2, 3, 7, 11, 19, 43, 67, 163}

print("=== the rotation-field law:  norm-N Eisenstein rhombus ⟹ rotation ∈ ℚ(√−(4N−1)) ===\n")
print(f"{'N':>4} {'Eis?':>5} {'4N−1':>6} {'sqfree':>7} {'rotation field':>16} {'escapes ℚ(√−3)?':>16} {'Heegner?':>9}")
norms = set(eisenstein_norms(60))
firstescape = None
rows = []
for N in range(1, 31):
    d = 4*N - 1
    sf = squarefree_part(d)
    is_eis = N in norms
    field = f"ℚ(√−{sf})"
    escapes = (sf != 3)
    heeg = (sf in HEEGNER)
    if is_eis and escapes and firstescape is None:
        firstescape = N
    if N <= 16 or (is_eis and escapes):
        rows.append((N, is_eis, d, sf, field, escapes, heeg))
for (N, is_eis, d, sf, field, escapes, heeg) in rows:
    mark = "  <-- MOSER" if (N == 3) else ""
    print(f"{N:>4} {('yes' if is_eis else 'no'):>5} {d:>6} {sf:>7} {field:>16} "
          f"{('YES' if escapes else 'no — stays Eis'):>16} {('yes' if heeg else 'no'):>9}{mark}")

print(f"\nsmallest representable Eisenstein norm that ESCAPES ℚ(√−3): N = {firstescape} "
      f"⟹ ℚ(√−{squarefree_part(4*firstescape-1)})  (= the Moser spindle, √−11)")
print(f"(N=2 would give the smaller √−7, but 2 ∉ Eisenstein norms {{1,3,4,7,9,12,13,...}})")

# the non-escaping norms: 4N-1 = 3·(odd square)
print("\nnon-escaping Eisenstein norms (rotation stays in ℚ(√−3), no new field):")
ne = [N for N in range(1, 200) if N in norms and squarefree_part(4*N-1) == 3]
print(f"  N ∈ {ne}   (4N−1 = 3·odd² : {[4*N-1 for N in ne]} = 3·{[ (4*N-1)//3 for N in ne]})")

# ---- verify the ℤ[ζ6] 3-coloring that caps lattice gadgets at χ≤3 ----
print("\n=== why χ≤3 inside ℤ[ζ6]: the explicit lattice 3-coloring c(a,b)=(a−b) mod 3 ===")
# units (norm-1 neighbours) of ℤ[ζ6] in (a,b) coords, point = a + b·ζ6:
units = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def col(a,b): return (a-b) % 3
ok = all(col(0,0) != col(u,v) for (u,v) in units)  # every neighbour differs in colour
diffs = sorted({(col(0,0)-col(u,v)) % 3 for (u,v) in units})
print(f"  the 6 unit neighbours = {units}")
print(f"  colour changes to neighbours (mod 3) = {diffs}  (0 absent ⟹ proper)  proper 3-colouring? {ok}")
print(f"  ⟹ the triangular lattice ℤ[ζ6] unit-distance graph is 3-colourable, so EVERY")
print(f"     unit-distance graph drawn inside it has χ ≤ 3. χ≥4 must leave the lattice — via (★).")

print("\n=== upshot ===")
print("  rotation field = ℚ(√−(4N−1)); the Moser spindle (N=3) is the minimal field-escape,")
print("  forced to √−11. The χ=3→χ=4 jump is exactly 'acquire a 2nd imaginary-quadratic field'.")
print("  Heegner −d appear early (N=3→−11, N=4→−15(no), N=9→−35(no)...); √−7,√−11,√−19,√−43,√−67")
print("  reachable as 4N−1 squarefree parts — a target list for the field-tower conjecture (HYP-2276).")
