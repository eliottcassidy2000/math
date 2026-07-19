#!/usr/bin/env python3
"""
boxeph-2026-07-19-S124 — the Kakeya-needle obstruction to M=3/38 (the depth-minimal gap value).

Owner: work the 3/38 residue system for achievability; think Kakeya needles.

3/38 is the unique numerator-3 value in the open gap (1/13, 2/25) (boxeph-S123). A family at
M=3/38 must (opus-S117 + S123): be COVERING, have a determinant-3 maximizing pair at s=38, and
place all 12 residues mod 38 in the safe band [3,35] with 3/38 the GLOBAL deepest hole.

KAKEYA LENS: each modulus q is a "needle direction"; the witness t=p/q probes a family's deepest
hole in that direction. M=3/38 requires blocking EVERY needle to depth <= 3/38, with the q=38
needle at exactly 3/38 the deepest. This script shows the obstruction is a NEEDLE-COVERING:

 (A) band-filled-mod-38 covering families are LOOSE -- the q=38 hole (3/38) is real but never the
     deepest; a deeper hole sits at a MEDIUM modulus q' in (13,38), giving M >= 2/25.
 (B) the mod-19 parity split (38=2*19): at t=m/38 an even speed 2u has ||2u m/38|| = ||u m/19||,
     a multiple of 1/19=2/38, so a band-satisfying even speed sits at >= 2/19 = 4/38 > 3/38 --
     the 3/38 hole is carried ENTIRELY by ODD speeds; and the mod-19 witness must be blocked to
     <= 1/19 (else a 2/19 hole beats 3/38).
 (C) NO single medium needle is universal: a family can BLOCK the mod-19 needle by containing 19
     (then w19=0), but such families are caught by a DIFFERENT medium needle. The UNION of
     medium-modulus needles covers every band-filled covering family -- you can evade any one, not
     all. (Adaptive needles, cf. codex-S16 "adaptive graphic rank".)

HONEST: not a proof. 3/38 is verified unachievable for all bases in [1,26] (kps-S12); the
unbounded-modulus escape tail (macmini-S36, HYP-4667) is untouched. This MAPS the obstruction.
"""
from fractions import Fraction as F
from math import gcd
import random


def fd(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(sp, t: F) -> F:
    return min(fd(F(v) * t) for v in sp)


def trueM_arg(sp):
    Ds = set()
    n = len(sp)
    for i in range(n):
        for j in range(i + 1, n):
            Ds.add(sp[i] + sp[j])
            Ds.add(abs(sp[i] - sp[j]))
    best, bt = F(0), F(0)
    for D in sorted(Ds):
        if D <= 0:
            continue
        for m in range(1, D):
            v = fmin(sp, F(m, D))
            if v > best:
                best, bt = v, F(m, D)
    return best, bt


def covering(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))


def mod_witness(sp, q):
    best = F(0)
    for n in range(1, q):
        best = max(best, fmin(sp, F(n, q)))
    return best


TH = F(3, 38)
print("=" * 74)
print("(A) band-filled-mod-38 covering families are LOOSE: q=38 hole dominated by")
print("    a deeper MEDIUM-modulus hole q' in (13,38).")
print("=" * 74)
cands = [
    [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35],
    [3, 5, 7, 8, 9, 10, 11, 12, 17, 21, 24, 35],
    [3, 7, 8, 9, 10, 11, 12, 15, 20, 22, 33, 35],
    [3, 5, 7, 9, 10, 11, 12, 14, 16, 24, 33, 35],
]
for C in cands:
    C = sorted(set(C))
    m38 = fmin(C, F(1, 38))
    M, t = trueM_arg(C)
    print(f"  C={C}")
    print(f"     hole@t=1/38 = {m38}=3/38 ; trueM = {M} ({float(M):.4f}) at q'={t.denominator} "
          f">2/25:{M >= F(2,25)}")

print()
print("=" * 74)
print("(B) mod-19 parity split (38 = 2*19): even speeds sit at >= 2/19 at the q=38 hole")
print("=" * 74)
print("  even v=2u: ||2u/38|| = ||u/19|| in {0, 1/19, 2/19, ...}; band (>=3/38) forces >= 2/19=4/38.")
for u in range(1, 7):
    v = 2 * u
    print(f"    v={v:>2}: ||v/38|| = {fd(F(v,38))}")
print("  => the 3/38 hole is carried only by ODD speeds; the mod-19 needle must be blocked to <=1/19.")

print()
print("=" * 74)
print("(C) NEEDLE-COVERING: no single medium needle is universal, but the union is")
print("=" * 74)
random.seed(124)
pool = list(range(3, 36))
tested = 0
beat_any_medium = 0
w19_beats = 0
has19_blocks = 0
for _ in range(8000):
    rest = random.sample([x for x in pool if x not in (3, 35)], 10)
    C = sorted({3, 35} | set(rest))
    if len(C) != 12 or not covering(C):
        continue
    # band-filled at m=1 by construction (subset of [3,35], 3 & 35 present)
    tested += 1
    M, t = trueM_arg(C)
    # some medium modulus 14..37 gives a deeper hole than 3/38?
    deep_medium = any(mod_witness(C, q) > TH for q in range(14, 38))
    if deep_medium:
        beat_any_medium += 1
    if mod_witness(C, 19) > TH:
        w19_beats += 1
    if 19 in C or 38 in C:
        has19_blocks += 1
print(f"  band-filled covering families tested: {tested}")
print(f"  with SOME medium needle q' in [14,37] deeper than 3/38 (=> M>3/38): {beat_any_medium}/{tested}")
print(f"  where the mod-19 needle alone beats 3/38: {w19_beats}/{tested}")
print(f"  that contain 19 (blocking the mod-19 needle): {has19_blocks}/{tested}  "
      f"-- these evade mod-19 but are caught by another needle")
print(f"  => union of medium needles covers ALL {tested} families: "
      f"{beat_any_medium == tested}")
print("\nDONE.")
