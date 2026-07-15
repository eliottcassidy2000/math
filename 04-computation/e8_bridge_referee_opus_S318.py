# opus-2026-07-15-S318 -- HYP-6945 / THM-868: THE E8 BRIDGE referee.
# (1) n=8: {d in Z^8 : all d odd, |d| <= 7, Sum d = 0}/2 subset E8 (D8+ test);
#     shell census: lattice points vs Landau-legal points per level x = Sum d^2;
#     x = 8 shell: near-regular = 70 of E8's 240 roots (the Sigma=0 half-coset).
# (2) Aut odd (=> solvable by Feit-Thompson) across all classes n = 4..7.
# (3) climb length = (n^3-n)/24 - [n even] n/8; increments = paired triangulars.
from math import comb
from fractions import Fraction
from collections import defaultdict
import sys, itertools
sys.path.insert(0, '04-computation')

# ---------- (1) the E8 bridge
print("(1) THE E8 BRIDGE at n = 8:")
def in_E8(v):   # v = tuple of 8 rationals; E8 = D8+ : all in Z with even sum,
    # or all in Z+1/2 with even sum
    from fractions import Fraction as F
    ints = all(x.denominator == 1 for x in v)
    halfs = all(x.denominator == 2 for x in v)
    if not (ints or halfs): return False
    return sum(v) % 2 == 0

n = 8
pts = []
for d in itertools.product(range(-7, 8, 2), repeat=7):
    last = -sum(d)
    if last % 2 == 1 and -7 <= last <= 7:
        pts.append(d + (last,))
print(f"   all-odd Sigma=0 vectors (|d|<=7): {len(pts)}")
all_in = all(in_E8(tuple(Fraction(x, 2) for x in v)) for v in pts)
print(f"   ALL half-vectors in E8 (D8+ membership): {all_in}")
shells_lat = defaultdict(int)
shells_srt = defaultdict(set)
for v in pts:
    x = sum(t*t for t in v)
    shells_lat[x] += 1
    shells_srt[x].add(tuple(sorted(v)))
def landau_ok(dvec):
    s = sorted((d + 7)//2 for d in dvec)
    return (sum(s) == comb(8, 2) and
            all(sum(s[:k+1]) >= comb(k+1, 2) for k in range(8)))
print(f"   x=8 shell: lattice points = {shells_lat[8]} (E8 roots in the "
      f"Sigma=0 half-coset: C(8,4) = {comb(8,4)}); E8 total roots = 240")
print("   per-level: x | lattice pts | sorted orbits | Landau-legal orbits")
for x in sorted(shells_lat)[:12]:
    legal = sum(1 for m in shells_srt[x] if landau_ok(m))
    print(f"      {x:4d} | {shells_lat[x]:6d} | {len(shells_srt[x]):4d} | {legal:4d}")

# ---------- (2) Aut odd
print("\n(2) Feit-Thompson face: all class |Aut| odd (n = 4..7):")
from smith_diagram_of_the_metagraph_opus_S307 import build
for nn in range(4, 8):
    B = build(nn)
    H_of, cls_of = B['H_of'], B['cls_of']
    fiber = defaultdict(int)
    for t in range(1 << (nn*(nn-1)//2 - (nn-1))): fiber[cls_of[t]] += 1
    auts = {c: H_of[c] // fiber[c] for c in fiber}
    allodd = all(a % 2 == 1 for a in auts.values())
    print(f"   n={nn}: classes={len(auts)}, all |Aut| odd: {allodd} "
          f"(=> all tournament symmetry groups SOLVABLE by Feit-Thompson)")

# ---------- (3) climb length
print("\n(3) climb length and paired-triangular increments:")
prev = None
for nn in range(3, 13):
    if nn % 2 == 1: climb = (nn**3 - nn)//24
    else: climb = (nn**3 - 4*nn)//24
    inc = climb - prev if prev is not None else None
    prev = climb
    print(f"   n={nn:2d}: climb={climb:3d}" + (f"  increment={inc}" if inc else ""))
print("   increments 1,3,3,6,6,10,10,15,15 = the triangulars, each twice")
