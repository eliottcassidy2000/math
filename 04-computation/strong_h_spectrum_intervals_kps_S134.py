#!/usr/bin/env python3
r"""
strong_h_spectrum_intervals_kps_S134.py
(kind-pasteur-2026-07-26-S134; companion to HYP-9029)

The strong-interval tiling law. Inputs: the canon exhaustive strong
H-spectra (monad s4/s5 outs for m<=8; the m=9 value dump regenerated
this session from monad's s6 iso-class generator, stored as
strong_H_spectrum_m9_values_kps_S134.out).  Computes the maximal
contiguous odd intervals and the junctions:

  strong(7) >= [65, 105]   (islands below: [25,61] with hole 63)
  strong(8) >= [69, 609]   (islands: 45, [49,53], [57,65])
  strong(9) >= [85, 2881]  (islands: 75, 81)

Junction law verified: c_{n+1} <= d_n at both junctions (69 <= 105,
85 <= 609), so the union covers EVERY odd in [65, 2881] by a strong
witness.  With the canon semigroup law (THM-1862/1936: H-spectrum =
multiplicative closure of strong values) and the machine-checked
base (products of m<=8 strong values give all odd <= 400 except
{7,21}), spectrum completeness (THM-1370's conjecture) REDUCES to
the tiling law: strong(n) contains [c_n, d_n] with c_{n+1} <= d_n
for all n >= 7.

Reproduction: python 04-computation/strong_h_spectrum_intervals_kps_S134.py
"""
import re

def fail(m):
    raise SystemExit("CHECK FAILED: " + m)

def intervals(v):
    out = []; start = prev = v[0]
    for x in v[1:]:
        if x == prev + 2:
            prev = x
        else:
            out.append((start, prev)); start = prev = x
    out.append((start, prev))
    return out

# m=7 from the exhaustive labeled run (also re-derived independently
# this session by metagraph engine cross-check)
v7 = [25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,65,67,
      69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,109,
      111,113,115,117,121,123,125,127,129,131,133,135,137,139,141,143,
      145,147,151,153,155,157,159,171,175,189]
s = open('05-knowledge/results/strong_H_spectrum_m8_isoclass_monad_s5.out',
         encoding='utf-8').read()
m8 = re.search(r'm=8 strong H-spectrum \(297 values\):\s*\[([^\]]+)\]', s)
v8 = sorted(int(x) for x in m8.group(1).split(','))
v9 = sorted(int(x) for x in
            open('05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out').read().split())
if len(v8) != 297 or len(v9) != 1482:
    fail("cardinalities")
big = {}
for name, v in (('7', v7), ('8', v8), ('9', v9)):
    iv = intervals(v)
    b = max(iv, key=lambda t: t[1] - t[0])
    big[name] = b
    print(f"strong({name}): |vals|={len(v)} giant run {b}; "
          f"islands below: {[t for t in iv if t[1] < b[0]]}")
if not (big['8'][0] <= big['7'][1] and big['9'][0] <= big['8'][1]):
    fail("junction law")
print(f"junctions: {big['8'][0]} <= {big['7'][1]} and "
      f"{big['9'][0]} <= {big['8'][1]}: PASS")
print(f"union covers every odd in [{big['7'][0]}, {big['9'][1]}]")
if big['7'][0] != 65 or big['9'][1] != 2881:
    fail("expected range")
print("ALL CHECKS PASSED")
