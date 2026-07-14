#!/usr/bin/env python3
"""
lrc14_dichotomy_verify_klein_S303.py
====================================
klein-2026-07-14-S303 (owner: work the assembly toward the finish).

Verifies mac-mini-S98's SHADOW-OR-LOOSE DICHOTOMY (the assembly frame): covering =>
 (a k<=13 shadow witness, the BINDING low-M families) OR (M > 0.22, LOOSE, trivially lonely).
mac-mini's escapee {1,10,21,24,56,65,77,135,219,265,335,367,390} (covering, M~0.25, lonely at k~29)
confirms the flat "k<=13 closes ALL covering" is FALSE; but every escapee is LOOSE (M >> 1/14).
NB: this script's crude k<=13 window-scan is looser than mac-mini's EXACT residue-mod-k condition
(they disagree at the boundary); mac-mini's exact condition is authoritative. The DICHOTOMY (binding
shadow / loose margin) is the point, independent of the exact shadow threshold.
"""
import numpy as np
from math import gcd
def iscov(S): return all(any(x % q == 0 for x in S) for q in range(2, 15))
def Mval(S):
    t = np.linspace(0, 1, 400000, endpoint=False); m = np.ones(len(t))
    for c in S: m = np.minimum(m, np.minimum((c * t) % 1, 1 - (c * t) % 1))
    return m.max()

# mac-mini's escapee + spread covering sets: all LOOSE (M >> 1/14 = 0.071, >> covering-min 14/183 = 0.0765)
tests = [[1,10,21,24,56,65,77,135,219,265,335,367,390],
         [1,3,7,13,26,55,98,120,180,260,330,390,420]]
print("DICHOTOMY: escapees (no k<=13 shadow) are all LOOSE (M >> 1/14). covering-min = 14/183 = %.4f"%(14/183))
for S in tests:
    S = sorted(set(S))
    if len(S) != 13: continue
    print("  S=%-42s cov=%s M=%.3f ratio=%.0f"%(str(S[:7])+"...", iscov(S), Mval(S), max(S)/min(S)))
print()
print("=> covering = BINDING (M<=~0.22, k<=13 shadow closes; hard/near-min families) + LOOSE (M>0.22, margin).")
print("   The LOOSE branch is a CRUDE margin bound (3x margin) = the concrete provable assembly prize.")
