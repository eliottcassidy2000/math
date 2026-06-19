#!/usr/bin/env python3
"""
ANGLE H part 3 — HARD HUNT for meas(S7(E)) > cap_k, focusing on the
DANGEROUS k=8 (margin cap_8 - meas_S7(consec_8) = 0.3815-0.3272 = 0.054, the
smallest).  If ANY integer E with 0 in E, |E|=8 has meas(S7(E)) > cap_8 =
2243/5880, the seven-sector reduction (THM-532/HYP-2603) has a hole.

Strategy:
  (a) Exhaustive larger box {0..B}, B up to 20 (8-subsets containing 0).
  (b) Targeted: APs (all give meas_S7 = consec by scale-invariance), perforated
      near-APs, multi-block, dilated, "fill all 7 sectors" greedy adversary.
  (c) Random large-spread integer E, many samples.
  (d) The theoretical MAX of meas(S7) over ALL real (not just integer) phase
      configs — if even the real sup is < cap_k we have margin; if the real
      sup exceeds cap but integers don't reach it, covering may matter.
"""
from fractions import Fraction as F
from itertools import combinations
import random, sys
sys.path.insert(0, "04-computation")
from lrc14_adversarial_chain_macmini_0618sH import meas_S7, M7

cap8 = F(2243,5880)
print(f"cap_8 = {cap8} = {float(cap8):.6f}")
print(f"M7(8) main term = {M7(8)} = {float(M7(8)):.6f}")
print(f"meas_S7(consec_8) = {meas_S7(list(range(8)))} = {float(meas_S7(list(range(8)))):.6f}")
print()

best = meas_S7(list(range(8))); argbest = tuple(range(8))
total_checked = 0

# (a) exhaustive box up to B
for B in [14, 16, 18, 20]:
    cnt = 0; exceed = 0; locbest = F(0); locarg=None
    for rest in combinations(range(1, B+1), 7):
        E = (0,)+rest
        v = meas_S7(list(E))
        cnt += 1
        if v > locbest:
            locbest = v; locarg = E
        if v > cap8:
            exceed += 1
        if v > best:
            best = v; argbest = E
    total_checked += cnt
    print(f"  box B={B}: checked {cnt:>7}  max_S7={float(locbest):.5f} at {locarg}"
          f"  (#>cap8 = {exceed})")

# (b) targeted structured families
print("\n[targeted families]")
fams = {
 "perforated {0..9}\\{1,7}": [0,2,3,4,5,6,8,9],
 "perforated {0..9}\\{3,6}": [0,1,2,4,5,7,8,9],
 "two-block 0-3 / 7-10":     [0,1,2,3,7,8,9,10],
 "two-block 0-3 / 20-23":    [0,1,2,3,20,21,22,23],
 "Sidon-ish":                [0,1,3,7,12,20,30,44],
 "geometric-ish":            [0,1,2,4,8,16,32,64],
 "spread-7 perfect":         [0,1,2,3,4,5,6,7],
 "dense-mod7 attempt":       [0,1,2,3,4,5,6,49],  # extra at 0 mod7
}
for name, E in fams.items():
    v = meas_S7(E)
    flag = "  <<< EXCEEDS cap8!" if v > cap8 else ""
    print(f"  {name:28} meas_S7={float(v):.5f}{flag}")
    if v > best: best=v; argbest=tuple(E)

# (c) random large-spread integer E
print("\n[random large-spread, 30000 samples up to spread 200]")
exceed_rand = 0; randbest=F(0); randarg=None
for _ in range(30000):
    spread = random.randint(10, 200)
    rest = random.sample(range(1, spread+1), 7)
    E = [0]+sorted(rest)
    v = meas_S7(E)
    if v > randbest: randbest=v; randarg=tuple(E)
    if v > cap8:
        exceed_rand += 1
        print(f"   EXCEEDS: E={E}  meas_S7={float(v):.5f}")
    if v > best: best=v; argbest=tuple(E)
print(f"  random max_S7 = {float(randbest):.5f} at {randarg}  (#>cap8={exceed_rand})")

print()
print("="*60)
print(f"GLOBAL max meas_S7 over ALL searched k=8 E: {float(best):.6f} = {best}")
print(f"  at E = {argbest}")
print(f"  cap_8 = {float(cap8):.6f};  margin = {float(cap8-best):.6f}")
print(f"  CONCLUSION: {'HOLE FOUND (S7>cap)' if best>cap8 else 'no integer E beats cap_8 — reduction holds in search'}")
print("="*60)
