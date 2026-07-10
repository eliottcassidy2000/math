from fractions import Fraction as F
import random
random.seed(96)
def lonely_at(S, tau):
    return all(14 * ((v * tau) % 1) >= 1 and 14 * (1 - (v * tau) % 1) >= 1 for v in S)
# CHARACTERIZATION: exists odd m: lonely at m/16  <=>  (no v = 0 mod 16) AND (odd speeds miss
# some unit +-class mod 16, classes {+-1,+-3,+-5,+-7})
viol = 0
for _ in range(300000):
    S = random.sample(range(1, 500), 13)
    has = any(lonely_at(S, F(m, 16)) for m in (1,3,5,7,9,11,13,15))
    no16 = all(v % 16 != 0 for v in S)
    ocl = {min(v % 16, 16 - v % 16) for v in S if v % 2 == 1}
    miss = len(ocl & {1,3,5,7}) < 4
    if has != (no16 and miss): viol += 1
print(f"depth-4 characterization: violations {viol} / 300k")
# even free pass: 2^a * odd mod 16 always in [2,14] for a=1,2,3
fp = all(min((2**a * u * m) % 16, 16 - (2**a * u * m) % 16) >= 2
         for a in (1,2,3) for u in range(1,16,2) for m in range(1,16,2))
print(f"even free-pass at depth 4 (all 2^a*odd*odd distances >= 2): {fp}")
# uniqueness of the sweet spot: depth 5 free pass FAILS
fp5 = all(min((2 * u * m) % 32, 32 - (2 * u * m) % 32) >= 3 for u in range(1,32,2) for m in range(1,32,2))
print(f"depth-5 free pass holds: {fp5}  (expected False -- 16 < 28 < 32 calibration)")
# DECIDE-FRACTION on covering sets
from math import gcd
from functools import reduce
def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
dec = tot = 0
blind_examples = []
while tot < 3000:
    S = sorted(random.sample(range(1, 300), 13))
    if not covering(S) or reduce(gcd, S) != 1: continue
    tot += 1
    no16 = all(v % 16 != 0 for v in S)
    ocl = {min(v % 16, 16 - v % 16) for v in S if v % 2 == 1}
    if no16 and len(ocl & {1,3,5,7}) < 4: dec += 1
    elif len(blind_examples) < 3: blind_examples.append(S)
print(f"covering sets DECIDED by depth 4 alone: {dec}/{tot} = {dec/tot:.1%}")
print("blind examples:", blind_examples[:2])
