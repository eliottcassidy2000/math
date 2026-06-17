import sys
sys.path.insert(0,'/tmp')
from lrc_base import *

# COVERING ENGINEERING: take a near-tight 12-speed base, then choose the 13th speed
# (a "plugger") to cover the surviving safe points. A large speed v has bands at k/v
# spaced 1/v apart, width 1/(7v); if 1/v < width can't help, but many bands can hit
# all 3 safe points 1/14,3/14,5/14 simultaneously only if v*{1,3,5}/14 all near integer+-<1/14.
# i.e. need v*1/14, v*3/14, v*5/14 all within 1/14 of an integer. 
# v*1 ≡ small mod 14, v*3, v*5 too. v≡0 mod 14 does it but then ||v tau||=0 trivially...
# but v≡0 mod14 means v multiple of 14 -> at tau=k/14, v*k/14 = integer*k -> norm 0 <1/14. 
# So a speed v=14 would KILL all three points 1/14,3/14,5/14! Let's test base 1..13 minus one + 14.
base12=list(range(1,13))  # 1..12
for plug in [13,14,28,42]:
    S=base12+[plug]
    if len(set(S))<13: continue
    m,at=M(S)
    print(f"1..12 + {plug}: M={m}={float(m):.6f} at {at} prim={primitive(S)}")

# The catch: adding 14 kills 1/14,3/14,5/14 but OPENS new safe points elsewhere
# (14 has bands only at k/14, leaving 1/13-type gaps). Let's see where M moves.
S=list(range(1,13))+[14]
m,at=M(S)
print(f"\n1..12,14: M={m}={float(m):.6f} at tau={at}")
print("  norms at that tau:", sorted(set(nrm(v*at) for v in S)))
# Try 1..13 + replace 13 by 14, plus try to also fix with the freed slot... only 13 speeds.
# Try base = 1..11 + {14, X} to kill 1/14,3/14,5/14 via 14 AND cover the new gap via X.
print("\nbase 1..11 + 14 + X (X plugs the gap 14 opens):")
best=(F(1),None)
b11=list(range(1,12))
for X in range(12,400):
    if X in b11+[14]: continue
    S=b11+[14,X]
    if len(set(S))<13: continue
    m,_=M(S)
    if m<best[0]: best=(m,tuple(sorted(S)))
print("  best:", best[0], float(best[0]), best[1], "prim=",primitive(list(best[1])) if best[1] else None)
