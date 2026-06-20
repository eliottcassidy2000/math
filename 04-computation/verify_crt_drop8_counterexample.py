#!/usr/bin/env python3
"""
Pin down the drop-8 counterexample to the 'all evens beat all odds' claim,
and check the halving-ratio claims precisely.
"""
from fractions import Fraction as F

def lonely(E):
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        for k in range(0,v+1):
            c=F(k,v); lo=c-w; hi=c+w
            if hi<=0 or lo>=1: continue
            bands.append((max(lo,F(0)),min(hi,F(1))))
    bands.sort(); merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]:
            if e>merged[-1][1]: merged[-1]=(merged[-1][0],e)
        else: merged.append((s,e))
    return F(1)-sum(e-s for s,e in merged)

# drop-8 is the deepest dyadic tower member (8=2^3). Claim says deleting it
# leaves consec minus 8; the cascade says +level3{8} takes L 0.045->0.
# So drop-8 L = exactly the level-2 value 0.04519290.
print("drop-8 L =", lonely([v for v in range(1,14) if v!=8]), "=", float(lonely([v for v in range(1,14) if v!=8])))
print("level-2 value (odd+lvl1+lvl2) =", lonely([1,3,5,7,9,11,13,2,6,10,4,12]),
      "=", float(lonely([1,3,5,7,9,11,13,2,6,10,4,12])))
print("These match BECAUSE drop-8 core = consec minus {8} = {odd}U{lvl1}U{lvl2}. Tautology.")
print()
# So drop-8 (an EVEN deletion) is NOT a good deletion: deleting the deepest tower
# element is BAD. The claim 'all 5 even deletions beat all 7 odd' contradicts this.
print("Odd deletions strictly better than drop-8:")
for e in [3,5,13]:
    print(f"  drop-{e}: L={float(lonely([v for v in range(1,14) if v!=e])):.8f} < drop-8 L=0.04519290")

print()
print("LEVEL-1 RATIO CLAIM: '+level1 ratio 0.490 ~ 1/2'")
L0=lonely([1,3,5,7,9,11,13])
L1=lonely([1,3,5,7,9,11,13,2,6,10])
print(f"  L0={L0}  L1={L1}  ratio={L1/L0}={float(L1/L0):.6f}")
print(f"  exactly 1/2? {L1/L0==F(1,2)}  -> it is {float(L1/L0):.4f}, NOT 1/2. The '~1/2' is approximate only.")
print(f"  level-2 ratio: {float(lonely([1,3,5,7,9,11,13,2,6,10,4,12])/L1):.4f}  (0.385, NOT 1/2)")
print(f"  level-3 ratio: 0  (8 closes everything)")
print("  => the halving cascade is NOT a clean factor of 1/2 per level. Only level1 is near 1/2.")
