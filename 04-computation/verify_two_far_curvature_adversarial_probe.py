import sys
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation')
from verify_two_far_curvature_independent_v2 import I_B, Phi_2
from fractions import Fraction as Fr
from math import gcd

B=(0,1,2,3,4,5,6,7); ph2=Phi_2(B)
def rd(u,v,lim):
    return min(abs(m*u+n*v) for m in range(-lim,lim+1) for n in range(-lim,lim+1) if (m,n)!=(0,0))

# Deterministic adversarial family: hide the relation just past the box.
# v/u ~ p/q with q in {8,9,10,11} (>7) and small u-grid. product(box7) blows up.
print("Adversarial v/u = p/q, q just > 7  (box<=7 misses relation):")
print(" u      v      relation         resdist7  resdistTrue  |dev|       prod7    prodTrue")
worst7=Fr(0)
for (p,q) in [(11,9),(13,9),(15,11),(13,8),(17,11)]:
    for u in [2000,4000,6000]:
        v=u*p//q
        if gcd(u,v)!=1: 
            v+=1
            if gcd(u,v)!=1: continue
        dev=abs(I_B(B,u,v)-ph2)
        r7=rd(u,v,7); rt=rd(u,v,15)
        p7=dev*r7; pt=dev*rt
        if p7>worst7: worst7=p7
        rel = f"{q}v-{p}u={q*v-p*u}"
        print(f"{u:6d} {v:6d}  {rel:15s}  {r7:7d}  {rt:9d}  {float(dev):.7f}  {float(p7):.5f}  {float(pt):.5f}")
print()
print("MAX prod7 in this family =", round(float(worst7),5), " (claimed bound ~0.01)")
