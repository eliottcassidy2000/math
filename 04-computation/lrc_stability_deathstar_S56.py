from fractions import Fraction as F
from math import gcd
import random
random.seed(9)
def M_exact(fam,Qcap=None):
    if Qcap is None: Qcap=2*max(fam)+2
    best=F(0)
    for q in range(2,Qcap+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best

# (1) j=7 SPREAD fast killers over a 6-core: is M(V)=M(core) (far, invisible)? => NOT a wall
print("(1) j=7 SPREAD fast killers over core -- M(V) = M(core)? (far => invisible, not a wall):")
for core,K in [(list(range(1,7)),[200,260,330,410,500,620,760]),
               (list(range(1,7)),[91,130,170,220,280,360,470])]:
    V=sorted(core+K); MV=M_exact(V,3*max(core)+4); MC=M_exact(core)
    print(f"   core{core}+7 killers: M(V)>={float(MV):.4f}  M(core)={float(MC):.4f}  1/14={1/14:.4f} => M(V)>1/14: {MV>F(1,14)}")

# (2) STABILITY LEMMA check: non-AP core V=W u {vmax}, M(V)<1/13 => vmax <= v2/(13*delta)
print("\n(2) Stability lemma: for M<1/13 families with NON-AP core, is vmax<=v2/(13 delta)?  (delta=M(W)-1/13)")
# construct non-AP-core M<1/13 families are RARE; instead verify the lemma's inequality direction on
# near-AP cores + comparable vmax by checking: if W non-AP (delta>0), M(V)<1/13 requires vmax large.
def M_lt_13(fam):
    Q=2*max(fam)+2
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 13*m>=q: return False
    return True
# take near-AP W (swap a couple), compute delta, then the lemma bound vmax<=v2/(13 delta); test if any comparable vmax gives M<1/13
tested=0; viol=0
for _ in range(3000):
    W=list(range(1,13)); 
    for _ in range(random.randint(1,2)):
        W[random.randrange(12)]=random.randint(1,30)
    W=sorted(set(W))
    if len(W)!=12: continue
    MW=M_exact(W); delta=MW-F(1,13)
    if delta<=0: continue   # W is AP-equivalent
    v2=max(W); bound=v2/(13*delta)   # lemma: any vmax giving M(V)<1/13 must be <= this
    # test a comparable vmax (<= bound and > v2): does it give M<1/13?
    for vmax in range(v2+1, min(int(bound)+2, v2+200)):
        V=sorted(W+[vmax])
        if len(set(V))!=13: continue
        if M_lt_13(V):
            tested+=1
            if vmax>float(bound)+1e-9: viol+=1; print(f"   LEMMA VIOLATED: W={W} vmax={vmax} bound={float(bound):.1f}")
            break
print(f"   near-AP families with M<1/13 for some comparable vmax: {tested}, lemma violations: {viol}")

# (3) Is there an M-GAP above 1/13 for non-AP 12-sets? (smallest delta over non-AP)
print("\n(3) M-gap: smallest delta=M(W)-1/13 over non-AP 12-sets (small => no gap => delicate):")
mind=F(1); arg=None
for _ in range(8000):
    W=sorted(random.sample(range(1,26),12))
    MW=M_exact(W); d=MW-F(1,13)
    if F(0)<d<mind: mind=d; arg=W
print(f"   smallest positive delta found: {float(mind):.5f} = {mind} at W={arg}  (1/13={1/13:.4f})")
