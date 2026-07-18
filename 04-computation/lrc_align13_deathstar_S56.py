from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
random.seed(2)
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(m,q)
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
# Test: covering-2..12 12-set with an optimal denominator D_W divisible by 13 => AP?
print("Do any NON-AP covering-2..12 12-sets have an optimal loneliness denom divisible by 13?")
found_nonAP=[]; ap_denoms=None
base=list(range(1,13)); pool=[v for v in range(1,40) if v%13 and v%14]
def test(W):
    global ap_denoms
    W=sorted(set(W))
    if len(W)!=12 or not covers(W,range(2,13)): return
    M,args=M_arg(W); denoms=sorted(set(q for a,q in args))
    has13=any(q%13==0 for q in denoms)
    if is_AP(W):
        if ap_denoms is None: ap_denoms=denoms
        return
    if has13: found_nonAP.append((W,denoms,float(M)))
# AP reference
test(base)
print(f"  AP {{1..12}}: optimal denominators = {ap_denoms}  (13-divisible: {[q for q in ap_denoms if q%13==0]})")
# search non-AP
cnt=0
for k in (1,2):
    for pos in combinations(range(12),k):
        combos=[tuple(random.choice(pool) for _ in range(k)) for _ in range(400)] if k==2 else [(x,) for x in pool]
        for nv in combos:
            W=base[:]
            for p,x in zip(pos,nv): W[p]=x
            cnt+=1; test(W)
print(f"  searched ~{cnt} near-AP covering cores; NON-AP with 13|D_W: {len(found_nonAP)}")
for W,d,m in found_nonAP[:8]: print(f"    {W}: denoms={d} M={m:.4f}")
if not found_nonAP:
    print("  => CONJECTURE HOLDS on sample: only the AP has a 13-divisible optimal denominator.")
    print("     (=> only the AP's good set is 182-aligned; the far element covers only the AP.)")
