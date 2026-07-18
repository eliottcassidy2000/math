from fractions import Fraction as F
from math import gcd
# meas(G_W) = measure of {t: min_{w in W} ||wt|| > 1/13}, exact via danger-arc complement.
def good_measure(W, lvl=F(1,13)):
    # danger of w at level lvl: intervals [(k-lvl)/w,(k+lvl)/w]; good = complement measure
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-lvl)/w; hi=(F(k)+lvl)/w
            ivs.append((max(lo,F(0)),min(hi,F(1))))
    ivs=[(a,b) for a,b in ivs if b>a]; ivs.sort()
    merged=[]
    for a,b in ivs:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    danger=sum(b-a for a,b in merged)
    return F(1)-danger
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
print("Measure necessary condition: M(V)<1/13 => meas(G_W) <= 2/13 = %.4f"%(2/13))
print("AP {1..12}: meas(G_W) =", float(good_measure(list(range(1,13)))), "(tight => 0)")
# how many non-AP valid cores have meas(G_W) > 2/13 (eliminated by MEASURE alone, no equidistribution)?
from itertools import combinations; import random; random.seed(1)
base=list(range(1,13)); pool=[v for v in range(1,36) if v%13 and v%14]
big=0; small=0; ex=[]
tested=0
for k in (1,2):
    for pos in combinations(range(12),k):
        combos=[tuple(random.choice(pool) for _ in range(k)) for _ in range(200)] if k==2 else [(x,) for x in pool]
        for nv in combos:
            W=base[:]
            for p,x in zip(pos,nv): W[p]=x
            W=sorted(set(W))
            if len(W)!=12 or is_AP(W) or not covers(W,range(2,13)) or any(v%13==0 or v%14==0 for v in W): continue
            tested+=1
            g=good_measure(W)
            if g>F(2,13): big+=1
            else:
                small+=1
                if len(ex)<6: ex.append((W,float(g),float(M_exact(W))))
print(f"valid non-AP cores tested: {tested}")
print(f"  meas(G_W) > 2/13 (ELIMINATED by measure alone, PROVED): {big}")
print(f"  meas(G_W) <= 2/13 (near-tight; need the discrepancy/equidistribution): {small}")
print("  near-tight examples (W, meas(G_W), M(W)):")
for W,g,m in ex: print(f"     {W}: meas={g:.4f} M(W)={m:.4f}")
