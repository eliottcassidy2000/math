# Wide adversarial sweep: confirm consecutive block is the UNIQUE argmax (up to dilation)
# over a large window at k=7, and spot-check k=8.
import itertools
from fractions import Fraction as F
from math import comb
def measS7(E):
    E=sorted(set(int(e) for e in E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len({int(((e*xm)%1)*7) for e in E})==7: total+=x1-x0
    return total
def M7(k): return sum(F((-1)**t*comb(6,t))*F(7-t,7)**(k-1) for t in range(7))
def corr(E): return measS7(E)-M7(len(E))
def is_block(E):
    d=sorted(E); return all(d[i+1]-d[i]==1 for i in range(len(d)-1))
# k=7 full window [0,13]
k=7; N=13
best=None; argmaxes=[]
nonblock_best=None
for combo in itertools.combinations(range(N+1),k):
    c=corr(list(combo))
    if best is None or c>best:
        best=c; argmaxes=[combo]
    elif c==best:
        argmaxes.append(combo)
    if not is_block(combo):
        if nonblock_best is None or c>nonblock_best[1]:
            nonblock_best=(combo,float(c))
print(f"k={k} window[0,{N}]: max corr={float(best):.6f}")
print(f"  #argmaxes={len(argmaxes)}; ALL consecutive blocks? {all(is_block(a) for a in argmaxes)}")
print(f"  argmaxes: {argmaxes}")
print(f"  best NON-block: {nonblock_best[0]} corr={nonblock_best[1]:.6f} deficit={float(best)-nonblock_best[1]:+.6f}")
# Is every argmax a single-step consecutive block (i.e. a translate of [0..k-1] or [1..k])?
# Note dilations leave [0,N] window if step>1 and span>N, so only step-1 blocks fit; that's expected.
