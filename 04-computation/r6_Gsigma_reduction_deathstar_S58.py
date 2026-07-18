from fractions import Fraction as F
# G(sigma) = largest gap in [0,1) minus union of 5 arcs width 1/7 centered at {0,s,2s,3s,4s} (mod 1)
def G(s):  # s a Fraction in [0,1)
    w=F(1,14)
    arcs=[]
    for m in range(5):
        c=(m*s)%1; lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: arcs.append((lo,hi))
        else: arcs.append((lo,F(1))); arcs.append((F(0),hi))
    arcs.sort()
    merged=[]
    for lo,hi in arcs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append((lo,hi))
    if not merged: return F(1)
    # circular gaps
    mg=F(0)
    for i in range(len(merged)):
        a=merged[i][1]; b=merged[(i+1)%len(merged)][0]+(1 if i+1==len(merged) else 0)
        if b-a>mg: mg=b-a
    return mg
# scan sigma to find where G(s) > 1/7, and the minimum of G, and value at s ~ 0.314
thr=F(1,7)
import numpy as np
N=20000
best_min=(F(1),None); vals=[]
for i in range(1,N):
    s=F(i,N)
    g=G(s); vals.append((s,g))
    if g<best_min[0]: best_min=(g,s)
print("min G(s) over grid =",float(best_min[0]),"at s=",float(best_min[1]),"  (1/7=%.5f)"%float(thr))
# where is G(s) > 1/7 ?  find bands
good=[float(s) for s,g in vals if g>thr]
if good:
    # merge into ranges
    good.sort(); bands=[]; a=good[0]; p=good[0]
    for x in good[1:]:
        if x-p>2.0/N+1e-9: bands.append((a,p)); a=x
        p=x
    bands.append((a,p))
    print("sigma-bands where G(s)>1/7:", [(round(a,4),round(b,4)) for a,b in bands])
    print("total measure of good sigma:", round(sum(b-a for a,b in bands),4))
else:
    print("NO sigma with G>1/7 (!)")
# value at the worst config's sigma
for sv in [F(2*15,98), F(2*9,56), F(157,500), F(157,500)]:
    print("  G(sigma=%.5f) = %.5f -> R_sharp~1/(7G)=%.4f"%(float(sv%1),float(G(sv%1)),1/(7*float(G(sv%1)))))
# exact worst config: t* is center of the max gap ~0.159; sigma=2t*
t_star=F(1585+1596,2*10000)  # approx from layout gap [0.158596,0.159592]
print("worst-config t* ~0.159, sigma=2t*~0.318: G=",float(G((2*t_star)%1)))
