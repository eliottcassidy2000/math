"""
Reconcile the upstream G_P measure constants and thr_k = 1 - min_{|P|=13-k} meas(G_P).
G_P = {x in [0,1): ||p x|| >= 1/14 for all p in P}, P subset {1..13}.
meas(G_P) computed exactly: for each p, the BAD set is {x: ||px||<1/14} = union of arcs
around k/p of half-width (1/14)/p. G_P is the complement intersection.
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor

def measGP(P):
    # complement of union over p of {||px||<1/14}
    # ||px||<1/14  <=>  exists integer k with |p x - k| < 1/14 <=> x in ((k-1/14)/p,(k+1/14)/p)
    bad=[]
    for p in P:
        p=int(p)
        for k in range(0,p+1):
            lo=(F(k)-F(1,14))/p; hi=(F(k)+F(1,14))/p
            lo=max(lo,F(0)); hi=min(hi,F(1))
            if lo<hi: bad.append((lo,hi))
    bad.sort(); tot=F(0); cur=cb=None
    for lo,hi in bad:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return F(1)-tot

# Full small-part P = {1..13}
full=list(range(1,14))
print("meas(G_{1..13}) =", measGP(full), float(measGP(full)))
# 10-element P from the prompt
P10=[1,2,3,5,7,8,9,11,12,13]
print("meas(G_P), P=", P10, "=", measGP(P10), float(measGP(P10)))
print("  is 14249/252252:", measGP(P10)==F(14249,252252))

# thr_k = 1 - min over |P|=13-k of meas(G_P)
print("\n=== thr_k = 1 - min_{|P|=13-k} meas(G_P), P subset {1..13} ===")
expected={8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7),13:F(0)}
for k in range(8,14):
    sz=13-k
    if sz==0:
        mn=measGP([])  # empty P => G_P = [0,1), meas 1
        argmin=()
    else:
        mn=None;argmin=None
        for P in combinations(range(1,14),sz):
            v=measGP(list(P))
            if mn is None or v<mn: mn=v; argmin=P
    thr=1-mn
    print(f"k={k} (|P|={sz}): min meas(G_P)={mn}={float(mn):.5f} at {argmin}  thr={thr}={float(thr):.5f}  matches:{thr==expected[k]}")
