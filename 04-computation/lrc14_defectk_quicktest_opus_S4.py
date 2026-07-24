import itertools
from fractions import Fraction as Fr
h=Fr(3,41)
def Lmax(V):
    segs=[]
    for v in V:
        for m in range(0,v+1):
            lo=(Fr(m)-h)/v; hi=(Fr(m)+h)/v
            if lo<0: segs.append((Fr(0),hi)); segs.append((lo+1,Fr(1)))
            elif hi>1: segs.append((lo,Fr(1))); segs.append((Fr(0),hi-1))
            else: segs.append((lo,hi))
    segs=[s for s in segs if s[1]>s[0]]; segs.sort()
    mg=[]
    for a,b in segs:
        if mg and a<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b))
        else: mg.append((a,b))
    best=Fr(0)
    for i in range(len(mg)):
        e=mg[i][1]; nxt=mg[i+1][0] if i+1<len(mg) else mg[0][0]+1
        if nxt-e>best: best=nxt-e
    return best
print("IMMEDIATE defect-k closure test (no scan):")
print("  klein lemma: sum_i 1/s_i >= B_k(C) = L_max(C)*(1-2kh)/(2h).")
print("  All far speeds >= 14 => sum_i 1/s_i <= k/14.")
print("  So if  min_C B_k(C) > k/14  then defect k is IMPOSSIBLE outright.\n")
for k in [2,3,4,5,6]:
    F=(1-2*k*h)/(2*h)
    if F<=0: print(f"  k={k}: lemma invalid (1-2kh<=0)"); continue
    worst=None
    for drop in itertools.combinations(range(1,14),k):
        C=[v for v in range(1,14) if v not in drop]
        B=Lmax(C)*F
        if worst is None or B<worst[0]: worst=(B,drop)
    B,drop=worst
    thr=Fr(k,14)
    ok = B>thr
    print(f"  k={k}: factor={F} min_C B_k={float(B):.4f} (worst drop {drop})  vs k/14={float(thr):.4f}"
          f"   -> {'CLOSED OUTRIGHT' if ok else 'not closed by this test'}")
