# LRC(14) realizability as NON-COVERABILITY: meas(union U_s) < 1 (safe set nonempty).
# U_s = {t: ||s t|| < 1/14} = s arcs of width 1/(7s). Sum of measures = 13/7 ~ 1.857 > 1.
# So covering is MEASURE-possible; LRC must come from FORCED OVERLAPS (additive structure).
from fractions import Fraction as Fr
from itertools import combinations
def unsafe_arcs(s, thr=Fr(1,14)):  # arcs (lo,hi) where ||s t|| < thr, in [0,1)
    arcs=[]
    for k in range(s):
        lo=Fr(k,s)-thr/s; hi=Fr(k,s)+thr/s
        arcs.append((lo,hi))
    return arcs
def meas_union(S, thr=Fr(1,14)):
    iv=[]
    for s in S:
        for lo,hi in unsafe_arcs(s,thr):
            iv.append((lo%1,hi))  # may wrap; normalize below
    # normalize wrapping: split intervals crossing 0/1
    norm=[]
    for lo,hi in iv:
        lo=lo%1
        hi_un=lo+(2*thr/ (hi-lo) * (hi-lo))  # keep width
    # redo cleanly: width = 2*thr/s
    iv=[]
    for s in S:
        w=Fr(2,14)/s  # = 1/(7s)
        for k in range(s):
            c=Fr(k,s); lo=(c-w/2)%1; hi=lo+w
            if hi<=1: iv.append((lo,hi))
            else: iv.append((lo,1)); iv.append((Fr(0),hi-1))
    iv.sort()
    tot=Fr(0); cur_lo,cur_hi=None,None
    for lo,hi in iv:
        if cur_hi is None: cur_lo,cur_hi=lo,hi
        elif lo<=cur_hi: cur_hi=max(cur_hi,hi)
        else: tot+=cur_hi-cur_lo; cur_lo,cur_hi=lo,hi
    if cur_hi is not None: tot+=cur_hi-cur_lo
    return tot
def safe_meas(S): return 1-meas_union(S)
def add_energy(S):
    from collections import Counter; c=Counter()
    for a in S:
        for b in S: c[a+b]+=1
    return sum(v*v for v in c.values())
sets={
 'consec {1..13}': list(range(1,14)),
 'AP d=2 {2,4..26}': [2*i for i in range(1,14)],
 'AP d=3': [3*i for i in range(1,14)],
 'Sidon-ish': [1,2,5,11,22,33,40,55,69,78,89,95,102][:13],
 'random A': [1,4,7,9,13,17,23,29,31,37,41,43,47],
 'powers-2-ish':[1,2,4,8,16,32,64,128,256,512,1024,2048,4096],
}
print(f"{'set':24} {'meas(safe)':>12} {'meas(union)':>12} {'A(E)':>7}")
for name,S in sets.items():
    sm=safe_meas(S); print(f"{name:24} {float(sm):>12.5f} {float(1-sm):>12.5f} {add_energy(S):>7}")
print("\nmeas(safe)>0 == LRC holds for that set. The MARGIN = how non-covering. Sum meas(U_s)=13/7=%.4f."%(13/7))
print("Realizability question: is meas(safe)>0 FORCED for all 13-sets? (the forced-overlap obstruction)")
