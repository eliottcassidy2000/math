# Validate the exact-tiler crux across small n (where LRC is PROVEN, n<=7).
# LRC(n): n-1 runners, gap 1/n. Unsafe U_s={||s t||<1/n} = s arcs width 2/(n s). meas(safe)=1-meas(union).
# CLAIM (S39/S41 cat-1 finish): the ONLY primitive exact-tilers (meas safe=0) are {1,...,n-1}.
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce
def meas_safe(S, n):
    iv=[]
    for s in S:
        w=Fr(2,n)/s
        for k in range(s):
            c=Fr(k,s); lo=(c-w/2)%1; hi=lo+w
            if hi<=1: iv.append((lo,hi))
            else: iv.append((lo,Fr(1))); iv.append((Fr(0),hi-1))
    iv.sort(); tot=Fr(0); clo=chi=None
    for lo,hi in iv:
        if chi is None: clo,chi=lo,hi
        elif lo<=chi: chi=max(chi,hi)
        else: tot+=chi-clo; clo,chi=lo,hi
    if chi is not None: tot+=chi-clo
    return 1-tot
for n in range(3,9):
    k=n-1; B=3*k  # search primitive (n-1)-sets, max speed <= B
    exact=[]; tot_prim=0
    for S in combinations(range(1,B+1), k):
        if reduce(gcd,S)!=1: continue
        tot_prim+=1
        if meas_safe(S,n)==0: exact.append(S)
    consec=tuple(range(1,n))
    others=[S for S in exact if S!=consec]
    print(f"n={n} ({k} runners, gap 1/{n}, B={B}): primitive sets={tot_prim}; exact-tilers={len(exact)}; "
          f"is consec there={consec in exact}; OTHERS={others[:4]}{'...' if len(others)>4 else ''}")
