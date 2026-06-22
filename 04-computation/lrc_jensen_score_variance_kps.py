"""
The Jensen/score-variance route to consec-extremality (kind-pasteur S31m).

The right tournament template (kps-S31l) is the H-MAXIMIZER proof (THM-027): the regular
tournament maximizes H by JENSEN on scores -- c_3 = C(k,3) - sum_i C(s_i,2), minimized
score-variance.  LRC analogue (HYP-2605): the winding tournament T(x) has scores s_i(x);
coverage tracks the cyclic content c_3(T(x)).  TEST: does the AP/consec MINIMIZE the
integrated score functional  Phi(E) = int_0^1 sum_i C(s_i(x),2) dx  (=> maximize c_3 =>
maximize coverage), the convex/Jensen extremality that keeps the sign (Angle F)?

Also answers mac-mini-S39: is p0 SCALING-invariant (singles out consec-multiples d*{1..13})
or TRANSLATION-invariant (general APs)?  [Score variance lives on the PERMUTOHEDRON of S_k =
the truncated octahedron for k=4 -- the convex arena where the AP is the balanced center.]
"""
import math
from itertools import combinations
from collections import Counter

def winding_scores(e, x):
    k=len(e); s=[0]*k
    for i in range(k):
        for j in range(k):
            if i==j: continue
            if 0 < ((e[i]-e[j])*x) % 1.0 < 0.5: s[i]+=1
    return s

def Phi(e, N=4000):   # integrated sum_i C(s_i,2)
    tot=0.0
    for n in range(1,N):
        s=winding_scores(e, n/N)
        tot += sum(si*(si-1)//2 for si in s)
    return tot/(N-1)

def p0(E):
    E=sorted(set(ee for ee in E if ee!=0)); bset={0.0,1.0}
    for e in E:
        for j in range(8):
            b=j/7.0; m=0
            while True:
                xv=(b+m)/e
                if xv>=1: break
                if xv>=0: bset.add(xv)
                m+=1
    B=sorted(bset); tot=0.0
    for lo,hi in zip(B,B[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int((ee*mid)%1*7) for ee in E)&set(range(1,7)))==6: tot+=hi-lo
    return tot

if __name__=="__main__":
    k=8
    print(f"k={k}: the Jensen functional Phi(E)=int sum_i C(s_i,2) (LOW = balanced = high c_3 = high coverage)")
    AP=list(range(k))
    cands = {
        "AP [0..k-1]": AP,
        "perturb [0..k-2,k]": list(range(k-1))+[k],
        "perturb [0..k-3,k-1,k]": list(range(k-2))+[k-1,k],
        "spread [0,1,2,4,8,16,32,64][:k]": [0,1,2,4,8,16,32,64][:k],
        "random [0,3,5,9,11,17,23,30]": [0,3,5,9,11,17,23,30][:k],
    }
    rows=[]
    for name,E in cands.items():
        rows.append((name, Phi(E), p0(E)))
    rows.sort(key=lambda r:r[1])
    print(f"  {'set':30s} {'Phi (Jensen)':>13s} {'p0 (cover)':>11s}")
    for name,ph,p in rows:
        print(f"  {name:30s} {ph:13.3f} {p:11.4f}")
    print(f"  => AP has the SMALLEST Phi? {rows[0][0].startswith('AP')};  and LARGEST p0?")
    print("\nmac-mini's question -- p0 scaling vs translation invariance:")
    base=list(range(k))
    print(f"  p0([0..{k-1}])      = {p0(base):.4f}")
    print(f"  p0(2*[0..{k-1}])    = {p0([2*x for x in base]):.4f}   (DILATION: equal => scaling-invariant)")
    print(f"  p0([1..{k}])        = {p0([x+1 for x in base]):.4f}   (TRANSLATION by 1)")
    print(f"  p0([2..{k+1}])      = {p0([x+2 for x in base]):.4f}   (TRANSLATION by 2)")
