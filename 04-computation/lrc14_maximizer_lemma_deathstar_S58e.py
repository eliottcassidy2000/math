"""Verifies the strict-interior MAXIMIZER LEMMA (death-star-S58e):
for M=val/q with 1/14<M<1/13 and maximizer t=a/q:
 (i) both band edges occupied: min residue=val, max=q-val;
 (ii) edge speeds satisfy q|(v*+w*), hence v*+w*=q;
 (iii) span = 11*val + g0 exactly, g0 = q-13*val in {1,...,val-1} (=> val>=2).
Also records the HARMONIC-COMPETITOR bound min_i|k r_i|_q <= val (all k), saturated by the AP."""
from fractions import Fraction as F
from math import gcd
def Mexact(V):
    mx=max(V);Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0:Qs.add(d)
    best=F(0);bq=None
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1:continue
            m=q
            for v in V:
                r=(v*a)%q;dd=r if r<q-r else q-r
                if dd<m:m=dd
            if m and(F(m,q)>best or(F(m,q)==best and(bq is None or m<bq[0]))):best=F(m,q);bq=(m,q,a)
    return best,bq
for V in [list(range(1,13))+[182],list(range(1,13))+[26],list(range(1,13))+[364],list(range(1,13))+[546]]:
    M,(val,q,a)=Mexact(V);R=sorted((v*a)%q for v in V)
    print("Vmax=%d M=%s val=%d q=%d: edges(val,q-val)=%s v*+w*=%d(=q?%s) span=%d=11val+g0=%d g0=%d(in[1,val-1]?%s)"
          %(V[-1],M,val,q,(R[0]==val and R[-1]==q-val),
            (min(v for v in V if (v*a)%q==val)+min(v for v in V if (v*a)%q==q-val)),
            "yes" if (min(v for v in V if (v*a)%q==val)+min(v for v in V if (v*a)%q==q-val))==q else "no",
            R[-1]-R[0],11*val+(q-13*val),q-13*val,1<=q-13*val<=val-1))
