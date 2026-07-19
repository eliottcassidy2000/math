"""Hunt for a SECOND interior small gap: strict-interior families (1/14<M<1/13, val>=2)
whose sorted residues have >=2 gaps < val (=> non-AP core). If none exist even WITHOUT
covering, the kernel is a pure maximizer statement. Logs covering status of any hit.
death-star-S58e."""
from fractions import Fraction as F
from math import gcd
import itertools, random
def Mexact(V):
    mx=max(V); Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            d=1
            while d*d<=s:
                if s%d==0:
                    if 2<=d<=2*mx: Qs.add(d)
                    if 2<=s//d<=2*mx: Qs.add(s//d)
                d+=1
    best=F(0); bq=None
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=q
            for v in V:
                r=(v*a)%q; dd=r if r<q-r else q-r
                if dd<m: m=dd
                if m==0: break
            if m and (F(m,q)>best or (F(m,q)==best and (bq is None or m<bq[0]))):
                best=F(m,q); bq=(m,q,a)
    return best,bq
def ngaps(V,bq):
    val,q,a=bq; R=sorted((v*a)%q for v in V)
    gaps=[R[i+1]-R[i] for i in range(len(R)-1)]
    return sum(1 for g in gaps if g<val),val,q,gaps
def cover14(V): return all(any(v%k==0 for v in V) for k in range(2,15))
LO,HI=F(1,14),F(1,13)
random.seed(7)
hits=0; tested=0
# curated far-from-AP cores + random cores; far elements over a range
curated=[
 [1,2,4,8,16,32,64,128,256,512,1024,2048],
 [1,2,3,5,8,13,21,34,55,89,144,233],
 [1,3,5,7,9,11,13,15,17,19,21,23],
 [1,2,3,4,5,6,8,10,12,14,16,18],
 [1,2,3,4,5,6,7,8,9,10,11,25],
 [2,3,5,7,11,13,17,19,23,29,31,37],
 [1,4,9,16,25,36,49,64,81,100,121,144],
]
def trial(C):
    global tested,hits
    mc=max(C)
    for w in list(range(mc+1, 13*mc+2, 1))[:600]:
        if w in C: continue
        V=sorted(C+[w])
        if len(set(V))<13: continue
        M,bq=Mexact(V); tested+=1
        if bq and LO<M<HI:
            ns,val,q,gaps=ngaps(V,bq)
            if ns>=2:
                hits+=1
                print("HIT V=%s M=%s val=%d q=%d cover14=%s #gaps<val=%d gaps=%s"%(V,M,val,q,cover14(V),ns,gaps),flush=True)
for C in curated: trial(C)
for _ in range(4000):
    C=sorted(random.sample(range(1,45),12))
    if C[0]!=1: C[0]=1; C=sorted(set(C))
    if len(C)==12: trial(C)
print("DONE: tested=%d strict-interior-with-2gaps hits=%d"%(tested,hits),flush=True)
