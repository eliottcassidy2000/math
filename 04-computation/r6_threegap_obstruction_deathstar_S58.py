from fractions import Fraction as F
# min-L config; analyze the finest-killer gap structure & what a three-gap/equidistribution bound gives
P=[1,2,5,8,9,11,12]; K=[316,322,326,328,330]
def danger(v):
    w=F(1,14*v);o=[]
    for j in range(v):
        c=F(j,v);lo=(c-w)%1;hi=(c+w)%1
        if lo<hi:o.append((lo,hi))
        else:o.append((lo,F(1)));o.append((F(0),hi))
    return o
def subf(S,a):
    for clo,chi in sorted(a):
        n=[]
        for lo,hi in S:
            if chi<=lo or clo>=hi:n.append((lo,hi));continue
            if clo>lo:n.append((lo,clo))
            if chi<hi:n.append((chi,hi))
        S=n
    return S
S0=[(F(0),F(1))]
for v in P:S0=subf(S0,danger(v))
# find core arc with the max gap, count k5-safe gaps in it
G=S0
for k in K: G=subf(G,danger(k))
L=max(hi-lo for lo,hi in G); k5=max(K)
# which core arc holds the max gap
best=max(((hi-lo,lo,hi) for lo,hi in G),key=lambda x:x[0])
gaplo,gaphi=best[1],best[2]
Acore=[(a,b) for a,b in S0 if a<=gaplo<=b][0]
W=Acore[1]-Acore[0]
print("min-L config core=%s K=%s"%(P,K))
print("L=%.7f = %.4f/(7k5)  (k5=%d) ; core arc width W=%.5f"%(float(L),float(L*7*k5),k5,float(W)))
print("k5-safe gaps in core arc: N5 = W*k5 ~ %.1f  (equidistribution in the 4-torus of the other killers needs N5 >> ??)"%(float(W)*k5))
print("full k5-safe gap = 6/(7k5) = %.6f ; observed L = %.6f = %.2f of a full gap"%(6/(7*k5),float(L),float(L)/(6/(7*k5))))
# heuristic three-gap/mean bound vs truth
Sig=sum(K)
print("crude mean bound (2/7)/Sigma_k = %.7f ; truth L=%.7f (truth/crude=%.1f x -- the gap three-gap must recover)"%(2/(7*Sig),float(L),float(L)/(2/(7*Sig))))
print("target 1/2331=%.7f ; 1/(7k5)=%.7f"%(1/2331,1/(7*k5)))
# how many core arcs are wide enough to even hold ONE full k5-gap (W*k5>=1)?
import itertools
def csafeF(P):
    S=[(F(0),F(1))]
    for v in P:S=subf(S,danger(v))
    return S
narrowcount=0;tot=0
for Pc in itertools.combinations(range(1,13),7):
    lo=13*max(Pc)+1
    for a,b in csafeF(list(Pc)):
        tot+=1
        if float(b-a)*332 < 4: narrowcount+=1  # <4 finest-killer gaps: equidist hopeless
print("core-safe arcs across all 792 cores: %d ; with <4 finest-killer gaps (W*332<4): %d (%.0f%%)"%(tot,narrowcount,100*narrowcount/tot))
