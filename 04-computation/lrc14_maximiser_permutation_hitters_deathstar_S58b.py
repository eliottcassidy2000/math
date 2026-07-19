from fractions import Fraction as F
import itertools
def bad_from_g(g):
    cuts=[]
    for x in g:
        h=F(7,6)*x
        a=max(h-F(1,6),F(0)); b=min(h,F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1-cur>L: L=1-cur
    return L<=F(1,6)
CENTRES=[tuple(F(x,4) for x in p) for p in itertools.permutations([1,2,3])]
def hits_and_sojourn(DS, N):
    # sojourn = fraction of u in [0,1) that is bad; also does it hit a centre exactly?
    badcount=0; hit=False; minsup=None
    for t in range(N):
        u=F(t,N)
        g=tuple((-d*u)%1 for d in DS)
        if bad_from_g(g): badcount+=1
        for c in CENTRES:
            dd=max(min(abs(gi-ci),1-abs(gi-ci)) for gi,ci in zip(g,c))
            if minsup is None or dd<minsup: minsup=dd
            if dd==0: hit=True
    return hit, float(badcount)/N, float(minsup)
N=5040
print("dir            hit?  sojourn(bad measure)   min-sup-dist   vs 2/21=%.4f"%(2/21))
for DS in [(1,2,3),(1,2,7),(1,6,3),(5,2,3),(1,2,11),(2,4,6),(2,4,14),(1,2,4),(2,4,7)]:
    h,s,m=hits_and_sojourn(DS,N)
    prop = (DS[1]==2*DS[0] and DS[2]==3*DS[0])
    flag=""
    if not prop and h: flag=" <== NON-PROPORTIONAL but HITS (kps missed)"
    if not prop and s>2/21: flag+=" *** EXCEEDS 2/21 ***"
    print("%-14s %-5s %.5f               %.5f      %s"%(str(DS),"YES" if h else "no",s,m,flag))
