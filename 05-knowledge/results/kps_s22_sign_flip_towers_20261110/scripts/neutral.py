"""neutral.py - neutral gadgets: pairs (X, P) with kappa_X(P) = Pi2(X) - 2 W_X(P) - 2 W3_X(P) = 0.
For such a gadget the offset transforms as r+ = r_X(P) * r for EVERY base state (M does not enter)."""
import sys, collections
from towers import tables
def adj(m, code):
    A=[[0]*m for _ in range(m)]; b=0
    for i in range(m):
        for j in range(i+1,m):
            if code>>b&1: A[i][j]=1
            else: A[j][i]=1
            b+=1
    return A
def pairmats(A):
    n=len(A); full=(1<<n)-1; en,st,hs=tables(A)
    Q2=[[0]*n for _ in range(n)]; Q3=[[0]*n for _ in range(n)]
    for X in range(1,full+1):
        rest=full^X
        ex=[(e,en[X][e]) for e in range(n) if X>>e&1 and en[X][e]]
        Y=rest
        while Y:
            Z=rest^Y; hz=hs[Z]
            for f in range(n):
                if Y>>f&1 and st[Y][f]:
                    c=st[Y][f]
                    for e,ce in ex:
                        if Z==0: Q2[e][f]+=ce*c
                        else: Q3[e][f]+=ce*c*hz
            Y=(Y-1)&rest
    pi2=sum(hs[X]*hs[full^X] for X in range(full+1))
    return hs[full],pi2,Q2,Q3
MMAX=int(sys.argv[1]) if len(sys.argv)>1 else 6
stats=collections.Counter(); ex={}
for m in range(1,MMAX+1):
    for line in open("cls_%d.txt"%m):
        code=int(line.split()[0]); A=adj(m,code); h,pi2,Q2,Q3=pairmats(A)
        for P in range(1,1<<m):
            vs=[i for i in range(m) if P>>i&1]
            w=sum(Q2[e][f] for e in vs for f in vs if e!=f); w3=sum(Q3[e][f] for e in vs for f in vs if e!=f)
            kap=pi2-2*w-2*w3; r=h-w
            if kap==0:
                key=(m,len(vs),r,h)
                stats[key]+=1; ex.setdefault(key,(code,vs))
print("neutral gadgets (kappa = 0) by (order, |P|, r_X, H(X)): count [first witness code, P]")
for k in sorted(stats): print(" ",k,stats[k],ex[k])
