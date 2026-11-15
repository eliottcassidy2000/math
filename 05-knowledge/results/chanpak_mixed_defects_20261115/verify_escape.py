from itertools import combinations
from collections import defaultdict

def vals(n):
    H=[1]*(1<<n)
    for x in range(1,1<<n):
        u,b=divmod(x,2)
        H[x]=H[u]+(0 if b else u+(u&1))
    return H

def witnesses(n,add):
    H=vals(n);by=defaultdict(list)
    for w,h in enumerate(H):by[h].append(w)
    for h in sorted(by):
        if any(h+j not in by for j in (2,4,6)):continue
        for i,j in (((2,4),) if add else combinations((2,4,6),2)):
            for a in by[h]:
                for b in by[h+i]:
                    x=a^b
                    for c in by[h+j]:
                        y=a^c
                        if not x&y and H[a^x^y]==h+12-i-j:
                            yield a,x,y

def exact(w,n):
    bit=lambda i:(w>>(n-1-i))&1
    pred=[0]*(n+1)
    for i in range(n):
        pred[i]=(1<<i)-1
        if bit(i):pred[i]|=1<<n
        else:pred[n]|=1<<i
    full=(1<<(n+1))-1
    dp=[[0]*(n+1) for _ in range(full+1)]
    for j in range(n+1):dp[1<<j][j]=1
    for S in range(1,full+1):
        if S.bit_count()==1:continue
        for j in range(n+1):
            if not (S>>j)&1:continue
            E=S^(1<<j);I=E&pred[j]
            while I:
                lsb=I&-I;i=lsb.bit_length()-1;I^=lsb
                dp[S][j]+=dp[E][i]
    return sum(dp[full])

def need(x):
    if not x:raise RuntimeError("interleaved reader escape failed")

for n0,w,x,y,counts,add in ((9,232,458,16,[213,215,219,217],False),
                          (11,1158,32,220,[721,723,725,727],True)):
    need([exact(z,n0) for z in (w,w^x,w^y,w^x^y)]==counts)
    need(not any(x&y for _ in (0,)))
    need(all(next(witnesses(n,add),None) is None for n in range(1,n0)))
    need((w,x,y) in list(witnesses(n0,add)))
print("8 independently counted physical tournaments; first spans 9 nonlinear and 11 additive")
