import sys, itertools, subprocess
from collections import defaultdict
def parse(line,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if line[k]=='1': out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def odd_cycles(out,n):
    cyc=[]
    for L in (3,5):
        if L>n: break
        for S in itertools.combinations(range(n),L):
            s0=S[0]
            for perm in itertools.permutations(S[1:]):
                seq=(s0,)+perm
                if all((out[seq[i]]>>seq[(i+1)%L])&1 for i in range(L)):
                    cyc.append(frozenset(S))
    return cyc
def strong(out,n):
    reach=1; fr=1
    while fr:
        nf=0
        for v in range(n):
            if (fr>>v)&1: nf|=out[v]
        nf&=~reach; reach|=nf; fr=nf
    if reach!=(1<<n)-1: return False
    inn=[0]*n
    for v in range(n):
        for u in range(n):
            if (out[v]>>u)&1: inn[u]|=1<<v
    reach=1; fr=1
    while fr:
        nf=0
        for v in range(n):
            if (fr>>v)&1: nf|=inn[v]
        nf&=~reach; reach|=nf; fr=nf
    return reach==(1<<n)-1
for n in (4,5,6):
    lines=subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split()
    occ=defaultdict(int)
    for line in lines:
        out=parse(line,n); cyc=odd_cycles(out,n)
        c3=sum(1 for c in cyc if len(c)==3); c5=sum(1 for c in cyc if len(c)==5)
        d=sum(1 for a,b in itertools.combinations([c for c in cyc if len(c)==3],2) if not (a&b))
        H=1+2*(c3+c5)+4*d
        occ[(strong(out,n),c3,c5,d,H)]+=1
    print(f"order {n}: (strong, c3, c5, d) -> H")
    for k in sorted(occ): print("   ",k[:4],"H =",k[4],"classes",occ[k])
    Hs=sorted(set(k[4] for k in occ)); print("   H values:",Hs)
