# Two-core refill: odd strong Y1 => Y2  ==  odd strong Z -> TT1 (same full slot polynomial)?  orders <= 10
from math import comb
def qv(F,T): return tuple(sum(F[j-1]*comb(t,j) for j in range(1,len(F)+1)) for t in range(1,T+1))
strong={}
for a in (3,5,7,9):
    strong[a]=[]
    for line in open('prof_%d.txt'%a):
        p=line.split()
        if int(p[2])==1: strong[a].append((p[1],tuple(map(int,p[7:]))))
for n in (6,8,10):
    T=n
    # t*Q_Z for Z odd strong of order n-1
    single={}
    for c,F in strong[n-1]:
        v=qv(F,T); key=tuple(t*x for t,x in zip(range(1,T+1),v)); single.setdefault(key,[]).append(c)
    hits=0
    for a in range(3,n,2):
        b=n-a
        if b<a or b<3: continue
        for c1,F1 in strong[a]:
            v1=qv(F1,T)
            for c2,F2 in strong[b]:
                v2=qv(F2,T); key=tuple(x*y for x,y in zip(v1,v2))
                if key in single:
                    hits+=1
                    if hits<=12: print('order',n,': Y%d %s => Y%d %s  ==  Z%d %s -> TT1   H=%d'%(a,c1,b,c2,n-1,single[key][:4],F1[0]*F2[0]))
    print('order',n,'two-core refills found:',hits)
