from functools import lru_cache

def cuts(rows):
    n=len(rows); full=(1<<n)-1
    h=[0]*(full+1); h[0]=1
    all_d=[]
    for b in range(n):
        d=[[0]*n for _ in h]; d[1<<b][b]=1
        for s in range(1,full+1):
            if not s>>b&1: continue
            for a,k in enumerate(d[s]):
                if not k: continue
                h[s]+=k
                t=rows[a] & (full^s)
                while t:
                    bit=t&-t; t-=bit
                    d[s|bit][bit.bit_length()-1]+=k
        all_d.append(d)
    w=[[sum(d[s][a]*h[full^s] for s in range(1,full+1) if not s.bit_count()%2) for d in all_d] for a in range(n)]
    k=[sum(w[a][b] for a in range(n) if mask>>a&1 for b in range(n) if not mask>>b&1) for mask in range(full+1)]
    return h[full],w,k

def extend(rows,x):
    n=len(rows)
    return tuple(row|((x>>a&1)<<n) for a,row in enumerate(rows))+(((1<<n)-1)^x,)

def direct(rows):
    n=len(rows); d=[[0]*n for _ in range(1<<n)]
    for a in range(n): d[1<<a][a]=1
    for s in range(1,1<<n):
        for a,k in enumerate(d[s]):
            if k:
                for b in range(n):
                    if rows[a]>>b&1 and not s>>b&1: d[s|1<<b][b]+=k
    return sum(d[-1])


from itertools import permutations
from math import isqrt

def prime(p):
    return p>=2 and all(p%q for q in range(2,isqrt(p)+1))

A=(102,252,248,241,225,192,128,1)
B=(16,5,17,7,74,31,175,63)
ac=tuple(255^row^(1<<i) for i,row in enumerate(A))
p=next(p for p in permutations(range(8))
       if all((B[i]>>j&1)==(ac[p[i]]>>p[j]&1) for i in range(8) for j in range(8)))
print('63 cores are converse under',p)
front={A:()}
pot=[]
while front:
    nxt={}; lane={}
    for a,witness in front.items():
        h,_,k=cuts(a)
        assert h==direct(a) and h<=301
        if h%7==0 and prime(h//7):
            lane.setdefault(h//7,witness)
        for x,c in enumerate(k[1:-1],1):
            if h+2*c<=301:
                assert c>0
                nxt[extend(a,x)]=witness+(x,)
    print('size, labelled extensions, prime quotient/masks',len(next(iter(front))),len(front),lane)
    pot.append(set(lane))
    front=nxt
assert pot==[set(),{19,23,37,43},{29,31,41,43}]
