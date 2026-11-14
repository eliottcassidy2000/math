from itertools import combinations
def k(a,b):
    d=abs(a-b)
    return (1<<(d-1))+(d==1)
def graph(n):
    v=[(a,b,s) for a,b in combinations(range(n),2) for s in (-1,1)]
    out=[0]*len(v)
    for i,(a,b,s) in enumerate(v):
        for j in range(i):
            c,d,t=v[j]
            if len({a,b,c,d})<4:continue
            if k(a,c)+s*k(b,c)+t*k(a,d)+s*t*k(b,d)==0:
                out[i]|=1<<j;out[j]|=1<<i
    return v,out
def maxcl(adj):
    best=()
    def go(cl,can):
        nonlocal best
        if len(cl)>len(best):best=cl
        while can:
            if len(cl)+can.bit_count()<=len(best):return
            x=can&-can;can^=x;i=x.bit_length()-1
            go((*cl,i),can&adj[i])
    go((),(1<<len(adj))-1)
    return best
for n in range(4,13):
    v,g=graph(n)
    for ix,(a,b,s) in enumerate(v):
        for j,(c,d,t) in enumerate(v[:ix]):
            if c<a:a1,b1,s1,c1,d1,t1=c,d,t,a,b,s
            else:a1,b1,s1,c1,d1,t1=a,b,s,c,d,t
            nested=a1<c1<d1<b1 and c1-a1==b1-d1
            predicted=nested and (s1!=t1 or s1==t1==-1 and d1-c1==c1-a1==1)
            assert bool(g[ix]>>j&1)==predicted
    q=maxcl(g)
    assert len(q)==(3 if n>=6 else 2)
    print(n,[v[i] for i in q])
assert k(0,2)-k(1,2)==0
print('PASS: all two-site interaction classifications and the separated carry exception')
