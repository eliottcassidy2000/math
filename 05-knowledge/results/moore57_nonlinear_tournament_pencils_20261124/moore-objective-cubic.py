from itertools import combinations
from collections import Counter

def tup(E,n=10):
    m=[0]*n
    for i,j in E:m[i]|=1<<j;m[j]|=1<<i
    return tuple(m)
def edges(r):
    return {(i,j) for i,j in combinations(range(len(r)),2) if r[i]>>j&1}
def defc(r):
    t=sum((r[i]&r[j]).bit_count() for i,j in edges(r))//3
    c=sum(u*(u-1)//2 for i,j in combinations(range(len(r)),2)
          if (u:=(r[i]&r[j]).bit_count())>1)//2
    return 6*t+8*c
def invariants(r):
    n=len(r)
    return defc(r),tuple(sorted(tuple(sorted((r[i]&r[j]).bit_count() for j in range(n) if i!=j and (r[i]>>j&1)==bit)) for i in range(n) for bit in (0,1)))
def iso(r,s):
    n=len(r);V=(1<<n)-1
    c=[V]*n;p={}
    def rec(usable):
        if not usable:return True
        i=min((i for i in range(n) if i not in p),key=lambda i:(c[i]&usable).bit_count())
        opts=c[i]&usable
        while opts:
            y=opts&-opts;opts^=y;j=y.bit_length()-1
            old=c[:];p[i]=j
            for k in range(n):
                if k not in p:c[k]&=s[j] if r[i]>>k&1 else (~s[j]&V)
            if rec(usable^y):return True
            c[:]=old;del p[i]
        return False
    return rec(V)
def neighbors(r):
    E=edges(r)
    for (a,b),(c,d) in combinations(E,2):
        if len({a,b,c,d})!=4:continue
        for x,y in ((c,d),(d,c)):
            u,v=tuple(sorted((a,x))),tuple(sorted((b,y)))
            if u not in E and v not in E:yield tup(E-{(a,b),(c,d)}|{u,v})
V=list(combinations(range(5),2))
r=tup({(i,j) for i,j in combinations(range(10),2) if not set(V[i])&set(V[j])})
queue=[r];B={invariants(r):[r]};receipt=[]
for g in queue:
    value=defc(g);lo=value;n=0
    for u in neighbors(g):
        n+=1
        key=invariants(u);lo=min(lo,key[0])
        if not any(iso(u,v) for v in B.get(key,())):
            B.setdefault(key,[]).append(u);queue.append(u)
    receipt.append((value,lo,n))
print('states',len(queue),'height sorted',sorted(receipt))
assert len(queue)==21 and all(not a or b<a for a,b,_ in receipt)
