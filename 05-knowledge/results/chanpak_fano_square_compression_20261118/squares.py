from collections import Counter
from math import isqrt

def paths(n):
    adj=[0]*n
    for u in range(n):
        for v in range(n):
            if u!=v and isqrt(u+v+2)**2==u+v+2: adj[u]|=1<<v
    full=(1<<n)-1
    def dfs(w,t):
        if not t:
            if w[0]<w[-1] or n==1: yield tuple(u+1 for u in w)
        else:
            avail=t&adj[w[-1]]
            while avail:
                b=avail&-avail;avail^=b
                yield from dfs(w+[b.bit_length()-1],t^b)
    for v in range(n): yield from dfs([v],full^(1<<v))

def face(w,v):
    i=w.index(v)
    if i in (0,len(w)-1): return True
    x=w[i-1]+w[i+1]
    return isqrt(x)**2==x

counts=[]
no_32=[]
for n in range(1,33):
    p=list(paths(n))
    closable=[w for w in p if len(w)>2 and isqrt(w[0]+w[-1])**2==w[0]+w[-1]]
    assert len(closable)%n==0
    counts.append((n,len(p),len(closable)//n))
    if n>=2:
        F=Counter(sum(face(w,v) for v in (n-1,n)) for w in p)
        if 2*F[1]==32: no_32.append(n)
        if n==32:
            print("oriented physical deletion at 31,32:",tuple(2*F[k] for k in range(3)))
    if n==31:
        print("G31 endpoint incidence 4,17:", sorted(Counter(
            (4 in (w[0],w[-1]),17 in (w[0],w[-1])) for w in p).items()))
    if n==32:
        print("G32 directions up to reversal:",Counter(
            ("end 4" if w[1 if w[0]==32 else -2]==4 else "end 17") if 32 in (w[0],w[-1])
            else "interior" for w in p))
        print("spanning cycle first opened at 32:",closable[0])
        print("newborn endpoint paths among 32 cuts:",sum(32 in (w[0],w[-1]) for w in closable))
print("support:",[n for n,h,c in counts if h])
print("cycles:",[(n,c) for n,h,c in counts if c])
print("counts from 25:",counts[24:])
print("orders 2–32 with 32 oriented one-face paths:",no_32)
