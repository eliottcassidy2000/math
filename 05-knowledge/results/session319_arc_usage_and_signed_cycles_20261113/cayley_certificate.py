from collections import Counter
V=[(a,b) for b in range(3) for a in range(7)]
S={(0,1),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(4,1),(4,2)}
def mul(x,y):
    a,b=x;c,d=y
    return ((a+pow(2,b,7)*c)%7,(b+d)%3)
out={u:{mul(u,s) for s in S} for u in V}
if any((v in out[u])==(u in out[v]) for u in V for v in V if u!=v):
    raise RuntimeError('not a tournament')
if any(any(mul(g,v) not in out[mul(g,u)] for v in out[u]) for g in V for u in V):
    raise RuntimeError('not Cayley')
seen={V[0]}
for _ in V:
    seen|=set().union(*(out[v] for v in seen))
if len(seen)!=21: raise RuntimeError('not strong')
def sig(out):
    Q=out[V[0]]
    d={v:len(Q & out[v]) for v in Q}
    special=[v for v in Q if d[v]==7]
    if len(special)!=1: raise RuntimeError('unique marker missing')
    return sorted(d[u] for u in Q & out[special[0]])
a=sig(out); b=sig({u:{v for v in V if u in out[v]} for u in V})
if (a,b)!=([3,4,4,4,5,5,5],[3,4,4,5,5,5,5]):
    raise RuntimeError('failed anti certificate')
print(a,b)
