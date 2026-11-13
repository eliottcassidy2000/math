from itertools import combinations, permutations, combinations_with_replacement
from collections import Counter, defaultdict
from fractions import Fraction

def rows(n,code):
    r=[0]*n
    for a,(i,j) in enumerate(combinations(range(n),2)):
        if not (code>>a&1):i,j=j,i
        r[i]|=1<<j
    return r

def orbit_scan(n):
    V=tuple(range(n)); pairs=tuple(combinations(V,2)); M=len(pairs)
    P=list(permutations(V)); pairidx={e:i for i,e in enumerate(pairs)}
    signatures=[]
    for p in P:
        plus=minus=0; edges=[]
        for u,v in zip(p,p[1:]):
            e=pairidx[min(u,v),max(u,v)]; edges.append(e)
            if u<v:plus|=1<<e
            else:minus|=1<<e
        signatures.append((plus,minus,edges))
    mask=(1<<M)-1
    seen={}
    for c in range(1<<M):
        pos=set(); wedges=set()
        for plus,minus,edges in signatures:
            if c&plus==plus and c&minus==0:
                pos.update(edges)
                wedges.update(tuple(sorted((a,b))) for a,b in zip(edges,edges[1:]))
        g={i:set()for i in pos}
        for i,j in wedges:g[i].add(j);g[j].add(i)
        if g:
            todo=[next(iter(g))]; visit=set(todo)
            while todo:
                for j in g[todo.pop()]-visit:visit.add(j);todo.append(j)
            if len(visit)!=len(g):raise ValueError((n,c,'disconnected'))
        key=frozenset(pos),frozenset(wedges); out=min(c,c^mask)
        if key in seen and seen[key]!=out:raise ValueError((n,c,'collision'))
        seen[key]=out
    if len(seen)!=(1<<(M-1)):raise ValueError((n,'incorrect image'))
    print('one-slot wedge images',n,len(seen))

def complete_slots(T,q):
    n=len(T); h=0
    edge=Counter(); start=Counter();end=Counter();wedge=Counter()
    for p in permutations(range(n)):
        must={i for i in range(1,n)if not T[p[i-1]]>>p[i]&1}
        for D in combinations_with_replacement(range(n+1),q-1):
            if not must<=set(D):continue
            h+=1; r=(0,)+D+(n,)
            for s,(i,j) in enumerate(zip(r,r[1:])):
                if i==j:continue
                start[s,p[i]]+=1;end[s,p[j-1]]+=1
                for a,b in zip(p[i:j],p[i+1:j]):edge[s,min(a,b),max(a,b)]+=1
                for a,b,c in zip(p[i:j],p[i+1:j],p[i+2:j]):
                    wedge[s,min(a,c),b,max(a,c)]+=1
    return h,edge,start,end,wedge

for n in range(2,6):orbit_scan(n)
# C3 block -> chain of two -> singleton -> C3 block.
A=[2|24,4|24,1|24,16|32,32,7]
B=[4|24,1|24,2|24,16|32,32,7]
for q in range(1,4):
    x,y=complete_slots(A,q),complete_slots(B,q)
    if x[:4]!=y[:4] or x[4]==y[4]:raise ValueError(q)
print('all first moments agree in three slot probes of literal cyclic inflation; wedges differ')
