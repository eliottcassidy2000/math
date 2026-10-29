"""Exact neighbourhood-recursion controls; standard library, no assert."""
import itertools, json
from functools import lru_cache

def need(x,msg):
    if not x: raise RuntimeError(msg)
def top(V,E):
    # Vertex belongs to a source SCC iff every predecessor reaches back.
    R={(x,x) for x in V}|set(E)
    for z in V:
        R|={(x,y) for x in V for y in V if (x,z) in R and (z,y) in R}
    return frozenset(x for x in V if all((x,y) in R for y in V if (y,x) in R))

def check_tournament(n,E):
    V=tuple(range(n))
    @lru_cache(None)
    def fixed(W):
        if not W:return frozenset()
        return top(W,{(x,y) for y in W for x in fixed(tuple(z for z in W if (z,y) in E))})
    @lru_cache(None)
    def depth(W):
        if len(W)<=1:return 0
        return 1+max([0]+[depth(tuple(z for z in W if (z,y) in E)) for y in W])
    @lru_cache(None)
    def iteration(W,k,seed):
        if not W:return frozenset()
        if k==0:
            cw=frozenset(x for x in W if all(x==y or (x,y) in E for y in W))
            return frozenset(W) if seed==0 or not cw else cw
        return top(W,{(x,y) for y in W for x in iteration(tuple(z for z in W if (z,y) in E),k-1,seed)})
    d=depth(V)
    need(d<=max(0,n-1),'depth bound')
    for s in (0,1):need(iteration(V,d,s)==fixed(V),'stabilization')
    return d

hist={}
for n in range(1,6):
    pairs=list(itertools.combinations(range(n),2));count=0
    for bits in range(1<<len(pairs)):
        E={(a,b) if bits>>i&1 else (b,a) for i,(a,b) in enumerate(pairs)}
        check_tournament(n,E);count+=1
    hist[n]=count
V=tuple(range(5))
E={(0,1),(2,0),(0,3),(4,0),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)}
a=(2,0,1,0,3); b=(4,0,1,2,3)
need(all((N[y],y) in E for N in (a,b) for y in V),'in-neighbour guard')
autos=[p for p in itertools.permutations(V) if {(p[x],p[y]) for x,y in E}==E]
need(len(autos)==1,'rigid host')
A=top(V,{(a[y],y) for y in V});B=top(V,{(b[y],y) for y in V})
need(A==frozenset((0,1,2)) and B==frozenset(V),'hostile choices')
print(json.dumps({'labelled_tournaments_checked':hist,'seeds':['all vertices','Condorcet winner if present, else all'],'hostile_host_edges':sorted(E),'host_automorphisms':len(autos),'singleton_neighbours_A':a,'singleton_neighbours_B':b,'winners_A':sorted(A),'winners_B':sorted(B)},indent=2,sort_keys=True))
