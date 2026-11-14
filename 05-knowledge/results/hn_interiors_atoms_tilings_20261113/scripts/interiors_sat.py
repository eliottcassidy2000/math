# m private interiors (default 8): interior-first SAT/CEGAR search for a forcing configuration in the paired-port
# selector architecture.  Exact planar necessary conditions as in six_interiors.py; colour groups of 2/4/6 points.
import sys, subprocess, itertools, networkx as nx
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
M=int(sys.argv[1]) if len(sys.argv)>1 else 8
NOGEO='--nogeo' in sys.argv
def admissible(H):
    if NOGEO: return True
    for Q in itertools.combinations(H.nodes,4):
        if all(H.has_edge(a,b) for a,b in itertools.combinations(Q,2)): return False
    for a,b in itertools.combinations(H.nodes,2):
        if len(set(H[a])&set(H[b]))>=3: return False
    return True
def list_colour(H,forb,ncol=6,kgroups=4):
    # find a proper colouring of H with colour(x) not in forb[x]; prefer free colours (>= kgroups)
    order=sorted(H.nodes,key=lambda v:-H.degree(v)); col={}
    pref=list(range(kgroups,ncol))+list(range(kgroups))
    def bt(i):
        if i==len(order): return True
        v=order[i]; used={col[u] for u in H[v] if u in col}
        for c in pref:
            if c in forb[v] or c in used: continue
            col[v]=c
            if bt(i+1): return True
            del col[v]
        return False
    return dict(col) if bt(0) else None
def search(H,sizes):
    k=len(sizes); n=M; pool=IDPool()
    pts=[(g,j) for g,s in enumerate(sizes) for j in range(s)]
    t={(p,x):pool.id(('t',p,x)) for p in range(len(pts)) for x in range(n)}
    a={(g,x):pool.id(('a',g,x)) for g in range(k) for x in range(n)}
    S=Solver(name='cadical153'); add=S.add_clause
    for g in range(k):
        for x in range(n):
            ps=[p for p,(gg,_) in enumerate(pts) if gg==g]
            add([-a[g,x]]+[t[p,x] for p in ps])
            for p in ps: add([-t[p,x],a[g,x]])
    nbr={x:set(H[x]) for x in H}
    cycles=[frozenset(c) for c in nx.simple_cycles(H) if len(c)>=3]
    cycles=list(set(cycles))
    P=len(pts)
    for p in (range(P) if not NOGEO else []):
        for x in range(n):                       # degree <= 2 inside a circle
            for tri in itertools.combinations(sorted(nbr[x]),3):
                add([-t[p,x]]+[-t[p,y] for y in tri])
            for tri in itertools.combinations(sorted(nbr[x]),3):   # interior x vs port p: <= 2 common nbrs
                add([-t[p,y] for y in tri])
        for C in cycles:
            if len(C)!=6: add([-t[p,c] for c in C])
            else:
                for y in range(n):
                    if y not in C: add([-t[p,c] for c in C]+[-t[p,y]])
    for p,q in (itertools.combinations(range(P),2) if not NOGEO else []):  # two port points share <= 2 contacts
        for tri in itertools.combinations(range(n),3):
            add([-t[p,y] for y in tri]+[-t[q,y] for y in tri])
    for x,y in (itertools.combinations(range(n),2) if not NOGEO else []):  # two interiors: common nbrs (H + port points) <= 2
        room=2-len(nbr[x]&nbr[y])
        if room<0: return 'bad-H'
        for grp in itertools.combinations(range(P),room+1):
            add([l for p in grp for l in (-t[p,x],-t[p,y])])
    for x in range(n):                          # criticality: deg_Q(x) >= 6
        need=6-H.degree(x)
        if need>k: return 'excluded'
        if need>0:
            for cl in CardEnc.atleast([a[g,x] for g in range(k)],bound=need,vpool=pool,encoding=EncType.seqcounter).clauses: add(cl)
    for g in range(k):                          # group degree >= 6
        need=7-k
        for cl in CardEnc.atleast([a[g,x] for x in range(n)],bound=need,vpool=pool,encoding=EncType.seqcounter).clauses: add(cl)
    it=0
    while S.solve():
        it+=1
        model=set(l for l in S.get_model() if l>0)
        forb={x:{g for g in range(k) if a[g,x] in model} for x in range(n)}
        phi=list_colour(H,forb,6,k)
        if phi is None:
            T=[[x for x in range(n) if t[p,x] in model] for p in range(P)]
            return ('SURVIVOR',T,it)
        cl=[a[phi[x],x] for x in range(n) if phi[x]<k]
        if not cl: return ('excluded',it)
        add(cl)
    return ('excluded',it)
cases={4:[(2,2,2,2)],3:[(2,2,2),(4,2,2)],2:[(2,2),(4,2),(4,4),(6,2)]}
tot=0; surv=[]; maxit=0
for k,sizelist in cases.items():
    out=subprocess.run(['geng','-q','-c',f'-d{6-k}',str(M)],capture_output=True,text=True).stdout.split()
    Hs=[nx.from_graph6_bytes(s.encode()) for s in out]
    Hs=[H for H in Hs if admissible(H)]
    for H in Hs:
        for sizes in sizelist:
            r=search(H,sizes); tot+=1
            if isinstance(r,tuple) and r[0]=='SURVIVOR': surv.append((nx.to_graph6_bytes(H,header=False).strip(),sizes,r[1])); print("SURVIVOR",surv[-1],flush=True)
            elif isinstance(r,tuple): maxit=max(maxit,r[1])
    print(f"k={k}: interior graphs {len(Hs)}, cases so far {tot}, survivors so far {len(surv)}, max CEGAR iterations {maxit}",flush=True)
print("TOTAL cases",tot,"survivors",len(surv))
