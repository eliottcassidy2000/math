#!/usr/bin/env python3
"""Exact, optimization-independent controls for shared-root diamond frameworks."""
from fractions import Fraction as F
from itertools import combinations, product
import json

# Q(sqrt(3),sqrt(11)); basis 1,sqrt(3),sqrt(11),sqrt(33).
Z = (F(0),)*4
ONE = (F(1),F(0),F(0),F(0))
def elt(a=0,b=0,c=0,d=0): return tuple(map(F,(a,b,c,d)))
def add(x,y): return tuple(a+b for a,b in zip(x,y))
def neg(x): return tuple(-a for a in x)
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    z=[F(0)]*4
    for i,a in enumerate(x):
        for j,b in enumerate(y):
            z[i^j] += a*b*(3 if i&j&1 else 1)*(11 if i&j&2 else 1)
    return tuple(z)
def norm2(p,q):
    d=[sub(x,y) for x,y in zip(p,q)]
    return add(mul(d[0],d[0]),mul(d[1],d[1]))
def demand(ok,msg):
    if not ok: raise RuntimeError(msg)
def coords(m):
    pts=[(Z,Z)]
    A,B=F(1),F(0)
    for k in range(m):
        pts += [(elt(b=A),elt(d=B)),
                (elt(b=A/2,c=-B/2),elt(a=A/2,d=B/2)),
                (elt(b=A/2,c=B/2),elt(a=-A/2,d=B/2))]
        A,B=(5*A-11*B)/6,(A+5*B)/6
    return pts

def graph(m):
    edges=set()
    for k in range(m):
        t,a,b=1+3*k,2+3*k,3+3*k
        edges |= {tuple(sorted(e)) for e in [(0,a),(0,b),(a,b),(a,t),(b,t)]}
        if k: edges.add((t-3,t))
    return 3*m+1,edges

def adjacency(n,E):
    adj=[0]*n
    for a,b in E: adj[a]|=1<<b;adj[b]|=1<<a
    return adj

def rank_mod(mat,p):
    a=[r[:] for r in mat]; rank=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]%p),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        inv=pow(a[rank][col],p-2,p)
        a[rank]=[v*inv%p for v in a[rank]]
        for i in range(rank+1,len(a)):
            c=a[i][col]%p
            if c: a[i]=[(x-c*y)%p for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a): break
    return rank

def primality(n): return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
P=next(p for p in range(100001,101000) if p%132==1 and primality(p))
S3=next(x for x in range(P) if x*x%P==3)
S11=next(x for x in range(P) if x*x%P==11)
def reduce_field(x):
    return sum((v.numerator*pow(v.denominator,-1,P)%P)*r for v,r in zip(x,(1,S3,S11,S3*S11)))%P

def matrix_rank(points,E):
    n=len(points); reduced=[tuple(map(reduce_field,p)) for p in points]; rows=[]
    for u,v in sorted(E):
        row=[0]*(2*n)
        for d in range(2):
            diff=(reduced[u][d]-reduced[v][d])%P
            row[2*u+d]=diff;row[2*v+d]=-diff%P
        rows.append(row)
    return rank_mod(rows,P)

def components(n,E,deleted=()):
    adj=adjacency(n,E); unseen=((1<<n)-1)&~sum(1<<v for v in deleted); sizes=[]
    while unseen:
        front=unseen&-unseen;seen=0
        while front:
            v=(front&-front).bit_length()-1;front &= front-1
            if seen>>v&1: continue
            seen|=1<<v;front|=adj[v]&unseen&~seen
        unseen &= ~seen;sizes.append(seen.bit_count())
    return sorted(sizes)

def hp_count(n,E,cycle=False):
    adj=adjacency(n,E); dp=[{} for _ in range(1<<n)]
    for v in ([0] if cycle else range(n)): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,count in dp[mask].items():
            cand=adj[v]&~mask
            while cand:
                b=cand&-cand;cand-=b;w=b.bit_length()-1
                nxt=dp[mask|b];nxt[w]=nxt.get(w,0)+count
    if cycle: return sum(c for v,c in dp[-1].items() if adj[v]&1)//2
    return sum(dp[-1].values())//2

def color_count(n,E,q):
    if q<=0:return 0
    adj=[set() for _ in range(n)]
    for u,v in E:adj[u].add(v);adj[v].add(u)
    colors=[-1]*n
    def rec(left):
        if not left:return 1
        v=max(left,key=lambda v:(len({colors[w] for w in adj[v] if colors[w]>=0}),len(adj[v])))
        forbidden={colors[w] for w in adj[v]};ans=0
        for c in range(q):
            if c not in forbidden:
                colors[v]=c;ans+=rec(left-{v})
        colors[v]=-1
        return ans
    return rec(set(range(n)))

def chromatic_transfer(m,q):
    A=(q-1)*(q-2);B=(q-2)*(q-3)
    x,y=A,(q-1)*B
    for _ in range(m-1): x,y=A*y,(q-1)*B*x+(q-2)*B*y
    return q*(x+y)

def linear_forest(n,E):
    adj=adjacency(n,E)
    return max((a.bit_count() for a in adj),default=0)<=2 and len(E)==n-len(components(n,E))

def bipartite(n,E):
    adj=adjacency(n,E);color={}
    for v in range(n):
        if v in color:continue
        color[v]=0;todo=[v]
        while todo:
            u=todo.pop();bits=adj[u]
            while bits:
                bit=bits&-bits;bits-=bit;w=bit.bit_length()-1
                if w in color:
                    if color[w]==color[u]:return False
                else:color[w]=1-color[u];todo.append(w)
    return True

def injects_integer_path(n,E):
    # Independently try all component placements within n integer addresses.
    # Unconstrained components can be translated apart afterwards.
    adj=adjacency(n,E);unseen=set(range(n))
    while unseen:
        root=min(unseen);component={root};todo=[root]
        while todo:
            v=todo.pop()
            for w in unseen-component:
                if adj[v]>>w&1:component.add(w);todo.append(w)
        unseen-=component;address={root:0}
        def rec():
            if len(address)==len(component):return True
            v=next(v for v in component-address.keys() if any(adj[v]>>w&1 for w in address))
            parents=[w for w in address if adj[v]>>w&1]
            candidates={address[parents[0]]-1,address[parents[0]]+1}
            for w in parents[1:]:candidates&={address[w]-1,address[w]+1}
            for x in candidates-set(address.values()):
                address[v]=x
                if rec():return True
                del address[v]
            return False
        if not rec():return False
    return True

def main():
    out={'field':'Q(sqrt(3),sqrt(11))','modular_rank_prime':P,'sqrt3_mod_p':S3,'sqrt11_mod_p':S11}
    coordrows=[];pair_checks=0
    for m in list(range(1,13))+[16,24]:
        pts=coords(m);n,E=graph(m)
        demand(len(set(pts))==n,'coordinate collision')
        actual=set()
        for i,j in combinations(range(n),2):
            pair_checks+=1
            if norm2(pts[i],pts[j])==ONE:actual.add((i,j))
        demand(actual==E,f'extra or missing edges at {m}')
        rank=matrix_rank(pts,E)
        demand(rank==6*m-1,'rigidity rank failure')
        if m>=2:
            demand(min(a.bit_count() for a in adjacency(n,E))==3,'degree')
            demand(all(len(components(n,E,[v]))==1 for v in range(n)),'cut vertex')
        if m>=5:demand(len(components(n,E,[0,4,10]))==5,'separator')
        coordrows.append({'diamonds':m,'vertices':n,'edges':len(E),'rigidity_rank':rank})
    out['exact_coordinate_rows']=coordrows;out['exact_squared_distance_checks']=pair_checks
    out['minimum_path_cover_certificates']=[]
    for m in range(1,25):
        n,E=graph(m);covers=[]
        for k in range(0,m,2):
            if k+1<m:covers.append([2+3*k,3+3*k,1+3*k,4+3*k,5+3*k,6+3*k])
            else:covers.append([2+3*k,3+3*k,1+3*k])
        original_count=len(covers)
        separator=[0]+[1+3*k for k in range(1,m-1,2)]
        lower_H=len(components(n,E,separator))-(len(separator)-1)
        demand(lower_H==original_count==(m+1)//2,'triangle path cover lower bound')
        if len(covers)>1:covers=[covers[0]+[0]+covers[1]]+covers[2:]
        else:covers[0]=[0]+covers[0]
        demand(sorted(v for path in covers for v in path)==list(range(n)),'cover partition')
        demand(all(tuple(sorted((u,v))) in E for path in covers for u,v in zip(path,path[1:])),'cover edge')
        demand(len(covers)==max(1,(m+1)//2-1),'path cover formula')
        out['minimum_path_cover_certificates'].append({'diamonds':m,'minimum_paths':len(covers),'cover':covers,'deleted_root_and_odd_tips':separator,'remaining_components':len(components(n,E,separator))})
    out['hamiltonian_rows']=[]
    for m in range(1,6):
        n,E=graph(m);hp=hp_count(n,E);hc=hp_count(n,E,True)
        demand(bool(hp)==(m<=4),'Hamilton path threshold')
        demand(bool(hc)==(m<=2),'Hamilton cycle threshold')
        out['hamiltonian_rows'].append({'diamonds':m,'unoriented_hamilton_paths':hp,'unoriented_hamilton_cycles':hc})
    out['color_rows']=[]
    for m in range(1,4):
        n,E=graph(m)
        for q in range(1,6):
            direct=color_count(n,E,q);formula=chromatic_transfer(m,q)
            demand(direct==formula,'chromatic transfer')
            out['color_rows'].append({'diamonds':m,'colors':q,'proper_colorings':direct})
    out['graph_classification_rows']=[]
    for n in range(1,6):
        pairs=list(combinations(range(n),2));lf=bp=total=0
        for mask in range(1<<len(pairs)):
            E={e for bit,e in enumerate(pairs) if mask>>bit&1};total+=1
            good=linear_forest(n,E);demand(good==injects_integer_path(n,E),'integer embedding classification')
            lf+=good;bp+=bipartite(n,E)
        out['graph_classification_rows'].append({'vertices':n,'all_graphs':total,'linear_forests':lf,'bipartite_graphs':bp})
    out['fresh_step_words']=[]
    for m in range(2,13):
        survivors=0
        for steps in product((-1,1),repeat=m-1):
            addresses=[0]
            for s in steps:addresses.append(addresses[-1]+s)
            survivors+=len(set(addresses))==m
        demand(survivors==2,'freshness pruning')
        out['fresh_step_words'].append({'diamonds':m,'local_words':2**(m-1),'injective_tip_words':survivors})
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
