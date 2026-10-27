#!/usr/bin/env python3
"""Exact finite controls for Session73 shared-port/rank-debt laws; no asserts."""
import itertools, json, random
from fractions import Fraction
from collections import Counter

def require(test, label):
    if not test: raise RuntimeError(label)

def bits(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def tournament(n, arcs):
    out=[0]*n
    for i,j in arcs: out[i] |= 1<<j
    require(all(((out[i]>>j)&1)+((out[j]>>i)&1)==1 for i,j in bits(n)), 'not a tournament')
    return tuple(out)
def from_mask(n,mask):
    return tournament(n, ((i,j) if mask>>k&1 else (j,i) for k,(i,j) in enumerate(bits(n))))
def path_tournament(n,mask):
    free=[(i,j) for i,j in bits(n) if j>i+1]
    arcs=[(i,i+1) for i in range(n-1)]
    arcs += [(i,j) if mask>>k&1 else (j,i) for k,(i,j) in enumerate(free)]
    return tournament(n,arcs)
def path_complement(out,P):
    n=len(out); pe={frozenset(e) for e in zip(P,P[1:])}
    require(all(out[i]>>j&1 for i,j in zip(P,P[1:])), 'illegal path')
    return tournament(n, ((i,j) if (bool(out[i]>>j&1) == (frozenset((i,j)) in pe)) else (j,i) for i,j in bits(n)))
def fas(out,weights=None):
    n=len(out)
    if weights is None:
        d=[0]*(1<<n)
        for S in range(1,1<<n):
            z=S; best=n*n
            while z:
                v=(z&-z).bit_length()-1; z &= z-1
                R=S^(1<<v); best=min(best,d[R]+(out[v]&R).bit_count())
            d[S]=best
        return d[-1]
    w=weights; d=[Fraction(0)]*(1<<n)
    for S in range(1,1<<n):
        z=S; best=None
        while z:
            v=(z&-z).bit_length()-1; z &= z-1
            R=S^(1<<v)
            c=d[R]+w[v]*sum((w[u] for u in range(n) if R>>u&1 and out[v]>>u&1),Fraction(0))
            if best is None or c<best: best=c
        d[S]=best
    return d[-1]
def rank_cost(out,order,weights=None):
    if weights is None: weights=[1]*len(out)
    return sum(weights[u]*weights[v] for i,u in enumerate(order) for v in order[i+1:] if out[v]>>u&1)
def brute_fas(out,weights=None): return min(rank_cost(out,p,weights) for p in itertools.permutations(range(len(out))))
def subst(q,children):
    offsets=[]; pos=0
    for t in children: offsets.append(pos); pos+=len(t)
    arcs=[]
    for i,t in enumerate(children):
        arcs.extend((offsets[i]+u,offsets[i]+v) for u in range(len(t)) for v in range(len(t)) if t[u]>>v&1)
    for i,j in bits(len(q)):
        for u in range(len(children[i])):
            for v in range(len(children[j])):
                x,y=offsets[i]+u,offsets[j]+v
                arcs.append((x,y) if q[i]>>j&1 else (y,x))
    return tournament(pos,arcs)
def paley(p):
    squares={x*x%p for x in range(1,p)}
    return tournament(p, ((i,j) if (j-i)%p in squares else (j,i) for i,j in bits(p)))
def triangles(out):
    tris=[]; cnt=Counter()
    for a,b,c in itertools.combinations(range(len(out)),3):
        z=[a,b,c]
        if all(sum(out[u]>>v&1 for v in z if u!=v)==1 for u in z):
            tris.append((a,b,c))
            for e in itertools.combinations(z,2):cnt[e]+=1
    return tris,cnt

def cycle_rank(n,edges):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]];x=parent[x]
        return x
    for x,y in edges:
        a,b=find(x),find(y)
        if a!=b:parent[a]=b
    c=len({find(x) for x in range(n)})
    return len(edges)-n+c,c

def exact_port_gluing(children_edges, child_ns, hinge_edges):
    """Compare literal beta to component-incidence multigraph beta."""
    offsets=[];n=0
    for sz in child_ns: offsets.append(n);n+=sz
    all_edges=[]; comps=[]; owner={}; child_beta=0
    for k,(sz,edges) in enumerate(zip(child_ns,children_edges)):
        beta,_=cycle_rank(sz,edges);child_beta+=beta
        adj=[set() for _ in range(sz)]
        for x,y in edges: adj[x].add(y);adj[y].add(x)
        seen=set()
        for x in range(sz):
            if x in seen:continue
            Q=[x];seen.add(x);C=[]
            for v in Q:
                C.append(v)
                for u in adj[v]:
                    if u not in seen:seen.add(u);Q.append(u)
            c=len(comps);comps.append(C)
            for v in C:owner[offsets[k]+v]=c
        all_edges.extend((offsets[k]+x,offsets[k]+y) for x,y in edges)
    all_edges+=hinge_edges
    beta,_=cycle_rank(n,all_edges)
    contracted=[(owner[x],owner[y]) for x,y in hinge_edges]
    bq,_=cycle_rank(len(comps),contracted)
    require(beta==child_beta+bq,'port component debt')
    return beta,bq,len(comps)

result={'status':'FINITE-EXACT; proofs in REPORT.md','seed':73073,'universes':{'complement':'all fixed-path presentations at n2..7; equivalent to all marked tournaments modulo path relabelling','weighted_substitution':'72 seeded arbitrary quotient/child instances at total order≤8 with positive rational vertex costs','port_gluing':'200 seeded simple child/cross-edge configurations, quotient multiplicities retained'}}
# Whole fixed-path universe, equivalent to all labelled marked tournaments up to relabelling P.
complement=[]
for n in range(2,8):
    width=(n-1)*(n-2)//2;count=1<<width;vals=[fas(path_tournament(n,m)) for m in range(count)]
    tariff=n-1-min(2,n//2)
    maximum=max(abs(vals[m]-vals[count-1-m]) for m in range(count))
    require(maximum<=tariff,'complement Lipschitz')
    witness=next(m for m in range(count) if abs(vals[m]-vals[count-1-m])==maximum)
    complement.append({'n':n,'marked_presentations':count,'bound':tariff,'max_difference':maximum,'witness_free_mask':witness})
result['complement_lipschitz']=complement
# Independent permutations DP cross-check through n5, all labelled tournaments through n4.
rng=random.Random(73073);checks=0
for n in range(2,5):
    for m in range(1<<(n*(n-1)//2)):
        t=from_mask(n,m);require(fas(t)==brute_fas(t),'DP/permutation mismatch');checks+=1
for _ in range(20):
    n=5;t=from_mask(n,rng.getrandbits(10));w=[Fraction(rng.randrange(1,8),rng.randrange(1,5)) for _ in range(n)]
    require(fas(t,w)==brute_fas(t,w),'weighted DP/permutation mismatch');checks+=1
result['independent_DP_checks']=checks
# Arbitrary uneven children and rational vertex weights; source law kept separate.
weighted=[]
for qn in range(2,5):
    for case in range(24):
        sizes=[rng.randrange(1,3) for _ in range(qn)]
        q=from_mask(qn,rng.getrandbits(qn*(qn-1)//2))
        children=[from_mask(n,rng.getrandbits(n*(n-1)//2)) for n in sizes]
        cw=[[Fraction(rng.randrange(1,8),rng.randrange(1,5)) for _ in range(n)] for n in sizes]
        w=[a for x in cw for a in x];A=[sum(x,Fraction(0)) for x in cw]
        actual=fas(subst(q,children),w);expected=sum((fas(t,a) for t,a in zip(children,cw)),Fraction(0))+fas(q,A)
        require(actual==expected,'weighted substitution')
        weighted.append(str(actual))
result['weighted_substitution']={'cases':len(weighted),'nonuniform_rational_weights':True,'cost_checksum':str(sum(map(Fraction,weighted),Fraction(0)))}
# Exact Paley counts and independent FAS DP; no all-word reachability conclusion.
paley_result=[]
for p in (3,7,11,19):
    t=paley(p);tri,inc=triangles(t)
    require(len(tri)==p*(p*p-1)//24,'Paley triangle total')
    require(set(inc.values())=={(p+1)//4},'Paley edge incidence')
    lb=(len(tri)+(p+1)//4-1)//((p+1)//4)
    row={'p':p,'triangles':len(tri),'edge_triangle_incidence':(p+1)//4,'FAS_lower':lb,'forest_one_complement_upper':2*p-4 if p>=4 else None}
    if p<=11:row['exact_FAS']=fas(t)
    if p>=11:require(lb>2*p-4,'Paley one-complement barrier')
    paley_result.append(row)
result['Paley']=paley_result
# Native inherited P7 predecessor, literal named path and forest witness.
p=tuple(range(7));t=paley(7);s=path_complement(t,p)
order=(0,3,1,6,4,2,5);pos={v:i for i,v in enumerate(order)}
d=[(i,j) for i,j in bits(7) if bool(s[i]>>j&1)!=(pos[i]<pos[j])]
beta,c=cycle_rank(7,d);require(beta==0 and c==1,'P7 predecessor is tree')
require(path_complement(s,p)==t,'P7 predecessor replay')
result['P7_predecessor']={'order':order,'exceptions':d,'beta':beta,'components':c}
# C3 substitution hierarchy: exact lower decomposes by first differing coordinate.
c3=tournament(3,[(0,1),(1,2),(2,0)]);t=(0,);towers=[]
for k in range(1,5):
    t=subst(c3,[t,t,t]);n=len(t);expected=n*(n-1)//6
    require(rank_cost(t,range(n))==expected,'tower rank upper')
    if n<=9:require(fas(t)==expected,'tower DP lower')
    debt,c=cycle_rank(n,[(i,j) for i,j in bits(n) if t[j]>>i&1])
    if k>=2:
        m=n//3;connected_order=[0]+list(range(m,3*m))+list(range(1,m))
        require(rank_cost(t,connected_order)==expected,'connected optimal tower rank cost')
        at={v:i for i,v in enumerate(connected_order)}
        connected_edges=[(i,j) for i,j in bits(n) if bool(t[i]>>j&1)!=(at[i]<at[j])]
        exact_debt,exact_components=cycle_rank(n,connected_edges)
        require(exact_components==1 and exact_debt==expected-n+1,'connected optimal tower debt')
    else:exact_debt=0
    towers.append({'k':k,'n':n,'minimum_cycle_debt_proved':exact_debt,'FAS_proved':expected,'natural_FAS_checked':expected,'beta_lower':max(0,expected-n+1),'natural_beta':debt,'natural_components':c,'DP_checked':n<=9})
result['C3_towers']=towers
# Inherited P7[C3^k] family: quotient optimum connected, hence exact cycle debt.
u=(0,);paley_towers=[]
for k in range(3):
    t=subst(paley(7),[u]*7);n=len(t);expected=n*(n-1)//6
    require(rank_cost(t,range(n))==expected,'P7/C3 tower edge certificate')
    debt,c=cycle_rank(n,[(i,j) for i,j in bits(n) if t[j]>>i&1])
    require(c==1 and debt==expected-n+1,'P7/C3 tower connected optimal debt')
    paley_towers.append({'k':k,'n':n,'FAS_proved':expected,'minimum_cycle_debt_proved':debt,'components':c})
    u=subst(c3,[u,u,u])
result['P7_C3_towers']=paley_towers
# Gluing exact component incidence law, preserving port identity/multiplicity.
gluing=[]
for _ in range(200):
    ns=[rng.randrange(1,5) for _ in range(rng.randrange(2,5))]
    es=[[e for e in bits(n) if rng.randrange(2)] for n in ns]
    starts=[];n=0;which=[]
    for i,sz in enumerate(ns):starts.append(n);which.extend([i]*sz);n+=sz
    he=[e for e in bits(n) if which[e[0]]!=which[e[1]] and rng.randrange(4)==0]
    gluing.append(exact_port_gluing(es,ns,he))
result['port_gluing']={'cases':len(gluing),'total_beta':sum(x[0] for x in gluing)}
# Hostile: merging distinct edges with same component endpoints destroys one cycle.
# Two child edges plus two cross edges form a 4-cycle; contracted interface has parallel edges.
exact=exact_port_gluing([[(0,1)],[(0,1)]],[2,2],[(0,2),(1,3)])
require(exact[0]==1 and exact[1]==1,'parallel interface cycle')
# Hostile: repeated named shared port in a tree doesn't independently create cycle.
shared=exact_port_gluing([[],[],[]],[1,1,1],[(0,1),(0,2)])
require(shared[0]==0,'shared-port star should not be cycle')
result['hostiles']={'contract_parallel_edges_beta':exact[0],'collapsed_simple_quotient_wrong_beta':0,'shared_port_star_beta':shared[0], 'unweighted_replacement_changes_value':fas(c3,[1,2,3])!=fas(c3)}
# Source-labelled deterministic operations preserve unequal masses and integer clocks.
mu={'a':Fraction(1,7),'b':Fraction(2,7),'c':Fraction(4,7)}
archive={s:{'mu':m,'clock':5,'ports':{'anchor':0,'cursor':2},'T':c3} for s,m in mu.items()}
for state in archive.values():state['T']=path_complement(state['T'],(0,1,2));state['clock']+=1
require({s:x['mu'] for s,x in archive.items()}==mu,'source masses changed')
require(all(x['clock']==6 and x['ports']=={'anchor':0,'cursor':2} for x in archive.values()),'clock or shared marks changed')
result['source_archive']={'masses':{s:str(m) for s,m in mu.items()},'integer_clock':6,'distinct_named_sources':len(archive)}
print(json.dumps(result,indent=2,sort_keys=True))
