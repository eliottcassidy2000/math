"""Independent walk-cancellation referee; no producer imports.

Prefix denominators are obtained from direct rational vertex ratios,
independently of accumulated edge products. Repeated vertices are retained.
"""
from fractions import Fraction
from itertools import product
from math import gcd,lcm
from functools import reduce
from pathlib import Path
import hashlib,json,sys

sys.stdout.reconfigure(newline='\n')
GATES=0
Q=91**6
H=Q//(42*177)

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)

def cgcd(xs):return reduce(gcd,xs,0)

def packet(walk):
    A=B=1
    for u,v in zip(walk,walk[1:]):
        d=gcd(u,v);A*=u//d;B*=v//d
    C=gcd(A,B)
    denoms=[Fraction(v,walk[0]).denominator for v in walk]
    L=lcm(*denoms);d=cgcd(walk);D=gcd(walk[0],walk[-1])
    need(L==walk[0]//d,'lcm of direct rational denominators is normalized initial coordinate')
    need(A%C==0 and L%(A//C)==0,'reduced endpoint denominator divides clearing lcm')
    J=L//(A//C)
    need(D==d*J,'literal endpoint gcd versus literal whole-walk gcd')
    need(A%L==0 and C%J==0,'exact divisibility J|C via L|A')
    need(Fraction(B,A)==Fraction(walk[-1],walk[0]),'edge-product orientation')
    need(d==cgcd(set(walk)),'repeated positions do not increase subset cardinality')
    return dict(A=A,B=B,C=C,L=L,J=J,collective_gcd=d,endpoint_gcd=D,
                entries=len(walk),distinct_vertices=len(set(walk)))

def atlas_edge(u,v):
    if u==v:return False
    n=(u+v)//gcd(u,v)
    if n>356:return False
    p=2
    while p*p<=n:
        if n%p==0:
            e=0
            while n%p==0:n//=p;e+=1
            if p%3!=2 or e>2:return False
        p+=1
    return n==1 or n%3==2

def main():
    here=Path(__file__).parent
    for name,h in {
      'overnight14_20260906_lrc_endpoint_walk.py':'e38ce143de236a35d089f0e7d4110fae5200ea819868e7b925a705f50428664a',
      'overnight14_20260906_lrc_endpoint_walk.out':'4178d510f86a4cdb4cca86afcb6405439e5179cddd924e48fffd1654a1729561'}.items():
        p=here/name
        if not p.exists() and name.endswith('.out'):p=here.parent/'05-knowledge'/'results'/name
        need(hashlib.sha256(p.read_bytes()).hexdigest()==h,'frozen producer input '+name)
    words=0
    for w in product(range(1,19),repeat=3):
        packet(w);words+=1
    adjacency={u:[v for v in range(1,15) if atlas_edge(u,v)] for u in range(1,15)}
    walks=[(u,) for u in range(1,15)];atlas_count=0;lengths=[]
    for edges in range(5):
        lengths.append(len(walks))
        for w in walks:
            packet(w);atlas_count+=1
        walks=[w+(v,) for w in walks for v in adjacency[w[-1]]]
    controls={}
    for name,w in [
      ('distinct_endpoint_hostile',(6,2,3)),
      ('closed_walk_distinct_endpoint_gate',(6,2,6)),
      ('same_vertices_extra_loop',(6,2,6,2,3)),
      ('three_vertex_strict_gain',(6,4,1)),
      ('five_vertex_strict_gain',(18,12,3,9,6)),
      ('seven_vertex_positive',(729,243,81,27,9,3,1))]:
        need(all(atlas_edge(u,v) for u,v in zip(w,w[1:])),'named control has actual atlas edges')
        result=packet(w)
        controls[name]=result
        for t in (1,2,5):
            scaled=packet(tuple(t*v for v in w))
            need(scaled['J']==result['J'] and scaled['C']==result['C'],'common physical scale preserves ratio packet')
            need(scaled['endpoint_gcd']//t==result['endpoint_gcd'],'physical endpoint gcd must be divided by t')
    need(controls['distinct_endpoint_hostile']['J']==3 and controls['distinct_endpoint_hostile']['collective_gcd']==1,
         'collective gcd alone fails even with distinct maximum endpoint')
    need(controls['closed_walk_distinct_endpoint_gate']['distinct_vertices']==2 and
         controls['closed_walk_distinct_endpoint_gate']['J']==3,'repeated endpoint does not supply a distinct pair')
    need(controls['same_vertices_extra_loop']['J']==3 and controls['same_vertices_extra_loop']['C']==9,
         'repeated old-vertex excursion inflates C while exact J stays fixed')
    need(controls['seven_vertex_positive']['distinct_vertices']==7 and
         controls['seven_vertex_positive']['J']==controls['seven_vertex_positive']['C']==1,
         'monotone seven-vertex positive packet')
    need(controls['three_vertex_strict_gain']['J']==1 and controls['three_vertex_strict_gain']['C']==2,
         'distinct-vertex strict gain is not solely a repeated-walk effect')
    need(controls['five_vertex_strict_gain']['J']==2<=H//11342250<controls['five_vertex_strict_gain']['C']==12,
         'exact five-vertex packet qualifies while content-only packet fails')
    core=(1,2,3,4,5,6,7,9,10,12,18);g=36*Q+1;row=core+(g,3*g)
    unseen=set(range(13));components=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(j for j in unseen-seen if atlas_edge(row[i],row[j]))
        unseen-=seen;components.append(seen)
    need({frozenset(c) for c in components}=={frozenset(range(11)),frozenset((11,12))},
         'strict-gain positive has the actual complete 11+2 decoder graph')
    need(len(set(row))==13 and cgcd(row)==1 and sum(row)<=Q*Q and g>2*Q*max(core),
         'strict-gain control has genuine finite-box equality by full crossing dominance')
    # Distinct endpoints and only two distinct vertices always give J=1,
    # irrespective of repetitions: the whole set equals the endpoint pair.
    for u,v in product(range(1,15),repeat=2):
        if u==v:continue
        p=packet((u,v,u,v,u,v))
        need(p['J']==1,'three distinct vertices are necessary for the distinct-endpoint hostile')
    caps={5:11342250,6:31950,7:90};thresholds={}
    for t in (1,2):
        thresholds[t]=[]
        for r,M in caps.items():
            threshold=t*H//M
            need(M*threshold<=t*H<M*(threshold+1),'exact inclusive threshold and first rejected J')
            thresholds[t].append(threshold)
    need(thresholds=={1:[6,2390,848756],2:[13,4781,1697513]},'all six inherited physical-to-primitive thresholds')
    print('STATUS: INDEPENDENT PASS; exact walk cancellation and conditional endpoint criterion')
    print('UNIVERSES:',words,'arbitrary three-entry words; actual atlas walk counts by edge length 0..4',lengths,'total',atlas_count)
    print('CONTROLS:',json.dumps(controls,sort_keys=True))
    print('THRESHOLDS:',json.dumps(thresholds,sort_keys=True))
    print('SCOPE: count distinct visited vertices; endpoint must differ from global maximum; no universal qualifying-walk claim')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps([controls,thresholds],sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
