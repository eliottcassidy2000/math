#!/usr/bin/env python3
"""Exact finite controls. No assert statements; normal and -O are identical."""
import itertools, json, math
from collections import Counter
from functools import lru_cache
from pathlib import Path

checks=Counter()
def require(ok,label):
    checks[label]+=1
    if not ok: raise RuntimeError(label)

def component_class(n,edges):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for a,b in edges:
        x,y=find(a),find(b)
        if x!=y:parent[x]=y
    nv=Counter(find(v) for v in range(n));ne=Counter(find(a) for a,b in edges)
    return all(ne[c]<v for c,v in nv.items()),all(ne[c]<=v for c,v in nv.items())

def independence(n,edges):
    adj=[0]*n
    for a,b in edges:adj[a]|=1<<b;adj[b]|=1<<a
    @lru_cache(None)
    def solve(S):
        if not S:return 0
        v=max((i for i in range(n) if S>>i&1),key=lambda i:(adj[i]&S).bit_count())
        rest=S&~(1<<v)
        return max(solve(rest),1+solve(rest&~adj[v]))
    return solve((1<<n)-1)

def good(n,red,blue):
    colors=[set(red),set(blue)]
    return all(not all(tuple(sorted(e)) in E for e in itertools.combinations(t,2)) for E in colors for t in itertools.combinations(range(n),3))

def capacities(n,red,blue):
    rf,rp=component_class(n,red);bf,bp=component_class(n,blue)
    hf,hp=component_class(n,list(red)+list(blue));a=independence(n,list(red)+list(blue))
    for applies,c,name in [(hf,2,'union_forest'),(hp,3,'union_pseudoforest'),((rf and bp) or (bf and rp),4,'mixed'),(rp and bp,5,'separate_pseudoforests')]:
        if applies:require(n<=c*a,name)
    return a,(hf,hp,rf,bf,rp,bp)

census=[]
for n in range(1,6):
    pairs=list(itertools.combinations(range(n),2));valid=0
    for code in itertools.product(range(3),repeat=len(pairs)):
        red=[e for e,c in zip(pairs,code) if c==0];blue=[e for e,c in zip(pairs,code) if c==1]
        if good(n,red,blue):
            valid+=1;capacities(n,red,blue)
    census.append({'n':n,'all_colorings':3**len(pairs),'no_red_or_blue_triangle':valid})

blocks={2:([(0,1)],[]),3:([(0,1),(1,2)],[(0,2)]),4:([(0,1),(1,2),(2,3)],[(0,2),(0,3),(1,3)]),5:([(0,1),(1,2),(2,3),(3,4),(0,4)],[(0,2),(0,3),(1,3),(1,4),(2,4)])}
sharp=[]
for c,(R,B) in blocks.items():
    for k in range(2,8):
        n=c*(k-1);red=[(a+c*j,b+c*j) for j in range(k-1) for a,b in R];blue=[(a+c*j,b+c*j) for j in range(k-1) for a,b in B]
        require(good(n,red,blue),'sharp_triangle_control')
        alpha,classes=capacities(n,red,blue);require(alpha==k-1,'sharp_alpha')
        require(classes[{2:0,3:1,4:2,5:4}[c]],'sharp_membership')
        if c==4:require(classes[3],'K4_second_forest')
        if c==5:require(classes[5],'K5_second_pseudoforest')
        sharp.append({'factor':c,'k':k,'order':n,'alpha':alpha})
for c in (3,4,5):
    R,B=blocks[c];require(c>(c-1)*independence(c,R+B),'lower_factor_hostile')

# Matchings and gluing: multiedges are intentionally retained.
def matchings(labels):
    if not labels:yield ();return
    a=labels[0]
    for i in range(1,len(labels)):
        b=labels[i]
        for rest in matchings(labels[1:i]+labels[i+1:]):yield ((a,b),)+rest

def loops(P,Q):
    parent=list(range(2*len(P)))
    def find(x):
        if parent[x]!=x:parent[x]=find(parent[x])
        return parent[x]
    for a,b in P+Q:parent[find(a)]=find(b)
    return len({find(v) for v in range(len(parent))})

def loops_permutation(P,Q):
    p={a:b for a,b in P for a,b in ((a,b),(b,a))};q={a:b for a,b in Q for a,b in ((a,b),(b,a))};unseen=set(p);cycles=0
    while unseen:
        x=next(iter(unseen));cycles+=1
        while x in unseen:unseen.remove(x);x=p[q[x]]
    require(cycles%2==0,'permutation_even_cycles')
    return cycles//2

loop_rows=[]
poly=[1]
for r in range(1,7):
    old=poly;poly=[0]*(len(old)+1)
    for i,v in enumerate(old):poly[i]+=2*(r-1)*v;poly[i+1]+=v
    Ms=list(matchings(tuple(range(2*r))));P=Ms[0]
    count=Counter(loops(P,Q) for Q in Ms)
    require([count[i] for i in range(r+1)]==poly,'matching_polynomial')
    require(count[1]==2**(r-1)*math.factorial(r-1),'one_loop_count')
    require(sum(count.values())==math.prod(range(1,2*r,2)),'double_factorial')
    if r<=4:
        for P in Ms:
            for Q in Ms:
                l=loops(P,Q);require(l==loops_permutation(P,Q),'two_loop_algorithms')
                require((l==r)==(P==Q),'matching_separation')
    loop_rows.append({'r':r,'coefficients_constant_first':poly[:],'matchings':len(Ms),'single_loop':count[1]})

# Fox color operations: the three quandle identities are the local move checks.
def star(a,b,p=3):return (2*b-a)%p
for a,b,c in itertools.product(range(3),repeat=3):
    require(star(a,a)==a,'Fox_R1');require(star(star(a,b),b)==a,'Fox_R2')
    require(star(star(a,b),c)==star(star(a,c),star(b,c)),'Fox_R3')
fox={}
for m in [0,1,2,3,4,5,6]:
    fixed=[]
    for a,b in itertools.product(range(3),repeat=2):
        x,y=a,b
        for _ in range(m):x,y=(2*x-y)%3,x
        if (x,y)==(a,b):fixed.append([a,b])
    fox[m]=len(fixed)
require(fox[1]==3 and fox[3]==9,'same_endpoint_permutation_distinct_knots')

# Recount the corrected 1729 endpoint channels, not just an old scalar.
q=11;QR={i*i%q for i in range(1,q)}
@lru_cache(None)
def ends(S,v):
    if S==(1<<0):return int(v==0)
    if v==0 or not(S>>v&1):return 0
    smaller=S&~(1<<v)
    return sum(ends(smaller,u) for u in range(q) if smaller>>u&1 and (v-u)%q in QR)
E=[ends((1<<q)-1,v) for v in range(q)]
require(set(E[v] for v in QR)=={628},'Paley_forward_channel')
require(set(E[v] for v in range(1,q) if v not in QR)=={1101},'Paley_reverse_channel')
H=q*sum(E);C=sum(E[v] for v in range(1,q) if (-v)%q in QR)
require((H,C)==(95095,5505),'Paley_H_C')
# Required arc count with a forward subset DP retaining its use bit.
D={(1<<s,s,False):1 for s in range(q)}
for length in range(1,q):
    N={}
    for (S,u,used),count in D.items():
        for v in range(q):
            if not(S>>v&1) and (v-u)%q in QR:
                key=(S|(1<<v),v,used or (u==0 and v==1));N[key]=N.get(key,0)+count
    D=N
arc_count=sum(v for (S,u,used),v in D.items() if used)
require(arc_count==17290,'Paley_arc_incidence_not_1729')
require(55*arc_count==(q-1)*H,'Paley_incidence_handshake')

out={'status':'PASS','census':census,'sharp_constructions':sharp,'loop_rows':loop_rows,'fox_3_coloring_counts':fox,'paley11':{'endpoint_row':E,'a':628,'b':1101,'s':1729,'H':H,'cycles':C,'specified_arc_incidence':arc_count},'checks':dict(sorted(checks.items()))}
print(json.dumps(out,indent=2,sort_keys=True))
