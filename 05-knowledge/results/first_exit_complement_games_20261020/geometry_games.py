#!/usr/bin/env python3
"""Exact complementary games and a nineteen-point hinged boundary gadget."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
    c=[F(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i^j]+=x*y*(-3 if i&j&1 else 1)*(-11 if i&j&2 else 1)
    return tuple(c)
def norm(a):return mul(a,(a[0],-a[1],-a[2],a[3]))

def colouring(n,edges,q):
    adj=[set() for _ in range(n)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    cs=[-1]*n
    def go():
        free=[i for i in range(n) if cs[i]<0]
        if not free:return cs[:]
        i=max(free,key=lambda v:(len({cs[u] for u in adj[v] if cs[u]>=0}),len(adj[v])))
        bad={cs[j] for j in adj[i]}
        for c in range(q):
            if c not in bad:
                cs[i]=c
                ans=go()
                if ans is not None:return ans
        cs[i]=-1
        return None
    return go()

def geometry():
    one=(F(1),F(0),F(0),F(0));zero=tuple(F(0) for _ in range(4))
    z=(F(1,2),F(1,2),F(0),F(0));r=[one]
    for _ in range(5):r.append(mul(r[-1],z))
    t=[add(r[2*i],r[2*i+1]) for i in range(3)]
    v=[add(x,neg(t[0])) for x in [zero]+r+t]
    rho=(F(5,6),F(0),F(1,6),F(0));q=[mul(rho,x) for x in v]
    need(norm(rho)==one,'rotation isometry')
    pts=v+q[:7]+q[8:]
    need(len(set(pts))==19,'distinct points')
    edges=[(i,j) for i,j in combinations(range(19),2) if norm(add(pts[i],neg(pts[j])))==one]
    cross=[(i,j) for i in range(10) for j in range(10)
           if v[i]!=q[j] and norm(add(v[i],neg(q[j])))==one]
    need(cross==[(0,0),(1,7),(2,7),(7,1),(7,2)],'complete cross ledger')
    need(len(edges)==37,'induced edge count')
    need(colouring(19,edges,3) is None,'three colour hostile')
    c4=colouring(19,edges,4);need(c4 is not None,'four colour witness')
    # The common tip is colour0. Count the remaining tips after the rim.
    freq=Counter();examples={}
    for rim in product(range(6),repeat=6):
        if any(rim[i]==rim[(i+1)%6] for i in range(6)):continue
        if 0 in rim[:2]:continue
        mask=sum(1<<c for c in range(6) if c not in rim)
        freq[mask]+=16
        examples.setdefault(mask,rim)
    need(sum(freq.values())==166720,'single-wing boundary count')
    hist=Counter();newblocked=0;localblocked=0
    for a,fa in freq.items():
        for b,fb in freq.items():
            k=a.bit_count()*b.bit_count()-(a&b).bit_count()
            # Independent literal colour-pair count.
            lit=sum(1 for x in range(6) for y in range(6)
                    if x!=y and a>>x&1 and b>>y&1)
            need(k==lit,'missing-palette kernel')
            hist[k]+=fa*fb
            if not a or not b:localblocked+=fa*fb
            elif not k:newblocked+=fa*fb
    need(all((1<<c) in freq for c in range(6)),'rank-six minor')
    return {'vertices':19,'induced_edges':edges,'cross_contacts_copy_labels':cross,
            'four_colouring':c4,'six_colour_boundary_count_fixed_common_tip':sum(hist.values()),
            'newly_blocked_despite_both_local_extensions':newblocked,
            'already_blocked_locally':localblocked,'extension_histogram':dict(sorted(hist.items())),
            'missing_masks_frequency':dict(sorted(freq.items())),
            'example_rims_missing_colour0':examples[1],
            'example_rims_missing_colour1':examples[2]}

def graph_info(n,mask,es):
    adj=[0]*n
    for k,(i,j) in enumerate(es):
        if mask>>k&1:adj[i]|=1<<j;adj[j]|=1<<i
    seen=1;old=0
    while old!=seen:
        old=seen
        for v in range(n):
            if seen>>v&1:seen|=adj[v]
    connected=seen==(1<<n)-1
    dp=[0]*(1<<n)
    for i in range(n):dp[1<<i]=1<<i
    for s in range(1,1<<n):
        for i in range(n):
            if dp[s]>>i&1:
                nxt=adj[i]&~s
                while nxt:
                    b=nxt&-nxt;nxt-=b;dp[s|b]|=b
    return bool(dp[-1]),connected

def games():
    rows=[]
    for n in range(1,7):
        es=list(combinations(range(n),2));full=(1<<len(es))-1
        info=[graph_info(n,m,es) for m in range(full+1)]
        counts=[[0,0],[0,0]];mins=[[None,None],[None,None]];minimal=[]
        for mask in range(full+1):
            for p in range(2):
                for q in range(2):
                    if not info[mask][p] and not info[full^mask][q]:
                        counts[p][q]+=1
                        e=mask.bit_count()
                        if mins[p][q] is None or e<mins[p][q]:mins[p][q]=e
            if n>=4 and mask.bit_count()==n-1 and not info[mask][0] and not info[full^mask][0]:
                deg=[sum(bool(mask>>k&1) for k,e in enumerate(es) if v in e) for v in range(n)]
                minimal.append(sorted(deg))
        need(counts[1][1]==0,'disconnected complement hostile')
        if n>=4:need(mins[0][0]==n-1 and mins[0][1]==n-1,'sharp deletion minimum')
        if n>=5:need(all(ds==[1]*(n-1)+[n-1] for ds in minimal),'unique star optimum')
        if n<=3:need(not any(sum(row) for row in counts),'small-order minimality')
        rows.append({'n':n,'graph_count':full+1,'winning_counts_H_W':counts,'minimum_edges_H_W':mins})
    connected_hostiles=[]
    for k in range(3,9):
        n=k+4;es=list(combinations(range(n),2))
        # u,v,w = 0,1,2; z=3; k independent twins=4,...,n-1.
        chosen={(0,1),(0,2),(1,2),(2,3)}|{(i,j) for i in (0,1) for j in range(4,n)}
        mask=sum(1<<j for j,e in enumerate(es) if e in chosen);full=(1<<len(es))-1
        need(graph_info(n,mask,es)==(False,True),'connected red obstruction')
        need(graph_info(n,full^mask,es)==(False,True),'connected blue obstruction')
        degrees=[sum(i in e for e in chosen) for i in range(n)]
        need(sum(d%2 for d in degrees)==(2 if k%2==0 else 4),'Euler parity')
        connected_hostiles.append({'n':n,'red_edges':sorted(chosen),'red_odd_degrees':sum(d%2 for d in degrees)})
    # Directed Euler edge-words equal Hamilton vertex-words of the line digraph.
    line_controls=0
    for n in range(2,4):
        arcs=[(a,b) for a in range(n) for b in range(n) if a!=b]
        for mask in range(1,1<<len(arcs)):
            chosen=[a for k,a in enumerate(arcs) if mask>>k&1]
            euler=sum(all(p[i][1]==p[i+1][0] for i in range(len(p)-1)) for p in permutations(chosen))
            line=sum(all(chosen[p[i]][1]==chosen[p[i+1]][0] for i in range(len(p)-1))
                     for p in permutations(range(len(chosen))))
            need(euler==line,'directed line-graph bijection');line_controls+=1
    torus_controls=0
    for p in range(2,15):
        for q in range(1,25):
            unseen=set(range(p));lengths=[]
            while unseen:
                x=min(unseen);start=x;k=0
                while x in unseen:
                    unseen.remove(x);k+=1;x=(x+q)%p
                need(x==start,'permutation cycles close');lengths.append(k)
            need(len(lengths)==gcd(p,q) and all(k==p//gcd(p,q) for k in lengths),
                 'torus component count')
            torus_controls+=1
    return {'universe':'all labelled simple graphs on1..6 vertices; H=Hamilton path,W=spanning walk',
            'rows':rows,'both_connected_family':connected_hostiles,'directed_line_graph_controls':line_controls,
            'torus_step_controls':torus_controls}

if __name__=='__main__':
    print(json.dumps({'status':'FINITE-EXACT geometry and graph controls',
                      'geometry':geometry(),'games':games()},indent=2,sort_keys=True))
