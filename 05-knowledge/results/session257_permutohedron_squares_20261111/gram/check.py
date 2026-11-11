#!/usr/bin/env python3
"""Exact finite controls for Session257 endpoint Gram and color-census bridge."""
import itertools as it, json, math
from pathlib import Path

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def sumgraph(n,targets=None):
    targets=set(targets or [k*k for k in range(2,math.isqrt(2*n)+1)])
    return {i:{t-i for t in targets if 1<=t-i<=n and t-i!=i} for i in range(1,n+1)}

def signature(g):
    """Unoriented spanning paths by unordered pair of distinct endpoints."""
    ans={}
    if len(g)<2:return ans
    n=len(g); full=set(g)
    def dfs(start,v,seen):
        if len(seen)==n:
            if start<v:ans[(start,v)]=ans.get((start,v),0)+1
            return
        for w in g[v]-seen:dfs(start,w,seen|{w})
    for start in sorted(g):dfs(start,start,{start})
    return ans

def dot(s,t):return sum(v*t.get(k,0) for k,v in s.items())
def cycle_count(g,seam):
    start=min(g); n=len(g); hist={}
    def dfs(path,seen):
        v=path[-1]
        if len(seen)==n:
            if start in g[v] and path[1]<path[-1]:
                e=list(zip(path,path[1:]+[start])); k=sum(frozenset(x) in seam for x in e)
                hist[k]=hist.get(k,0)+1
            return
        for w in g[v]-seen:dfs(path+[w],seen|{w})
    dfs([start],{start});return hist

def glue(g,h):
    # Node integers encode side and original label, retaining common port labels.
    G={2*v:set(2*w for w in ns) for v,ns in g.items()}
    G.update({2*v+1:set(2*w+1 for w in ns) for v,ns in h.items()})
    seam=set()
    for v in g.keys()&h.keys():
        G[2*v].add(2*v+1);G[2*v+1].add(2*v);seam.add(frozenset((2*v,2*v+1)))
    return G,seam

def reflection_path(s,t):
    A,B=s*s,t*t; g=sumgraph(B-1,[A,B,A+B]); ends=[v for v in g if len(g[v])==1]
    check(len(ends)==2 and all(len(g[v]) in [1,2] for v in g), "Reflection graph degree condition failed")
    p=[min(ends)]; prev=None
    while True:
        nxt=g[p[-1]]-({prev} if prev is not None else set())
        if not nxt:break
        prev=p[-1];p.append(next(iter(nxt)))
    check(len(p)==B-1 and len(set(p))==B-1, "Reflection path is not spanning and simple")
    colors=[(A,B,A+B).index(a+b) for a,b in zip(p,p[1:])]
    return p,colors

def decode(c):
    a,b,c=c; e=b-a-c
    check(e in [0,1], "Invalid census parity coordinate")
    A,B=2*a+1+e,2*b+2-e
    s,t,u=map(math.isqrt,(A,B,A+B))
    check(s*s==A and t*t==B and u*u==A+B and math.gcd(s,t)==1, "Census does not decode a primitive triple")
    odd=s if s%2 else t
    m,n=math.isqrt((u+odd)//2),math.isqrt((u-odd)//2)
    check(m*m+n*n==u and m*m-n*n==odd and m>n>0, "Spinor recovery failed")
    return s,t,u,m,n

def embed(colors,targets,n):
    witnesses=[]
    for x in range(1,n+1):
        p=[x]
        for c in colors:p.append(targets[c]-p[-1])
        if set(p)==set(range(1,n+1)) and len(p)==n:witnesses.append(p)
    return witnesses

def tournament_h(m,n):
    # Sorted-spinor carrier independent subset-DP.
    x,y=m-n,n; size=x+y+1; adj=[set() for _ in range(size)]
    for i in range(size-1):
        for j in range(i+1,size-1):adj[i].add(j)
        if i<x:adj[size-1].add(i)
        else:adj[i].add(size-1)
    dp={ (1<<i,i):1 for i in range(size)}
    for mask in range(1,1<<size):
        for v in range(size):
            q=dp.get((mask,v),0)
            if q:
                for w in adj[v]:
                    if not mask>>w&1:dp[(mask|1<<w,w)]=dp.get((mask|1<<w,w),0)+q
    return sum(dp.get(((1<<size)-1,v),0) for v in range(size))

def main():
    receipt={"status":"FINITE-EXACT controls; infinite statements proved in PROOF.md"}
    qs={n:signature(sumgraph(n)) for n in [12,15,16,17,18,23,24,25]}
    receipt['square_endpoint_signatures']={str(n):{f'{a},{b}':v for (a,b),v in s.items()} for n,s in qs.items()}
    receipt['square_gram_orders']=list(qs)
    receipt['square_gram']=[[dot(s,t) for t in qs.values()] for s in qs.values()]
    # All simple graphs on 3 labelled ports: independent cycle recount.
    graphs=[]
    es=list(it.combinations(range(1,4),2))
    for mask in range(8):
        g={v:set() for v in range(1,4)}
        for k,(a,b) in enumerate(es):
            if mask>>k&1:g[a].add(b);g[b].add(a)
        graphs.append(g)
    for g,h in it.product(graphs,repeat=2):
        G,seam=glue(g,h)
        check(cycle_count(G,seam).get(2,0)==dot(signature(g),signature(h)), "Two-seam independent recount mismatch")
    receipt['two_seam_checks']=64
    # Four-seam matching block: no Hamiltonian paths inside either half.
    matchings=[[(1,2),(3,4)],[(1,3),(2,4)],[(1,4),(2,3)]]
    block=[]
    for a in matchings:
        row=[]
        for b in matchings:
            gs=[]
            for mm in [a,b]:
                g={v:set() for v in range(1,5)}
                for u,v in mm:g[u].add(v);g[v].add(u)
                gs.append(g)
            G,seam=glue(*gs);row.append(cycle_count(G,seam).get(4,0))
        block.append(row)
    check(block==[[0,1,1],[1,0,1],[1,1,0]], "Four-seam matrix mismatch")
    arithmetic_ports=[[1,3,7,9],[1,7,3,9],[1,7,9,3]]
    for portmap,mm in zip(arithmetic_ports,matchings):
        actual=[(i+1,j+1) for i,j in it.combinations(range(4),2) if math.isqrt(portmap[i]+portmap[j])**2==portmap[i]+portmap[j]]
        check(actual==mm, "Square embedding matching mismatch")
    receipt['four_seam_square_embeddings']=arithmetic_ports
    receipt['four_seam_matrix']=block;receipt['negative_test_vector']=[1,-1,0];receipt['negative_quadratic_value']=-2
    triples=[]
    for m in range(2,13):
        for n in range(1,m):
            if math.gcd(m,n)>1 or (m-n)%2==0:continue
            s,t=sorted((m*m-n*n,2*m*n));p,col=reflection_path(s,t)
            c=[col.count(k) for k in range(3)]
            check(decode(c)==(s,t,m*m+n*n,m,n), "Color census decoder round trip failed")
            if t<=15: check(len(embed(col,[s*s,t*t,s*s+t*t],t*t-1))==1, "Small reflection itinerary embedding mismatch")
            H=2+(2**(m-n)-1)*(2**n-1)
            if m<=10:check(tournament_h(m,n)==H, "Independent tournament DP mismatch")
            triples.append({'spinor':[m,n],'triple':[s,t,m*m+n*n],'vertices':t*t-1,'census':c,'tournament_H':H})
    receipt['primitive_spinor_range']='2 <= m <= 12, 1 <= n < m, gcd=1 and opposite parity'
    receipt['primitive_controls']=triples
    p,col=reflection_path(3,4);bad=None
    for i,j in it.combinations(range(len(col)),2):
        if col[i]==col[j]:continue
        cc=col[:];cc[i],cc[j]=cc[j],cc[i]
        if not embed(cc,[9,16,25],15):bad={'swap_positions_zero_based':[i,j],'colors':cc};break
    check(bad is not None, "Equal-census hostile not found")
    receipt['census_not_itinerary_hostile']={'valid_path':p,'valid_colors':col,'invalid_equal_census':bad}
    # Keep the independent color roles explicitly distinct from user's 2026-10-23 edge partition.
    receipt['color_scope']='Colors here mean ordered square targets, not the recovered user three-color decomposition.'
    Path(__file__).with_name('results.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'checks':'PASS','primitive_controls':len(triples),'square_path_totals':{n:sum(s.values()) for n,s in qs.items()},'four_seam_matrix':block},indent=2))
if __name__=='__main__':main()
