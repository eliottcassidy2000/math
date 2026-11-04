#!/usr/bin/env python3
"""Dependency-free exact certificate. Run: python3 verify.py > RESULTS.json"""
import itertools as it, json, math, random
from collections import Counter

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def matchings(xs, planar=True):
    if not xs:
        yield ();return
    a=xs[0]
    for j in range(1,len(xs)):
        b=xs[j]
        for rest in matchings(xs[1:j]+xs[j+1:],False):
            p=tuple(sorted(((a,b),)+rest))
            if not planar or all(not(a<c<b<d or c<a<d<b) for (a,b),(c,d) in it.combinations(p,2)):
                yield p

def states(k, planar=True):
    for d in it.product(range(3),repeat=k):
        ones=tuple(i for i,x in enumerate(d) if x==1)
        if ones and len(ones)%2==0:
            for p in matchings(ones,planar):yield (d,p)

def connected(vertices,edges):
    adj={v:set() for v in vertices}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    seen=set(); todo=[next(iter(adj))]
    while todo:
        v=todo.pop()
        if v in seen:continue
        seen.add(v);todo.extend(adj[v]-seen)
    return len(seen)==len(adj)

def compatible(s,t):
    d,p=s;e,q=t
    if any(x+y!=2 for x,y in zip(d,e)):return False
    ones=[i for i,x in enumerate(d) if x==1]
    return connected(ones,p+q)

def profile(vertices,edges,boundary):
    vertices=tuple(sorted(vertices));edges=tuple(sorted(edges));bs=set(boundary);inter=set(vertices)-bs
    out=Counter()
    for mask in range(1<<len(edges)):
        es=[e for j,e in enumerate(edges) if mask>>j&1]
        deg=Counter(v for e in es for v in e)
        if any(deg[v]!=2 for v in inter) or any(deg[v]>2 for v in bs):continue
        adj={v:[] for v in vertices}
        for a,b in es:adj[a].append(b);adj[b].append(a)
        unseen={v for v in vertices if deg[v]};pair=[];valid=True
        while unseen:
            v=min(unseen);todo=[v];comp=set()
            while todo:
                w=todo.pop()
                if w in comp:continue
                comp.add(w);todo.extend(adj[w])
            unseen-=comp
            ends=sorted(boundary.index(w) for w in comp if deg[w]==1)
            if len(ends)!=2:valid=False;break
            pair.append(tuple(ends))
        if valid and pair:out[(tuple(deg[v] for v in boundary),tuple(sorted(pair)))]+=1
    return out

def probe(s,offset=1000):
    """Realize a unit profile: put degree-2 ports on endpoint collar chains, then subdivide all edges."""
    d,p=s;k=len(d);D=[i for i,x in enumerate(d) if x==1];raw=[];tips={}
    for z,u in enumerate(D):
        v=D[(z+1)%len(D)];chain=[u];j=(u+1)%k
        while j!=v:
            if d[j]==2:chain.append(j)
            j=(j+1)%k
        raw.extend(zip(chain,chain[1:]));tips[u]=chain[-1]
    raw.extend((tips[a],tips[b]) for a,b in p)
    edges=[]
    for j,(a,b) in enumerate(raw):edges.extend([(min(a,offset+j),max(a,offset+j)),(min(b,offset+j),max(b,offset+j))])
    return set(range(k))|set(range(offset,offset+len(raw))),set(edges)

def cycles(vertices,edges):
    adj={v:set() for v in vertices}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    start=min(vertices);ans=0
    def walk(path,left):
        nonlocal ans
        if not left:
            if start in adj[path[-1]] and path[1]<path[-1]:ans+=1
            return
        for v in sorted(adj[path[-1]]&left):walk(path+[v],left-{v})
    if len(vertices)>=3:walk([start],set(vertices)-{start})
    return ans

def dot(a,b):return sum(x*y for s,x in a.items() for t,y in b.items() if compatible(s,t))

def wheel(k,offset,own_boundary=True,stack=1):
    vs=set(range(k))|{offset};es={tuple(sorted((i,offset))) for i in range(k)}
    if own_boundary:es|={tuple(sorted((i,(i+1)%k))) for i in range(k)}
    # Extra vertices lie in the triangular sector bounded by 0,1,offset.
    old=offset
    for j in range(stack):
        w=offset+j+1;vs.add(w)
        es|={tuple(sorted((w,u))) for u in (0,1,old)};old=w
    return vs,es

def rank_mod2(mat):
    rows=[sum((x&1)<<j for j,x in enumerate(row)) for row in mat];r=0
    for c in range(len(mat[0])):
        q=next((i for i in range(r,len(rows)) if rows[i]>>c&1),None)
        if q is None:continue
        rows[r],rows[q]=rows[q],rows[r]
        for j in range(len(rows)):
            if j!=r and rows[j]>>c&1:rows[j]^=rows[r]
        r+=1
    return r

def det(mat):
    from fractions import Fraction
    a=[[Fraction(x) for x in row] for row in mat];out=Fraction(1)
    for j in range(len(a)):
        q=next((i for i in range(j,len(a)) if a[i][j]),None)
        if q is None:return 0
        if q!=j:a[q],a[j]=a[j],a[q];out=-out
        p=a[j][j];out*=p
        for i in range(j+1,len(a)):
            z=a[i][j]/p
            for c in range(j,len(a)):a[i][c]-=z*a[j][c]
    return int(out)

def sq(n):return {v:{w for w in range(1,n+1) if w!=v and math.isqrt(v+w)**2==v+w} for v in range(1,n+1)}

def hp_witness(G):
    n=len(G)
    if n==1:return [next(iter(G))]
    if any(not G[v] for v in G):return None
    leaves=[v for v in G if len(G[v])==1]
    if len(leaves)>2:return None
    starts=leaves[:1] if leaves else sorted(G)
    def rec(path,left):
        if not left:return path
        for w in sorted(G[path[-1]]&left,key=lambda x:len(G[x]&left)):
            ans=rec(path+[w],left-{w})
            if ans:return ans
    for start in starts:
        ans=rec([start],set(G)-{start})
        if ans:return ans
    return None

def main():
    report={'status':'FINITE-EXACT; no external libraries','state_counts':[]}
    allst={k:list(states(k)) for k in range(2,7)}
    probe_count=0
    for k in range(2,7):
        ss=allst[k];rowdegs=[sum(compatible(s,t) for t in ss) for s in ss]
        row={'k':k,'ambient':math.comb(2*k+2,k+1)//(k+2),'active':len(ss),'completion_multiplicities':dict(sorted(Counter(rowdegs).items())),'unrestricted_active':len(list(states(k,False)))}
        need(row['active']==row['ambient']-2**k,'Catalan identity')
        if k<=5:need(set(rowdegs)=={1},'unique mate')
        report['state_counts'].append(row)
        for s in ss:
            vs,es=probe(s)
            need(profile(vs,es,tuple(range(k)))==Counter({s:1}),'unit probe')
            probe_count+=1
    report['unit_profiles_verified']=probe_count
    six=list(matchings(tuple(range(6))));mat=[[int(connected(range(6),a+b)) for b in six] for a in six]
    need(det(mat)==-2 and rank_mod2(mat)==4,'six-port block')
    report['six_port']={'matchings':six,'matrix':mat,'determinant':det(mat),'rank_mod2':rank_mod2(mat)}
    # Direct physical graph tests, independent vertex-order enumerator.
    rng=random.Random(190);tests=[]
    for k in (3,4,5):
        for variant in range(7):
            av,ae=wheel(k,100,True,1);bv,be=wheel(k,200,False,1)
            if variant:
                ae={e for e in sorted(ae) if rng.randrange(5)};be={e for e in sorted(be) if rng.randrange(5)}
            for deleted in [None]+sorted(av|bv):
                va=av-{deleted};vb=bv-{deleted};ea={e for e in ae if deleted not in e};eb={e for e in be if deleted not in e};B=tuple(x for x in range(k) if x!=deleted)
                pa=profile(va,ea,B);pb=profile(vb,eb,B)
                inferred=dot(pa,pb);direct=cycles(va|vb,ea|eb)
                need(inferred==direct,'physical gluing count')
                tests.append({'k':k,'variant':variant,'deleted':deleted,'cycles':direct,'A':{'vertices':sorted(va),'edges':sorted(ea)},'B':{'vertices':sorted(vb),'edges':sorted(eb)}})
    report['physical_gluing_cases']=tests
    # Four-port degree-only hostile: same endpoints and degrees; different pairing.
    pairings=list(matchings(tuple(range(4))));a=((1,)*4,pairings[0]);b=((1,)*4,pairings[1])
    av,ae=probe(a,100);bv,be=probe(a,200);cv,ce=probe(b,300)
    need(cycles(av|bv,ae|be)==0 and cycles(av|cv,ae|ce)==1,'degree hostile')
    report['degree_only_hostile']={'matching_A':a,'matching_B':b,'same_matchings_cycle_count':0,'different_matchings_cycle_count':1}
    # Low-degree witness for every interval graph before the known cycle threshold.
    report['square_low_degree']=[{'n':n,'vertex':min(sq(n),key=lambda v:(len(sq(n)[v]),v)),'minimum_degree':min(map(len,sq(n).values()))} for n in range(2,32)]
    need(all(r['minimum_degree']<=2 for r in report['square_low_degree']),'square degree bound')
    deletion=[]
    for n in (18,19,20,21,22,24):
        G=sq(n);need(hp_witness(G) is None,'claimed failure')
        yes={};no=[]
        for v in G:
            H={w:nb-{v} for w,nb in G.items() if w!=v};p=hp_witness(H)
            if p:
                need(len(p)==n-1 and len(set(p))==n-1 and all(b in H[a] for a,b in zip(p,p[1:])),'path witness')
                yes[v]=p
            else:no.append(v)
        deletion.append({'n':n,'restoring_deletions':yes,'nonrestoring_deletions':no})
    report['square_deletion_profiles']=deletion
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()
