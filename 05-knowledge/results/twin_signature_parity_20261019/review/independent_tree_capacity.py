from itertools import permutations,product
from math import factorial,comb
from collections import defaultdict
import json

checks=0
def check(x,msg):
    global checks
    checks+=1
    if not x:raise RuntimeError(msg)

def adjacency(n,edges):
    a=[[] for _ in range(n)]
    for u,v in edges:a[u].append(v);a[v].append(u)
    return a

def rooted_shape(adj,v,parent=-1):
    return '('+''.join(sorted(rooted_shape(adj,w,v) for w in adj[v] if w!=parent))+')'

def canon(n,edges):
    a=adjacency(n,edges)
    return min(rooted_shape(a,v) for v in range(n))

def trees_to(nmax):
    rows={1:{'()':()}}
    for n in range(2,nmax+1):
        row={}
        for ed in rows[n-1].values():
            for p in range(n-1):
                ee=ed+((p,n-1),);row[canon(n,ee)]=ee
        rows[n]=row
    return rows

def automorphisms(n,edges):
    a=adjacency(n,edges);classes=defaultdict(list)
    for v in range(n):classes[len(a[v])].append(v)
    classes=list(classes.values());es={tuple(sorted(e)) for e in edges};out=[]
    for blocks in product(*(list(permutations(cls)) for cls in classes)):
        p=list(range(n))
        for cls,im in zip(classes,blocks):
            for x,y in zip(cls,im):p[x]=y
        if all(tuple(sorted((p[x],p[y]))) in es for x,y in edges):out.append(tuple(p))
    return out

def fixed_masks(n,edges,autos):
    ed=sorted(tuple(sorted(e)) for e in edges);lookup={e:i for i,e in enumerate(ed)}
    out=[]
    for p in autos:
        maps=[lookup[tuple(sorted((p[x],p[y])))] for x,y in ed]
        flips=sum((p[x]>p[y])<<maps[i] for i,(x,y) in enumerate(ed))
        fix=[]
        for mask in range(1<<(n-1)):
            image=flips
            for i,j in enumerate(maps):image^=((mask>>i)&1)<<j
            if image==mask:fix.append(mask)
        out.append(set(fix))
    return out

def capacity(adj,v,parent=-1):
    groups=defaultdict(list)
    for w in adj[v]:
        if w!=parent:groups[rooted_shape(adj,w,v)].append(w)
    ans=1
    for ws in groups.values():
        choices=2*capacity(adj,ws[0],v)
        ans*=comb(choices,len(ws)) if choices>=len(ws) else 0
    return ans

def centers(adj):
    # Use minimum eccentricity rather than leaf stripping.
    ecc=[]
    for v in range(len(adj)):
        ds={v:0};todo=[v]
        for x in todo:
            for y in adj[x]:
                if y not in ds:ds[y]=ds[x]+1;todo.append(y)
        ecc.append(max(ds.values()))
    return [i for i,x in enumerate(ecc) if x==min(ecc)]

def verify_tree(n,edges):
    adj=adjacency(n,edges);aut=automorphisms(n,edges);fix=fixed_masks(n,edges,aut)
    nonid=[f for p,f in zip(aut,fix) if p!=tuple(range(n))]
    asym=sum(not any(mask in f for f in nonid) for mask in range(1<<(n-1)))
    check(asym%len(aut)==0,'free orientation action')
    exact=asym//len(aut)
    for v in range(n):
        ids=[i for i,p in enumerate(aut) if p[v]==v]
        bad=[fix[i] for i in ids if aut[i]!=tuple(range(n))]
        rooted_asym=sum(not any(mask in f for f in bad) for mask in range(1<<(n-1)))
        check(rooted_asym%len(ids)==0,'rooted orbit divisibility')
        check(rooted_asym//len(ids)==capacity(adj,v),'rooted capacity against direct permutations')
    cc=centers(adj)
    if len(cc)==1:pred=capacity(adj,cc[0])
    else:
        x,y=cc
        pred=capacity(adj,x,y)*capacity(adj,y,x)
        if rooted_shape(adj,x,y)!=rooted_shape(adj,y,x):pred*=2
    check(pred==exact,'unrooted center factor')
    # Every nontrivial oriented automorphism group must have even order.
    for mask in range(1<<(n-1)):
        size=sum(mask in f for f in fix)
        check(size==1 or size%2==0,'oriented tree symmetry has even order')
    return {'automorphisms':len(aut),'asymmetric_oriented_classes':exact}

def valuation(x):
    ans=0
    while x:x//=2;ans+=x
    return ans

def main():
    trees=trees_to(8);expected=[1,1,1,2,3,6,11,23];rows=[]
    for n,tt in trees.items():
        check(len(tt)==expected[n-1],'unlabeled tree count')
        results=[verify_tree(n,e) for e in tt.values()]
        rows.append({'n':n,'trees':len(tt),'all_orientations_automorphic_trees':sum(r['asymmetric_oriented_classes']==0 for r in results)})
    witness=((0,1),(0,2),(0,3),(0,4),(0,5),(1,6),(2,7),(6,8))
    ww=verify_tree(9,witness)
    check(ww=={'automorphisms':6,'asymmetric_oriented_classes':0},'nine witness')
    for k in range(1,9):
        dp=[0]+[10**9]*64
        for _ in range(k):
            dp=[min(dp[s-t]+valuation(t) for t in range(s+1)) for s in range(65)]
        for m,got in enumerate(dp):
            pred=0;j=1
            while (1<<j)<=m+k:
                pred+=max(0,(m+k+(1<<j)-1)//(1<<j)-k);j+=1
            check(got==pred,'independent factorial-budget optimization')
    print(json.dumps({'status':'PASS','checks':checks,'tree_census':rows,'nine_vertex_witness':ww,'factorial_budget_k_range':[1,8],'factorial_budget_m_range':[0,64],'scope':'Independent permutation-group verification of every orientation on every unlabeled tree through eight vertices and the specific nine-vertex witness; no embedding-parity minimality claim from capacity alone.'},indent=2))

if __name__=='__main__':main()
