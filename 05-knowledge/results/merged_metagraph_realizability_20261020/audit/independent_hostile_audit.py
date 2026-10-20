#!/usr/bin/env python3
"""Independent hostile check by full labelled-orbit closure under adjacent swaps.
No imports from the packet under audit. All integer computations are exact.
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def orbit_partition(n, merge_converse):
    pairs=list(combinations(range(n),2)); ix={p:k for k,p in enumerate(pairs)}
    size=1<<len(pairs); parent=list(range(size))
    def root(a):
        while a!=parent[a]:
            parent[a]=parent[parent[a]]; a=parent[a]
        return a
    def union(a,b):
        a,b=root(a),root(b)
        if a!=b: parent[max(a,b)]=min(a,b)
    for c in range(size):
        for v in range(n-1):
            p=list(range(n)); p[v],p[v+1]=p[v+1],p[v]
            d=0
            for bit,(i,j) in enumerate(pairs):
                a,b=p[i],p[j]; z=(c>>ix[min(a,b),max(a,b)])&1
                d|=(z^(a>b))<<bit
            union(c,d)
        if merge_converse: union(c,c^(size-1))
    return [root(c) for c in range(size)]

def arc(c,a,b,pairs):
    ix={p:k for k,p in enumerate(pairs)}
    return ((c>>ix[min(a,b),max(a,b)])&1)^(a>b)

def main():
    n=5; pairs=list(combinations(range(n),2)); free=[p for p in pairs if p[1]>p[0]+1]
    ix={p:k for k,p in enumerate(free)}; N=1<<len(free)
    ordinary=orbit_partition(n,False); merged=orbit_partition(n,True)
    assert len(set(ordinary))==12 and len(set(merged))==10
    codes=[sum((1 if b==a+1 else (x>>ix[a,b])&1)<<k for k,(a,b) in enumerate(pairs)) for x in range(N)]
    good=[merged[c] for c in codes]; bad=good.copy(); bad[4],bad[26]=bad[26],bad[4]
    reps=sorted(set(good)); reflection=[sum(((x>>k)&1)<<ix[n-1-b,n-1-a] for k,(a,b) in enumerate(free)) for x in range(N)]
    assert reps==[0,2,4,5,8,10,12,40,41,76]
    assert good[4]==good[9]==12 and good[26]==41
    mass=Counter(good); assert mass==Counter(bad)
    assert [mass[u] for u in reps]==[1,2,10,1,9,9,11,5,13,3]
    def matrices(q):
        b=Counter(); k=Counter()
        for x,u in enumerate(q):
            (b if reflection[x]==x else k)[u,q[x^(N-1)]]+=1
        return b,k
    B,K=matrices(good); Bb,Kb=matrices(bad); assert K==Kb
    for u in reps:
        assert sum(B[u,v] for v in reps)==sum(Bb[u,v] for v in reps)
        assert B[u,u]==Bb[u,u]
    def profile(c): return tuple(sorted(sum(arc(c,a,b,pairs) for b in range(n) if b!=a) for a in range(n)))
    def fold(s):return min(tuple(sorted(s)),tuple(sorted(n-1-a for a in s)))
    def c3(s):return 10-sum(a*(a-1)//2 for a in s)
    def rows(mat):
        ans={u:Counter() for u in reps}
        for (u,v),w in mat.items():
            if w:ans[u][c3(profile(v))-c3(profile(u)),fold(profile(v))]+=w
        return ans
    assert rows(B)==rows(Bb) and rows(K)==rows(Kb)
    # Independently reconstruct HP/Aut endpoint inventories from each representative.
    hp_rows=[{},{}]
    for u in reps:
        s=[sum(arc(u,a,b,pairs) for b in range(n) if b!=a) for a in range(n)]
        aut=0; anti=0; raw=[Counter(),Counter()]
        for p in permutations(range(n)):
            transformed=sum(arc(u,p[a],p[b],pairs)<<k for k,(a,b) in enumerate(pairs))
            aut+=transformed==u; anti+=transformed==(u^1023)
            if not all(arc(u,p[i],p[i+1],pairs) for i in range(n-1)):continue
            tile=sum(arc(u,p[a],p[b],pairs)<<k for k,(a,b) in enumerate(free))
            target=[n-1-a for a in s]; target[p[0]]+=1; target[p[-1]]-=1
            raw[0 if reflection[tile]==tile else 1][s[p[0]]-s[p[-1]]-1,fold(target)]+=1
        multiplier=1 if anti else 2
        for col in (0,1):
            assert all(v*multiplier%aut==0 for v in raw[col].values())
            hp_rows[col][u]=Counter({k:v*multiplier//aut for k,v in raw[col].items()})
    assert hp_rows==[rows(Bb),rows(Kb)]
    # No graph isomorphism survives, even allowing arbitrary node relabelling.
    groups={m:[u for u in reps if mass[u]==m] for m in set(mass.values())}
    num=0; isos=0
    for image in product(*(list(permutations(v)) for v in groups.values())):
        p={u:v for us,vs in zip(groups.values(),image) for u,v in zip(us,vs)}; num+=1
        isos+=all(B[u,v]==Bb[p[u],p[v]] and K[u,v]==Kb[p[u],p[v]] for u in reps for v in reps)
    assert num==4 and isos==0
    # Exact Walsh coefficients by independent in-place butterfly transform.
    charges=[]
    for q in (good,bad):
        z={}
        for u in reps:
            a=[int(v==u) for v in q]; width=1
            while width<N:
                for start in range(0,N,2*width):
                    for j in range(start,start+width):
                        s,t=a[j],a[j+width]; a[j],a[j+width]=s+t,s-t
                width*=2
            z[u]=a
            assert all(sum(a[s]*a[s^t] for s in range(N))==N*a[t] for t in range(N))
        assert all(sum(z[u][s] for u in reps)==(N if s==0 else 0) for s in range(N))
        assert all(q[x]==q[reflection[x]] for x in range(N))
        charges.append([z[u][-1] for u in reps])
    assert charges[0]==charges[1]==[1,-2,2,-1,-5,7,3,5,-11,1]
    word=(0,4,1,2,3)
    assert all(arc(codes[4],word[i],word[i+1],pairs) for i in range(n-1))
    normalized=sum(arc(codes[4],word[a],word[b],pairs)<<k for k,(a,b) in enumerate(free))
    assert normalized==9 and bad[4]!=bad[9]
    report={'status':'PASS','method':'all 1024 labelled tournaments, orbit closure under adjacent vertex swaps plus converse','ordinary_classes':len(set(ordinary)),'merged_classes':len(set(merged)),'mass_preserving_node_permutations':num,'coloured_graph_isomorphisms':isos,'endpoint_profiles_against_HP_Aut':'PASS','boolean_Walsh_convolutions':'PASS','reflection_invariance':'PASS','top_charges':charges[0],'broken_transport':{'source':4,'path':word,'normalized':normalized}}
    (HERE/'independent_hostile_audit.json').write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2))

if __name__=='__main__':main()
