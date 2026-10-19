"""Independent source-permutation and transitive-host controls for cherry extensions."""
from independent_tree_capacity import automorphisms, canon
from itertools import permutations
import json

checks=0
def check(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def arcs(edges,mask):
    return {(v,u) if mask>>i&1 else (u,v) for i,(u,v) in enumerate(edges)}

def swap_preserves(aa,u,v):
    def p(x):return v if x==u else u if x==v else x
    return {(p(x),p(y)) for x,y in aa}==aa

def linear_extensions(n,aa):
    pred=[0]*n
    for u,v in aa:pred[v]|=1<<u
    dp=[0]*(1<<n);dp[0]=1
    for mask in range(1<<n):
        if not dp[mask]:continue
        for v in range(n):
            if not (mask>>v&1) and not (pred[v]&~mask):dp[mask|(1<<v)]+=dp[mask]
    return dp[-1]

def orient_auto(aa,autos):
    return sum({(p[u],p[v]) for u,v in aa}==aa for p in autos)

def component_sizes(n,edges):
    adj=[[] for _ in range(n)]
    for u,v in edges:adj[u].append(v);adj[v].append(u)
    seen=set();sizes=[]
    for u in range(n):
        if u in seen:continue
        q=[u];seen.add(u)
        for v in q:
            for w in adj[v]:
                if w not in seen:seen.add(w);q.append(w)
        sizes.append(len(q))
    return sorted(sizes)

def main():
    two=((0,1),(0,2),(0,3),(0,4),(1,5),(1,6),(2,7),(7,8))
    auto=automorphisms(9,two);check(len(auto)==4,'two-cherry automorphism group')
    asym=direct=reduced=0
    for mask in range(256):
        aa=arcs(two,mask);g=orient_auto(aa,auto);asym+=g==1
        if swap_preserves(aa,3,4) or swap_preserves(aa,5,6):
            direct+=1;check(g%2==0,'direct even symmetry')
        else:
            reduced+=1
            e=next(e for e in aa if set(e)=={0,3});rev=(aa-{e})|{e[::-1]}
            check(swap_preserves(rev,3,4),'first reverse cancels')
            deleted=aa-{e};f=next(f for f in deleted if set(f)=={1,5});rev2=(deleted-{f})|{f[::-1]}
            check(swap_preserves(rev2,5,6),'second reverse cancels')
            check(swap_preserves(deleted-{f},3,5),'two isolated vertices cancel')
        check(linear_extensions(9,aa)%2==0,'two-cherry transitive parity')
    check((direct,reduced,asym)==(192,64,64),'two-cherry state counts')
    concrete=arcs(two,40)
    check(linear_extensions(9,concrete)==260,'named transitive exact count')
    check(orient_auto(concrete,auto)==1,'named directed self-host exact count')
    d8=((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7))
    profiles={}
    for i,e in enumerate(d8):
        prof=tuple(component_sizes(9,d8[:i]+d8[i+1:]));profiles[str(e)]=prof
        check(prof in ((1,1,7),(1,2,6),(1,4,4)),'D8 derivative profile')
    d8_counts=[linear_extensions(9,arcs(d8,mask)) for mask in range(128)]
    check(all(c%2==1 for c in d8_counts),'all D8 plus isolate transitive parities')
    ext=[]
    for attach,existing_leaf in ((0,3),(2,6)):
        ee=d8+((attach,8),);autos=automorphisms(9,ee)
        check(len(autos)==2,'extended D8 automorphism group')
        odd=asymmetric=0
        for mask in range(256):
            aa=arcs(ee,mask);same=swap_preserves(aa,existing_leaf,8)
            c=linear_extensions(9,aa);check(c%2==int(not same),'D8 leaf extension parity')
            odd+=c%2;asymmetric+=orient_auto(aa,autos)==1
        check(odd==128 and asymmetric==128,'D8 orientation classes')
        ext.append({'attachment':attach,'automorphisms':len(autos),'odd_labelled_orientations':odd,'odd_orientation_classes':odd//len(autos),'asymmetric_labelled_orientations':asymmetric})
    check(canon(9,d8+((0,8),))!=canon(9,d8+((2,8),)),'D8 extensions distinct')
    print(json.dumps({'status':'PASS','checks':checks,'two_cherry':{'direct':direct,'two_deletion':reduced,'asymmetric_labelled':asym,'asymmetric_classes':asym//len(auto),'named_transitive_count':260,'named_selfhost_count':1},'D8_derivative_profiles':profiles,'D8_plus_isolate_odd_orientations':128,'D8_extensions':ext,'scope':'Independent source-side deletion certificates, direct permutation automorphism counts, and exact transitive-host linear extensions. Universal tournament claims rest on the accompanying mathematical review, not on sampling.'},indent=2))

if __name__=='__main__':main()
