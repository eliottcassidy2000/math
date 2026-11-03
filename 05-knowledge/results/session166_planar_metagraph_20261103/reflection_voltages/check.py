"""Independent finite controls for reflection-voltage classification.
No third-party packages. Every check survives python -O.
"""
from collections import Counter, defaultdict
from itertools import permutations, product
from fractions import Fraction
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament_data(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    chords = [(i,j) for i in range(n) for j in range(i+2,n)]
    idx = {e:i for i,e in enumerate(pairs)}
    ci = {e:i for i,e in enumerate(chords)}
    perm_maps=[]
    for p in permutations(range(n)):
        source=[idx[tuple(sorted((p[i],p[j])))] for i,j in pairs]
        reversal=sum((p[i]>p[j])<<k for k,(i,j) in enumerate(pairs))
        perm_maps.append((source,reversal))
    def relabel(mask, source, rev):
        answer=rev
        for k,j in enumerate(source): answer ^= ((mask>>j)&1)<<k
        return answer
    cache={}
    def canon(mask):
        if mask not in cache:
            orbit={relabel(mask,a,b) for a,b in perm_maps}
            c=min(orbit)
            for x in orbit: cache[x]=c
        return cache[mask]
    m=len(chords); allmask=(1<<m)-1
    C=[]
    for t in range(1<<m):
        full=sum(((t>>k)&1)<<idx[e] for k,e in enumerate(chords))
        C.append(canon(full))
    classes=sorted(set(C))
    op={c:canon(c ^ ((1<<len(pairs))-1)) for c in classes}
    sidx=[ci[n-1-j,n-1-i] for i,j in chords]
    S=[sum(((t>>k)&1)<<j for k,j in enumerate(sidx)) for t in range(1<<m)]
    require(all(C[S[t]]==op[C[t]] for t in range(1<<m)), 'converse action')
    require(all(S[t]!=t^allmask for t in range(1<<m)), 'SF has no fixed flag')
    edges=[(t,t^allmask) for t in range(1<<(m-1))]
    fixed={c for c in classes if c==op[c]}
    nodes=sorted({min(c,op[c]) for c in classes})
    free=[u for u in nodes if u not in fixed]
    blue=[]; black=[]; seen=set(); types=Counter(); loop_witnesses={}
    for a,b in edges:
        if S[a]==a:
            blue.append((C[a],C[b]))
            require(C[a] in fixed and C[b] in fixed,'blue endpoints fixed')
            continue
        key=tuple(sorted((a,b,S[a],S[b])))
        if key in seen:continue
        seen.add(key)
        u,v=min(C[a],op[C[a]]),min(C[b],op[C[b]])
        sign=None if u in fixed or v in fixed else (1 if (C[a]==u)==(C[b]==v) else -1)
        black.append((u,v,sign))
        if u in fixed and v in fixed:typ='SC-SC loop' if u==v else 'SC-SC cross'
        elif u in fixed or v in fixed:typ='SC-NS'
        elif u==v:typ='NS loop untwisted' if sign==1 else 'NS loop twisted'
        else:typ='NS-NS parallel' if sign==1 else 'NS-NS crossed'
        types[typ]+=1
        if u==v and sign is not None and str(sign) not in loop_witnesses:
            loop_witnesses[str(sign)]={'node':u,'converse_class':op[u],'tilings':[a,b,S[a],S[b]],'class_endpoints':[C[a],C[b],C[S[a]],C[S[b]]]}
    # Incidence matrices use two on every loop diagonal.
    def matrix(edge_list):
        A=Counter()
        for u,v in edge_list:A[u,v]+=1;A[v,u]+=1
        return +A
    actual_B=matrix([(C[a],C[b]) for a,b in edges if S[a]==a])
    actual_K=matrix([(C[a],C[b]) for a,b in edges if S[a]!=a])
    lifted=[]
    for u,v,sign in black:
        for eps in (0,1):
            cu=u if u in fixed or eps==0 else op[u]
            ev=eps if sign is None or sign==1 else 1-eps
            cv=v if v in fixed or ev==0 else op[v]
            lifted.append((cu,cv))
    require(matrix(lifted)==actual_K,'black reconstruction including loops')
    require(matrix(blue)==actual_B,'blue reconstruction')
    # Raw black-orbit count, inherited color/fiber parity, quotient degree.
    d=(n-1)**2//4
    require(len(blue)==2**(d-1),'blue count')
    require(len(black)==((1<<m)-(1<<d))//4,'black orbit count')
    Q=matrix([(u,v) for u,v,s in black])
    require(all(sum(Q[u,v] for v in nodes)%2==1 for u in free),'free base degrees odd')
    require(all(sum(actual_B[c,v] for v in classes)%2==1 for c in fixed),'blue fixed degrees odd')
    # Exact invariant/anti-invariant decomposition without square roots.
    # E columns are fixed indicators / sums, R columns are sheet differences.
    A=actual_B+actual_K
    M=Counter()
    for (c,e),value in A.items():M[min(c,op[c]),min(e,op[e])]+=value
    signed=Counter()
    for u,v,s in black:
        if s is not None:signed[u,v]+=s;signed[v,u]+=s
    for c in classes:
        u=min(c,op[c]);size=1 if u in fixed else 2
        for v in nodes:
            direct=sum(A[c,e] for e in classes if min(e,op[e])==v)
            require(Fraction(direct)==Fraction(M[u,v],size),'invariant block A E = E D^-1 M')
        for v in free:
            direct=A[c,v]-A[c,op[v]]
            predicted=0 if u in fixed else (1 if c==u else -1)*signed[u,v]
            require(direct==predicted,'anti-invariant block A R = R A_signed')
    # Colored-word endpoint identities and character traces, exact through length 4.
    word_checks=0; trace_checks=0
    if n<=6:
        def multiply(P,R,labels):
            rows=defaultdict(list)
            for (u,v),x in R.items():
                if x:rows[u].append((v,x))
            out=Counter()
            for (u,k),x in P.items():
                if x:
                    for v,y in rows[k]:out[u,v]+=x*y
            return out
        products={():Counter({(c,c):1 for c in classes})}
        signed_power=Counter({(u,u):1 for u in free})
        for length in range(1,5):
            signed_power=multiply(signed_power,signed,free)
            for word in product('BK',repeat=length):
                W=multiply(products[word[:-1]],actual_B if word[-1]=='B' else actual_K,classes)
                products[word]=W
                if 'B' in word:
                    require(all(W[c,d]==W[c,op[d]]==W[op[c],d] for c in classes for d in classes),'blue-word endpoint bijection')
                    word_checks+=1
                else:
                    difference=sum(W[c,c]-W[c,op[c]] for c in classes)
                    require(difference==2*sum(signed_power[u,u] for u in free),'black signed trace character')
                    trace_checks+=1
    # Independently calculate signed cycle rank and canonical gauge normalization.
    ff=[(u,v,s) for u,v,s in black if s is not None]
    adj=defaultdict(list)
    for i,(u,v,s) in enumerate(ff):adj[u].append((v,i,s));adj[v].append((u,i,s))
    gauge={};tree=set();components=0
    for root in free:
        if root in gauge:continue
        components+=1;gauge[root]=1;stack=[root]
        while stack:
            u=stack.pop()
            for v,i,s in adj[u]:
                if v not in gauge:gauge[v]=gauge[u]*s;tree.add(i);stack.append(v)
    normalized=[s*gauge[u]*gauge[v] for u,v,s in ff]
    require(all(normalized[i]==1 for i in tree),'tree signs normalize positive')
    rank=len(ff)-len(free)+components
    require(rank==len(ff)-len(tree),'cycle rank')
    # Gauge-invariant holonomy against every possible gauge for n<=5.
    gauge_checked=0
    if len(free)<=10:
        tree_pairs=[ff[i][:2] for i in tree]
        for gs in product((1,-1),repeat=len(free)):
            g=dict(zip(free,gs));changed=[s*g[u]*g[v] for u,v,s in ff]
            h={}
            for root in free:
                if root in h:continue
                h[root]=1;stack=[root]
                while stack:
                    u=stack.pop()
                    for v,i,_ in adj[u]:
                        if i in tree and v not in h:h[v]=h[u]*changed[i];stack.append(v)
            require([s*h[u]*h[v] for (u,v,_),s in zip(ff,changed)]==normalized,'gauge normal form invariant')
            gauge_checked+=1
    report={'n':n,'classes':len(classes),'SC':len(fixed),'NS_nodes':len(free),'blue_lines':len(blue),'black_orbits':len(black),'black_orbit_types':dict(sorted(types.items())),'NS_subgraph_edges':len(ff),'NS_subgraph_components':components,'NS_cycle_rank':rank,'normalized_negative_chords':sum(s<0 for s in normalized),'gauges_checked':gauge_checked,'NS_loop_witnesses':loop_witnesses,'blue_word_matrix_checks':word_checks,'black_trace_checks':trace_checks}
    return report


def abstract_gauge_controls():
    # Loop, multiedge, disconnected and isolated cases; all signatures,
    # all gauges; quotient has exactly 2^(E-V+c) switching classes.
    examples=[(1,[]),(1,[(0,0)]),(2,[(0,1),(0,1)]),(3,[(0,1),(1,2),(2,0)]),(4,[(0,1),(0,1),(2,2)])]
    output=[]
    for n,edges in examples:
        representatives=set(); stabilizer_sizes=set()
        for signs in product((1,-1),repeat=len(edges)):
            orbit=set()
            for gs in product((1,-1),repeat=n):orbit.add(tuple(s*gs[u]*gs[v] for (u,v),s in zip(edges,signs)))
            representatives.add(min(orbit))
            stabilizer_sizes.add(sum(all(s*gs[u]*gs[v]==s for (u,v),s in zip(edges,signs)) for gs in product((1,-1),repeat=n)))
        parent=list(range(n))
        def find(u):
            while parent[u]!=u:u=parent[u]
            return u
        for u,v in edges:parent[find(u)]=find(v)
        comp=len({find(u) for u in range(n)});rank=len(edges)-n+comp
        require(len(representatives)==1<<rank,'all signature orbit count')
        require(stabilizer_sizes=={1<<comp},'all signature gauge stabilizer counts')
        output.append({'vertices':n,'edges':edges,'cycle_rank':rank,'gauge_classes':len(representatives)})
    return output


if __name__=='__main__':
    report={'tournaments':[tournament_data(n) for n in (3,4,5,6,7)],'abstract_gauge_controls':abstract_gauge_controls(),'status':'PASS'}
    print(json.dumps(report,indent=2,sort_keys=True))
