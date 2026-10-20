#!/usr/bin/env python3
"""Standalone surface audit. F9 elements are tuples, not producer encodings."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
import json

CHECKS=Counter()
def need(ok,tag):
    CHECKS[tag]+=1
    if not ok:raise RuntimeError(tag)

def components(vertices,adj):
    unseen=set(vertices);answer=[]
    while unseen:
        stack=[unseen.pop()];part=set(stack)
        while stack:
            for w in adj[stack.pop()]:
                if w in unseen:
                    unseen.remove(w);part.add(w);stack.append(w)
        answer.append(frozenset(part))
    return answer

def rank_mod3(columns,rows):
    pivots={}
    for c in columns:
        v={i:x%3 for i,x in c.items() if x%3}
        while v:
            k=min(v)
            if k not in pivots:
                factor=pow(v[k],-1,3)
                pivots[k]={i:x*factor%3 for i,x in v.items()};break
            factor=v[k]
            for i,x in pivots[k].items():
                y=(v.get(i,0)-factor*x)%3
                if y:v[i]=y
                elif i in v:del v[i]
    return len(pivots)

def audit(p,r=1):
    q=p**r;els=list(product(range(p),repeat=r));zero=(0,)*r;one=(1,)+(0,)*(r-1)
    add=lambda x,y:tuple((a+b)%p for a,b in zip(x,y))
    neg=lambda x:tuple(-a%p for a in x)
    sub=lambda x,y:add(x,neg(y))
    def mul(x,y):
        if r==1:return ((x[0]*y[0])%p,)
        a,b=x;c,d=y
        return ((a*c-b*d)%p,(a*d+b*c)%p)
    va=lambda v,w:(add(v[0],w[0]),add(v[1],w[1]))
    vn=lambda v:(neg(v[0]),neg(v[1]))
    scale=lambda a,v:(mul(a,v[0]),mul(a,v[1]))
    canon=lambda v:min(v,vn(v))
    det=lambda v,w:sub(mul(v[0],w[1]),mul(v[1],w[0]))
    V=sorted({canon(v) for v in product(els,repeat=2) if v!=(zero,zero)})
    vs={v:i for i,v in enumerate(V)}
    E=[(i,j) for i,j in combinations(range(len(V)),2) if det(V[i],V[j]) in (one,neg(one))]
    eset=set(E);ei={e:i for i,e in enumerate(E)}
    adj={i:set() for i in range(len(V))}
    for i,j in E:adj[i].add(j);adj[j].add(i)
    triangles={tuple(sorted((i,j,k))) for i,j in E for k in adj[i]&adj[j]}
    F=sorted(triangles);fe=defaultdict(set)
    for f,t in enumerate(F):
        for e in combinations(t,2):fe[e].add(f)
    need((len(V),len(E),len(F))==((q*q-1)//2,q*(q*q-1)//4,q*(q*q-1)//6),'f_vector')
    need(len(components(range(len(V)),adj))==1,'graph_connected')
    need(all(len(adj[i])==q for i in adj),'graph_degree_q')
    need(all(len(f)==2 for f in fe.values()),'two_faces_per_edge')

    # Independent characteristic-p link calculation at every vertex.
    link_lengths=[]
    for i in adj:
        la={j:adj[j]&adj[i] for j in adj[i]}
        need(all(len(x)==2 for x in la.values()),'link_degree_two')
        sizes=sorted(map(len,components(la,la)));link_lengths.append(sizes)
        need(sizes==[p]*(q//p),'link_components_are_p_cycles')

    # Derive orientation signs on each sorted face, then inspect boundaries.
    face_signs=[];boundary=Counter()
    for t in F:
        cyclic_orders=set()
        reps=0
        for order in permutations(t):
            for signs in product((one,neg(one)),repeat=3):
                a,b,c=[scale(s,V[i]) for s,i in zip(signs,order)]
                if va(va(a,b),c)==(zero,zero) and det(a,b)==one:
                    cyclic_orders.add(order);reps+=1
        need(reps==6 and len(cyclic_orders)==3,'unique_triangle_orientation')
        order=min(cyclic_orders)
        need(cyclic_orders=={order,order[1:]+order[:1],order[2:]+order[:2]},'orientation_is_cyclic')
        sign=(-1)**sum(a>b for i,a in enumerate(order) for b in order[i+1:])
        face_signs.append(sign)
        for a,b in zip(order,order[1:]+order[:1]):boundary[tuple(sorted((a,b)))]+=1 if a<b else -1
    need(all(boundary[e]==0 for e in E),'coherent_orientation_cancels_every_edge')

    # Integral boundary signs evaluated mod 3, independent of the genus formula.
    d1=[{i:-1,j:1} for i,j in E]
    d2=[{ei[(b,c)]:1,ei[(a,c)]:-1,ei[(a,b)]:1} for a,b,c in F]
    rank1=rank_mod3(d1,len(V));rank2=rank_mod3(d2,len(E))
    betti=[len(V)-rank1,len(E)-rank1-rank2,len(F)-rank2]
    genus=1+(p*p-1)*(p-6)//24
    nc=q*(q*q-1)//(p*(p*p-1))
    need(betti[0]==1 and betti[2]==nc,'homology_component_fundamental_classes')
    need(betti[1]==2*nc*genus+nc*(p*p-1)//2-len(V)-nc+1,'homology_after_vertex_identifications')

    # Dual-face connectivity first, without knowing subfield planes.
    dual={i:set() for i in range(len(F))}
    for fs in fe.values():
        a,b=sorted(fs);dual[a].add(b);dual[b].add(a)
    dualparts=components(dual,dual)
    need(len(dualparts)==nc,'normalization_component_number')
    part_vertices=[]
    for fs in dualparts:
        verts=frozenset(i for f in fs for i in F[f])
        edges={e for f in fs for e in combinations(F[f],2)}
        need((len(verts),len(edges),len(fs))==((p*p-1)//2,p*(p*p-1)//4,p*(p*p-1)//6),'normalization_surface_f_vector')
        part_vertices.append(verts)
    need(set(Counter(v for vs0 in part_vertices for v in vs0).values())=={q//p},'branches_through_each_vertex')

    # Matrix cosets are evaluated for F9, including the side of the coset.
    coset_summary=None
    if q==9:
        def mm(m,n):
            a,b,c,d=m;e,f,g,h=n
            return (add(mul(a,e),mul(b,g)),add(mul(a,f),mul(b,h)),add(mul(c,e),mul(d,g)),add(mul(c,f),mul(d,h)))
        G={m for m in product(els,repeat=4) if sub(mul(m[0],m[3]),mul(m[1],m[2]))==one}
        H={m for m in G if all(x[1]==0 for x in m)}
        small=[(a,0) for a in range(p)]
        def plane(m):
            a,b,c,d=m;u=(a,c);v=(b,d)
            return frozenset(vs[canon(va(scale(x,u),scale(y,v)))] for x,y in product(small,repeat=2) if (x,y)!=(zero,zero))
        groups=defaultdict(set)
        for m in G:groups[plane(m)].add(m)
        need(len(G)==720 and len(H)==24 and len(groups)==30,'F9_SL2_counts')
        need(set(groups)==set(part_vertices),'F9_cosets_equal_dual_surface_components')
        for W,group in groups.items():
            m=min(group)
            need({mm(m,h) for h in H}==group,'F9_each_plane_is_one_left_coset_gH')
        for W,Z in combinations(groups,2):
            need(not any(a in W&Z and b in W&Z for a,b in E),'F9_distinct_planes_share_no_edges')
        need(betti==[1,51,30],'F9_singular_complex_homology')
        # Exactly 30 K4s: filling clique dimension three would change the object.
        k4s={tuple(sorted((a,b,c,d))) for a,b,c in F for d in adj[a]&adj[b]&adj[c]}
        need(len(k4s)==30 and set(map(frozenset,k4s))==set(groups),'F9_thirty_tetrahedral_boundaries_not_filled')
        coset_summary={'SL2_F9':len(G),'SL2_F3':len(H),'left_cosets':len(groups),'K4_cliques_not_filled':len(k4s)}
    if q==3:
        need((len(V),len(E),len(F))==(4,6,4) and betti==[1,0,1],'p3_boundary_not_clique_complex')
    if r==1:need(betti==[1,2*genus,1],'prime_surface_homology')
    return {'field_order':q,'characteristic':p,'f_vector':[len(V),len(E),len(F)],'genus_per_normalized_component':genus,'normalization_components':nc,'vertex_branches':q//p,'betti_mod3':betti,'link_lengths':link_lengths[0],'coset_summary':coset_summary}

runs=[audit(p) for p in [3,5,7,11]]+[audit(3,2)]
print(json.dumps({'status':'PASS','checks':dict(CHECKS),'total_checks':sum(CHECKS.values()),'runs':runs},indent=2,sort_keys=True))
