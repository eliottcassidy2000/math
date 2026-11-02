"""Exact Collatz growth, lattice capacity, phase and palette controls."""
from collections import Counter,deque
from fractions import Fraction as F
from itertools import permutations,product,combinations
from pathlib import Path
from math import gcd
import hashlib,importlib.util,json
from chart_checks import main as chart_checks
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def C(n,p=3):return n//2 if n%2==0 else p*n+1
def U(n,p=3):
    x=p*n+1;a=0
    while x%2==0:x//=2;a+=1
    return x,a
def preds(n,p=3):
    out=[2*n]
    if (n-1)%p==0 and (n-1)//p>0 and ((n-1)//p)%2:out.append((n-1)//p)
    return out
def sf5_children(n):
    A=[2,4,6] if n%3==1 else [3,5,7]
    out=[((2**a*n-1)//3,a) for a in A if (2**a*n-1)%9!=0]
    need(len(out)==2 and all(u>n and u%2 and u%3 and U(u)==(n,a) for u,a in out),'SF5 children')
    return out
def honey_neighbors(v):
    a,b,s=v
    return [(a,b,1),(a-1,b,1),(a,b-1,1)] if s==0 else [(a,b,0),(a+1,b,0),(a,b+1,0)]
def triangular_neighbors(v):
    a,b=v
    return [(a+x,b+y) for x,y in [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]]
def ball_counts(root,neighbors,maxr):
    seen={root};front={root};out=[1]
    for r in range(1,maxr+1):
        nxt={w for v in front for w in neighbors(v)}-seen
        seen|=nxt;front=nxt;out.append(len(seen))
    return out
def apply_affine(sig,v,x):
    return tuple(x[sig[i]]^v[i] for i in range(3))

def joint_port_checks():
    # Exact native expansion, compared to independent path unions and pair meetings.
    nodes={():5};parent={};edge_clock={};levels=[[()]]
    for h in range(1,9):
        layer=[]
        for address in levels[-1]:
            for bit,(u,a) in enumerate(sf5_children(nodes[address])):
                child=address+(bit,);nodes[child]=u;parent[child]=address
                edge_clock[child]=a+1;layer.append(child)
        levels.append(layer)
    paths={};positions={}
    for address,n in nodes.items():
        path=[n]
        while path[-1]!=5:
            need(len(path)<=8*len(address),'native path must terminate by the stored clock')
            path.append(C(path[-1]))
        paths[address]=path;positions[address]={x:i for i,x in enumerate(path)}
    pairs=0;shifted=0;witness=None
    for u,v in combinations(nodes,2):
        common=0
        while common<min(len(u),len(v)) and u[common]==v[common]:common+=1
        q=u[:common];parent_label=nodes[q]
        tu=sum(edge_clock[u[:i]] for i in range(common+1,len(u)+1))
        tv=sum(edge_clock[v[:i]] for i in range(common+1,len(v)+1))
        if common==min(len(u),len(v)):
            shift=0
        else:
            shift=min(edge_clock[u[:common+1]],edge_clock[v[:common+1]])-1
        expected=2**shift*parent_label
        first=next(x for x in paths[u] if x in positions[v])
        iu=positions[u][first];iv=positions[v][first]
        need(first==expected and (iu,iv)==(tu-shift,tv-shift),'native joint-port meeting formula')
        need(iu+iv==tu+tv-2*shift,'native tree metric')
        pairs+=1;shifted+=bool(shift)
        if nodes[u]==13 and nodes[v]==53:
            witness={'sources':[13,53],'odd_lca':5,'native_meeting':first,'meeting_clocks':[iu,iv],'distance':iu+iv}
    def polynomial(n,h):
        if h==0:return Counter({0:1})
        children=sf5_children(n);amax=max(a for _,a in children)
        out=Counter({i:1 for i in range(amax+1)})
        for u,a in children:
            for d,c in polynomial(u,h-1).items():out[d+a+1]+=c
        return out
    rows=[]
    for h in range(9):
        selected=[a for a in nodes if len(a)<=h]
        union={x for a in selected for x in paths[a]}
        literal={}
        for a in selected:
            for i,x in enumerate(paths[a]):
                distance=len(paths[a])-1-i
                need(x not in literal or literal[x]==distance,'unique native first-hit clock')
                literal[x]=distance
        predicted=polynomial(5,h)
        need(Counter(literal.values())==predicted,'shared-tail weighted recurrence')
        size_formula=len(selected)+sum(max(a for _,a in sf5_children(nodes[q])) for q in selected if len(q)<h)
        naive=1+sum(edge_clock[q] for q in selected if q)
        duplicate=sum(min(a for _,a in sf5_children(nodes[q])) for q in selected if len(q)<h)
        need(len(union)==size_formula==naive-duplicate,'exact gluing correction')
        rows.append({'depth':h,'odd_nodes':len(selected),'native_nodes':len(union),'naive_subdivision_nodes':naive,'duplicate_even_states':duplicate,'max_native_clock':max(predicted)})
    need(witness is not None and witness['native_meeting']==40,'13/53 meeting witness')
    # Tail overlap is tested from literal labelled edge sets, then by determinant.
    exponents=[3,5,9,11,15,17]
    tails=[{(2**i*5,2**(i-1)*5) for i in range(1,a+1)} for a in exponents]
    K=[[len(s&t) for t in tails] for s in tails]
    need(K==[[min(a,b) for b in exponents] for a in exponents],'tail Gram kernel')
    M=[[F(x) for x in row] for row in K];det=F(1)
    for i in range(len(M)):
        pivot=M[i][i];need(pivot!=0,'tail Gram pivot');det*=pivot
        for j in range(i+1,len(M)):
            ratio=M[j][i]/pivot
            for k in range(i,len(M)):M[j][k]-=ratio*M[i][k]
    predicted=exponents[0]
    for a,b in zip(exponents,exponents[1:]):predicted*=b-a
    need(det==predicted,'exact Gram determinant')
    leaf_ranks={}
    for h,level in enumerate(levels):
        leaf_paths=[paths[a] for a in level]
        edges=[set(zip(path,path[1:])) for path in leaf_paths]
        if h==0:
            need([[len(s&t) for t in edges] for s in edges]==[[0]],'stopped root Gram rank is zero')
            leaf_ranks[h]=0
        else:
            for i,path in enumerate(leaf_paths):
                private=(path[0],path[1])
                need(sum(private in s for s in edges)==1,'private leaf edge gives full Gram rank')
            leaf_ranks[h]=len(leaf_paths)
    alias=[]
    for source in (13,853):
        pth=[source]
        while pth[-1]!=5:pth.append(C(pth[-1]))
        fixed=[53,160,80,40,20,10,5]
        meet=next(x for x in pth if x in fixed)
        alias.append({'source':source,'a':U(source)[1],'meeting_with53':meet,'clocks':[pth.index(meet),fixed.index(meet)]})
    need([x['meeting_with53'] for x in alias]==[40,160] and alias[0]['a']%6==alias[1]['a']%6,'periodic clock phase loses fixed-port response')
    return {'pairs_checked':pairs,'pairs_with_earlier_even_meeting':shifted,'witness':witness,'expansion_rows':rows,
        'tail_gram_exponents':exponents,'tail_gram_determinant':int(det),'leaf_gram_rank':len(leaf_paths),'leaf_gram_ranks_by_depth':leaf_ranks,'periodic_phase_hostile':alias}
def main():
    joint=joint_port_checks()
    # Independent inverse generation and forward replay of every literal label.
    layers=[[5]];seen={5};depth={5:0}
    for r in range(1,27):
        new=sorted({u for n in layers[-1] for u in preds(n)})
        need(not (seen&set(new)),'root5 tree has no repeated vertex')
        for u in new:
            need(C(u) in set(layers[-1]),'native edge replay')
            x=u
            for t in range(r):
                need(x!=5,'first-hit clock');x=C(x)
            need(x==5,'whole native path replay');depth[u]=r
        seen.update(new);layers.append(new)
    counts=[];total=0
    for layer in layers:total+=len(layer);counts.append(total)
    need(counts[21]==704 and 1+3*21*22//2==694,'honeycomb finite obstruction')
    need(all(counts[r]<=1+3*r*(r+1)//2 for r in range(21)),'first radial count obstruction at21 for this root')
    need(counts[26]==2263 and 1+3*26*27==2107,'triangular finite obstruction')
    need(all(counts[r]<=1+3*r*(r+1) for r in range(26)),'first triangular radial obstruction')
    H=ball_counts((0,0,0),honey_neighbors,105);T=ball_counts((0,0),triangular_neighbors,105)
    need(all(H[r]==1+3*r*(r+1)//2 for r in range(106)),'literal honeycomb balls')
    need(all(T[r]==1+3*r*(r+1) for r in range(106)),'literal triangular balls')
    # Recovered source-labelled tree, rooted at5, up to depth13.
    selected=[(5,(),0)];all_nodes={5};sf5_rows=[]
    for h in range(1,14):
        nxt=[]
        for n,address,clock in selected:
            for bit,(u,a) in enumerate(sf5_children(n)):
                need(u not in all_nodes,'unique inverse address');all_nodes.add(u)
                x=u
                for _ in range(a+1):x=C(x)
                need(x==n,'macro/native clock');nxt.append((u,address+(bit,),clock+a+1))
        selected=nxt;need(len(selected)==2**h,'binary level count')
        sf5_rows.append({'depth':h,'level_count':len(selected),'cumulative':len(all_nodes),'max_raw_clock':max(t for _,_,t in selected)})
    need(len(all_nodes)==16383 and H[104]==16381,'all-size proof finite witness')
    # Smaller geometric rank is explicit only on this selected proper family.
    need(U(7)[0]>7,'known convergent integer outside monotone tree');x=7
    for _ in range(100):
        if x==1:break
        x=C(x)
    need(x==1,'excluded integer is nevertheless convergent')
    # A rank-decreasing labelled binary drawing, not a fixed-scale lattice embedding.
    coordinates={};drawing_checks=0
    for h in range(8):
        for j in range(2**h):
            coordinates[h,j]=(F(2*j+1,2**(h+1)),1-F(1,2**h))
            if h:
                x,y=coordinates[h,j];u,v=coordinates[h-1,j//2]
                need((x-u)**2+(y-v)**2==F(5,2**(2*h+2)),'dyadic drawing edge scale');drawing_checks+=1
    need(len(set(coordinates.values()))==len(coordinates),'drawing injectivity')
    # Fixed mod3^k charts fail under the same legal inverse branch a=2.
    phase=[]
    for k in range(2,9):
        q=19;other=q+2*3**k;a=2
        need(q%2==1 and q%3==1,'phase seed')
        u=(2**a*q-1)//3;v=(2**a*other-1)//3
        need(q%3**k==other%3**k and u%3**k!=v%3**k,'lost3-adic digit')
        need(any(w==u for w,b in sf5_children(q)) and any(w==v for w,b in sf5_children(other)),'same legal branch')
        phase.append({'k':k,'roots':[q,other],'outputs':[u,v],'modulus':3**k})
    # Generalization at p5: four increasing unit predecessors in one order25 block.
    p5=[]
    for n in range(1,202,2):
        if n%5==0:continue
        cand=[(a,(2**a*n-1)//5) for a in range(3,23) if (2**a*n-1)%5==0 and (2**a*n-1)%25!=0]
        need(len(cand)==4 and all(u>n and U(u,5)==(n,a) and u%5 for a,u in cand),'5n+1 branching hostile')
        p5.append([n,[a for a,u in cand]])
    need([U(x,5)[0] for x in [13,33,83]]==[33,83,13],'nonprincipal5n+1 cycle')
    # Exhaust all48 palette affine maps. Recovered holonomy theorem instantiated.
    histogram=Counter();parity_pass=0;false_pass=0
    frames=list(product(range(2),repeat=3))
    for sig in permutations(range(3)):
        cycles=3 if sig==(0,1,2) else (2 if any(sig[i]==i for i in range(3)) else 1)
        image={tuple(x[i]^x[sig[i]] for i in range(3)) for x in frames}
        for v in frames:
            fixed=[x for x in frames if apply_affine(sig,v,x)==x]
            predicted=2**cycles if v in image else 0
            need(len(fixed)==predicted,'affine fixed-frame classification');histogram[len(fixed)]+=1
            if sum(v)%2==0:parity_pass+=1;false_pass+=not fixed
    need(dict(histogram)=={8:1,0:33,4:6,2:8} and parity_pass==24 and false_pass==9,'complete palette census')
    # All labelled total maps on up to5 vertices: pseudoforest and component colouring count.
    maps_checked=0
    for n in range(1,6):
        for f in product(range(n),repeat=n):
            E={tuple(sorted((i,f[i]))) for i in range(n) if f[i]!=i}
            todo=set(range(n));poly=1
            while todo:
                root=min(todo);part={root};stack=[root]
                while stack:
                    u=stack.pop()
                    for e in E:
                        if u in e:
                            w=e[0] if e[1]==u else e[1]
                            if w not in part:part.add(w);stack.append(w)
                todo-=part;edges=[e for e in E if e[0] in part]
                need(len(edges)<=len(part),'functional pseudoforest')
                if len(edges)==len(part)-1:poly*=3*2**(len(part)-1)
                else:
                    deg=Counter(x for e in edges for x in e);remaining=set(part);queue=[x for x in part if deg[x]==1]
                    while queue:
                        u=queue.pop();remaining.remove(u)
                        for e in edges:
                            if u in e:
                                w=e[0] if e[1]==u else e[1]
                                if w in remaining:
                                    deg[w]-=1
                                    if deg[w]==1:queue.append(w)
                    ell=len(remaining);need(ell>=3,'simple cycle length');poly*=2**(len(part)-ell)*(2**ell+(-1)**ell*2)
            literal=sum(all(color[u]!=color[v] for u,v in E) for color in product(range(3),repeat=n))
            need(literal==poly>0,'chromatic product on simple functional graph');maps_checked+=1
    baseline=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('hostile_controls',baseline);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    hostile=mod.hostile_checks()
    raw='\n'.join(f'{n},{d}' for n,d in sorted(depth.items()))+'\n'
    (HERE/'inverse_root5_depth26.tsv').write_text('integer\tfirst_hit_steps\n'+raw.replace(',','\t'),encoding='utf-8')
    return {'status':'PASS','inverse_ball_counts':counts,'inverse_level_counts':[len(x) for x in layers],
        'root5_certificate_sha256':hashlib.sha256(raw.encode()).hexdigest(),'sf5_rows':sf5_rows,
        'honeycomb_obstruction':{'radius':21,'sources':704,'capacity':694},'triangular_obstruction':{'radius':26,'sources':2263,'capacity':2107},
        'lattice_ball_checks':212,'drawing_checks':drawing_checks,'phase_hostiles':phase,'p5_parents_checked':len(p5),
        'p5_cycle':[13,33,83],'palette_histogram':dict(histogram),'parity_accepts':parity_pass,'parity_false_accepts':false_pass,
        'functional_maps_checked':maps_checked,'hostile81_27':hostile,'joint_ports':joint,'finite_horizon_charts':chart_checks()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
