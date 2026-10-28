"""Exact finite bridge tests. Standard library; every guard survives -O."""
from itertools import combinations, permutations, product
from collections import Counter
from math import gcd, isqrt, prod
import json

BANK=(4,2,5,3,7,13,11)
EDGES={4:(0,7),2:(0,3),5:(2,3),3:(0,0),7:(0,2),13:(2,7),11:(0,3)}

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def factors(n):
    ans={};p=2
    while p*p<=n:
        while n%p==0: ans[p]=ans.get(p,0)+1;n//=p
        p+=1
    if n>1:ans[n]=1
    return ans

def sigma(n):
    return prod((p**(e+1)-1)//(p-1) for p,e in factors(n).items())

def square(n):return isqrt(n)**2==n
def signature(n):return frozenset(p for p,e in factors(n).items() if e%2)
def selected(mask):return tuple(x for i,x in enumerate(BANK) if mask>>i&1)
def legal(xs):return all(gcd(a,b)==1 for a,b in combinations(xs,2))

def arithmetic():
    rows=[];zero=[];states=[]
    for mask in range(128):
        xs=selected(mask);n=prod(xs);formal=prod(sigma(x) for x in xs)
        degrees={v:0 for v in (0,2,3,7)}
        for x in xs:
            a,b=EDGES[x];degrees[a]+=1;degrees[b]+=1
        need(square(formal)==all(d%2==0 for d in degrees.values()),'incidence law')
        if legal(xs):
            states.append(mask);need(sigma(n)==formal,'multiplicative guard')
        if square(formal):
            zero.append(mask)
            rows.append({'labels':xs,'n':n,'formal_sigma':formal,'legal':legal(xs),'sigma':sigma(n)})
    circuits=[selected(m) for m in zero if m and not any(k and k!=m and k&m==k for k in zero)]
    actual=sorted(r['n'] for r in rows if r['legal'])
    need(actual==[1,3,22,66,70,210,364,385,1092,1155,2860,8580],'complete bank list')
    need(len(states)==96 and len(zero)==16,'bank counts')
    # Literal guarded append automaton, observing legality and square sigma.
    nums={m:prod(selected(m)) for m in states}
    out={m:int(square(sigma(nums[m]))) for m in states}
    def step(m,i):
        ok=gcd(nums[m],BANK[i])==1
        return ok,m|(1<<i) if ok else m
    cls=out.copy();rounds=[]
    while True:
        keys={};new={}
        for m in states:
            key=(out[m],tuple((step(m,i)[0],cls[step(m,i)[1]]) for i in range(7)))
            new[m]=keys.setdefault(key,len(keys))
        rounds.append(len(keys))
        same=all((cls[x]==cls[y])==(new[x]==new[y]) for x in states for y in states)
        cls=new
        if same:break
    groups={}
    for m,c in cls.items():groups.setdefault(c,[]).append(nums[m])
    merged=sorted(sorted(g) for g in groups.values() if len(g)>1)
    need(len(groups)==86 and len(merged)==10,'minimal automaton')
    def FUS(n):
        f=factors(n);r=len(f)
        if n==1:return True
        F=prod(e+1 for e in f.values())-2
        S=2**r-1-int(all(e==1 for e in f.values()))
        U=r-int(len(f)==1 and next(iter(f.values()))==1)
        return F==U+S
    def refined_count(extra):
        colors={m:(out[m],extra(nums[m])) for m in states}
        while True:
            keys={};new={}
            for m in states:
                key=(out[m],extra(nums[m]),tuple((step(m,i)[0],colors[step(m,i)[1]]) for i in range(7)))
                new[m]=keys.setdefault(key,len(keys))
            same=all((colors[x]==colors[y])==(new[x]==new[y]) for x in states for y in states)
            colors=new
            if same:return len(keys)
    extra_counts={'with_FUS':refined_count(FUS),'with_tau':refined_count(lambda n:prod(e+1 for e in factors(n).values()))}
    need(extra_counts=={'with_FUS':87,'with_tau':96},'observer hierarchy')
    # A complete bounded proof of a global minimum, below a supplied witness.
    first_shape=[n for n in range(2,365) if sorted(factors(n).values())==[1,1,2] and square(sigma(n))]
    need(first_shape==[364],'least p2qr square sigma')
    # Common legal context separating the forward-equivalent 4 and 13.
    contexts=[m for m in range(1,206) if gcd(m,52)==1 and square(sigma(4*m))!=square(sigma(13*m))]
    need(contexts==[205],'least common square context')
    fus=[]
    for n in actual:
        if n==1:continue
        f=factors(n);r=len(f);F=prod(e+1 for e in f.values())-2
        S=2**r-1-int(all(e==1 for e in f.values()));U=r-int(len(f)==1 and next(iter(f.values()))==1)
        if F==U+S:fus.append(n)
    need(fus==[3,364],'FUS bank intersection')
    return {'bank':BANK,'graph_edges':EDGES,'square_signature_rows':rows,'minimal_circuits':circuits,
            'legal_states':96,'formal_kernel_size':16,'square_sources':actual,'refinement_class_counts':rounds,
            'minimal_guarded_states':86,'merged_pairs':merged,'least_p2qr_square_sigma':364,
            'additional_observer_class_counts':extra_counts,
            'least_common_context':{'m':205,'sigma_m':sigma(205),'sigma_820':sigma(820),'sigma_2665':sigma(2665)},
            'FUS_square_bank_sources':fus}

def height_and_trees():
    e={2:5,3:7,5:3,7:4,11:1,13:2};N=prod(p**a for p,a in e.items())
    def height(b):
        count=0;n=N
        while n%b==0:n//=b;count+=1
        return count
    cycle=[4,2,5,3,7]
    need([height(x) for x in cycle]==cycle[1:]+cycle[:1] and height(13)==2,'height wiring')
    divs=[prod(p**a for p,a in zip(e,exp)) for exp in product(*(range(a+1) for a in e.values()))]
    hist={}
    for d in divs:
        if d>1:hist[height(d)]=hist.get(height(d),0)+1
    need(hist=={1:5616,2:120,3:16,4:4,5:2,7:1},'height census')
    mats=(((1,-2,2),(2,-1,2),(2,-2,3)),((1,2,2),(2,1,2),(2,2,3)),((-1,2,2),(-2,1,2),(-2,2,3)))
    def mul(M,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in M)
    nodes=[('',(3,4,5),(1,1,1))];allnodes=[]
    for depth in range(6):
        allnodes+=nodes
        nxt=[]
        for w,p,q in nodes:
            need(p[0]**2+p[1]**2==p[2]**2,'Pythagorean shell')
            need(q[0]**2+q[1]**2==q[2]**2+1,'norm one shell')
            for i,M in enumerate(mats):nxt.append((w+'ABC'[i],mul(M,p),mul(M,q)))
        nodes=nxt
    need(len({p for _,p,_ in allnodes})==364 and len({q for _,_,q in allnodes})==364,'address collisions')
    x,y=364,420;t=isqrt(x+y);a,b,c=(5,5,7)
    u=a*a*(x+y)-x;v=b*b*(x+y)-y
    need((u+x,x+y,y+v,u+v)==((a*t)**2,t*t,(b*t)**2,(c*t)**2),'absorber')
    return {'N':N,'cycle':cycle,'extra_arrow':[13,2],'exact_height_counts':hist,
            'branch_words_through_depth':5,'paired_tree_nodes':len(allnodes),
            'absorber_path':[u,x,y,v],'absorber_edge_square_roots':[a*t,t,b*t],'host_edge_square_root':c*t}

def graph_controls():
    norm=lambda a,b:tuple(sorted((a,b)))
    outer=[norm(i,(i+1)%5) for i in range(5)]
    inner=[norm(i+5,(i+2)%5+5) for i in range(5)]
    spokes=[(i,i+5) for i in range(5)];edges=outer+inner+spokes
    def cycles_of(es,vertices):
        adj={v:set() for v in vertices}
        for a,b in es:adj[a].add(b);adj[b].add(a)
        if any(len(v)!=2 for v in adj.values()):return None
        unseen=set(vertices);sizes=[]
        while unseen:
            stack=[next(iter(unseen))];comp=set()
            while stack:
                x=stack.pop()
                if x in comp:continue
                comp.add(x);stack.extend(adj[x]-comp)
            unseen-=comp;sizes.append(len(comp))
        return sorted(sizes)
    direct=set()
    for chosen in combinations(edges,10):
        shape=cycles_of(chosen,range(10))
        if shape is not None:
            need(shape==[5,5],'Petersen components');direct.add(frozenset(chosen))
    lifted=set()
    for mask in range(1024):
        chosen=[e for i,e in enumerate(outer+inner) if mask>>i&1]
        O=[0]*5;I=[0]*5
        for a,b in chosen:
            row=O if a<5 else I;row[a%5]+=1;row[b%5]+=1
        if all(O[i]==I[i] and O[i] in (1,2) for i in range(5)):
            lifted.add(frozenset(chosen+[spokes[i] for i in range(5) if O[i]==1]))
    need(lifted==direct and len(direct)==6,'full quotient lift')
    left=(2,8,9);right=(3,6,7);K={norm(a,b) for a in left for b in right}
    marks={norm(2,6),norm(3,9),norm(7,8)}
    hc=set()
    vs=sorted(left+right);root=vs[0]
    for tail in permutations(vs[1:]):
        w=(root,)+tail;es=frozenset(norm(w[i],w[(i+1)%6]) for i in range(6))
        if es<=K:hc.add(es)
    need(len(hc)==6 and sum(marks<=c for c in hc)==2,'suppression marks')
    need(sorted(abs(a-b) for a in (0,1,2) for b in (3,6,9))==list(range(1,10)),'K33 graceful')
    need(all(r*(5-r)!=5 for r in range(6)),'K5 graceful obstruction')
    return {'Petersen_2factors':6,'component_sizes':[5,5],'K5_subsets_checked':1024,
            'K33_Hamilton_cycles':6,'marked_cycles_lifting':2,'K33_graceful_labels':[[0,1,2],[3,6,9]]}

def fano_extension():
    # The bit coordinates are valuations at 2,3,5 of sigma, modulo two.
    by_vector={1:7,2:2,3:5,4:19,5:89,6:59,7:29}
    for v,p in by_vector.items():
        need(factors(p)=={p:1},'prime realization')
        need(signature(sigma(p))==frozenset(q for i,q in enumerate((2,3,5)) if v>>i&1),'Fano column')
        need(not any(factors(q)=={q:1} and signature(sigma(q))==signature(sigma(p)) for q in range(2,p)), 'least prime in specified square class')
    words=[];weight_counts={};cosets={s:[] for s in range(8)}
    for mask in range(128):
        vs=[v for v in range(1,8) if mask>>(v-1)&1];syndrome=0
        for v in vs:syndrome^=v
        n=prod(by_vector[v] for v in vs)
        cosets[syndrome].append(frozenset(vs))
        need(square(sigma(n))==(syndrome==0),'Fano arithmetic code')
        if syndrome==0:
            words.append(vs);weight_counts[len(vs)]=weight_counts.get(len(vs),0)+1
    need(weight_counts=={0:1,3:7,4:7,7:1},'Hamming enumerator')
    marginal_tests=0
    for r in range(4):
        for coordinates in combinations(range(1,8),r):
            laws=[Counter(tuple(int(v in w) for v in coordinates) for w in cosets[s]) for s in (0,1)]
            need(laws[0]==laws[1] and set(laws[0].values())=={2**(4-r)},'three-wise marginal hostile')
            marginal_tests+=1
    for s in (0,1):need(Counter(len(w)%2 for w in cosets[s])=={0:8,1:8},'same sign marginals')
    lines={frozenset(w) for w in words if len(w)==3}
    autos=sum({frozenset(perm[v-1] for v in line) for line in lines}==lines for perm in permutations(range(1,8)))
    need(autos==168,'Fano automorphisms')
    # F8 = F2[alpha]/(alpha^3+alpha+1), basis (1,alpha,alpha^2).
    def mul(x,y):
        z=0
        while y:
            if y&1:z^=x
            y>>=1;x<<=1
            if x&8:x^=11
        return z
    powers=[1]
    for _ in range(6):powers.append(mul(powers[-1],2))
    need(powers==[1,2,4,3,6,7,5],'phase convention')
    logs={v:i for i,v in enumerate(powers)}
    def arc(x,y):return (logs[y]-logs[x])%7 in (1,2,4)
    hp=sum(all(arc(w[i],w[i+1]) for i in range(6)) for w in permutations(range(1,8)))
    ta=sum(all(arc(x,y)==arc(w[x-1],w[y-1]) for x in range(1,8) for y in range(1,8) if x!=y) for w in permutations(range(1,8)))
    need((hp,ta)==(189,21),'marked Paley control')
    frob={by_vector[v]:by_vector[mul(v,v)] for v in range(1,8)}
    normal=(3,5,7)
    need([mul(v,v) for v in normal]==[5,7,3],'normal basis orbit')
    values=[]
    for bits in product((0,1),repeat=3):
        value=0;rotated=0
        for bit,v in zip(bits,normal):
            if bit:value^=v
        for bit,v in zip((bits[2],bits[0],bits[1]),normal):
            if bit:rotated^=v
        need(mul(value,value)==rotated,'Frobenius is normal-coordinate rotation');values.append(value)
    need(len(set(values))==8,'normal basis independence')
    return {'prime_by_nonzero_vector':by_vector,'square_subset_weight_enumerator':weight_counts,
            'square_triples':sorted(sorted(by_vector[v] for v in line) for line in lines),
            'Fano_automorphism_order':autos,'chosen_phase_prime_order':[by_vector[v] for v in powers],
            'Frobenius_on_prime_labels':frob,'chosen_Paley_Hamilton_paths':hp,'chosen_Paley_automorphism_order':ta,
            'normal_basis_vectors':normal,'normal_basis_prime_labels':[by_vector[v] for v in normal],
            'uniform_coset_hostile':{'syndromes':[0,1],'sources_per_law':16,'coordinate_subsets_compared':marginal_tests,
                'largest_marginal_order':3,'even_odd_counts_each':[8,8],'square_probabilities':[1,0]},
            'weight_transport_control':{'source':70,'image':7*19*89,'source_sigma':sigma(70),'image_sigma':sigma(7*19*89)}}

if __name__=='__main__':
    print(json.dumps({'arithmetic':arithmetic(),'height_and_Berggren':height_and_trees(),'graph_controls':graph_controls(),'Fano_extension':fano_extension()},indent=2,sort_keys=True))
