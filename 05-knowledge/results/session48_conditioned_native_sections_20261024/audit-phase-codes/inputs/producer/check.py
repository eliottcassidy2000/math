#!/usr/bin/env python3
"""Exact native affine phase-code controls. Standard library only; no asserts."""
from collections import Counter, deque
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def add(a,b): return tuple((x+y)%3 for x,y in zip(a,b))
def sub(a,b): return tuple((x-y)%3 for x,y in zip(a,b))
def dot(a,b): return sum(x*y for x,y in zip(a,b))%3
def wt(v): return sum(x!=0 for x in v)
def dist(a,b): return sum(x!=y for x,y in zip(a,b))

def rref(rows,n):
    a=[list(x) for x in rows]; piv=[]; p=0
    for j in range(n):
        k=next((i for i in range(p,len(a)) if a[i][j]%3),None)
        if k is None: continue
        a[p],a[k]=a[k],a[p]
        c=pow(a[p][j]%3,-1,3); a[p]=[(x*c)%3 for x in a[p]]
        for i in range(len(a)):
            if i!=p:
                c=a[i][j]%3; a[i]=[(x-c*y)%3 for x,y in zip(a[i],a[p])]
        piv.append(j);p+=1
    return tuple(tuple(x) for x in a[:p]),tuple(piv)

def rank(rows,n): return len(rref(rows,n)[0])
def span(rows,n):
    return {tuple(sum(c*v[i] for c,v in zip(cs,rows))%3 for i in range(n)) for cs in product(range(3),repeat=len(rows))}

def solve(A,b,n):
    rows=[list(a)+[c%3] for a,c in zip(A,b)]
    rr,piv=rref(rows,n)
    # rref only returns pivot rows, so independently detect incompatible zero rows.
    if rank(rows,n+1)>len(rr): return None
    free=[j for j in range(n) if j not in piv]
    origin=[0]*n
    for row,p in zip(rr,piv):origin[p]=row[-1]
    basis=[]
    for j in free:
        v=[0]*n;v[j]=1
        for row,p in zip(rr,piv):v[p]=-row[j]%3
        basis.append(tuple(v))
    return tuple(origin),tuple(basis)

def all_spaces(n):
    for d in range(n+1):
        for piv in combinations(range(n),d):
            free=[(i,j) for i,p in enumerate(piv) for j in range(p+1,n) if j not in piv]
            for vals in product(range(3),repeat=len(free)):
                rows=[[0]*n for _ in piv]
                for i,p in enumerate(piv): rows[i][p]=1
                for (i,j),x in zip(free,vals):rows[i][j]=x
                yield tuple(tuple(row) for row in rows)

def low_weight_basis(A,n,B):
    """Select original low-weight vectors, not possibly heavy row-reduced rows."""
    chosen=[];examined=0
    for size in range(1,min(B,n)+1):
        for support in combinations(range(n),size):
            for values in product((1,2),repeat=size):
                v=[0]*n
                for i,z in zip(support,values):v[i]=z
                v=tuple(v);examined+=1
                if all(dot(a,v)==0 for a in A) and rank(chosen+[v],n)>len(chosen):chosen.append(v)
    return tuple(chosen),examined

def basis_route(basis,x,y):
    delta=sub(y,x);d=len(basis);n=len(x)
    transposed=tuple(tuple(basis[j][i] for j in range(d)) for i in range(n))
    solved=solve(transposed,delta,d)
    if solved is None:return None
    coefficients=solved[0];route=[x]
    for c,v in zip(coefficients,basis):
        if c:route.append(add(route[-1],tuple(c*z%3 for z in v)))
    return route

def components(states,B):
    todo=set(states); out=[]
    while todo:
        root=min(todo); todo.remove(root); found={root}; queue=[root]
        while queue:
            x=queue.pop()
            fresh={y for y in todo if dist(x,y)<=B}
            todo-=fresh;found|=fresh;queue.extend(fresh)
        out.append(found)
    return out

P=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2))
G=(4,0,2,3,1,5,6,7)

def source(r):
    core=[[0]*8 for _ in range(8)]
    for z,(i,j) in enumerate(combinations(range(8),2)):
        core[i][j]=(29558785>>(27-z))&1;core[j][i]=1-core[i][j]
    q=8*r
    return tuple(tuple(0 if i==j else core[i%8][j%8] if i<q and j<q and i//8==j//8 else int(i<j) if i<q and j<q else int(i<q) if (i<q)!=(j<q) else int((j-i)%3==1) for j in range(q+3)) for i in range(q+3))

def path(r,key,k):
    return tuple(8*j+x for j in range(r) for x in P[(key>>j)&1])+tuple(8*r+(k+t)%3 for t in range(3))

def partner(a,p):
    pos={u:i for i,u in enumerate(p)}
    return tuple(tuple(0 if u==v else a[u][v] if abs(pos[u]-pos[v])==1 else a[v][u] for v in range(len(a))) for u in range(len(a)))

def word(a,p):return tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i+2,len(p)))
def transpose(a):return tuple(zip(*a))
def iso(a,b,p):return all(a[u][v]==b[p[u]][p[v]] for u in range(len(a)) for v in range(len(a)))

def native_checks():
    reports=[]
    for r in range(1,5):
        a=source(r);q=8*r;n=q+3
        pp={(key,k):path(r,key,k) for key in range(1<<r) for k in range(3)}
        targets={z:partner(a,p) for z,p in pp.items()}
        words={z:word(a,p) for z,p in pp.items()}
        forward=set(words.values());backward={word(transpose(a),p[::-1]) for p in pp.values()}
        reservoirs=forward|backward|{tuple(1-x for x in w) for w in forward|backward}
        need(len(forward)==1<<r,'distinct prefix keys are distinct actual words')
        need(len(reservoirs)==4*(1<<r),'four-word normalized reservoirs are disjoint')
        for key in range(1<<r):
            need(len({words[key,k] for k in range(3)})==1,'all three kernel lifts share one resource')
            prefix=pp[key,0][:-3]
            actual={prefix+t for t in permutations(range(q,q+3)) if all(a[u][v] for u,v in zip(prefix+t,(prefix+t)[1:]))}
            need(actual=={pp[key,k] for k in range(3)},'complete prefix-pinned record menu')
        for (key,k),p in pp.items():
            need(sorted(p)==list(range(n)) and all(a[u][v] for u,v in zip(p,p[1:])),'literal full Hamiltonian path')
            need(p[q]-q==k,'linear guard phase is actual clock-q owner')
            rp=p[::-1]
            need(rp[2]-q==k and rp[0]-q==(k+2)%3,'authenticated converse phase extractors')
            need(partner(transpose(a),rp)==transpose(targets[key,k]),'actual converse partner square')
            for delta in (1,2):
                perm=tuple(u if u<q else q+(u-q+delta)%3 for u in range(n))
                new=(key,(k+delta)%3)
                need(iso(a,a,perm),'actual suffix rotation source automorphism')
                need(tuple(perm[u] for u in p)==pp[new],'actual owner transport to changed phase')
                need(iso(targets[key,k],targets[new],perm),'actual kernel partner transport')
                need(dist(p,pp[new])==3,'one changed phase costs exactly three owner clocks')
            for j in range(r):
                g=G if not(key>>j)&1 else tuple(G.index(u) for u in range(8))
                perm=tuple(8*j+g[u%8] if u//8==j else u for u in range(n))
                need(iso(targets[key,k],targets[key^(1<<j),k],perm),'endpoint-marked partner gluing across all prefix keys')
        reflection=tuple(u if u<q else q+(-(u-q))%3 for u in range(n))
        need(not iso(a,a,reflection),'phase reflection is not a same-source C3 automorphism')
        reports.append(dict(prefix_blocks=r,source_order=n,full_records=len(pp),resources=1<<r,reservoir_words=len(reservoirs)))
    seed=tuple(row[:8] for row in source(1)[:8]);degrees=[sum(row) for row in seed]
    autos=[p for p in permutations(range(8)) if all(degrees[i]==degrees[p[i]] for i in range(8)) and iso(seed,seed,p)]
    need(autos==[tuple(range(8))],'rigid seed checked against all eight-factorial permutations')
    for start in range(8):
        todo=[start];seen={start}
        while todo:
            u=todo.pop()
            fresh={v for v in range(8) if seed[u][v] and v not in seen};seen|=fresh;todo.extend(fresh)
        need(len(seen)==8,'seed is a strong component')
    return dict(native=reports,seed_automorphisms=len(autos))

def no_return_dag(C,x,y,B=None):
    """Exact minimum bottleneck, or minimum event count under B, plus route."""
    delta=sub(y,x);support=[i for i,z in enumerate(delta) if z]
    nodes={}
    for stage in product(range(3),repeat=len(support)):
        z=list(x)
        for i,t in zip(support,stage):z[i]=(x[i]+(0,-delta[i],delta[i])[t])%3
        if tuple(z) in C:nodes[stage]=tuple(z)
    stages=sorted(nodes,key=lambda t:(sum(t),t));start=(0,)*len(support);end=(2,)*len(support)
    best={start:0};prev={}
    for a in stages:
        if a not in best:continue
        for b in stages:
            if a==b or any(u>v for u,v in zip(a,b)):continue
            cost=dist(nodes[a],nodes[b])
            if B is not None and cost>B:continue
            value=best[a]+1 if B is not None else max(best[a],cost)
            if b not in best or value<best[b]:best[b]=value;prev[b]=a
    if end not in best:return None,[]
    route=[];z=end
    while z!=start:route.append(nodes[z]);z=prev[z]
    route.append(x);route.reverse()
    return best[end],route

def no_return_bfs(C,x,y,B):
    """Independent literal visited-record-mask search."""
    n=len(x);start=(x,tuple(1<<v for v in x));seen={start};todo=deque([(start,0)])
    while todo:
        (z,masks),length=todo.popleft()
        if z==y:return length
        for w in C:
            changed=[i for i in range(n) if z[i]!=w[i]]
            if not 0<len(changed)<=B or any(masks[i]>>w[i]&1 for i in changed):continue
            nxt=(w,tuple(mask|(1<<v) for mask,v in zip(masks,w)))
            if nxt not in seen:seen.add(nxt);todo.append((nxt,length+1))
    return None

def affine_checks():
    rows=[];nr_cases=0
    for n in range(5):
        ambient=set(product(range(3),repeat=n));spaces=cosets=graphs=0
        r=max(1,(n-1).bit_length());pp=[[path(r,i,k) for k in range(3)] for i in range(n)]
        for G0 in all_spaces(n):
            C=span(G0,n);d=len(G0);spaces+=1
            dual=[z for z in ambient if all(dot(g,z)==0 for g in G0)]
            A=rref(dual,n)[0]
            need(rank(A,n)==n-d,'annihilator has complementary rank')
            model=solve(A,(0,)*len(A),n)
            need(model is not None and span(model[1],n)==C,'Gaussian kernel equals independently generated code')
            enum=Counter(wt(v) for v in C)
            for B in range(n+1):
                short=[v for v in C if 0<wt(v)<=B];rho=rank(short,n)
                cc=components(C,B)
                original_basis,examined=low_weight_basis(A,n,B)
                from math import comb
                need(examined==sum(comb(n,j)*2**j for j in range(1,B+1)),'constructive algorithm enumerates only bounded-support phase vectors')
                need(len(original_basis)==rho and all(wt(v)<=B for v in original_basis),'routing basis retains original low-weight generators')
                for delta in C:
                    route=basis_route(original_basis,(0,)*n,delta)
                    need((route is not None)==(delta in span(original_basis,n)),'constructive endpoint membership and coefficient recovery')
                    if route is not None:
                        need(route[-1]==delta and len(route)-1<=rho,'at most rank basis-generator events')
                        need(all(z in C for z in route) and all(dist(x,y)<=B for x,y in zip(route,route[1:])),'actual affine state space keeps every basis route legal')
                        need(sum(3*dist(x,y) for x,y in zip(route,route[1:]))<=3*B*rho,'basis-route owner-clock work bound')
                pivots=rref(G0,n)[1];eigenvalues=[]
                for character in product(range(3),repeat=d):
                    counts=Counter(sum(character[j]*v[pivots[j]] for j in range(d))%3 for v in short)
                    need(counts[1]==counts[2],'general budget generator character sum is an exact real integer')
                    eigenvalue=len(short)-counts[0]+counts[1]
                    need(eigenvalue==3*counts[1] and eigenvalue>=0,'all ternary budget Laplacian eigenvalues are multiples of three')
                    eigenvalues.append(eigenvalue)
                need(eigenvalues.count(0)==3**(d-rho),'general character zero multiplicity equals component count')
                need(sum(eigenvalues)==3**d*len(short),'general budget Laplacian exact trace')
                need(len(cc)==3**(d-rho),'exact low-weight-span component count')
                need({len(c) for c in cc}=={3**rho},'all components have the predicted affine size')
                need({frozenset(sub(y,min(c)) for y in c) for c in cc}=={frozenset(span(rref(short,n)[0],n))},'each graph component is a low-weight-span coset')
                graphs+=1
            unseen=set(ambient)
            while unseen:
                base=min(unseen);actual_coset={add(base,c) for c in C};unseen-=actual_coset;cosets+=1
                b=tuple(dot(a,base) for a in A)
                literal={ks for ks in ambient if all(sum(a[i]*(pp[i][ks[i]][8*r]-8*r) for i in range(n))%3==t for a,t in zip(A,b))}
                need(literal==actual_coset,'literal clock-owner linear predicates give exactly affine code')
                solved=solve(A,b,n)
                need(solved is not None and {add(solved[0],v) for v in span(solved[1],n)}==literal,'Gaussian affine model matches literal records')
                need(Counter(dist(base,k) for k in literal)==enum,'affine distance enumerator is kernel weight enumerator')
                for ks in literal:
                    reflected=tuple(pp[i][ks[i]][::-1][0]-8*r for i in range(n))
                    need(all(dot(a,reflected)==(t+2*sum(a))%3 for a,t in zip(A,b)),'converse head-extractor guard is A ell=b+2 A1')
            if n<=3:
                zero=(0,)*n
                for y in C:
                    best,route=no_return_dag(C,zero,y)
                    for B in range(n+1):
                        count,rte=no_return_dag(C,zero,y,B)
                        literal_count=no_return_bfs(C,zero,y,B)
                        need(count==literal_count,'no-return DAG agrees with independent visited-record search')
                        need((count is not None)==(best<=B),'minimum no-return bottleneck matches all budget decisions')
                        nr_cases+=1
        rows.append(dict(slots=n,linear_codes=spaces,affine_fibres=cosets,budget_graphs=graphs))
    # Explicit contradictory systems, duplicate rows and a nonzero affine right side.
    for A,b,n in [(((0,0),),(1,),2),(((1,0),(1,0)),(0,1),2),(((1,2),(2,1)),(1,1),2)]:
        need(solve(A,b,n) is None,'inconsistent affine system rejected, never negative exponent counted')
    return dict(catalogues=rows,no_return_budget_endpoint_tests=nr_cases)

def mode_code(L):
    return {tuple(x for j in range(L) for x in (a[j],(a[j]+m)%3))+(m,) for a in product(range(3),repeat=L) for m in range(3)}

def mode_generators(L):
    same=[];bridges=[]
    for j in range(L):
        for delta in (1,2):same.append(tuple(delta if i==j else 0 for i in range(L))+(0,))
    for delta in (1,2):
        bridges.extend(tuple(a)+(delta,) for a in product((0,(-delta)%3),repeat=L))
    return same+bridges

def mode_phase(z):
    a=z[:-1];m=z[-1]
    return tuple(x for x0 in a for x in (x0,(x0+m)%3))+(m,)

def mode_checks():
    reports=[]
    for L in range(1,6):
        C=mode_code(L);n=2*L+1;gs=mode_generators(L)
        need(len(C)==3**(L+1),'all three helper phases retained')
        enum=Counter(wt(x) for x in C)
        # W=(1+2z^2)^L+2z(2z+z^2)^L.
        expected=Counter()
        from math import comb
        for j in range(L+1):
            expected[2*j]+=comb(L,j)*2**j
            expected[L+1+j]+=2*comb(L,j)*2**(L-j)
        need(enum==expected,'all-mode exact weight enumerator')
        need(min(wt(x) for x in C if x[-1])==L+1,'first mode-changing packet costs exactly L+1 slots')
        need({len(c) for c in components(C,2)}==({9} if L==1 else {3**L}),'budget2 mode-fibre sizes')
        need(len(components(C,2))==(1 if L==1 else 3),'budget2 has exactly three modes for L>=2')
        need(len(gs)==2*L+2**(L+1) and len(set(gs))==len(gs),'selected phase graph degree')
        need(all(mode_phase(g) in C for g in gs),'selected generators are actual accepted phase differences')
        need({wt(mode_phase(g)) for g in gs}=={2,L+1},'selected generator exact native slot costs')
        full={x for x in C if 0<wt(x)<=L+1}
        need(({mode_phase(g) for g in gs}==full)==(L<=2),'selected graph differs from full budget graph exactly when extra moves enter')
        spectrum=Counter()
        for r in product(range(3),repeat=L):
            s=wt(r)
            for t in range(3):
                chi=r+(t,)
                counts=Counter(dot(chi,g) for g in gs)
                need(counts[1]==counts[2],'real character sum proved by inversion pairing')
                exact_adjacency=counts[0]-counts[1]
                cross=(2**(L-s+1)*(-1)**s if (t+sum(r))%3==0 else -2**(L-s)*(-1)**s)
                formula=3*s+2**(L+1)-cross
                need(len(gs)-exact_adjacency==formula,'exact integer Fourier eigenvalue agrees with formula')
                spectrum[formula]+=1
        need(sum(spectrum.values())==len(C) and spectrum[0]==1,'Fourier spectrum covers every character, connected')
        gap=min(v for v in spectrum if v)
        need(gap==3+3*2**(L-1),'all-size spectral-gap formula on stated finite range')
        need(spectrum[gap]==(6 if L==1 else 4*L),'spectral-gap multiplicity')
        reports.append(dict(L=L,resources=n,states=len(C),weight_enumerator=dict(sorted(enum.items())),selected_degree=len(gs),full_budget_Lplus1_degree=len(full),spectral_gap=gap,spectrum=dict(sorted(spectrum.items()))))
    L=2;C=mode_code(L);n=5;r=3;q=24;a=source(r)
    pp={(i,k):path(r,i,k) for i in range(n) for k in range(3)}
    ww={(i,k):word(a,p) for (i,k),p in pp.items()}
    accepted=[]
    for ks in product(range(3),repeat=n):
        owners=[pp[i,ks[i]][q]-q for i in range(n)]
        if all((owners[2*j+1]-owners[2*j]-owners[-1])%3==0 for j in range(L)):accepted.append(ks)
    need(set(accepted)==C,'27-state example checked directly against actual suffix owners')
    for ks in C:
        need(len({ww[i,k] for i,k in enumerate(ks)})==n,'all five distinct actual normalized resources occupied')
        for kt in C:
            clocks=sum(dist(pp[i,ks[i]],pp[i,kt[i]]) for i in range(n))
            need(clocks==3*dist(ks,kt),'actual owner-clock packet work is three times slot support')
    need(dict(reports[1]['spectrum'])=={0:1,9:8,12:6,15:12},'27-state Laplacian spectrum')
    # Full adjacency matrix validation without numerical eigensolvers: exact traces and annihilating polynomial.
    states=sorted(C);N=len(states)
    lap=[[(-1 if i!=j and dist(x,y)<=3 else 12 if i==j else 0) for j,y in enumerate(states)] for i,x in enumerate(states)]
    def mul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
    identity=[[int(i==j) for j in range(N)] for i in range(N)]
    poly=identity
    for lam in (0,9,12,15):poly=mul(poly,[[lap[i][j]-(lam if i==j else 0) for j in range(N)] for i in range(N)])
    need(not any(any(row) for row in poly),'literal 27-by-27 Laplacian satisfies exact split polynomial')
    power=identity
    traces=[]
    for exponent in range(4):
        actual=sum(power[i][i] for i in range(N));expected=sum(m*lam**exponent for lam,m in reports[1]['spectrum'].items())
        need(actual==expected,'exact Laplacian traces fix stated multiplicities');traces.append(actual);power=mul(power,lap)
    from fractions import Fraction
    normalized=[]
    for L in range(1,65):
        value=Fraction(3+3*2**(L-1),2*L+2**(L+1))
        need(value>=Fraction(27,40) and (value==Fraction(27,40))==(L==4),'normalized Laplacian gap sharp lower bound')
        if L<=8:normalized.append(dict(L=L,numerator=value.numerator,denominator=value.denominator))
    return dict(families=reports,L2_direct_matrix_trace_powers_0_to_3=traces,normalized_gaps=normalized)

def no_return_control():
    n=5;v=(1,0,0,1,1);w=(0,1,1,0,1);C=span((v,w),n);zero=(0,)*n;target=sub(v,w)
    low={x for x in C if 0<wt(x)<=3}
    need(low=={v,w,tuple(-x%3 for x in v),tuple(-x%3 for x in w)},'only signed v,w are budget3 generators')
    need(len(components(C,3))==1 and len(components(C,2))==9,'ordinary connectivity threshold three')
    need(all(x[-1]!=0 for x in low) and target[-1]==0,'every budget3 departure changes an endpoint-fixed coordinate')
    best,rte=no_return_dag(C,zero,target)
    need(best==4 and no_return_bfs(C,zero,target,3) is None,'sharp full-record no-return budget four')
    need(no_return_bfs(C,zero,target,4)==1 and rte==[zero,target],'direct endpoint packet attains no-return optimum')
    A=((0,1,2,0,0),(1,0,0,2,0),(2,2,0,0,1))
    need({x for x in product(range(3),repeat=n) if all(dot(a,x)==0 for a in A)}==C,'explicit native linear relation matrix for no-return control')
    r=3;a=source(r)
    ordinary=[zero,v,target]
    for z in ordinary:
        need(len({word(a,path(r,i,k)) for i,k in enumerate(z)})==5,'ordinary return route keeps five resources distinct')
    need([z[-1] for z in ordinary]==[0,1,0],'budget3 route visibly revisits helper full record')
    need(all(dist(x,y)<=3 for x,y in zip(ordinary,ordinary[1:])),'ordinary two-packet route budget three')
    return dict(resources=5,source_order=27,states=9,matrix=A,basis=[v,w],endpoint=target,ordinary_minimum_budget=3,no_return_minimum_budget=4,ordinary_route=ordinary,no_return_route=rte)

def coordinate_checks():
    n=3;A=((1,2,0),(0,1,1));b=(2,1);ambient=list(product(range(3),repeat=n));cases=0
    for p in permutations(range(n)):
        for scales in product((1,2),repeat=n):
            for shift in product(range(3),repeat=n):
                def convert(z):return tuple((scales[i]*z[p[i]]+shift[i])%3 for i in range(n))
                newA=tuple(tuple(sum(a[i]*scales[i] for i in range(n) if p[i]==j)%3 for j in range(n)) for a in A)
                newb=tuple((t-dot(a,shift))%3 for a,t in zip(A,b))
                need({z for z in ambient if all(dot(a,convert(z))==t for a,t in zip(A,b))}=={z for z in ambient if all(dot(a,z)==t for a,t in zip(newA,newb))},'affine monomial reparameterization transforms guard exactly')
                need(all(dist(convert((0,)*n),convert(z))==wt(z) for z in ambient),'monomial coordinate change preserves literal support cost')
                cases+=1
    need(wt((1,0))!=wt((1,1)),'general invertible linear maps need not preserve slot weight')
    return dict(affine_monomial_coordinate_cases=cases,native_claim='Only actual matrix/path transports, not arbitrary coordinate naming, certify native symmetries.')

def raw_owner_translation_controls():
    cases=[]
    for r,a,beta in ((1,(1,0),0),(2,(1,0),0),(1,(1,1),1),(3,(1,0),0)):
        q=8*r;shift=q*sum(a)%3;good=[];bad=[]
        for k in product(range(3),repeat=2):
            owners=tuple(path(r,i,t)[q] for i,t in enumerate(k))
            phases=tuple((o-q)%3 for o in owners)
            head_phases=tuple((path(r,i,t)[::-1][0]-q)%3 for i,t in enumerate(k))
            need(phases==k,'canonical phase extractor subtracts q from raw owners')
            need((dot(a,phases)==beta)==(dot(a,owners)==(beta+shift)%3),'raw-owner RHS transports by q A1')
            need((dot(a,phases)==beta)==(dot(a,head_phases)==(beta+2*sum(a))%3),'converse-head phase RHS transports separately by 2 A1')
            if dot(a,phases)==beta:good.append(k)
            if dot(a,owners)==beta:bad.append(k)
        need((good!=bad)==(shift!=0),'same-RHS raw-owner equation changes the accepted fibre when shift is nonzero')
        cases.append(dict(r=r,row=a,canonical_rhs=beta,raw_rhs=(beta+shift)%3,converse_head_phase_rhs=(beta+2*sum(a))%3,canonical_states=good,unshifted_raw_states=bad))
    need(cases[0]['canonical_states']==[(0,0),(0,1),(0,2)] and cases[0]['unshifted_raw_states']==[(1,0),(1,1),(1,2)],'r1 sum-one hostile control')
    need(cases[1]['raw_rhs']==1 and cases[1]['converse_head_phase_rhs']==2,'r2 distinguishes raw-owner from converse-head shifts')
    return dict(formula='A o=b+(8r)A1 modulo3 in the declared suffix numbering',cases=cases)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    receipt=dict(schema='native-affine-phase-codes-v1',native=native_checks(),affine=affine_checks(),modes=mode_checks(),no_return=no_return_control(),coordinates=coordinate_checks(),raw_owner_coordinates=raw_owner_translation_controls())
    receipt['checks_before_receipt_comparison']=CHECKS
    canonical=json.dumps(receipt,sort_keys=True,indent=2)+'\n'
    if args.write:(HERE/'receipt.json').write_text(canonical)
    else:need(json.loads((HERE/'receipt.json').read_text())==json.loads(canonical),'frozen deterministic receipt')
    print('PASS native affine phase codes;',CHECKS,'explicit checks; all outputs exact integers')
    print('27-state spectrum: 0^1, 9^8, 12^6, 15^12; no-return control: budget3 connected, budget4 sharp')
    print('Receipt SHA256:',hashlib.sha256(canonical.encode()).hexdigest())

if __name__=='__main__':main()
