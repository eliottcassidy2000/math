#!/usr/bin/env python3
"""Independent finite audit. No producer imports, assertions, or third-party modules."""
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import argparse, json

HERE=Path(__file__).resolve().parent
CHECKS=0

def require(value, message):
    global CHECKS
    CHECKS+=1
    if not value: raise RuntimeError(message)

def plus(x,y): return tuple((a+b)%3 for a,b in zip(x,y))
def times(c,x): return tuple(c*a%3 for a in x)
def minus(x,y): return plus(x,times(2,y))
def inner(x,y): return sum(a*b for a,b in zip(x,y))%3
def weight(x): return sum(a!=0 for a in x)
def distance(x,y): return sum(a!=b for a,b in zip(x,y))

def extend(space,v):
    return frozenset(plus(x,times(c,v)) for x in space for c in range(3))

def original_basis(vectors,n):
    basis=[]; current=frozenset({(0,)*n})
    for v in vectors:
        if v not in current:
            basis.append(v);current=extend(current,v)
    return tuple(basis),current

def all_subspaces(n):
    """BFS adjoining actual ambient vectors, not an RREF catalogue."""
    ambient=tuple(product(range(3),repeat=n))
    zero=frozenset({(0,)*n}); seen={zero}; queue=[zero]
    for space in queue:
        for v in ambient:
            if v not in space:
                larger=extend(space,v)
                if larger not in seen:seen.add(larger);queue.append(larger)
    return sorted(seen,key=lambda s:(len(s),sorted(s)))

def coefficients(basis,n):
    return {tuple(sum(c*v[i] for c,v in zip(cs,basis))%3 for i in range(n)):cs
            for cs in product(range(3),repeat=len(basis))}

def components(points,B):
    remaining=set(points); output=[]
    while remaining:
        seen={min(remaining)};queue=list(seen);remaining-=seen
        for x in queue:
            fresh={y for y in remaining if distance(x,y)<=B}
            remaining-=fresh;seen|=fresh;queue.extend(fresh)
        output.append(frozenset(seen))
    return output

def source_matrix(r):
    core=[[0]*8 for _ in range(8)]
    bits=f'{29558785:028b}'
    for bit,(i,j) in zip(bits,combinations(range(8),2)):
        core[i][j]=int(bit);core[j][i]=1-int(bit)
    blocks=[range(8*j,8*j+8) for j in range(r)]+[range(8*r,8*r+3)]
    n=8*r+3; a=[[0]*n for _ in range(n)]
    for bi,block in enumerate(blocks):
        for u,v in permutations(block,2):
            a[u][v]=core[u%8][v%8] if bi<r else int((v-u)%3==1)
        for later in blocks[bi+1:]:
            for u in block:
                for v in later:a[u][v]=1
    return tuple(map(tuple,a))

PREFIXES=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2))
GLUE=(4,0,2,3,1,5,6,7)

def record(r,key,phase):
    return tuple(8*j+u for j in range(r) for u in PREFIXES[(key>>j)&1])+tuple(8*r+(phase+t)%3 for t in range(3))

def partner(a,p):
    keep=set(zip(p,p[1:]))
    return tuple(tuple(0 if i==j else a[i][j] if (i,j) in keep or (j,i) in keep else a[j][i] for j in range(len(a))) for i in range(len(a)))

def opposite(a):return tuple(zip(*a))
def full_word(a,p):return tuple(a[u][v] for i,u in enumerate(p) for v in p[i+1:])
def maps(a,b,f):return all(a[i][j]==b[f[i]][f[j]] for i in range(len(a)) for j in range(len(a)))

def native_checks():
    output=[]
    for r in range(1,6):
        a=source_matrix(r);q=8*r;reservoirs={};allpaths={}
        for key in range(2**r):
            paths=[record(r,key,k) for k in range(3)];allpaths[key]=paths
            prefixes=paths[0][:-3]
            complete={prefixes+t for t in permutations(range(q,q+3)) if all(a[u][v] for u,v in zip(prefixes+t,(prefixes+t)[1:]))}
            require(complete==set(paths),'complete pinned suffix menu')
            actual=set()
            for k,p in enumerate(paths):
                require(len(set(p))==q+3 and all(a[u][v] for u,v in zip(p,p[1:])),'literal source Hamiltonian path')
                b=partner(a,p)
                actual|={full_word(a,p),full_word(b,p),full_word(opposite(a),p[::-1]),full_word(opposite(b),p[::-1])}
                require(partner(opposite(a),p[::-1])==opposite(b),'converse K square')
                require((p[::-1][2]-q,p[::-1][0]-q)==(k,(k+2)%3),'two converse phase extractors')
                for d in (1,2):
                    f=tuple(range(q))+tuple(q+(t+d)%3 for t in range(3))
                    moved=paths[(k+d)%3]
                    require(maps(a,a,f) and tuple(f[u] for u in p)==moved,'actual suffix source transport')
                    require(maps(b,partner(a,moved),f),'actual suffix partner transport')
                    require(distance(p,moved)==3,'exact three-clock support')
                for j in range(r):
                    g=GLUE if not (key>>j)&1 else tuple(GLUE.index(i) for i in range(8))
                    f=tuple(8*j+g[i-8*j] if 8*j<=i<8*j+8 else i for i in range(q+3))
                    target=record(r,key^(1<<j),k)
                    require(maps(b,partner(a,target),f),'endpoint-marked prefix partner transport')
                    require(f[p[8*j]]==target[8*j] and f[p[8*j+7]]==target[8*j+7],'glue preserves both actual block endpoints')
            require(len(actual)==4,'three lifts have one four-word reservoir')
            reservoirs[key]=actual
        require(len(set().union(*reservoirs.values()))==4*2**r,'distinct normalized resources are disjoint')
        reflection=tuple(range(q))+tuple(q+(-t)%3 for t in range(3))
        require(not maps(a,a,reflection),'arbitrary reflection is not native same-source symmetry')
        require(reservoirs[0]=={full_word(a,record(r,0,k)) for k in range(3)}|{full_word(partner(a,record(r,0,k)),record(r,0,k)) for k in range(3)}|{full_word(opposite(a),record(r,0,k)[::-1]) for k in range(3)}|{full_word(opposite(partner(a,record(r,0,k))),record(r,0,k)[::-1]) for k in range(3)},'duplicate key phases cannot create capacity')
        output.append({'r':r,'order':q+3,'resources':2**r,'records':3*2**r,'full_words':4*2**r})
    seed=tuple(row[:8] for row in source_matrix(1)[:8]);degrees=list(map(sum,seed))
    autos=[p for p in permutations(range(8)) if all(degrees[i]==degrees[p[i]] for i in range(8)) and maps(seed,seed,p)]
    require(autos==[tuple(range(8))],'seed rigid by all 8! permutations')
    # Hostile named partners: same isomorphism class does not give a common K host.
    a=source_matrix(1);p=record(1,0,0);q=record(1,1,0)
    require(p[:2]==(3,0) and q[:2]==(3,1),'actual conflicting observer successors')
    require(partner(a,p)[3][0]!=partner(a,q)[3][0],'equal target class is not same named partner')
    return {'carriers':output,'seed_automorphisms':len(autos)}

def no_return_distances(C,B):
    """Literal full-record histories. Each coordinate starts at phase zero."""
    n=len(next(iter(C)));zero=(0,)*n
    initial=(zero,(1,)*n);seen={initial};queue=deque([(initial,0)]);best={zero:0}
    while queue:
        (x,history),steps=queue.popleft()
        for y in C:
            changed=[i for i in range(n) if x[i]!=y[i]]
            if not 0<len(changed)<=B or any(history[i]&(1<<y[i]) for i in changed):continue
            state=(y,tuple(mask|(1<<phase) for mask,phase in zip(history,y)))
            if state in seen:continue
            seen.add(state);queue.append((state,steps+1));best.setdefault(y,steps+1)
    return best,len(seen)

def monotone_distance(C,y,B):
    zero=(0,)*len(y);stages={}
    for x in C:
        if any(t==0 and s!=0 for s,t in zip(x,y)):continue
        stages[x]=tuple(0 if s==0 else 2 if s==t else 1 for s,t in zip(x,y))
    order=sorted(stages,key=lambda x:sum(stages[x]));best={zero:0}
    for x in order:
        if x not in best:continue
        for z in order:
            if x==z or distance(x,z)>B:continue
            if all(a<=b for a,b in zip(stages[x],stages[z])):
                best[z]=min(best.get(z,10**9),best[x]+1)
    return best.get(y)

def integer_spectrum(C,basis,B):
    coords=coefficients(basis,len(next(iter(C))))
    moves=[g for g in C if 0<weight(g)<=B];spectrum=Counter()
    for chi in product(range(3),repeat=len(basis)):
        residues=Counter(inner(chi,coords[g]) for g in moves)
        require(residues[1]==residues[2],'inverse generator pairing makes spectrum integral')
        eigenvalue=len(moves)-residues[0]+residues[1]
        require(eigenvalue>=0 and eigenvalue%3==0,'all full-budget eigenvalues nonnegative multiples of three')
        spectrum[eigenvalue]+=1
    N=len(C);degree=len(moves);move_set=set(moves)
    triangles=sum(times(2,plus(g,h)) in move_set for g in moves for h in moves)
    direct=[N,N*degree,N*(degree*degree+degree),N*(degree**3+3*degree*degree-triangles)]
    require([sum(mult*lam**j for lam,mult in spectrum.items()) for j in range(4)]==direct,'independent graph trace moments 0..3')
    if N<=27:
        # Every Fourier column is checked against the complete literal adjacency matrix.
        for chi in product(range(3),repeat=len(basis)):
            residues=Counter(inner(chi,coords[g]) for g in moves)
            lam=degree-residues[0]+residues[1]
            character={x:((1,0),(0,1),(-1,-1))[inner(chi,cs)] for x,cs in coords.items()}
            for x in C:
                neighbours=[y for y in C if 0<distance(x,y)<=B]
                lhs=tuple(degree*character[x][j]-sum(character[y][j] for y in neighbours) for j in range(2))
                require(lhs==tuple(lam*a for a in character[x]),'literal matrix times exact cyclotomic character')
    return spectrum

def affine_code_checks():
    output=[];nr_endpoints=0;affine_graphs=0;spectral_graphs=0
    for n in range(5):
        ambient=set(product(range(3),repeat=n));spaces=all_subspaces(n);coset_count=0;graphs=0
        for C in spaces:
            basis,regenerated=original_basis(sorted(C),n);require(regenerated==C,'independent closure basis')
            d=len(basis);dual=[u for u in ambient if all(inner(u,v)==0 for v in C)]
            A,_=original_basis(sorted(dual),n)
            require(len(A)+d==n,'dual dimension')
            require({v for v in ambient if all(inner(a,v)==0 for a in A)}==set(C),'explicit literal linear kernel')
            for B in range(n+1):
                selected=[];current=frozenset({(0,)*n});candidate_count=0
                for size in range(1,B+1):
                    for support in combinations(range(n),size):
                        for values in product((1,2),repeat=size):
                            v=[0]*n
                            for i,t in zip(support,values):v[i]=t
                            v=tuple(v);candidate_count+=1
                            if v in C and v not in current:selected.append(v);current=extend(current,v)
                require(candidate_count==sum(comb(n,j)*2**j for j in range(1,B+1)),'fixed-budget candidate bound')
                require(all(0<weight(v)<=B for v in selected),'retain ORIGINAL low-weight basis')
                cc=components(C,B);rho=len(selected)
                require(len(cc)==3**(d-rho) and all(len(c)==3**rho for c in cc),'component count and sizes')
                require(all(frozenset(minus(x,min(c)) for x in c)==current for c in cc),'each component is exact low-weight-span coset')
                for endpoint,cs in coefficients(selected,n).items():
                    path=[(0,)*n]
                    for c,v in zip(cs,selected):
                        if c:path.append(plus(path[-1],times(c,v)))
                    require(path[-1]==endpoint and len(path)-1<=rho,'rank-many original-generator route')
                    require(all(x in C for x in path) and all(distance(x,y)<=B for x,y in zip(path,path[1:])),'every intermediate route state legal')
                    require(sum(3*distance(x,y) for x,y in zip(path,path[1:]))<=3*B*rho,'owner-clock total work bound')
                sp=integer_spectrum(C,basis,B)
                require(sp[0]==len(cc),'Fourier zero multiplicity equals components')
                spectral_graphs+=1
                literal,length=no_return_distances(C,B)
                for y in C:
                    require(literal.get(y)==monotone_distance(C,y,B),'literal visited-record history equals ternary stage DAG')
                    nr_endpoints+=1
                graphs+=1
            unseen=set(ambient)
            while unseen:
                origin=min(unseen);F={plus(origin,c) for c in C};unseen-=F;coset_count+=1
                b=tuple(inner(a,origin) for a in A)
                # Actual owner extraction at clock 8r, with distinct prefix keys.
                r=max(1,(n-1).bit_length());q=8*r
                extracted={ks for ks in ambient if all(sum(a[i]*(record(r,i,ks[i])[q]-q) for i in range(n))%3==t for a,t in zip(A,b))}
                require(extracted==F,'affine equations equal literal record-owner predicates')
                require(Counter(distance(origin,x) for x in F)==Counter(weight(c) for c in C),'based affine distance enumerator')
                for B in range(n+1):
                    actual=components(F,B)
                    short=[c for c in C if 0<weight(c)<=B];bb,subspace=original_basis(sorted(short),n)
                    require(len(actual)==3**(d-len(bb)),'affine component count including singleton and zero budget')
                    require(all(frozenset(minus(x,min(part)) for x in part)==subspace for part in actual),'literal affine graph components')
                    affine_graphs+=1
                for x in F:
                    head=tuple(record(r,i,k)[::-1][0]-q for i,k in enumerate(x))
                    require(all(inner(a,head)==(t+2*sum(a))%3 for a,t in zip(A,b)),'converse head changes right side by 2 A1')
        output.append({'slots':n,'subspaces':len(spaces),'affine_fibres':coset_count,'kernel_budget_graphs':graphs})
    # Explicit basis-cost hostile case: reducing the chosen rows creates weight four.
    u=(1,1,1,0,0);v=(1,0,0,1,1);heavy=minus(u,v)
    require(weight(u)==weight(v)==3 and weight(heavy)==4,'row reduction need not retain packet budget')
    require((heavy[0],heavy[1])==(0,1),'heavy row is the second reduced pivot row')
    return {'catalogues':output,'affine_budget_graphs':affine_graphs,'fourier_graphs':spectral_graphs,'no_return_endpoint_budget_tests':nr_endpoints,'rref_cost_hostile':{'original':[u,v],'reduced':[v,heavy]}}

def arbitrary_system_checks():
    tested=empty=0
    for n in range(4):
        ambient=tuple(product(range(3),repeat=n))
        for rows in range(3):
            for flat in product(range(3),repeat=n*rows):
                A=tuple(tuple(flat[i*n:(i+1)*n]) for i in range(rows))
                C={x for x in ambient if all(inner(a,x)==0 for a in A)}
                for b in product(range(3),repeat=rows):
                    F={x for x in ambient if tuple(inner(a,x) for a in A)==b};tested+=1
                    if not F:empty+=1;continue
                    x=min(F)
                    require(F=={plus(x,c) for c in C},'arbitrary rectangular/redundant affine system translation')
    require({x for x in product(range(3),repeat=2) if inner((0,0),x)==1}==set(),'zero-row contradiction empty')
    return {'rectangular_systems':tested,'inconsistent_systems':empty}

def mode_embedding(z):
    return tuple(v for a in z[:-1] for v in (a,(a+z[-1])%3))+(z[-1],)

def mode_checks():
    output=[]
    for L in range(1,7):
        params=tuple(product(range(3),repeat=L+1));C={mode_embedding(z) for z in params}
        shifts=[]
        for j in range(L):
            for delta in (1,2):shifts.append(tuple(delta if i==j else 0 for i in range(L))+(0,))
        for delta in (1,2):
            shifts.extend(tuple(v)+(delta,) for v in product((0,(-delta)%3),repeat=L))
        require(len(C)==3**(L+1),'all three helper phases, exact count')
        enum=Counter(weight(x) for x in C);formula=Counter()
        for j in range(L+1):
            formula[2*j]+=comb(L,j)*2**j
            formula[L+1+j]+=2*comb(L,j)*2**(L-j)
        require(enum==formula,'mode code exact weight polynomial')
        require(min(weight(x) for x in C if x[-1])==L+1,'sharp first mode bridge')
        selected={mode_embedding(z) for z in shifts};full={x for x in C if 0<weight(x)<=L+1}
        degree=2*L+2**(L+1)
        require(len(selected)==degree and selected<=full,'designed graph native legal generators')
        require((selected==full)==(L<=2),'designed is not full-budget for L>=3')
        spectrum=Counter()
        for chi in params:
            k=weight(chi[:-1]);h=(chi[-1]+sum(chi[:-1]))%3
            cross=(-1)**k*2**(L-k)*(2 if h==0 else -1)
            lam=3*k+2**(L+1)-cross
            counts=Counter(inner(chi,g) for g in shifts)
            require(counts[1]==counts[2] and degree-counts[0]+counts[1]==lam,'exact designed Fourier character sum')
            spectrum[lam]+=1
        require(spectrum[0]==1 and min(x for x in spectrum if x)==3+3*2**(L-1),'designed spectral gap')
        require(spectrum[3+3*2**(L-1)]==(6 if L==1 else 4*L),'gap eigenspace multiplicity')
        if L<=3:
            basis,_=original_basis(sorted(C),2*L+1)
            for B in (2,L,L+1):
                cc=components(C,B)
                require(len(cc)==(1 if B>=L+1 else len(C) if B<2 else 3),'mode components at relevant budgets')
            if L==2:
                require(all({(x[i],x[j]) for x in C}==set(product(range(3),repeat=2)) for i,j in combinations(range(5),2)),'all binary projections are full, so binary closure has 243 not 27 states')
                spec=integer_spectrum(C,basis,3)
                require(spec==Counter({0:1,9:8,12:6,15:12}),'literal full-budget 27-state matrix spectrum')
        output.append({'L':L,'states':len(C),'resources':2*L+1,'selected_degree':degree,'full_budget_degree':len(full),'gap':min(x for x in spectrum if x),'spectrum':dict(sorted(spectrum.items())),'weight_enumerator':dict(sorted(enum.items()))})
    for L in range(1,257):
        ratio=Fraction(3+3*2**(L-1),2*L+2**(L+1))
        require(ratio>=Fraction(27,40) and (ratio==Fraction(27,40))==(L==4),'sharp normalized gap finite control')
        require(40*(3+3*2**(L-1))-27*(2*L+2**(L+1))==6*(2**L-9*L+20),'all-size gap comparison identity')
    return {'families':output,'normalized_bound':[27,40],'normalized_equality_L':4,'L2_pairwise_projection_closure_states':243}

def no_return_hostile():
    v=(1,0,0,1,1);w=(0,1,1,0,1);C=extend(extend(frozenset({(0,)*5}),v),w)
    target=minus(v,w);zero=(0,)*5
    require({x for x in C if 0<weight(x)<=3}=={v,w,times(2,v),times(2,w)},'budget3 only signed v,w')
    require(len(components(C,2))==9 and len(components(C,3))==1,'ordinary threshold exactly three')
    ordinary=[zero,v,target]
    require(all(distance(x,y)==3 for x,y in zip(ordinary,ordinary[1:])),'two-event ordinary route')
    require([x[4] for x in ordinary]==[0,1,0],'fixed endpoint coordinate returns to first full record')
    lower,_=no_return_distances(C,3);upper,_=no_return_distances(C,4)
    require(target not in lower and upper[target]==1 and weight(target)==4,'exact no-return threshold four')
    A=((0,1,2,0,0),(1,0,0,2,0),(2,2,0,0,1))
    require({x for x in product(range(3),repeat=5) if all(inner(a,x)==0 for a in A)}==set(C),'explicit guard matrix exact code')
    a=source_matrix(3)
    for x in C:
        reservoirs=[]
        for i,k in enumerate(x):
            p=record(3,i,k);b=partner(a,p)
            reservoirs.append({full_word(a,p),full_word(b,p),full_word(opposite(a),p[::-1]),full_word(opposite(b),p[::-1])})
        require(len(set().union(*reservoirs))==20,'all nine states occupy five distinct normalized units')
    return {'matrix':A,'basis':[v,w],'endpoint':target,'ordinary_route':ordinary,'no_return_route':[zero,target],'ordinary_budget':3,'no_return_budget':4}

def coordinate_checks():
    n=3;A=((1,2,0),(0,1,1));b=(2,1);ambient=tuple(product(range(3),repeat=n));count=0
    for p in permutations(range(n)):
        for slopes in product((1,2),repeat=n):
            for shift in ambient:
                images={z:tuple((slopes[i]*z[p[i]]+shift[i])%3 for i in range(n)) for z in ambient}
                transformedA=tuple(tuple(next(a[i]*slopes[i]%3 for i in range(n) if p[i]==j) for j in range(n)) for a in A)
                transformedb=tuple((t-inner(a,shift))%3 for a,t in zip(A,b))
                require({z for z,x in images.items() if tuple(inner(a,x) for a in A)==b}=={z for z in ambient if tuple(inner(a,z) for a in transformedA)==transformedb},'affine monomial guard formula')
                require(all(distance(images[x],images[y])==distance(x,y) for x in ambient for y in ambient),'all-pair affine monomial distance preservation')
                count+=1
    require(weight((1,0))!=weight((1,1)),'invertible shear need not preserve support')
    return {'all_pair_monomial_cases':count,'general_linear_support_preservation':'refuted by (x,y)->(x,x+y)'}

def raw_owner_coordinate_checks():
    """Canonical k=owner(q)-q is not the raw owner residue; check both charts."""
    cases=0;nonzero=0
    for r in range(1,6):
        q=8*r;ambient=tuple(product(range(3),repeat=2))
        raw={ks:tuple(record(r,i,k)[q] for i,k in enumerate(ks)) for ks in ambient}
        head={ks:tuple(record(r,i,k)[::-1][0]-q for i,k in enumerate(ks)) for ks in ambient}
        require(all(raw[ks]==tuple(q+k for k in ks) for ks in ambient),'literal raw owners equal q times one plus canonical phase')
        for a in product(range(3),repeat=2):
            offset=q*sum(a)%3
            for beta in range(3):
                canonical={ks for ks in ambient if inner(a,ks)==beta}
                literal={ks for ks in ambient if inner(a,tuple(o-q for o in raw[ks]))==beta}
                transported={ks for ks in ambient if inner(a,raw[ks])==(beta+offset)%3}
                unshifted={ks for ks in ambient if inner(a,raw[ks])==beta}
                require(canonical==literal==transported,'raw-owner equation requires beta plus q A1')
                require((canonical!=unshifted)==(offset!=0),'same-RHS raw-owner equation fails exactly when its nonzero translation matters')
                require({ks for ks in ambient if inner(a,head[ks])==(beta+2*sum(a))%3}==canonical,'converse head phase separately requires beta plus 2 A1')
                cases+=1;nonzero+=int(offset!=0)
    r=1;q=8;a=(1,0);beta=0
    good={ks for ks in product(range(3),repeat=2) if inner(a,ks)==beta}
    bad={ks for ks in product(range(3),repeat=2) if inner(a,tuple(record(r,i,k)[q] for i,k in enumerate(ks)))==beta}
    require({ks[0] for ks in good}=={0} and {ks[0] for ks in bad}=={1},'r1 sum-one hostile control rejects raw same-RHS conflation')
    require((16*sum(a)%3,2*sum(a)%3)==(1,2),'r2 raw-owner and converse-head shifts are distinct')
    return {'systems':cases,'nonzero_raw_translation_systems':nonzero,'raw_owner_formula':'A o = beta + (8 r) A1 modulo 3','converse_head_phase_formula':'A ell = beta + 2 A1 modulo 3','hostile':{'r':1,'q':8,'row':a,'canonical_rhs':beta,'raw_owner_rhs':2,'correct_first_phase':0,'unshifted_raw_first_phase':1},'distinct_shifts_at_r2':{'raw_owner':1,'converse_head_phase':2}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    data={'schema':'independent-native-affine-phase-audit-v1','native':native_checks(),'affine':affine_code_checks(),'rectangular':arbitrary_system_checks(),'modes':mode_checks(),'no_return':no_return_hostile(),'coordinates':coordinate_checks(),'raw_owner_coordinates':raw_owner_coordinate_checks()}
    data['explicit_checks']=CHECKS
    text=json.dumps(data,sort_keys=True,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')
    print('PASS independent native phase-code audit:',CHECKS,'checks; SHA256',sha256(text.encode()).hexdigest())

if __name__=='__main__':main()
