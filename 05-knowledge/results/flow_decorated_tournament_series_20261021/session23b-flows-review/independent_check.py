#!/usr/bin/env python3
"""Independent exact controls. Standard library only; no import of producer code."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib, json

CHECKS=0
def check(ok, msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def arcs(n,m):
    return tuple((v,u) if m>>k&1 else (u,v) for k,(u,v) in enumerate(combinations(range(n),2)))

def mask(n, aa):
    aa=set(aa)
    return sum((1<<k) for k,(u,v) in enumerate(combinations(range(n),2)) if (v,u) in aa)

def relabel(aa,p,anti=False):
    return tuple((p[v],p[u]) if anti else (p[u],p[v]) for u,v in aa)

def is_strong(n,aa):
    # Independent directed-cut criterion.
    return all(any(u in S and v not in S for u,v in aa) for k in range(1,n) for S in map(set,combinations(range(n),k)))

def cycles(n,aa):
    aset=set(aa); out=[]
    for k in range(3,n+1):
        for V in combinations(range(n),k):
            for tail in permutations(V[1:]):
                seq=(V[0],)+tail
                cycle=tuple(zip(seq,seq[1:]+seq[:1]))
                if all(e in aset for e in cycle): out.append(tuple(int(e in cycle) for e in aa))
    return out

def bounded_vectors(d,b):
    if d==0:
        yield ()
    else:
        for x in range(b+1):
            for tail in bounded_vectors(d-1,b-x): yield (x,)+tail

def tree_solve_flows(n,aa,b):
    # Solve the incidence equations on the star at vertex 0. Non-tree arc
    # weights are enumerated directly; no cycle decomposition or source cone map.
    free=[i for i,e in enumerate(aa) if 0 not in e]
    star=[i for i,e in enumerate(aa) if 0 in e]
    out=set()
    for vv in bounded_vectors(len(free),b):
        ff=[0]*len(aa); div=[0]*n
        for i,x in zip(free,vv):
            ff[i]=x;u,v=aa[i];div[u]+=x;div[v]-=x
        for i in star:
            u,v=aa[i]
            ff[i]=div[v] if u==0 else -div[u]
        if min(ff,default=0)<0 or sum(ff)>b:continue
        out.add(tuple(ff))
    return out

def edge_map(src,dst,p,anti=False):
    inds={e:i for i,e in enumerate(dst)}
    return tuple(inds[e] for e in relabel(src,p,anti))

def push(f,em):
    out=[0]*len(f)
    for i,j in enumerate(em):out[j]=f[i]
    return tuple(out)

# K4: direct parity kernel and orientation census, without producer routines.
E=tuple(combinations(range(4),2))
Z=[]
for bits in product((0,1),repeat=6):
    if all(sum(bits[i] for i,e in enumerate(E) if v in e)%2==0 for v in range(4)):
        Z.append(sum(x<<i for i,x in enumerate(bits)))
Z=set(Z); points=Z-{0}
faces={p for p in points if p.bit_count()==3};quads=points-faces
lines={frozenset((a,b,a^b)) for a,b in combinations(points,2)}
check(len(Z)==8 and len(faces)==4 and len(quads)==3 and len(lines)==7,'K4 binary kernel')
F4=list(permutations(range(4)))
induced=[]
for p in F4:
    mp={x:sum(1<<E.index(tuple(sorted((p[u],p[v])))) for i,(u,v) in enumerate(E) if x>>i&1) for x in Z}
    induced.append(tuple(mp[x] for x in sorted(Z)))
check(len(set(induced))==24,'faithful S4 cycle action')
# Abstract incidence automorphisms, not the basis enumeration used in producer.
incidence_aut=[]
pts=sorted(points)
for p in permutations(pts):
    mp=dict(zip(pts,p))
    if {frozenset(mp[x] for x in L) for L in lines}==lines:incidence_aut.append(mp)
check(len(incidence_aut)==168,'Fano incidence group')
check(sum({m[x] for x in faces}==faces for m in incidence_aut)==24,'partition stabilizer')
profiles=Counter();fibers=defaultdict(list)
for m in range(64):
    aa=arcs(4,m); cc=cycles(4,aa)
    pp=tuple(sorted(sum(x<<i for i,x in enumerate(c)) for c in cc))
    prof=tuple(sorted(sum(c) for c in cc));profiles[prof]+=1
    if is_strong(4,aa):
        check(prof==(3,3,4),'all strong K4 cycle lengths')
        check(pp[0]^pp[1]^pp[2]!=0,'noncollinear strong supports')
        fibers[pp].append(m)
check(profiles==Counter({():24,(3,):16,(3,3,4):24}),'K4 tournament profile')
check(len(fibers)==12 and all(len(v)==2 and v[0]^v[1]==63 for v in fibers.values()),'only converse support ambiguity')

T4=((0,1),(1,2),(2,0),(1,3),(3,0),(2,3))
flows4=tree_solve_flows(4,T4,24)
coords4={tuple((a+b+c,a+c,a,b,b+c,c)) for a,b,c in product(range(9),repeat=3) if 3*a+3*b+4*c<=24}
check(flows4==coords4,'T4 star-kernel equals free unimodular cone')
# The integer generator matrix has a literal 3x3 identity coordinate minor.
r4=cycles(4,T4)
check(len(r4)==3,'three T4 primitive rays')
identity_rows=[tuple(r[i] for r in r4) for i in (2,3,5)]
check(sorted(identity_rows)==[(0,0,1),(0,1,0),(1,0,0)],'T4 lattice-unit coordinate minor')
a4=edge_map(T4,T4,(1,0,3,2),True)
check(all(push(f,a4)==(f[0],f[4],f[3],f[2],f[1],f[5]) for f in flows4),'T4 anti action')
raw4=Counter(map(sum,flows4));fix4=Counter(sum(f) for f in flows4 if push(f,a4)==f)
check(fix4==Counter(6*a+4*c for a,c in product(range(7),repeat=2) if 6*a+4*c<=24),'T4 anti fixed series')

# R5 direct incidence-kernel enumeration to weight 24; check parametrization,
# twelve rays, the actual anti-fixed sum, ordinary/anti Burnside and orbit types.
R5=tuple((i,(i+d)%5) for i in range(5) for d in (1,2))
RBOUND=24
rr=cycles(5,R5)
check(Counter(map(sum,rr))==Counter({3:5,4:5,5:2}),'R5 cycle ray count')
rcoords=[]
for r in rr:
    y=r[1::2];tt={r[2*i]+y[i]+y[(i-1)%5] for i in range(5)}
    check(len(tt)==1,'ray pentagon parameter')
    rcoords.append((next(iter(tt)),y))
expected={(1,(0,)*5),(2,(1,)*5)}
for i in range(5):
    expected.add((1,tuple(int(j==i) for j in range(5))))
    expected.add((1,tuple(int(j in {i,(i+2)%5}) for j in range(5))))
check(set(rcoords)==expected,'all twelve pentagon extreme rays')
rf=tree_solve_flows(5,R5,RBOUND)
for f in rf:
    y=f[1::2];tt={f[2*i]+y[i]+y[(i-1)%5] for i in range(5)}
    check(len(tt)==1 and sum(f)==5*next(iter(tt))-sum(y),'all R5 flows pentagon cone')
g5=[p for p in permutations(range(5)) if set(relabel(R5,p))==set(R5)]
a5=[p for p in permutations(range(5)) if set(relabel(R5,p,True))==set(R5)]
check(len(g5)==len(a5)==5,'R5 symmetry group sizes')
gm=[edge_map(R5,R5,p) for p in g5]; am=[edge_map(R5,R5,p,True) for p in a5]
can=lambda f:min(push(f,g) for g in gm)
reps=set(map(can,rf)); fixedrep={f for f in reps if can(push(f,am[0]))==f}
raw5=Counter(map(sum,rf));orb5=Counter(map(sum,reps));fix5=Counter(map(sum,fixedrep))
ordinary_burn=Counter();anti_burn=Counter()
for f in rf:
    ordinary_burn[sum(f)]+=sum(push(f,g)==f for g in gm)
    anti_burn[sum(f)]+=sum(push(f,a)==f for a in am)
check(all(ordinary_burn[w]==5*orb5[w] and anti_burn[w]==5*fix5[w] for w in range(RBOUND+1)),'R5 coset Burnside')
check(all(5*orb5[w]==raw5[w]+(4*(w//5+1) if w%5==0 else 0) for w in range(RBOUND+1)),'R5 ordinary formula 17')
# Bounds: under t>=max(...), all y<=t and 2sum(y)<=5t;
# hence weight>=5t/2, so t<=floor(2RBOUND/5).
paramfixed=Counter()
for t in range(2*RBOUND//5+1):
    for a,b,c in product(range(t+1),repeat=3):
        if t>=max(a+b,2*b,a+c):
            w=5*t-2*a-2*b-c
            if w<=RBOUND:paramfixed[w]+=1
check(paramfixed==fix5,'R5 explicit anti-fixed parameter sum 18')
neg=edge_map(R5,R5,tuple((-i)%5 for i in range(5)),True)
check(Counter(sum(f) for f in rf if push(f,neg)==f)==paramfixed,'R5 chosen reflection fixed points equal 18')
# A cone point need not have unique primitive-cycle coordinates in R5.
sums={};relation=None
for i,r in enumerate(rr):
    for j,s in enumerate(rr[i:],i):
        v=tuple(a+b for a,b in zip(r,s))
        if v in sums and sums[v]!=(i,j):relation=(sums[v],(i,j));break
        sums[v]=(i,j)
    if relation:break
check(relation is not None,'R5 primitive cycle decomposition nonunique hostile')

# All strong types through order 5, including genuinely non-self-converse atoms.
# Rebuild isomorphism types from labelled adjacency, and derive h,b directly from
# Hamilton words and positional reversals. No inherited atom table is imported.
types=[]
for n in (1,3,4,5):
    pp=list(permutations(range(n)));todo=set(range(1<<(n*(n-1)//2)))
    while todo:
        m=min(todo);aa=arcs(n,m)
        orbit={mask(n,relabel(aa,p)) for p in pp};todo-=orbit
        if not is_strong(n,aa):continue
        G=[p for p in pp if set(relabel(aa,p))==set(aa)]
        A=[p for p in pp if set(relabel(aa,p,True))==set(aa)]
        hp=[p for p in pp if all((u,v) in aa for u,v in zip(p,p[1:]))]
        blues=0
        for p in hp:
            reverse=[0]*n
            for u,v in zip(p,p[::-1]):reverse[u]=v
            blues+=tuple(reverse) in A
        check(len(hp)%len(G)==0 and blues%len(G)==0,'Hamilton orbit division')
        types.append(dict(n=n,m=m,aa=aa,G=G,A=A,h=len(hp)//len(G),b=blues//len(G)))
# Every strong type through order five happens to be self-converse. Add an
# actual order-six converse pair, detected by a non-palindromic score multiset.
n=6; pp=list(permutations(range(n)))
for m in range(1<<15):
    aa=arcs(n,m)
    deg=sorted(sum(u==v for u,w in aa) for v in range(n))
    if deg!=sorted(n-1-d for d in deg) and is_strong(n,aa):break
else:raise RuntimeError('missing order-six non-self-converse witness')
for source in (aa,tuple((v,u) for u,v in aa)):
    m=min(mask(n,relabel(source,p)) for p in pp);aa6=arcs(n,m)
    G=[p for p in pp if set(relabel(aa6,p))==set(aa6)]
    AA=[p for p in pp if set(relabel(aa6,p,True))==set(aa6)]
    hp=[p for p in pp if all((u,v) in aa6 for u,v in zip(p,p[1:]))]
    check(not AA and len(hp)%len(G)==0,'genuine non-self-converse order-six atom')
    types.append(dict(n=n,m=m,aa=aa6,G=G,A=AA,h=len(hp)//len(G),b=0))
N,Q=12,12
for i,T in enumerate(types):
    target=next((j,p) for j,S in enumerate(types) if S['n']==T['n'] for p in permutations(range(T['n'])) if set(relabel(T['aa'],p,True))==set(S['aa']))
    T['op'],T['anti']=target
    T['gm']=[edge_map(T['aa'],T['aa'],p) for p in T['G']]
    fs=tree_solve_flows(T['n'],T['aa'],Q)
    T['flows']=set(min(push(f,g) for g in T['gm']) for f in fs)
check(any(i!=T['op'] for i,T in enumerate(types)),'non-self-converse strong atoms included')
letters=[];opp={};atom=[];fixedatom=[]
for i,T in enumerate(types):
    a=Counter();k=Counter();S=types[T['op']]
    anti=edge_map(T['aa'],S['aa'],T['anti'],True)
    for f in sorted(T['flows']):
        l=(i,f);o=(T['op'],min(push(push(f,anti),g) for g in S['gm']))
        opp[l]=o;letters.append((l,T['n'],sum(f)));a[T['n'],sum(f)]+=1
        if l==o:k[T['n'],sum(f)]+=1
    atom.append(a);fixedatom.append(k)
check(all(opp[opp[l]]==l for l in opp),'converse well-defined involution on all atom flow orbits')

# Explicit words and literal converse orbits, grouped by underlying h, degree.
ordinary=defaultdict(Counter);actualfixed=defaultdict(Counter);undersc=defaultdict(Counter)
bydegree=defaultdict(set);ordinary_all=Counter();sccword_count=0

def process(word,n,w):
    global sccword_count
    sccword_count+=1
    ts=tuple(l[0] for l in word); rev=tuple(opp[l] for l in word[::-1])
    h=1
    for i in ts:h*=types[i]['h']
    sc=ts==tuple(types[i]['op'] for i in ts[::-1])
    ordinary[h][n,w]+=1;ordinary_all[n,w]+=1
    if sc:undersc[h][n,w]+=1
    if word==rev:actualfixed[h][n,w]+=1
    if sc:
        x=1
        for i in ts[:len(ts)//2]:x*=types[i]['h']
        b=types[ts[len(ts)//2]]['b'] if len(ts)%2 else 1
        check((h-x*b)%2==0,'underlying self-converse even deficit')
        degree=(h-x*b)//2
    else:degree=h
    bydegree[degree,n,w].add(min(word,rev))
    for l,d,u in letters:
        if n+d<=N and w+u<=Q:process(word+(l,),n+d,w+u)
process((),0,0)

ONE=Counter({(0,0):1})
def plus(*ss):
    r=Counter()
    for s in ss:r.update(s)
    return +r

def times(a,b):
    r=Counter()
    for (n,w),v in a.items():
        for (m,u),x in b.items():
            if n+m<=N and w+u<=Q:r[n+m,w+u]+=v*x
    return +r

def renewal(a):
    out=ONE.copy();power=ONE
    for _ in range(N):power=times(power,a);out.update(power)
    return +out

def square_sub(a):return Counter({(2*n,2*w):v for (n,w),v in a.items() if 2*n<=N and 2*w<=Q})

def pair_same_type(a):
    r=Counter()
    for (n,w),v in a.items():
        for (m,u),x in a.items():
            check(n==m,'one underlying type per independent square')
            if n+m<=N and w+u<=Q:r[n+m,w+u]+=v*x
    return +r

neutral=plus(*(a for a,T in zip(atom,types) if T['h']==1))
A=renewal(neutral)
pairatom=[pair_same_type(a) for a in atom]
AI=renewal(plus(*(a for a,T in zip(pairatom,types) if T['h']==1)))
P=defaultdict(Counter);PI=defaultdict(Counter)
for a,p,T in zip(atom,pairatom,types):
    if T['h']>1:P[T['h']].update(a);PI[T['h']].update(p)

# Equation 9 and equation 11 via finite multiplicative factorization.
from functools import lru_cache
@lru_cache(None)
def factor_tuples(h):
    if h==1:return ((),)
    return tuple((j,)+tail for j in sorted(P) if h%j==0 for tail in factor_tuples(h//j))

def tuple_series(h,carrier,polys):
    r=Counter()
    for seq in factor_tuples(h):
        term=carrier
        for j in seq:term=times(times(term,polys[j]),carrier)
        r.update(term)
    return +r
F={h:tuple_series(h,A,P) for h in ordinary}
G={h:tuple_series(h,AI,PI) for h in ordinary}
check(all(F[h]==ordinary[h] for h in ordinary),'all fixed underlying h ordinary finite tuple formulas')
J=defaultdict(Counter);U=defaultdict(Counter);L=defaultdict(Counter);M=defaultdict(Counter)
centers=[(1,1,ONE,ONE)]+[(T['h'],T['b'],atom[i],fixedatom[i]) for i,T in enumerate(types) if T['op']==i]
for x in ordinary:
    for a,b,ca,ck in centers:
        h=x*x*a;d=x*(x*a-b)
        J[h].update(times(square_sub(F[x]),ck));U[h].update(times(G[x],ca))
        L[d].update(times(G[x],ca));M[d].update(times(square_sub(F[x]),ck))
check(all(+J[h]==actualfixed[h] and +U[h]==undersc[h] for h in ordinary),'fixed and independently decorated center formulas 10 and 12')
check(all(not s or h in ordinary for h,s in J.items()),'no missing h slice in self-converse formula')
ks={k for k,n,w in bydegree}
for k in ks:
    numerator=plus(L[2*k],M[2*k])
    if k%2:
        numerator.update(F.get(k,{}));numerator.subtract(U[k]);numerator=+numerator
    expected=Counter({(n,w):len(s) for (d,n,w),s in bydegree.items() if d==k})
    check(all(v%2==0 for v in numerator.values()),'degree Burnside integrality')
    check(Counter({p:v//2 for p,v in numerator.items()})==expected,'correct fixed underlying degree formula 14')
# Degree-zero exact alphabet and independent versus mirror neutral carrier.
zero_centers=[(ca,ck) for a,b,ca,ck in centers if a==b]
check({T['n'] for T in types if T['op']==types.index(T) and T['h']==T['b']}=={1,3,5},'pure-blue strong atoms through5')
B0numerator=plus(times(AI,plus(*(ca for ca,ck in zero_centers))),times(square_sub(A),plus(*(ck for ca,ck in zero_centers))))
check(B0numerator==plus(L[0],M[0]),'degree-zero formula 15 including R5')
neutralall=times(AI,plus(ONE,neutral));neutralfixed=times(square_sub(A),plus(ONE,neutral))
check(neutralall[6,3]==2 and neutralfixed[6,3]==0,'independent vs mirror neutral hostile')
# q=0 reduction: every underlying ordinary atom has exactly one zero flow,
# and fixed zero flow occurs precisely for underlying self-converse atoms.
check(all(atom[i][T['n'],0]==1 and fixedatom[i][T['n'],0]==(T['op']==i) for i,T in enumerate(types)),'zero-flow forgetful specialization')

base=Path(__file__).parent.parent/'session23b-flows'
report={
 'status':'PASS','checks':CHECKS,
 'source_sha256':{p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in ('PROOF.md','check.py','normal.json','optimized.json')},
 'k4':{'orientation_profiles':{str(k):v for k,v in profiles.items()},'fano_incidence_automorphisms':len(incidence_aut),'T4_total_use_bound':24,'T4_flows':len(flows4)},
 'R5':{'direct_incidence_enumeration_bound':RBOUND,'raw_flows':len(rf),'flow_orbits':len(reps),'self_converse_orbits':len(fixedrep),'raw_coefficients':[raw5[w] for w in range(RBOUND+1)],'ordinary_coefficients':[orb5[w] for w in range(RBOUND+1)],'self_converse_coefficients':[fix5[w] for w in range(RBOUND+1)],'nonunique_primitive_cycle_pair_indices':relation},
 'all_strong_atoms_through5':[{'n':T['n'],'mask':T['m'],'h':T['h'],'b':T['b'],'Aut_order':len(T['G']),'converse_type':T['op'],'self_converse':i==T['op']} for i,T in enumerate(types)],
 'word_controls':{'vertex_bound':N,'use_bound':Q,'decorated_words':sccword_count,'ordinary_h_slices':sorted(ordinary),'underlying_degree_slices':sorted(ks),'neutral_hostile_all':neutralall[6,3],'neutral_hostile_fixed':neutralfixed[6,3]},
 'scope':'Exact finite controls and proof audit, not a census above the stated bounds. Formulas are conditional on the named inherited all-order finite-alphabet and center theorems.'}
print(json.dumps(report,sort_keys=True,indent=2))
