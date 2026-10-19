#!/usr/bin/env python3
"""Bounded exact audit of distinct E-families and their typed connections.
Python stdlib + SymPy; no network, random sampling, or Python assertions.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
from math import factorial, comb, isqrt
import json
import sympy as sp

checks = 0
def require(condition, message):
    global checks
    checks += 1
    if not condition: raise RuntimeError(message)

def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    if not a: return 0
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None: continue
        a[r], a[p] = a[p], a[r]
        v = a[r][c]; a[r] = [x/v for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                v = a[i][c]; a[i] = [x-v*y for x,y in zip(a[i],a[r])]
        r += 1
    return r

def adjacency(n, es):
    a = [set() for _ in range(n)]
    for u,v in es: a[u].add(v); a[v].add(u)
    return a

def connected_bipartite(n, es):
    a=adjacency(n,es); col={0:0}; stack=[0]; bip=True
    while stack:
        u=stack.pop()
        for v in a[u]:
            if v in col:
                if col[v]==col[u]: bip=False
            else: col[v]=1-col[u]; stack.append(v)
    return len(col)==n,bip

def signless_rank(n, es):
    return rank([[int(v==u or v==w) for v in range(n)] for u,w in es])

def endpoint_count(n):
    arc=lambda a,b: ((a<b) and (a,b)!=(0,n-1)) or (a,b)==(n-1,0)
    paths=[p for p in permutations(range(n)) if all(arc(a,b) for a,b in zip(p,p[1:]))]
    expected={tuple(range(n))}
    for mask in range(1<<(n-2)):
        s=[i+1 for i in range(n-2) if mask>>i&1]
        t=[i+1 for i in range(n-2) if not mask>>i&1]
        expected.add(tuple(s+[n-1,0]+t))
    require(set(paths)==expected, f'endpoint address n{n}')
    return len(paths)

# Curve classes a H + sum b_i E_i, H^2=1, E_i^2=-1.
def dot(v,w): return v[0]*w[0]-sum(a*b for a,b in zip(v[1:],w[1:]))
def curves(r):
    out=[]
    for i in range(r):
        v=[0]*(r+1); v[i+1]=1; out.append((f'E{i}',tuple(v)))
    for i,j in combinations(range(r),2):
        v=[0]*(r+1); v[0]=1; v[i+1]=v[j+1]=-1; out.append((f'L{i}{j}',tuple(v)))
    if r==5: out.append(('Q',tuple([2]+[-1]*r)))
    return out

def curve_graph(r):
    cs=curves(r); K=tuple([-3]+[1]*r)
    for name,v in cs: require(dot(v,v)==dot(v,K)==-1,f'curve {r} {name}')
    es=[]
    for i,j in combinations(range(len(cs)),2):
        z=dot(cs[i][1],cs[j][1]); require(z in (0,1),f'intersection {r} {i} {j}')
        if z: es.append((i,j))
    return cs,es

out={}
out['endpoint_reversal_paths']={n:endpoint_count(n) for n in range(3,9)}
require(out['endpoint_reversal_paths'][3]==3 and out['endpoint_reversal_paths'][4]==5 and out['endpoint_reversal_paths'][5]==9,'E345 counts')

# Explicitly different scalar observers on the same set.
S=tuple(range(1,6)); rep=Counter(a+b for a in S for b in S)
energy={k:sum(v**k for v in rep.values()) for k in range(2,6)}
schur=sum(a+b in S for a in S for b in S)
esym={k:sum(sp.prod(t) for t in combinations(S,k)) for k in range(3,6)}
require(energy=={2:85,3:325,4:1333,5:5725},'energy values')
require(schur==10 and esym=={3:225,4:274,5:120},'scalar taxonomy')
out['same_set_distinct_observers']={'set':S,'Schur_E3':schur,'energy':energy,'elementary_symmetric':{k:int(v) for k,v in esym.items()}}

# Exact rational scattering coefficients after Euler's even-zeta formula.
scattering={}
for r,expected in [(3,F(2835,8)),(4,F(23625,8)),(5,F(3274425,128))]:
    gamma_factor=F(factorial(2*r-2),4**(r-1)*factorial(r-1)**2) # gamma ratio / pi
    even_coefficient=F(abs(sp.bernoulli(2*r))*2**(2*r-1),factorial(2*r))
    actual=gamma_factor/even_coefficient
    require(actual==expected,f'Eisenstein coefficient r{r}')
    scattering[r]={'zeta_weight':2*r-1,'coefficient_times_pi_power':str(actual),'pi_denominator_power':2*r-1}
out['real_analytic_Eisenstein_scattering']=scattering
x,t=sp.symbols('x t')
out['Euler_polynomials_3_4_5']={n:str(sp.euler(n,x)) for n in (3,4,5)}
for n in (3,4,5):
    require(sp.expand(sp.euler(n,x+1)+sp.euler(n,x)-2*x**n)==0,f'Euler translation n{n}')
    require(sp.expand(sp.euler(n,1-x)-(-1)**n*sp.euler(n,x))==0,f'Euler reflection n{n}')
out['Euler_numbers_sech_3_4_5']={n:int(sp.euler(n)) for n in (3,4,5)}
require(out['Euler_numbers_sech_3_4_5']=={3:0,4:5,5:0},'Euler parity')

# Extended E3,E4,E5 del Pezzo root convention, with the graph data kept distinct.
geo={}
for r in (3,4,5):
    cs,es=curve_graph(r); a=adjacency(len(cs),es); K=tuple([-3]+[1]*r)
    # Cauchy gives (9-r)*a^2 <= 2r for roots: |a|<=1, |b_i|<=1.
    # Thus this finite scan is complete, not a conjectural height cutoff.
    roots=[v for v in product(range(-2,3),repeat=r+1) if dot(v,K)==0 and dot(v,v)==-2]
    # For (-1)-classes, Cauchy gives (1-3a)^2 <= r*(a^2+1).
    # For r=3,4 this forces a in {0,1}; for r=5 it allows {0,1,2}.
    all_curves=[v for v in product(range(-2,3),repeat=r+1) if dot(v,K)==dot(v,v)==-1]
    require(set(all_curves)=={v for _,v in cs},f'complete minus-one class census r{r}')
    require(len(roots)=={3:8,4:20,5:40}[r],f'root finite census r{r}')
    require(len(cs)=={3:6,4:10,5:16}[r],f'curve count r{r}')
    require(set(map(len,a))=={{3:2,4:3,5:5}[r]},f'curve degree r{r}')
    require(connected_bipartite(len(cs),es)[0],f'curve connected r{r}')
    geo[r]={'root_count':len(roots),'curve_vertices':len(cs),'curve_edges':len(es),'degree':len(a[0])}
    if r>3:
        removed=r-1 # marked exceptional E_(r-1)
        keep=[i for i in range(len(cs)) if i!=removed and i not in a[removed]]
        previous=curves(r-1)
        projected=[(cs[i][0],cs[i][1][:-1]) for i in keep]
        require(set(projected)==set(previous),f'disjoint-curve deletion r{r}')
        require(all(dot(cs[i][1],cs[j][1])==dot(cs[i][1][:-1],cs[j][1][:-1]) for i,j in combinations(keep,2)),f'intersection transport r{r}')
    if r==4:
        pairs=[]
        for name,v in cs:
            if name.startswith('E'): pairs.append(frozenset((int(name[1:]),4)))
            else: pairs.append(frozenset(set(range(4))-set(map(int,name[1:]))))
        require(len(set(pairs))==10,'Petersen pair bijection')
        require(all((not (pairs[i]&pairs[j]))==(j in a[i]) for i,j in combinations(range(10),2)),'Petersen incidence')
        def pair_lattice(v):
            return tuple(v[0]+v[i+1] for i in range(4))+(2*v[0]+sum(v[1:]),)
        def h(u,v): return F(sum(u)*sum(v),4)-sum(a*b for a,b in zip(u,v))
        basis=[tuple(int(i==j) for i in range(5)) for j in range(5)]
        require(abs(sp.Matrix.hstack(*(sp.Matrix(pair_lattice(v)) for v in basis)).det())==2,'pair lattice index two')
        require(all(h(pair_lattice(u),pair_lattice(v))==dot(u,v) for u in basis for v in basis),'Picard pair lattice isometry')
        require(pair_lattice(K)==(-2,)*5,'canonical class in pair lattice')
        A4={tuple(int(k==i)-int(k==j) for k in range(5)) for i in range(5) for j in range(5) if i!=j}
        require({pair_lattice(v) for v in roots}==A4,'E4 roots equal A4 differences')
        require(all(pair_lattice(v)==tuple(int(i in pair) for i in range(5)) for (_,v),pair in zip(cs,pairs)),'curve and pair maps commute')
out['del_pezzo_E3_E4_E5']=geo

# Every connected simple labelled graph through n=5: exact signless rank.
census=0
for n in range(2,6):
    all_es=list(combinations(range(n),2))
    for mask in range(1<<len(all_es)):
        es=[e for i,e in enumerate(all_es) if mask>>i&1]
        conn,bip=connected_bipartite(n,es)
        if not conn: continue
        census+=1
        require(signless_rank(n,es)==n-int(bip),f'signless rank n{n} mask{mask}')
out['connected_graph_signless_rank_census']=census
require(census==771,'connected graph universe')
petverts=list(combinations(range(5),2)); pet=[(i,j) for i,j in combinations(range(10),2) if not set(petverts[i])&set(petverts[j])]
K5=list(combinations(range(5),2)); K33=[(i,j) for i in range(3) for j in range(3,6)]
Q=lambda n:[(a-1,b-1) for a,b in combinations(range(1,n+1),2) if isqrt(a+b)**2==a+b]
targets=[('Petersen',10,pet),('K5',5,K5),('K33',6,K33),('Q14',14,Q(14)),('Q15',15,Q(15))]
out['named_signless_ranks']={name:signless_rank(n,es) for name,n,es in targets}
require(out['named_signless_ranks']=={'Petersen':10,'K5':5,'K33':5,'Q14':13,'Q15':15},'named rank values')

# Two quadratic-cycle controls; derivative and secant products are different observers.
c=F(-29,16); cyc=[F(-7,4),F(5,4),F(-1,4)]
require(all(cyc[(i+1)%3]==cyc[i]**2+c for i in range(3)),'rational orbit')
sums=[cyc[i]+cyc[(i+1)%3] for i in range(3)]
require(sp.prod(sums)==1 and sp.prod(2*a for a in cyc)==F(35,8),'distinct multipliers')
require(all((sums[i]-sums[(i+1)%3]+sums[(i+2)%3])/2==cyc[i] for i in range(3)),'odd-cycle decoder')
P=x**3+x**2/2-sp.Rational(9,4)*x-sp.Rational(1,8)
f=x*x-sp.Rational(7,4); f2=sp.expand(f.subs(x,f)); f3=sp.expand(f.subs(x,f2))
require(sp.expand(f3-x-(f-x)*P**2)==0,'parabolic factorization')
require(sp.rem((x+f)*(f+f2)*(f2+x)-1,P,x)==0,'parabolic secant product')
require(sp.rem(8*x*f*f2-1,P,x)==0,'parabolic derivative product')
out['quadratic_cycles']={'rational_sums':list(map(str,sums)),'rational_secant_product':'1','rational_multiplier':'35/8','parabolic_polynomial':str(P),'parabolic_secant_product':'1','parabolic_multiplier':'1'}

# These test the rational prefactors of cited period identities, not integrals.
periods={}
for n in (4,5,6):
    wheel=comb(2*n-2,n-1)
    zigzag=F(4,n)*wheel*(1-F(1-(-1)**n,2**(2*n-3)))
    periods[n]={'zeta_weight':2*n-3,'wheel_coefficient':wheel,'zigzag_coefficient':str(zigzag)}
require([periods[n]['wheel_coefficient'] for n in (4,5,6)]==[20,70,252],'wheel coefficients')
require([periods[n]['zigzag_coefficient'] for n in (4,5,6)]==['20','441/8','168'],'zigzag coefficients')
out['cited_graph_period_prefactors']=periods

# Eulerian metagraphs: graph isomorphism-class adjacency retains a basis.
def even_meta(n, basis_edges):
    edges=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(edges)}
    ps=list(permutations(range(n)))
    def canon(mask):
        return min(sum(1<<idx[tuple(sorted((p[u],p[v])))] for k,(u,v) in enumerate(edges) if mask>>k&1) for p in ps)
    def path(u,v):
        a=adjacency(n,basis_edges); queue=[(u,[])]; seen={u}
        for w,p in queue:
            if w==v:return p
            for z in a[w]:
                if z not in seen: seen.add(z);queue.append((z,p+[tuple(sorted((w,z)))]))
        raise RuntimeError('tree disconnected')
    basis=[]
    for e in edges:
        if e not in basis_edges: basis.append(sum(1<<idx[f] for f in path(*e)+[e]))
    states=[]
    for b in range(1<<len(basis)):
        g=0
        for i,v in enumerate(basis):
            if b>>i&1:g^=v
        require(all(sum(bool(g>>idx[e]&1) for e in edges if u in e)%2==0 for u in range(n)),'even parity')
        states.append(g)
    classes={g:canon(g) for g in states}; meta=set()
    for g in states:
        for b in basis:
            u,v=classes[g],classes[g^b]
            if u!=v:meta.add(tuple(sorted((u,v))))
    return len(set(classes.values())),len(meta)
out['Eulerian_metagraphs']={}
for n in (3,4,5):
    star=even_meta(n,[(0,i) for i in range(1,n)])
    path=even_meta(n,[(i,i+1) for i in range(n-1)])
    require(star[0]==path[0]=={3:2,4:3,5:7}[n],f'even classes n{n}')
    out['Eulerian_metagraphs'][n]={'star':star,'path_envelope':path}
require(out['Eulerian_metagraphs'][4]=={'star':(3,2),'path_envelope':(3,3)},'MISTAKE495 hostile')
out['checks_passed']=checks
out['status']='PASS: exact finite controls; no new irrationality, graceful existence, or global period claim'
print(json.dumps(out,indent=2,sort_keys=True))
