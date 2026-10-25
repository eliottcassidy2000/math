#!/usr/bin/env python3
"""Exact controls for root slices, genus-two quotients and the Fano/pants codec.
Only Python3 + SymPy. No network, repository, hidden state, or assert statements.
"""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from math import comb
import json
from pathlib import Path
import sympy as sp

CHECKS = 0
def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)

def bitpair(a,b): return (a & b).bit_count() % 2
def W(a):
    require(a.bit_count()%2==0, ('odd mask',a))
    return a if a.bit_count()<=2 else a^63

def pmask(pair): return sum(1<<(i-1) for i in pair)
def matchings(S):
    if not S:
        yield ()
    else:
        a=S[0]
        for b in S[1:]:
            for M in matchings(tuple(x for x in S if x not in (a,b))):
                yield ((a,b),)+M

def permute(a,p): return W(sum(1<<p[i] for i in range(6) if a>>i&1))

# Character dimensions, independently compared with complete-intersection coefficients.
rows=[]
for n in range(4,10):
    even=[m for m in range(1,1<<n) if m.bit_count()%2==0]
    curve_dims=[max(0,m.bit_count()//2-1) for m in even]
    surf_dims=[comb(m.bit_count()//2-1,2) if m.bit_count()>=6 else 0 for m in even]
    g=1+(1<<(n-3))*(n-4)
    pg=sum((-1)**j*comb(n-3,j)*comb(2*n-7-2*j,n-1)
           for j in range((n-6)//2+1)) if n>=6 else 0
    require(sum(curve_dims)==g,('curve dimensions',n))
    require(sum(surf_dims)==pg,('surface dimensions',n))
    rows.append({'n':n,'curve_genus':g,'surface_pg':pg,
                 'positive_curve_factor_dimensions':dict(sorted(Counter(d for d in curve_dims if d).items())),
                 'size_six_characters':sum(m.bit_count()==6 for m in even)})

# The graph of lifted pants pieces: n−2 base pieces, each lift a sphere with six holes.
lifted=[]
for n in (6,7):
    d=1<<(n-1); inertia=[1<<i for i in range(n-1)]+[d-1]
    def label(S):
        a=0
        for i in S:a ^= inertia[i]
        return a
    partitions=[({0},{1},set(range(2,n)))]
    partitions += [(set(range(k+1)),{k+1},set(range(k+2,n))) for k in range(1,n-3)]
    partitions += [(set(range(n-2)),{n-2},{n-1})]
    require(len(partitions)==n-2,('base pants',n))
    subgroups=[]
    for parts in partitions:
        a,b,c=map(label,parts)
        require(a and b and c and a^b^c==0 and len({a,b,c})==3,('rank2 cuff labels',n))
        subgroups.append((0,a,b,c))
    coset=lambda x,H:min(x^h for h in H)
    verts={(v,coset(x,H)) for v,H in enumerate(subgroups) for x in range(d)}
    adjacency=defaultdict(set); edges=[]
    for r in range(2,n-1):
        g=label(range(r)); H0,H1=subgroups[r-2:r]
        for x in range(d):
            if x> x^g:continue
            a=(r-2,coset(x,H0));b=(r-1,coset(x,H1))
            require(a in verts and b in verts,('edge endpoint',n))
            edges.append((a,b));adjacency[a].add(b);adjacency[b].add(a)
    seen=set();stack=[next(iter(verts))]
    while stack:
        x=stack.pop()
        if x in seen:continue
        seen.add(x);stack.extend(adjacency[x]-seen)
    require(seen==verts,('lift connected',n))
    require(len(verts)==(n-2)*(1<<(n-3)),('lift vertex formula',n))
    require(len(edges)==(n-3)*(1<<(n-2)),('lift edge formula',n))
    require(len(edges)-len(verts)+1==1+(1<<(n-3))*(n-4),('pants genus',n))
    lifted.append({'n':n,'vertices':len(verts),'edges':len(edges),'cycle_rank':len(edges)-len(verts)+1})

# Genus-two 2-torsion: all 16 classes and the full S6 action.
V=[0]+[pmask(p) for p in combinations(range(1,7),2)]
require(len(set(W(a) for a in range(64) if a.bit_count()%2==0))==16,'W dimension4')
for a,b,c in product(V,repeat=3):
    require(bitpair(W(a^b),c)==(bitpair(a,c)^bitpair(b,c)),'pairing bilinear')
for a in V:
    require(bitpair(a,a)==0,'pairing alternating')
    require(a==0 or any(bitpair(a,b) for b in V),'pairing nondegenerate')
actions=set()
for p in permutations(range(6)):
    A=tuple(permute(a,p) for a in V)
    actions.add(A)
    require(all(bitpair(a,b)==bitpair(permute(a,p),permute(b,p)) for a,b in product(V,repeat=2)),'S6 symplectic')
require(len(actions)==720,'faithful S6')
# Count symplectic bases, independently proving |Sp4(F2)|=720.
symp_bases=[]
for a,b in product(V[1:],repeat=2):
    if not bitpair(a,b):continue
    perp=[v for v in V[1:] if not bitpair(v,a) and not bitpair(v,b)]
    for c,d0 in product(perp,repeat=2):
        if bitpair(c,d0):symp_bases.append((a,b,c,d0))
require(len(symp_bases)==720,'Sp4 order')

# Every branch pairing is a maximal isotropic subgroup and a theta-pants codec.
Ms=list(matchings(tuple(range(1,7))))
require(len(Ms)==15,'15 matchings')
Lspaces=set();codec_rows=[]
for M in Ms:
    a1,a2,a3=map(pmask,M)
    require(W(a1^a2)==a3 and not bitpair(a1,a2),'matching Lagrangian')
    Lspaces.add(tuple(sorted((0,a1,a2,a3))))
    b1,b2=next((b1,b2) for b1,b2 in product(V,repeat=2)
               if (bitpair(a1,b1),bitpair(a2,b1),bitpair(a1,b2),bitpair(a2,b2),bitpair(b1,b2))==(1,0,0,1,0))
    fibers=defaultdict(list)
    for v in V:
        eps=tuple(bitpair(v,a) for a in (a1,a2,a3))
        tau=(bitpair(v,b1),bitpair(v,b2))
        require(eps[0]^eps[1]^eps[2]==0,'pants compatibility')
        dec=0
        for bit,a in zip((tau[0],eps[0],tau[1],eps[1]),(a1,b1,a2,b2)):
            if bit:dec^=a
        require(W(dec)==v,'pants exact decoder')
        fibers[eps].append(v)
    require(sorted(map(len,fibers.values()))==[4]*4,'four seam choices')
    require(set(fibers[(0,0,0)])=={0,a1,a2,a3},'seam kernel equals matching subgroup')
    codec_rows.append({'matching':M,'dual_basis_masks':[b1,b2]})
require(len(Lspaces)==15,'different matching kernels')
all_L=set()
for a,b in combinations(V[1:],2):
    if not bitpair(a,b):all_L.add(tuple(sorted((0,a,b,W(a^b)))))
require(Lspaces==all_L,'all Lagrangians covered')

# Old Paley/Fano matching, including its nontrivial C3 action.
fano=((1,3),(2,6),(4,5));F={0,*map(pmask,fano)}
sigma=tuple((2*(i+1)%7)-1 for i in range(6))
require({permute(v,sigma) for v in F}==F,'Fano C3 invariant kernel')
require([v for v in V if permute(v,sigma)==v]==[0],'C3 fixed-free on nonzero two torsion')
# Affine Borel21 equivariance of all seven marked pairings.
P=lambda a:frozenset(frozenset(((a+q)%7,(a+3*q)%7)) for q in (1,2,4))
for scale,shift,a in product((1,2,4),range(7),range(7)):
    moved=frozenset(frozenset((scale*x+shift)%7 for x in pair) for pair in P(a))
    require(moved==P((scale*a+shift)%7),'Borel21 marked equivariance')

# Actual hourglass slice [E:x:z]=[s:t:2s+5t].
u=sp.symbols('u')
L=[-1-6*u,sp.Integer(1),1+u,1-u,3+5*u,-1-5*u,3+6*u]
branch=[sp.oo if sp.degree(f,u)==0 else sp.solve(f,u)[0] for f in L]
require(len(set(branch))==7,'seven distinct branch points')
all_nine_branch=branch+[sp.Rational(-1,4),sp.Rational(-3,4)]
require(len(set(all_nine_branch))==9,'slice avoids all nine-line intersections')
richelot=[]
for a in range(7):
    pairs=[((a+q)%7,(a+3*q)%7) for q in (1,2,4)]
    Fs=[sp.expand(L[i]*L[j]) for i,j in pairs]
    det=sp.Matrix([[f.coeff(u,k) for k in range(3)] for f in Fs]).det()
    require(det!=0,('Richelot determinant',a))
    Gs=[sp.expand(sp.diff(Fs[(i+1)%3],u)*Fs[(i+2)%3]-Fs[(i+1)%3]*sp.diff(Fs[(i+2)%3],u)) for i in range(3)]
    target=sp.expand(sp.prod(Gs))
    require(sp.degree(target,u)==6 and sp.degree(sp.gcd(target,sp.diff(target,u)),u)==0,('smooth target',a))
    source=sp.expand(sp.prod(Fs))
    require(source==sp.expand(sp.prod(L[i] for i in range(7) if i!=a)),('factor map',a))
    require(sp.degree(source,u) in (5,6) and sp.degree(sp.gcd(source,sp.diff(source,u)),u)==0,('smooth source',a))
    richelot.append({'omitted_label':a,'pairs':pairs,'quadratics':[str(f) for f in Fs],
                    'determinant':int(det),'brackets':[str(f) for f in Gs],
                    'geometric_target_polynomial':str(target)})
require([r['determinant'] for r in richelot]==[-162,474,-537,-105,-60,-129,-68],'pinned determinants')
# Strong negative controls: a smooth source can have a split target;
# a label action can fail on the fixed branch configuration.
Fs=[u*u-c for c in (1,4,9)]
require(sp.degree(sp.gcd(sp.prod(Fs),sp.diff(sp.prod(Fs),u)),u)==0,'split control smooth source')
require(sp.Matrix([[f.coeff(u,k) for k in range(3)] for f in Fs]).det()==0,'split control determinant zero')
M=(-25*u-19)/(25*u+15)
require(sp.limit(M,u,sp.oo)==branch[2] and M.subs(u,branch[2])==branch[4], 'Mobius first two conditions')
require((25*u+15).subs(u,branch[4])==0,'Mobius third condition')
require(M.subs(u,branch[3])==sp.Rational(-11,10) and M.subs(u,branch[3])!=branch[6],'not the Paley C3 geometry')

# Finite field trace checks: direct root fibers versus all quotient fibers.
# Arithmetic normalization/twists of Richelot targets are NOT tested or asserted.
finite=[]
coeff=[(int(f.coeff(u,0)),int(f.coeff(u,1))) for f in L]
for p in list(sp.primerange(5,102)):
    points=[(1,t) for t in range(p)]+[(0,1)]
    roots=[[j for j,(s,t) in enumerate(points) if (a*s+b*t)%p==0] for a,b in coeff]
    if len({r[0] for r in roots})<7:continue
    nr=[sum(x*x%p==v for x in range(p)) for v in range(p)]
    chi=lambda v:0 if v%p==0 else (1 if pow(v%p,(p-1)//2,p)==1 else -1)
    for inds in (tuple(range(1,7)),tuple(range(7))):
        n=len(inds); total=0; traces=Counter()
        masks=[m for m in range(1,1<<n) if m.bit_count()%2==0]
        for s,t in points:
            values=[(coeff[i][0]*s+coeff[i][1]*t)%p for i in inds]
            numer=sum(sp.prod(nr[(c*v)%p] for v in values) for c in range(1,p))
            require(numer%(p-1)==0,'projective fiber integrality')
            total+=int(numer//(p-1))
            for m in masks:
                z=1
                for j,v in enumerate(values):
                    if m>>j&1:z=z*v%p
                traces[m]+=chi(z)
        require(total-(p+1)==sum(traces.values()),('trace sum',p,n))
        require(all(traces[m]==0 for m in masks if m.bit_count()==2),('rational quotients',p,n))
        finite.append({'p':p,'n':n,'root_curve_points':total,
                       'trace_sum_by_branch_size':dict(sorted(Counter({r:sum(t for m,t in traces.items() if m.bit_count()==r) for r in range(2,n+1,2)}).items()))})

out={'status':'PASS','checks':CHECKS,'character_dimensions':rows,'lifted_pants':lifted,
     'two_torsion':{'classes':16,'symplectic_permutation_actions':len(actions),'symplectic_bases':len(symp_bases),'lagrangians':15,'pants_codecs':codec_rows},
     'slice_branch_points':[str(x) for x in branch],'richelot':richelot,'finite_field_checks':finite,
     'hostiles':{'smooth_source_split_target':True,'Paley_C3_not_automorphism_of_displayed_slice':True},
     'scope':'Exact finite controls; proofs in THEOREM.md. No integer point, global colouring, Collatz or fixed-surface PSL2(7) claim.'}
print(json.dumps(out,indent=2,sort_keys=True))
