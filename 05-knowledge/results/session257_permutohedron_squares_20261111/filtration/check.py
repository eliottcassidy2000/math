#!/usr/bin/env python3
"""Finite exact controls; runtime checks survive python -O."""
import itertools as it, math, json
from fractions import Fraction as F

def ck(c, msg):
    if not c: raise RuntimeError(msg)
def compose(a,b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a): return tuple(a.index(i) for i in range(len(a)))
def cosets(G,H):
    remaining=set(G); ans=[]
    while remaining:
        g=min(remaining); c=frozenset(compose(g,h) for h in H)
        ans.append(c); remaining-=c
    return ans
G=list(it.permutations(range(4))); ident=tuple(range(4))
t=(1,0,2,3); d=(1,0,3,2)
marks=[]
for h in [t,d]:
    X=cosets(G,[ident,h]); lookup={g:i for i,c in enumerate(X) for g in c}
    action=[lookup[compose(t,min(c))] for c in X]
    f=sum(i==j for i,j in enumerate(action)); k=(len(X)+f)//2
    marks.append({'fixed_ports':f,'orbit_count':k,'fixed_section_dimension':k-1})
ck([a['fixed_ports'] for a in marks]==[2,0], 'S4 hostile marks')
# Exact Q12 graph, no inherited component list used in discovery.
adj={v:[] for v in range(1,13)}
for a,b in it.combinations(adj,2):
    if math.isqrt(a+b)**2==a+b: adj[a].append(b); adj[b].append(a)
remaining=set(adj); paths=[]
while remaining:
    seed=min(remaining); stack=[seed]; component=set()
    while stack:
        v=stack.pop()
        if v in component: continue
        component.add(v); stack+=adj[v]
    ends=sorted(v for v in component if len(adj[v])==1)
    ck(len(ends)==2 and all(len(adj[v]) in [1,2] for v in component),'path structure')
    path=[ends[0]]
    while len(path)<len(component): path.append(next(v for v in adj[path[-1]] if v not in path))
    paths.append(tuple(path)); remaining-=component
paths.sort(key=len); ck([len(p) for p in paths]==[3,4,5],'sizes')
endpoints={v for p in paths for v in [p[0],p[-1]]}
actions=[]
for bits in it.product([0,1],repeat=3):
    action={}
    for bit,p in zip(bits,paths): action.update(zip(p,p[::-1] if bit else p))
    actions.append(action)
orbits=[]; remaining=set(adj)
while remaining:
    v=min(remaining); orbit={g[v] for g in actions}; orbits.append(orbit); remaining-=orbit
fixed={v for v in adj if all(g[v]==v for g in actions)}
ck(fixed=={3,7} and not fixed&endpoints,'fixed ports versus endpoints')
ck(sorted(map(len,orbits))==[1,1,2,2,2,2,2],'orbit sizes')
words=[]
for order in it.permutations(range(3)):
    for bits in it.product([0,1],repeat=3):
        word=tuple(v for i in order for v in (paths[i][::-1] if bits[i] else paths[i]))
        words.append(word)
ck(len(set(words))==48,'word bank')
for w in words:
    ck(len({tuple(g[v] for v in w) for g in actions})==8,'free word action')
# Generic fixed-section vertices: independent averaging of all permutation
# vertices versus predicted interval-average vertices. Projection can contain
# interior points, so compare support maxima for all small integer coefficients.
orbit_sizes=[1,2,2]; n=sum(orbit_sizes); blocks=[]; s=0
for m in orbit_sizes: blocks.append(tuple(range(s,s+m))); s+=m
projected={tuple(sum(F(v[i]) for i in B)/len(B) for B in blocks) for v in it.permutations(range(1,n+1))}
pred=[]
for order in it.permutations(range(len(blocks))):
    x=[None]*len(blocks); s=0
    for j in order: x[j]=F(2*s+orbit_sizes[j]+1,2); s+=orbit_sizes[j]
    pred.append(tuple(x))
ck(len(set(pred))==6,'section vertex prediction')
for c in it.product(range(-2,3),repeat=3):
    objective=lambda x:sum(c[i]*orbit_sizes[i]*x[i] for i in range(3))
    ck(max(map(objective,projected))==max(map(objective,pred)),'section support')
v=(5,7,3,4,1,2); w=(4,7,3,4,1,2)
images={}; edges=set()
for perm in it.permutations(range(6)):
    old=tuple(v[i] for i in perm); new=tuple(w[i] for i in perm)
    images.setdefault(new,[]).append(old)
ck(len(images)==360 and all(len(x)==2 for x in images.values()),'wall fibers')
for pair in images.values():
    a,b=pair; differences=[i for i in range(6) if a[i]!=b[i]]
    ck(len(differences)==2 and {a[i] for i in differences}=={4,5},'collapsed adjacent ranks')
# Small functional graph characteristic polynomials independently via sympy.
import sympy as sp
states=[13,2,5,3,7]
cycle={13:2,2:5,5:3,3:7,7:13}; normalized={**cycle,2:13}
polys=[]
for mapping in [cycle,normalized]:
    mat=sp.zeros(5)
    for i,q in enumerate(states): mat[states.index(mapping[q]),i]=1
    z=sp.Symbol('z'); polys.append(str(mat.charpoly(z).as_expr()))
ck(polys==['z**5 - 1','z**5 - z**3'],'functional graph characteristic polynomials')
# Exact common carrier: reversal-orbit fixed section -> Young fixed hexagon.
u=list(map(len,orbits))
groups=[[j for j,O in enumerate(orbits) if O<=set(path)] for path in paths]
weights=list(map(len,paths)); coarse_vertices=set()
for order in it.permutations(range(3)):
    y=[0]*3; L=0
    for i in order: y[i]=weights[i]*L; L+=weights[i]
    coarse_vertices.add(tuple(y))
images_of_fine=set()
for order in it.permutations(range(7)):
    r=[F(0)]*7; z=[0]*7; L=0
    for j in order:
        r[j]=F(2*L+u[j]+1,2)
        z[j]=u[j]*(r[j]-F(u[j]+1,2))
        ck(z[j]==u[j]*L,'common carrier affine identity')
        L+=u[j]
    by_projection=tuple(sum(z[a] for a in group)-sum(u[a]*u[b] for a,b in it.combinations(group,2)) for group in groups)
    by_average=[]
    for i,group in enumerate(groups):
        coarse_r=sum(u[a]*r[a] for a in group)/weights[i]
        by_average.append(weights[i]*(coarse_r-F(weights[i]+1,2)))
    ck(by_projection==tuple(by_average),'Reynolds equals refinement projection')
    images_of_fine.add(by_projection)
ck(coarse_vertices<=images_of_fine,'hexagon vertices attained')
ck(len(images_of_fine-coarse_vertices)>0,'some fine vertices project to nonvertices')
for c in it.product(range(-2,3),repeat=3):
    objective=lambda x:sum(c[i]*x[i] for i in range(3))
    ck(max(map(objective,images_of_fine))==max(map(objective,coarse_vertices)),'coarse projected support')
# Signed Paley action on the Hamilton-word subset of Pi_7.
residues={1,2,4}
hamilton={q for q in it.permutations(range(7)) if all((q[j+1]-q[j])%7 in residues for j in range(6))}
ck(len(hamilton)==189,'Paley inherited positive control')
for a in range(1,7):
    for b in range(7):
        transformed=set()
        for q in hamilton:
            r=tuple((a*x+b)%7 for x in q)
            if a not in residues: r=r[::-1]
            transformed.add(r)
        ck(transformed==hamilton,'signed Paley permutation action')
plain_nonsquare={tuple((-x)%7 for x in q) for q in hamilton}
ck(not plain_nonsquare&hamilton,'unsigned nonsquare hostile')
print(json.dumps({'status':'FINITE-EXACT PASS','S4_hostile_marks':marks,
    'Q12_paths':paths,'Q12_fixed_vertices':sorted(fixed),'Q12_fixed_endpoints':sorted(fixed&endpoints),
    'Q12_orbit_sizes':sorted(map(len,orbits)),'Q12_component_traversal_words':len(words),
    'small_section_projected_points':len(projected),'small_section_vertices':len(pred),
    'small_section_objectives_checked':125,'valuation_original_vertices':720,
    'valuation_normalized_vertices':len(images),'collapsed_edges':len(images),
    'cycle_characteristic_polynomials':polys,'Q12_Young_hexagon_vertices':len(coarse_vertices),'Q12_distinct_projected_fine_vertices':len(images_of_fine),'Q12_projected_nonvertices':len(images_of_fine-coarse_vertices),'Paley_Hamilton_words':len(hamilton),'Paley_signed_maps_checked':42,'Paley_unsigned_nonsquare_overlap':len(plain_nonsquare&hamilton)},indent=2))
