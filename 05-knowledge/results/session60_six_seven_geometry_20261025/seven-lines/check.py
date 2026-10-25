#!/usr/bin/env python3
"""Exact seven-line, binary-code and quotient-incidence certificate; stdlib only."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import json

W=((1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1))
NAMES=('E','E+x','E-x','E+z','E-z','E+x+z','E-x-z')
T0=((0,1,2),(0,3,4),(0,5,6))
T3=tuple(sorted(T0+((1,3,6),(2,4,5))))

def need(x,msg):
    if not x: raise RuntimeError(msg)
def xor(vs):
    a=0
    for v in vs:a^=v
    return a
def det(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def rank(vs,p=0):
    a=[[v%p for v in row] if p else list(map(Fraction,row)) for row in vs]
    r=0
    for j in range(3):
        k=next((k for k in range(r,len(a)) if a[k][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        b=pow(a[r][j],-1,p) if p else 1/a[r][j]
        a[r]=[(v*b)%p if p else v*b for v in a[r]]
        for k in range(len(a)):
            if k!=r:
                b=a[k][j]
                a[k]=[(x-b*y)%p if p else x-b*y for x,y in zip(a[k],a[r])]
        r+=1
    return r

def span_binary(vs):
    out={0}
    for v in vs:out|={x^v for x in out}
    return frozenset(out)
def characteristic(indices,p):
    ans=Counter()
    for mask in range(1<<len(indices)):
        subset=[W[indices[i]] for i in range(len(indices)) if mask>>i&1]
        ans[3-rank(subset,p)]+=(-1)**len(subset)
    return {str(k):v for k,v in sorted(ans.items()) if v}
def coset(v,sub):return min(v^h for h in sub)
def plane_kernel(labels):
    hamming={mask for mask in range(128) if xor(labels[i] for i in range(7) if mask>>i&1)==0}
    even={mask for mask in hamming if mask.bit_count()%2==0}
    simplex={sum((((a&labels[i]).bit_count()%2)<<i) for i in range(7)) for a in range(8)}
    need(even==simplex,'even Hamming = simplex')
    need(Counter(m.bit_count() for m in hamming)=={0:1,3:7,4:7,7:1},'Hamming enumerator')
    need(Counter(m.bit_count() for m in simplex)=={0:1,4:7},'simplex enumerator')
    need(all((a&b).bit_count()%2==0 for a in simplex for b in simplex),'Lagrangian orthogonality')
    return tuple(sorted(simplex))

triples=list(combinations(range(7),3))
dets={str(t):det(*(W[i] for i in t)) for t in triples}
need(Counter(map(abs,dets.values()))=={0:3,1:18,2:12,3:2},'determinant histogram')
matroids={}
for p in (0,2,3,5,7,11):
    dep=[t for t in triples if rank([W[i] for i in t],p)<3]
    if p in (0,5,7,11):need(tuple(dep)==T0,'generic dependent triples')
    if p==3:need(tuple(dep)==T3,'char3 dependent triples')
    matroids[str(p)]={'dependent_triples':dep,'characteristic_polynomial':characteristic(list(range(7)),p),
                     'deletions':[{'label':NAMES[j],'triple_count':sum(j not in t for t in dep),
                                   'characteristic_polynomial':characteristic([i for i in range(7) if i!=j],p)} for j in range(7)]}
    # Independent complete rank method: determinant of all available minors.
    for mask in range(128):
        ids=[i for i in range(7) if mask>>i&1]
        r=0 if not ids else 1
        if any(any((W[i][a]*W[j][b]-W[i][b]*W[j][a])%p if p else W[i][a]*W[j][b]-W[i][b]*W[j][a]
                   for a,b in combinations(range(3),2)) for i,j in combinations(ids,2)):r=2
        if any(dets[str(t)]%p if p else dets[str(t)] for t in combinations(ids,3)):r=3
        need(r==rank([W[i] for i in ids],p),'independent rank check')
need(Counter(tuple(v%2 for v in row) for row in W)=={(1,0,0):1,(1,1,0):2,(1,0,1):2,(1,1,1):2},'char2 collapse')
# Direct affine-complement check in small fields.
for p in (2,3,5,7,11):
    count=sum(all(sum(a*b for a,b in zip(row,v))%p for row in W) for v in product(range(p),repeat=3))
    chi=matroids[str(p)]['characteristic_polynomial']
    need(count==sum(c*p**int(d) for d,c in chi.items()),'complement count')

# Canonical seven characters: unique dependence is the full seven-set.
canonical=[127^(1<<i) for i in range(7)]
need(len(span_binary(canonical))==64,'canonical rank6')
need([m for m in range(1,128) if xor(canonical[i] for i in range(7) if m>>i&1)==0]==[127],'canonical sole relation')

# Enumeration by 5040 bijections / GL(3,2), then independent Fano block systems.
quotients={}
planes=set()
for labels in permutations(range(1,8)):
    kernel=plane_kernel(labels)
    fano=tuple(t for t in triples if xor(labels[i] for i in t)==0)
    planes.add(fano)
    if kernel in quotients:
        need(quotients[kernel]['fano_lines']==fano,'kernel determines plane')
        quotients[kernel]['representatives']+=1
        continue
    d=sum(t in fano for t in T0)
    # Residual monodromy on normalization of the preimage of E.
    a=labels[0]
    mon=span_binary([a]+[labels[b]^labels[c] for _,b,c in T0])
    comps=sorted({coset(g,mon) for g in range(8)})
    edges=[]
    for _,b,c in T0:
        local=span_binary([a,labels[b],labels[c]])
        singular_points=sorted({coset(g,local) for g in range(8)})
        for point in singular_points:
            branches=[coset(point,mon),coset(point^labels[b],mon)]
            edges.append({'triple':[0,b,c], 'type':'1/4(1,1)' if labels[b]^labels[c]==a else 'A1',
                          'components':branches})
    need(len(comps)=={0:1,1:2,3:4}[d],'central normalization components')
    need(Counter(e['type'] for e in edges)=={k:v for k,v in {'A1':3-d,'1/4(1,1)':2*d}.items() if v},'local singularity count')
    # Adjunction computes strict rational central components' self intersections.
    self_int=[]
    degree=4//len(comps)
    for c0 in comps:
        bad_branches=sum(e['components'].count(c0) for e in edges if e['type']=='1/4(1,1)')
        canonical_intersection=Fraction(degree,2)-Fraction(bad_branches,2)
        self_int.append(-2-canonical_intersection)
    need(self_int==[{0:-4,1:-2,3:-1}[d]]*len(comps),'central self intersection')
    quotients[kernel]={'labels':labels,'kernel':kernel,'fano_lines':fano,'dependent_crossings':d,
                       'dependent_char3_crossings':sum(t in fano for t in T3),'representatives':1,
                       'central_components':comps,'central_degree_each':degree,
                       'central_self_intersection_each':[str(x) for x in self_int],'singularity_incidence':edges}
need(len(quotients)==len(planes)==30,'thirty quotients')
need(all(v['representatives']==168 for v in quotients.values()),'free GL action')
need(Counter(v['dependent_crossings'] for v in quotients.values())=={0:16,1:12,3:2},'16/12/2 census')
need(Counter(v['dependent_char3_crossings'] for v in quotients.values())=={0:12,1:10,2:6,4:2},'char3 incidence census')
# All order8 constant-weight4 subgroups of even signs, independently enumerated.
weight4=[m for m in range(128) if m.bit_count()==4]
simplex_spaces={span_binary(v) for v in combinations(weight4,3) if len(span_binary(v))==8 and all(m==0 or m.bit_count()==4 for m in span_binary(v))}
need(simplex_spaces=={frozenset(k) for k in quotients},'independent simplex kernel enumeration')
# Explicit fixed-locus test on zero coordinate supports, independent of labels.
for k,v in quotients.items():
    fixed=[t for t in T0 if (127^sum(1<<i for i in t)) in k]
    need(len(fixed)==v['dependent_crossings'],'geometric fixed locus test')
    need(not any(c in k for c in canonical),'no root-deletion involution in Fano kernel')

# Realized D12 arrangement action, supplied independently by symmetry lane.
R=(0,5,6,2,1,3,4)
S=(0,3,4,1,2,5,6)
def compose(p,q):return tuple(p[q[i]] for i in range(7))
def permute_kernel(p,k):
    return tuple(sorted(sum(((m>>i)&1)<<p[i] for i in range(7)) for m in k))
group={tuple(range(7))}
queue=list(group)
while queue:
    p=queue.pop()
    for h in (R,S):
        z=compose(p,h)
        if z not in group:group.add(z);queue.append(z)
need(len(group)==12,'D12 order')
remaining=set(quotients)
orbits=[]
while remaining:
    k=min(remaining)
    orb={permute_kernel(g,k) for g in group}
    need(orb<=set(quotients),'D12 action closes on kernels')
    remaining-=orb
    d=quotients[k]['dependent_crossings']
    need(all(quotients[l]['dependent_crossings']==d for l in orb),'D12 preserves crossing type')
    stabilizer=sorted(g for g in group if permute_kernel(g,k)==k)
    need(len(orb)*len(stabilizer)==12,'orbit stabilizer')
    orbits.append({'dependent_crossings':d,'size':len(orb),'stabilizer_order':len(stabilizer),
                   'representative_labels':quotients[k]['labels'],'representative_kernel':k,
                   'stabilizer_permutations':stabilizer,'kernels':sorted(orb)})
need(sorted((v['dependent_crossings'],v['size']) for v in orbits)==[(0,4),(0,12),(1,6),(1,6),(3,2)],'D12 five orbits')

out={'status':'FINITE-EXACT plus elementary proof in RESULT.md','universe':{'integer_branch_rows':W,'labels':NAMES,
     'rank_tests':6*128,'label_bijections':5040,'binary_constant_weight4_generator_triples':len(list(combinations(weight4,3)))},
     'determinants':dets,'determinant_histogram':dict(Counter(map(abs,dets.values()))),'matroids':matroids,
     'canonical_character_rank':6,'canonical_only_nonempty_relation_mask':127,
     'Fano_quotient_count':30,'generic_crossing_histogram':dict(Counter(v['dependent_crossings'] for v in quotients.values())),
     'char3_crossing_histogram':dict(Counter(v['dependent_char3_crossings'] for v in quotients.values())),
     'D12_generators':{'R':R,'S':S},'D12_orbits':orbits,
     'quotients':sorted(quotients.values(),key=lambda q:(q['dependent_crossings'],q['labels']))}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS: 768 full-subset ranks by two methods; 35 determinants; small-field complements.')
print('PASS: 5040 label bijections / 168 = 30 Fano kernels, independently constant-weight enumerated.')
print('PASS: generic dependent crossings 0:16, 1:12, 3:2; characteristic3 0:12, 1:10, 2:6, 4:2.')
print('PASS: canonical rank6, deletion involutions, central component degrees, incidence and adjunction census.')
print('PASS: realized D12 orbits: d0 sizes4+12; d1 sizes6+6; d3 size2.')
