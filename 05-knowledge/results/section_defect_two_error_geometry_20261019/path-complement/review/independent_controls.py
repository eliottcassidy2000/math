#!/usr/bin/env python3
"""Independent dense-grid and full-coordinate support controls; no imports from source."""
from collections import Counter
from pathlib import Path
import hashlib
import json

p=29
ds=((1,0),(0,1),(3,9))
checks=Counter()
def require(condition,label):
    checks[label]+=1
    if not condition: raise RuntimeError(label)
def index(x,y):return x%p+p*(y%p)
def shifted(i,d,k=1):return index(i%p+k*d[0],i//p+k*d[1])
def determinant(a,b):return (a[0]*b[1]-a[1]*b[0])%p
def label(i,d):return determinant(d,(i%p,i//p))
def difference(f,d):return [f[shifted(i,d)]-f[i] for i in range(p*p)]
def three(f):
    for d in ds:f=difference(f,d)
    return f
def path(origin,e,k):
    f=[0]*(p*p)
    for j in range(1,k+1):f[shifted(origin,e,j)]=1
    return f
def support_certificate(paths,ei,component):
    d=ds[component]
    axis=ds[ei] if ei!=component else ds[(ei+1)%3]
    changes=[difference(f,d) for f in paths]
    supports=[{i for i,v in enumerate(f) if v} for f in changes]
    counts=Counter(label(i,axis) for i in supports[0]|supports[1])
    labels=[{label(i,axis) for i in s} for s in supports]
    return {'columns':len(counts),'maximum_positions':max(counts.values()),
            'shared_columns':len(labels[0]&labels[1])}

records=[];repairs=0;origins=set()
for ei,e in enumerate(ds):
    quotient=[determinant(e,d) for di,d in enumerate(ds) if di!=ei]
    require(quotient[0] and quotient[1] and quotient[0] not in (quotient[1],-quotient[1]%p),
            'quotients_distinct_up_to_sign')
    for di,d in enumerate(ds):
        if di==ei:continue
        for sign in (-1,1):
            for phase in range(p):
                origin=shifted(shifted(0,d,sign),e,phase)
                require((ei,origin) not in origins,'unique_reduced_origin')
                origins.add((ei,origin))
                paths=[path(0,e,14),path(origin,e,14)]
                before=[support_certificate(paths,ei,j) for j in range(3)]
                require([j for j in range(3) if j!=ei and before[j]['shared_columns']] == [di],
                        'only_one_shared_component')
                crowded=any(row['maximum_positions']>27 for row in before)
                require(crowded == (phase in (14,15)), 'exact_crowded_phases')
                if crowded:
                    repairs+=1
                    whole=[int(label(i,e)==0) for i in range(p*p)]
                    revised=[a-b for a,b in zip(paths[0],whole)]
                    require(Counter(revised)=={0:826,-1:15},'negative_fifteen_point_complement')
                    require(difference(revised,e)==difference(paths[0],e),'unchanged_integer_dipole')
                    require(three(revised)==three(paths[0]),'unchanged_integer_mixed_cube')
                    paths[0]=revised
                after=[support_certificate(paths,ei,j) for j in range(3)]
                for j in range(3):
                    require(after[j]['columns']<=4 and after[j]['maximum_positions']<=27,
                            'all_repaired_column_bounds')
                    if crowded and j!=ei:require(after[j]['maximum_positions']<=15,'repair_stronger_bound_fifteen')
                    for f in paths:
                        sums=[0]*p
                        for i,value in enumerate(difference(f,ds[j])):sums[label(i,ds[j])]+=value
                        require(sums==[0]*p,'all_integer_cycle_sums_zero')
                records.append({'exceptional_direction':ei,'possibly_shared_component':di,
                                'line_shift_sign':sign,'relative_phase':phase,'repair_needed':crowded,
                                'before':before,'after':after})
require(len(records)==348 and repairs==24,'complete_348_cases_and_24_repairs')

source=Path(__file__).resolve().parent.parent
submitted=json.loads((source/'path-complement-exact-results.json').read_text())
require(records==submitted['certificates'],'all_submitted_before_after_records_match')

# Extend the finite check to every ordered pair of short same-direction
# dipoles on distinct cycles, up to translation of the first origin.
# Coordinates are (s,t) in basis (e,f).  Each correction column is a 29-bit
# cyclic set. This representation is independent of the dense point method.
full=(1<<p)-1
masks={(s,k):sum(1<<((s+j)%p) for j in range(1,k+1)) for s in range(p) for k in range(1,15)}
def rotate(mask,offset):
    offset%=p
    return ((mask<<offset)|(mask>>(p-offset)))&full
def columns_for(mask1,mask2,t,a,b):
    columns={}
    for height,mask in ((0,mask1),((-b)%p,rotate(mask1,-a)),(t,mask2),((t-b)%p,rotate(mask2,-a))):
        columns[height]=columns.get(height,0)|mask
    return columns

all_count=all_repairs=0
full_repair_cases=[]
for ei,e in enumerate(ds):
    f=ds[(ei+1)%3]
    inv=pow(determinant(e,f),-1,p)
    coordinates=[]
    for di,d in enumerate(ds):
        if di==ei:continue
        a=determinant(d,f)*inv%p
        b=determinant(e,d)*inv%p
        require(b!=0,'nonparallel_grid_coordinate')
        coordinates.append((di,a,b))
    for t in range(1,p):
        shared=[di for di,a,b in coordinates if t in (b,-b%p)]
        require(len(shared)<=1,'complete_quotient_at_most_one_shared_component')
        for s in range(p):
            for k1 in range(1,15):
                for k2 in range(1,15):
                    all_count+=1
                    m1,m2=masks[0,k1],masks[s,k2]
                    before=[columns_for(m1,m2,t,a,b) for di,a,b in coordinates]
                    crowded=any(mask.bit_count()>27 for cols in before for mask in cols.values())
                    if crowded:
                        all_repairs+=1
                        require(k1==k2==14,'only_two_length_fourteen_paths_can_crowd')
                        m1^=full
                        full_repair_cases.append([ei,s,t,k1,k2])
                    after=[columns_for(m1,m2,t,a,b) for di,a,b in coordinates]
                    require(all(len(cols)<=4 and max(mask.bit_count() for mask in cols.values())<=27 for cols in after),
                            'all_short_same_direction_geometries_normalize')
                    if crowded:
                        require(all(max(mask.bit_count() for mask in cols.values())<=15 for cols in after),
                                'all_complement_repairs_have_capacity_fifteen')
require(all_count==3*28*29*14*14==477456 and all_repairs==24,'complete_477456_geometry_universe')

# Deliberate sign error: replacing S by +indicator(complement) reverses
# the endpoint dipole, so the exact-preservation check must reject it.
f=path(0,ds[0],14)
positive=[int(label(i,ds[0])==0)-f[i] for i in range(p*p)]
require(difference(positive,ds[0])!=difference(f,ds[0]),'wrong_complement_sign_rejected')

receipt={'verdict':'PASS','reduced_geometries':len(records),'reduced_repairs':repairs,
         'complete_short_same_direction_geometries':all_count,'complete_repairs':all_repairs,
         'complete_repair_cases_basis_e_next_direction':full_repair_cases,
         'frozen_readme_sha256':hashlib.sha256((source/'README.md').read_bytes()).hexdigest(),
         'certificates':records,'checks':dict(sorted(checks.items())),
         'scope':'Normalization is automatic for same-direction dipoles on distinct cycles; disjoint cube supports and homogeneous dimension four still required for exclusion.'}
Path(__file__).with_name('independent-results.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ('certificates','complete_repair_cases_basis_e_next_direction')},indent=2))
