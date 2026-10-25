#!/usr/bin/env python3
"""Exact controls for projective sign characters and canonical root-cover forms.

Standard library only. The proof packet supplies all-order geometry; these
checks compare independent finite formulas and code/character operations.
"""
from collections import Counter
from itertools import combinations
from math import comb
import json

checks=0


def check(test,label):
    global checks
    checks+=1
    if not test:
        raise RuntimeError(label)


def parity(x):
    return x.bit_count()%2


def rank2(vectors):
    pivots={}
    for v in vectors:
        while v:
            p=v.bit_length()-1
            if p in pivots:
                v^=pivots[p]
            else:
                pivots[p]=v
                break
    return len(pivots)


def h0_plane(degree):
    return comb(degree+2,2) if degree>=0 else 0


rows=[]
for n in range(4,13):
    full=(1<<n)-1
    # Direct finite-cover eigensheaves.
    sheaf={s:h0_plane(s.bit_count()//2-3)
           for s in range(1<<n) if not parity(s)}
    # Independent parity-free monomial module over three degree-two variables.
    module=Counter()
    for roots in range(1<<n):
        leftover=n-6-roots.bit_count()
        if leftover>=0 and leftover%2==0:
            module[full^roots]+=h0_plane(leftover//2)
    check({s:d for s,d in sheaf.items() if d}==dict(module),
          'canonical character twist versus eigensheaf')
    pg=sum(sheaf.values())
    # Third path: coefficient of complete-intersection Hilbert series.
    m=n-6
    hilbert=sum((-1)**j*comb(n-3,j)*comb(m-2*j+n-1,n-1)
                for j in range(n-2) if m-2*j>=0)
    check(pg==hilbert,'Hilbert coefficient agrees')
    weights=Counter()
    for s,d in sheaf.items():
        if d:
            weights[s.bit_count()]+=d
    rows.append({'n':n,'pg':pg,'canonical_dimensions_by_branch_weight':dict(weights)})
    # All characters are well defined modulo simultaneous projective sign.
    for s in sheaf:
        check(parity(s&full)==0,'projective sign character')
    # Chart basis: even subsets biject with all subsets away from anchor j.
    for j in range(n):
        image={s&~(1<<j) for s in sheaf}
        check(len(image)==1<<(n-1),'normalized root chart basis')

check([row['pg'] for row in rows[:6]]==[0,0,1,7,31,111],
      'inherited first genera')

n=7
full=127
canonical=[full^(1<<i) for i in range(7)]
check(rank2(canonical)==6,'seven canonical characters have rank six')
relations=[]
for selection in range(1<<7):
    value=0
    for i in range(7):
        if selection>>i&1:
            value^=canonical[i]
    if value==0:
        relations.append(selection)
check(relations==[0,127],'unique nonempty seven-term relation')
for triple in combinations(canonical,3):
    check(triple[0]^triple[1]^triple[2]!=0,'not a Fano-line character triple')

# Label branch lines by the seven nonzero binary three-vectors.
labels=range(1,8)
simplex={sum(parity(a&v)<<i for i,v in enumerate(labels))
         for a in range(8)}
check(len(simplex)==8 and rank2(simplex)==3,'simplex code dimension')
check(Counter(v.bit_count() for v in simplex)==Counter({0:1,4:7}),
      'simplex code weights')
hamming={v for v in range(128) if all(not parity(v&s) for s in simplex)}
check(len(hamming)==16 and rank2(hamming)==4,'dual Hamming code')
check(Counter(v.bit_count() for v in hamming)==Counter({0:1,3:7,4:7,7:1}),
      'Hamming code weights')
check(full in hamming,'projective diagonal lies in Hamming kernel')
kernel={min(v,v^full) for v in hamming}
check(len(kernel)==8,'projective free candidate kernel order')
check(all(s not in canonical for s in simplex if s),'Campedelli has no canonical invariant')
restrictions={tuple((-1)**parity(s&h) for h in sorted(kernel)) for s in canonical}
check(len(restrictions)==7 and tuple(1 for h in kernel) not in restrictions,
      'seven canonical lines restrict to seven distinct nontrivial characters')
for h in kernel:
    trace=sum((-1)**parity(s&h) for s in canonical)
    check(trace==(7 if h==0 else -1),'augmentation representation trace')

# A Fano triple is exactly a dependent triple of nonzero branch labels.
fano=[]
for I in combinations(range(7),3):
    mask=sum(1<<i for i in I)
    dependent=labels[I[0]]^labels[I[1]]^labels[I[2]]==0
    check((mask in hamming)==dependent,'Fano incidence equals fixed-point kernel')
    if dependent:
        fano.append(list(I))
check(len(fano)==7,'seven Fano lines')

# Exterior-cube invariants are another exact character operation. Full deck
# has no three-term relation, but restriction to H recovers Fano incidence.
full_exterior=[]
kernel_exterior=[]
for I in combinations(range(7),3):
    char=canonical[I[0]]^canonical[I[1]]^canonical[I[2]]
    if char==0:
        full_exterior.append(list(I))
    if all(not parity(char&h) for h in kernel):
        kernel_exterior.append(list(I))
check(full_exterior==[],'full deck has no invariant alternating three-form')
check(len(kernel_exterior)==7,'order-eight kernel has seven exterior-cube invariants')
check(kernel_exterior==fano,'restricted exterior triples are exactly Fano lines')
check({v for v in hamming if not parity(v)}==simplex,
      'unique even representatives of projective kernel form simplex code')

print(json.dumps({'status':'FINITE-EXACT controls for proved character formulas',
 'checks':checks,'passed':True,
 'universe':'n=4..12 even masks and canonical parities; one labelled Fano code and all35triples',
 'genera':rows,
 'seven_canonical_rank':6,
 'simplex_code_weight_enumerator':{'0':1,'4':7},
 'hamming_code_weight_enumerator':{'0':1,'3':7,'4':7,'7':1},
 'fano_lines_on_zero_based_labels':fano,
 'exterior_cube_invariant_dimensions':{'full_deck':0,'order_eight_kernel':7},
 'scope':'No incidence census, surface isomorphism, rational-point, F1, or RH claim'},
 indent=2,sort_keys=True))
