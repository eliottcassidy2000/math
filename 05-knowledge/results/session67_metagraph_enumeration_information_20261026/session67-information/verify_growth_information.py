#!/usr/bin/env python3
"""Independent six-vertex extension hostile and F2 Laplacian guardrail."""
from fractions import Fraction
from itertools import combinations,permutations
from pathlib import Path
import json
from verify_information import require
ROOT=Path(__file__).resolve().parent

def c3blowup(weights):
    labels=[i for i,k in enumerate(weights) for _ in range(k)];n=len(labels)
    return [[int(i<j) if labels[i]==labels[j] else int((labels[j]-labels[i])%3==1) for j in range(n)] for i in range(n)]

def aut(A):
    n=len(A);scores=list(map(sum,A))
    return sum(all(scores[i]==scores[p[i]] for i in range(n)) and all(A[i][j]==A[p[i]][p[j]] for i,j in combinations(range(n),2)) for p in permutations(range(n)))

def core_order(A):
    n=len(A);modules=[]
    for mask in range(1,1<<n):
        ins=[i for i in range(n) if mask>>i&1];out=[i for i in range(n) if not mask>>i&1]
        if any(len({A[v][u] for u in ins})>1 for v in out):continue
        if any(A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i] for i,j,k in combinations(ins,3)):continue
        modules.append(mask)
    return sum(not any(m!=r and m&r==m for r in modules) for m in modules)

def rank2(rows):
    basis={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i in basis:x^=basis[i]
            else:basis[i]=x;break
    return len(basis)

rows=[]
for weights in [(3,1,1),(1,2,2)]:
    A=c3blowup(weights);n=len(A);stabilizer=aut(A);parentcore=core_order(A)
    require(stabilizer==1 and parentcore==3,'parent core/rigidity')
    selected=[]
    for b in range(1<<n):
        B=[row+[(b>>i)&1] for i,row in enumerate(A)]+[[1-((b>>i)&1) for i in range(n)]+[0]]
        if core_order(B)==3:selected.append({'extension_bits':b,'child_automorphisms':aut(B)})
    total=sum(x['child_automorphisms'] for x in selected)
    kval=Fraction(total,(n+1)*stabilizer)
    rows.append({'weights':weights,'aut_parent':stabilizer,'core_order_parent':parentcore,'selected_extension_count':len(selected),'symmetry_weighted_extension_sum':total,'K_target_core3':str(kval),'selected_extensions':selected})
require([r['K_target_core3'] for r in rows]==['4/3','2'],'exact weighted-growth collision')
# Independently check integer formulas over an unbounded theorem's finite controls.
family=[]
for k in range(2,21):
    a=Fraction(3*k+2,3*k);b=Fraction(3*k+2*(k+1),3*k)
    require(b-a==Fraction(2,3),'infinite-family arithmetic')
    family.append({'k':k,'n':3*k-1,'K_A':str(a),'K_B':str(b)})
# Incidence matrix rows (one per vertex) of K_{2,3}; Laplacian over F2.
edges=[(i,j) for i in range(2) for j in range(2,5)]
D=[sum(1<<e for e,(u,v) in enumerate(edges) if i in (u,v)) for i in range(5)]
L=[sum(((D[i]&D[j]).bit_count()%2)<<j for j in range(5)) for i in range(5)]
require(rank2(D)==4 and rank2(L)==2,'F2 Laplacian hostile')
out={'status':'PASS','extension_universe':'All 32 labelled one-vertex extensions of each displayed five-vertex tournament; automorphism and core computations by independent direct permutations/subsets.','rows':rows,'family_arithmetic_controls':family,'F2_laplacian_guardrail':{'graph':'K_2,3','incidence_rank':4,'laplacian_rank':2,'vertices':5}}
(ROOT/'growth_information_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
