#!/usr/bin/env python3
"""Exact controls: triangular lattice, K7 quotient, and retained winding.
Only Python's standard library is used. No source-paper proof is checked.
"""
from collections import Counter
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import json

S = ((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))
q = lambda v: (v[0]+3*v[1])%7
norm = lambda v: v[0]*v[0]+v[0]*v[1]+v[1]*v[1]
color = lambda v: (v[0]-v[1])%3
add = lambda x,y: (x[0]+y[0],x[1]+y[1])
check_count = 0

def need(p, message):
    global check_count
    check_count += 1
    if not p:
        raise ValueError(message)

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def mv(A,v):
    return [dot(row,v) for row in A]

need(set(map(q,S))==set(range(1,7)), 'six quotient directions')
need(all(norm(s)==1 and color(s)!=0 for s in S), 'unit steps and proper color')
periods=((-3,1),(1,2))
need(abs(periods[0][0]*periods[1][1]-periods[0][1]*periods[1][0])==7,'index seven')
need(all(q(p)==0 for p in periods),'kernel generators')
need(all(norm(p)==7 for p in periods),'equal shortest period candidates')
need(all(color(p)!=0 for p in periods),'ordinary three-coloring is nonperiodic')
A=[[int(i!=j) for j in range(7)] for i in range(7)]
need(mv(A,[1]*7)==[6]*7,'top eigenvector')
for j in range(1,7):
    v=[int(i==j)-int(i==0) for i in range(7)]
    need(mv(A,v)==[-x for x in v], 'six independent negative eigenvectors')

# Literal lifting: quotient edge label is the unique Euclidean unit step.
step_for={q(s):s for s in S}
triangle_winding=Counter()
triangles=[]
for i,j,k in combinations(range(7),3):
    winding=(0,0)
    for a,b in ((i,j),(j,k),(k,i)):
        winding=add(winding,step_for[(b-a)%7])
    need(q(winding)==0,'closed quotient triangle has period displacement')
    triangle_winding['zero' if winding==(0,0) else 'nonzero']+=1
    triangles.append({'vertices':[i,j,k],'winding':list(winding)})
need(triangle_winding=={'zero':14,'nonzero':21}, 'triangle census')
need(next(t['winding'] for t in triangles if t['vertices']==[0,1,2])==[3,-1], 'winding hostile')

# The true lift keeps displacement. The quotient remembers only q(displacement).
walks=Counter({(0,0):1})
moments=[]
for n in range(9):
    plane=walks[(0,0)]
    quotient=sum(count for displacement,count in walks.items() if q(displacement)==0)
    expected=(6**n+6*(-1)**n)//7
    need(quotient==expected, 'spectral quotient moment')
    need(quotient>=plane, 'forgetting displacement only adds returns')
    moments.append({'length':n,'lattice_returns_per_vertex':plane,'quotient_returns_per_vertex':quotient,'winding_returns':quotient-plane})
    nxt=Counter()
    for p,cnt in walks.items():
        for s in S:
            nxt[add(p,s)]+=cnt
    walks=nxt
need(moments[3]['lattice_returns_per_vertex']==12 and moments[3]['quotient_returns_per_vertex']==30, 'third moment split')

# Laurent identity: lambda(s,t)+3 = |1+e^(is)+e^(it)|^2.
base=((0,0),(1,0),(0,1))
laurent=Counter((a[0]-b[0],a[1]-b[1]) for a,b in product(base,repeat=2))
need(laurent==Counter({(0,0):3,**{s:1 for s in S}}),'exact spectral positivity identity')

# Finite strips/patches are genuine subgraphs and inherit the 3-coloring.
patch_records=[]
for size in range(1,10):
    points=set(product(range(-size,size+1),repeat=2))
    edges={tuple(sorted((p,add(p,s)))) for p in points for s in S if add(p,s) in points}
    need(all(color(x)!=color(y) for x,y in edges),'finite patch proper coloring')
    patch_records.append({'radius':size,'vertices':len(points),'edges':len(edges)})

# Seam gluing example: arbitrary two half-plane restrictions need not agree
# across a seam. The concrete unit pair has valid singleton restrictions.
p=(Fraction(-1,2),Fraction(0)); r=(Fraction(1,2),Fraction(0))
need(sum((p[i]-r[i])**2 for i in range(2))==1,'cross-seam unit pair')

receipt={'status':'FINITE-EXACT','checks':check_count,'unit_steps':[list(s) for s in S],
 'period_basis':[list(p) for p in periods],'quotient':'K7','quotient_spectrum':{'6':1,'-1':6},
 'triangular_lattice_chromatic_number':3,'quotient_chromatic_number':7,
 'triangle_winding_census':dict(triangle_winding),'third_trace':{'quotient':210,'lifted_zero_winding':84},
 'return_moments':moments,'patches':patch_records,'triangles':triangles,
 'scope':'Exact quotient and lift controls, not a six-chromatic Euclidean graph or Lean build.'}

# Two orientations of the SAME lattice and K7 quotient. Full enumeration is
# independent of the subset-DP path count below.
from itertools import permutations
from math import factorial

def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(A))) for j in range(len(B[0]))] for i in range(len(A))]

def charpoly(A):
    n=len(A); power=[[int(i==j) for j in range(n)] for i in range(n)]; traces=[0]
    for r in range(1,n+1):
        power=matmul(power,A); traces.append(sum(power[i][i] for i in range(n)))
    coeff=[Fraction(1)]
    for m in range(1,n+1):
        coeff.append(-sum(coeff[m-i]*traces[i] for i in range(1,m+1))/m)
    need(all(x.denominator==1 for x in coeff),'integral characteristic polynomial')
    return [int(x) for x in coeff]

def hamiltonian_dp(A):
    n=len(A); dp={(1<<v,v):1 for v in range(n)}
    for mask in range(1,1<<n):
        for v in range(n):
            count=dp.get((mask,v),0)
            for w in range(n):
                if count and not(mask>>w&1) and A[v][w]:
                    key=(mask|1<<w,w); dp[key]=dp.get(key,0)+count
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

orientations={
 'carousel':{'steps':((1,0),(0,1),(-1,1)),'residues':(1,2,3),'paths':175,'cycles':17},
 'paley':{'steps':((1,0),(0,-1),(-1,1)),'residues':(1,2,4),'paths':189,'cycles':24}
}
for name,d in orientations.items():
    Sdir=d['steps']; residue=set(map(q,Sdir))
    need(residue==set(d['residues']),'orientation quotient')
    M=[[int((j-i)%7 in residue) for j in range(7)] for i in range(7)]
    need(all(sum(row)==3 for row in M),'regular outdegree')
    need(all(M[i][j]+M[j][i]==int(i!=j) for i in range(7) for j in range(7)),'tournament')
    paths=sum(all(M[p[i]][p[i+1]] for i in range(6)) for p in permutations(range(7)))
    cycles=sum(all(M[p[i]][p[(i+1)%7]] for i in range(7)) for p in permutations(range(7)) if p[0]==0)
    need(paths==d['paths']==hamiltonian_dp(M),'independent Hamiltonian path counts')
    need(cycles==d['cycles'],'Hamiltonian cycle count modulo start')
    d['characteristic_polynomial_descending']=charpoly(M)
    d['moments']=[]
    walk=Counter({(0,0):1}); power=[[int(i==j) for j in range(7)] for i in range(7)]
    for n in range(13):
        lifted=walk[(0,0)]
        quotient=sum(v for p,v in walk.items() if q(p)==0)
        if name=='carousel':
            expected=int(n==0)
        else:
            expected=factorial(n)//factorial(n//3)**3 if n%3==0 else 0
        need(lifted==expected,'all-length formula through 12')
        need(7*quotient==sum(power[i][i] for i in range(7)),'independent matrix trace')
        d['moments'].append({'length':n,'zero_winding_per_vertex':lifted,'quotient_per_vertex':quotient})
        nxt=Counter()
        for p,count in walk.items():
            for s in Sdir:nxt[add(p,s)]+=count
        walk=nxt; power=matmul(power,M)
    d['first_spurious_return_length']=next(m['length'] for m in d['moments'] if m['quotient_per_vertex']>m['zero_winding_per_vertex'])
    d['steps']=[list(s) for s in d['steps']]
need(all(s[0]+2*s[1]>0 for s in orientations['carousel']['steps']),'acyclic height witness')
need(tuple(map(sum,zip(*orientations['paley']['steps'])))==(0,0),'Paley local triangle')
receipt['orientations']=orientations; receipt['checks']=check_count

# Paley spectral recurrence: (t-3)(t^2+t+2)^3.
def poly_mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
p=[-3,1]
for _ in range(3):p=poly_mul(p,[2,1,1])
need(list(reversed(p))==orientations['paley']['characteristic_polynomial_descending'],'Paley spectral factorization')
u=[2,-1]
for n in range(2,13):u.append(-u[-1]-2*u[-2])
for n,row in enumerate(orientations['paley']['moments']):
    need(3**n+3*u[n]==7*row['quotient_per_vertex'],'Paley exact recurrence versus walk DP')
receipt['paley_recurrence']={'u_initial':[2,-1],'rule':'u[n] = -u[n-1] - 2*u[n-2]',
    'quotient_returns':'(3^n + 3*u[n])/7','u_through_12':u}
receipt['checks']=check_count
Path(__file__).with_name('exact_checks.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ('triangles','patches')},indent=2))
