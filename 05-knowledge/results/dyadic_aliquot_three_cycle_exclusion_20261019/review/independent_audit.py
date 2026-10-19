#!/usr/bin/env python3
"""Independent exact controls for the dyadic three-cycle proof.

No author checker is imported. Direct rational Gaussian elimination provides
an independent numerical oracle. SymPy is used solely for polynomial identities.
Finite controls are regression evidence, not the all-exponent proof.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from hashlib import sha256
import json
from pathlib import Path
import sympy as S

counts = {}

def require(condition, label):
    if not condition:
        raise RuntimeError(label)

def identity(left, right, label):
    require(S.cancel(left-right) == 0, label)
    counts['symbolic_identities'] = counts.get('symbolic_identities', 0) + 1

def gaussian(cores):
    x,y,z=cores
    rows=[[Q(-(x-1)),Q(y),Q(0),Q(2*x-1)],
          [Q(0),Q(-(y-1)),Q(z),Q(2*y-1)],
          [Q(x),Q(0),Q(-(z-1)),Q(2*z-1)]]
    for col in range(3):
        pivot=next(i for i in range(col,3) if rows[i][col])
        rows[col],rows[pivot]=rows[pivot],rows[col]
        q=rows[col][col]
        rows[col]=[t/q for t in rows[col]]
        for i in range(3):
            if i != col:
                q=rows[i][col]
                rows[i]=[a-q*b for a,b in zip(rows[i],rows[col])]
    return tuple(row[3] for row in rows)

def determinant(x,y,z):
    return x*y*z-(x-1)*(y-1)*(z-1)

def last_numerator(x,y,z):
    # Derived by substituting the third edge into the first and then the second.
    return (y-1)*(x-1)*(2*z-1)+x*(y-1)*(2*x-1)+x*y*(2*y-1)

def v2(n):
    require(n>0,'valuation positivity')
    return (n & -n).bit_length()-1

c,d,e=S.symbols('c d e')
D=determinant(c,d,e)
M=S.Matrix([[-(c-1),d,0],[0,-(d-1),e],[c,0,-(e-1)]])
b=S.Matrix([2*c-1,2*d-1,2*e-1])
raw_solution=M.inv()*b
identity(M.det(),D,'matrix determinant')
for i,(x,y,z) in enumerate(((d,e,c),(e,c,d),(c,d,e))):
    identity(raw_solution[i],last_numerator(x,y,z)/D,'independent matrix solution')
Nf=last_numerator(c,d,e)
Nr=last_numerator(d,c,e)
expanded=2*c*c*d-2*c*c+2*c*d*d+2*c*d*e-3*c*d-2*c*e+2*c-2*d*e+d+2*e-1
identity(Nf,expanded,'forward printed numerator')
identity(Nr-Nf,(c-d)*(2*(c+d)-1),'orientation difference')
low=S.expand(Nr-(c-1)*D)
high=S.expand((3*c-1)*D-Nf)
identity(low.coeff(e),(c-1)*(d-c-1),'lower slope')
identity(low.subs(e,2*d),d*(4*(c-1)*d-c*(c+1)+3)+c*(c-1),'lower boundary')
identity((4*(c-1)*d-c*(c+1)+3).subs(d,2*c),7*c*c-9*c+3,'lower bracket')
identity(high.coeff(e),3*c*c+c*d-2*c+d-1,'upper slope')
identity(high.subs(e,2*d),c*((7*c-5)*d-c+2)+2*d*(d-1),'upper boundary')
identity(((7*c-5)*d-c+2).subs(d,2*c),14*c*c-11*c+2,'upper bracket')
A=2*c*c-c+d-1
identity(Nf-(2*c-1)*D,c*(2*d*d-1)-e*A,'forward forced equation')
identity(Nr-(2*c-1)*D,2*(c-1)*(c+d*d)+d-e*A,'reverse forced equation')
identity((2*(c-1)*(c+d*d)+d).subs(d,2*c),2*c*c*(4*c-3),'reverse adjacent RHS')
identity(A.subs(d,2*c),2*c*c+c-1,'reverse adjacent odd factor')
identity((2*c*c+c-1)-(4*c-3),2*c*c-3*c+2,'reverse residual')

C,T=S.symbols('C T')
H=determinant(C,C,T)
pp,qq,rr=[S.cancel(t.subs({c:C,d:C,e:T})) for t in raw_solution]
k=S.cancel((pp+1)/C)
identity(H,C*C+2*C*T-2*C-T+1,'repeated determinant')
identity(qq,(C-1)*k+1,'repeated first edge')
identity(k,(4*C*T-3*C+2*T*T-3*T+2)/H,'repeated k formula')
identity((k-2)*H,(T-C)*(2*(T+C)-1),'repeated k difference')
identity(rr,(4*C**3+2*C*C*T-5*C*C-4*C*T+3*C+2*T-1)/H,'repeated r formula')
identity((rr-(2*C-1))*H,(C-T)*(2*C*C-1),'repeated upper')
identity((rr-(C-2))*H,3*C**3-C*C+C*T-2*C+1,'repeated lower')
identity((rr-(C-1))*H,C*C*(3*C-2)-(C-1)*T,'repeated forced equation')

# Check every ordered triple, including all cyclic rotations and repeated cores,
# by exact Gaussian elimination rather than the printed formulas.
integer_solutions=[]
for exponents in product(range(1,21),repeat=3):
    cores=tuple(1<<a for a in exponents)
    sol=gaussian(cores)
    x,y,z=cores
    den=determinant(x,y,z)
    require(den>0 and den%2==1,'positive odd determinant')
    for i,(u,v,w) in enumerate(((y,z,x),(z,x,y),(x,y,z))):
        require(sol[i]==Q(last_numerator(u,v,w),den),'Gaussian formula comparison')
    flags=[v.denominator==1 for v in sol]
    require(len(set(flags))==1,'integrality propagation')
    require(all(v>0 for v in sol),'positive formal variables')
    if all(flags):
        require(all(v.numerator%2 for v in sol),'odd formal variables')
        integer_solutions.append({'exponents':exponents,'variables':[int(v) for v in sol]})
        expected=(len(set(exponents))==1 or exponents in ((1,1,4),(1,4,1),(4,1,1)))
        require(expected,'unexpected integral word')
    if len(set(exponents))==3:
        maximum_index=cores.index(max(cores))
        require(sol[maximum_index].denominator>1,'largest variable integral')
counts['gaussian_all_ordered_exponent_triples_1_to_20']=20**3
require(len(integer_solutions)==23,'fixed-point and ghost inventory')

for aa,bb,ee in combinations(range(1,65),3):
    x,y,z=(1<<i for i in (aa,bb,ee))
    den=determinant(x,y,z)
    for u,v in ((x,y),(y,x)):
        num=last_numerator(u,v,z)
        require((x-1)*den < num < (3*x-1)*den,'global trap finite control')
        require(den%x==1 and num%x==x-1,'retained congruence')
        require(num%den!=0,'unexpected largest-variable integrality')
counts['sorted_distinct_exponents_1_to_64_both_orientations']=2*len(list(combinations(range(1,65),3)))

sparse=(1,2,3,17,128,129,512,1024,4096)
for exponents in product(sparse,repeat=3):
    cores=tuple(1<<i for i in exponents)
    sol=gaussian(cores)
    if len(set(cores))==3:
        require(sol[cores.index(max(cores))].denominator>1,'large-gap nonintegrality')
    else:
        require(all(v.denominator==1 for v in sol)==(len(set(cores))==1),'large-gap repeated classification')
counts['sparse_extreme_gaussian_triples_through_exponent_4096']=len(sparse)**3

for a,b in combinations(range(1,129),2):
    x,y=1<<a,1<<b
    rhs=2*(x-1)*(x+y*y)+y
    require(v2(2*x*x-x+y-1)==0,'forced coefficient odd')
    require(v2(x*(2*y*y-1))==a,'forward valuation')
    if b>a+1:
        require(v2(rhs)==a+1,'reverse separated valuation')
    else:
        require(v2(rhs)==2*a+1,'reverse adjacent valuation')
        require(2*x*x+x-1 != 4*x-3,'reverse adjacent odd factors unequal')
counts['valuation_pair_controls_1_to_128']=len(list(combinations(range(1,129),2)))

# Hostile controls ensure we are retaining the actual power-of-two and prime guards.
require(gaussian((2,2,16))==(Q(23),Q(13),Q(1)),'ghost formal cycle')
s16=sum(i for i in range(1,16) if 16%i==0)
require(s16==15 and s16!=46,'ghost is not an aliquot cycle')
require(gaussian((2,44,158))==(Q(393),Q(9),Q(3)),'even nondyadic hostile')
require(all(v%2==0 for v in (2,44,158)),'even nondyadic control')
counts['native_and_nondyadic_hostile_controls']=2

report={'status':'PASS','scope':'Independent audit of the exact proof; finite controls do not supply its universal quantifiers.',
        'counts':counts,'integer_controls':integer_solutions,
        'hostile_control':{'cores':[2,44,158],'variables':[393,9,3],
                           'lesson':'All-even cores and the size/congruence trap alone do not imply exclusion; exact powers of two are essential.'},
        'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(report,indent=2,sort_keys=True))
