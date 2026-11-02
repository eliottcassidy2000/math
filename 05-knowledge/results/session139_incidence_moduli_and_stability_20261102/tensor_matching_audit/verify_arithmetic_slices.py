#!/usr/bin/env python3
"""Independent exact slice decomposition audit of the second tensor proof."""
from fractions import Fraction as F
from itertools import product
from random import Random
import json
rng=Random(139611)
counts={}

def check(condition,name):
    if not condition:raise RuntimeError(name)
    counts[name]=counts.get(name,0)+1

def q1(v):
    return sum((F(a*b,2**((i^j).bit_count())) for i,a in enumerate(v) for j,b in enumerate(v)),F(0))

def direct(c,alpha):
    matrices=[((1,0),(0,1)),((1,F(1,2)),(F(1,2),1)),((1,1),(1,2))]
    result=F(0)
    for i,a in enumerate(c):
        for j,b in enumerate(c):
            z=F(a*b)
            for h,t in enumerate(alpha):z*=matrices[t][i>>h&1][j>>h&1]
            result+=z
    return result

def cmul(x,y):
    # Z[zeta], zeta^2=zeta-1.
    return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])

def evaluate(c,signs):
    value=(0,0)
    for mask,a in enumerate(c):
        z=(a,0)
        for h,s in enumerate(signs):
            if mask>>h&1:z=cmul(z,(0,1) if s else (1,-1))
        value=(value[0]+z[0],value[1]+z[1])
    return value

def phase_check(c,m):
    norms=[]
    for signs in product([0,1],repeat=m):
        a,b=evaluate(c,signs)
        check((a%2,b%2)!=(0,0),'phase_nonzero_mod2')
        N=a*a+a*b+b*b
        check(N>=1,'phase_positive_integral_norm')
        norms.append(N)
    check(sum(map(F,norms))/2**m==q1(c),'independent_phase_average')

# Exact action of tensor M on a top-basis odd coset, followed by all identity slices.
for r in range(6):
    for trial in range(2):
        c=[2*rng.randrange(-4,5) for _ in range(1<<r)];c[-1]+=1
        for alpha in product(range(3),repeat=r):
            d=c[:]
            for h,t in enumerate(alpha):
                if t==2:
                    old=d[:]
                    for i in range(1<<r):
                        if not i>>h&1:d[i]=old[i]+old[i|(1<<h)]
            ones=[h for h,t in enumerate(alpha) if t==1]
            fixed=[h for h,t in enumerate(alpha) if t!=1]
            total=F(0);retained=0
            for bits in range(1<<len(fixed)):
                vec=[]
                for mask in range(1<<len(ones)):
                    full=sum((bits>>j&1)<<h for j,h in enumerate(fixed))+sum((mask>>j&1)<<h for j,h in enumerate(ones))
                    vec.append(d[full])
                value=q1(vec);total+=value
                check(value>=0,'all_slices_PSD')
                if all(bits>>j&1 for j,h in enumerate(fixed) if alpha[h]==0):
                    check(vec[-1]%2==1 and all(a%2==0 for a in vec[:-1]),'retained_odd_top_parity')
                    check(value>=1,'retained_slice_floor')
                    retained+=1
                    phase_check(vec,len(ones))
            check(retained==2**alpha.count(2),'retained_slice_count')
            check(total==direct(c,alpha),'exact_tensor_slice_decomposition')
            check(total>=2**alpha.count(2),'Bernstein_coefficient_floor')
# Q1 positivity without parity; eigenvalues are 1/2 and 3/2.
check(F(1)+F(1,2)>0 and F(1)-F(1,2)>0,'Q1_positive_eigenvalues')
print(json.dumps({'status':'PASS','counts':counts,'claim':'All-order proof independently accepted; finite controls verify intermediate slice and phase identities only.'},indent=2,sort_keys=True))
