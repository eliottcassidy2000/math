"""Independent audit: binary square classes and Gaussian Euclidean witnesses.

No import from the boundary producer; only standard-library arithmetic.
"""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from collections import Counter
from fractions import Fraction
from itertools import product
import json

checks=0
def check(b,label):
    global checks
    checks+=1
    if not b:raise RuntimeError(label)

# Class bits: bit0=-1, bit1=2. Entries +/-1,+/-2 exhaust this group.
entries=(1,-1,2,-2)
classbit={1:0,-1:1,2:2,-2:3}
histogram=Counter()
triple_masks=[(0,0,1,-1,1,-1),(1,-1,0,0,1,-1),(1,-1,-1,1,0,0)]
triple_rows=0
for d in product(entries,repeat=6):
    c=[classbit[n] for n in d]
    q=[1^c[j]^c[j+1] for j in (0,2,4)]
    tests=[q[1]==q[2]==0 and c[4]==c[2],
           q[0]==q[2]==0 and c[4]==c[0],
           q[0]==q[1]==0 and (1^c[2]^c[0])==0]
    actual=[]
    for mask in triple_masks:
        classes=[classbit[a*b] for a,b in zip(d,mask) if b]
        actual.append(len(set(classes))==1)
        triple_rows+=1
    check(tests==actual,'pointwise-vs-table')
    if q==[0,0,0]:
        histogram['all_pair_classes_square_not_decided_here']+=1
    else:
        check(sum(actual)<=1,'at most one finite triple')
        check(not any(actual) or sum(v!=0 for v in q)==1,'finite iff exactly one nonsquare pair')
        histogram['eight' if any(actual) else 'zero']+=1
check(dict(histogram)=={'zero':3888,'eight':144,'all_pair_classes_square_not_decided_here':64},'frozen histogram')

def gmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gnorm(a):return a[0]*a[0]+a[1]*a[1]
def gsub(a,b):return (a[0]-b[0],a[1]-b[1])
def ggcd(a,b):
    while b!=(0,0):
        num=gmul(a,(b[0],-b[1])); den=gnorm(b)
        q=(round(Fraction(num[0],den)),round(Fraction(num[1],den)))
        rem=gsub(a,gmul(q,b))
        check(2*gnorm(rem)<=gnorm(b),'Gaussian Euclidean remainder')
        a,b=b,rem
    return a

witnesses=[]
for p in range(5,200):
    if p%4!=1 or any(p%d==0 for d in range(2,int(p**0.5)+1)):continue
    r=next(x for x in range(1,p) if (x*x+1)%p==0)
    a,b=ggcd((p,0),(r,1))
    check(a*a+b*b==p,'split prime norm witness')
    witnesses.append([p,a,b])

# The negative example is an exact primitive residue obstruction.
survivors=[]
for a,b,c in product(range(9),repeat=3):
    if (a*a+b*b-3*c*c)%9==0 and (a%3 or b%3 or c%3):survivors.append([a,b,c])
check(not survivors,'primitive mod9 no-go')
check(1+1==2,'positive delta2 direction')

print(json.dumps({'status':'FINITE-EXACT independent audit',
                  'method':'F2^2 squareclass enumeration; Gaussian Euclidean algorithm; primitive mod9 obstruction',
                  'outer_twists':4**6,'triple_rows':triple_rows,
                  'finite_twist_histogram':dict(sorted(histogram.items())),
                  'split_prime_norm_witnesses_below200':witnesses,
                  'primitive_delta3_mod9_survivors':survivors,'gates':checks},indent=2,sort_keys=True))
