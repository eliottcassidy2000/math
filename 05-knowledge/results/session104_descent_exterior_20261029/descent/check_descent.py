#!/usr/bin/env python3
"""Exact finite controls; standard library; no removable assertions."""
from fractions import Fraction as F
import json

checks = 0

def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

def integer(x, low=0):
    if type(x) is not int or x < low:
        raise ValueError('integer domain')

def step(n, b=3):
    integer(n, 1); integer(b, 3)
    if n % 2 != 1 or b % 2 != 1:
        raise ValueError('odd domain')
    z, a = b*n+1, 0
    while z % 2 == 0:
        z //= 2; a += 1
    return z, a

def potential(n, q):
    integer(n, 1); integer(q, 1)
    if n % 2 != 1:
        raise ValueError('odd source')
    return F(n+1) * F(3,2)**q

def budget(p):
    if not isinstance(p, F) or p < 3:
        raise ValueError('active rational potential')
    k, z = 0, p
    while z >= 3:
        z *= F(27,28); k += 1
    require(z < 3 and z/F(27,28) >= 3, 'minimal budget')
    return k

def rank(n, q):
    integer(q)
    if q == 0:
        integer(n,1)
        if n % 2 != 1:
            raise ValueError('odd terminal')
        return 0
    return 3*budget(potential(n,q))+q

def controller_edge(n,q):
    integer(q,1)
    m,a = step(n)
    return m, q+(2 if a%3 == 0 else -1), a

def word_data(b, word):
    integer(b,3)
    if b%2 != 1:
        raise ValueError('odd multiplier')
    C,D,J=0,0,0
    for a in word:
        integer(a,1)
        C=b*C+2**D; D+=a; J+=1
    return J,D,C

# Controller pointwise ranks, arbitrary queue sidecar.
edge_count=0
for n in range(1,602,2):
    for q in range(1,6):
        m,s,a=controller_edge(n,q)
        old=potential(n,q)
        if s:
            new=potential(m,s)
            require(new <= old*(F(27,28) if a%3==0 else 1), 'potential edge')
        require(rank(m,s)<rank(n,q), 'integer rank edge')
        edge_count+=1
max_steps=0
for n0 in range(1,202,2):
    for q0 in range(1,5):
        n,q=n0,q0
        initial=rank(n,q); length=0
        while q:
            if length >= initial:
                raise RuntimeError('rank bound exceeded')
            n,q,a=controller_edge(n,q); length+=1
        require(length<=initial,'return step bound')
        max_steps=max(max_steps,length)

require([rank(n,q) for n,q in [(13,1),(5,3),(1,2),(1,1),(1,0)]] == [163,162,38,4,0], 'displayed rank example')

# Independent literal trajectory versus affine recurrence and payment recurrence.
affine_count=0
for b in (3,5):
    for source in range(1,202,2):
        n=source; J=D=C=L=0
        for _ in range(16):
            m,a=step(n,b)
            Lnext=b*L+2**D*((2**a-b)*source-1)
            C=b*C+2**D; D+=a; J+=1
            require(2**D*m==b**J*source+C,'affine identity')
            L=(2**D-b**J)*source-C
            require(L==Lnext==2**D*(source-m),'payment recurrence')
            modulus=2**(D+1)
            residue=((2**D-C)*pow(b**J,-1,modulus))%modulus
            require(source%modulus==residue,'exact source cylinder')
            require((L>0)==(m<source),'original-source payment iff')
            n=m; affine_count+=1

# Actual 5 cycle, strict rank impossibility, contracting slope and zero payment.
cycle=[]; n=13
for _ in range(3):
    cycle.append(n); n,a=step(n,5)
require(cycle==[13,33,83] and n==13,'5 cycle')
J,D,C=word_data(5,[1,1,5])
require((J,D,C)==(3,7,39),'5 carry')
require(F(5**J,2**D)<1 and (2**D-5**J)*13-C==0,'slope not payment')
n,q=83,1
for expected in [(13,3),(33,2),(83,1)]:
    n,a=step(n,5); q+=2 if a%5==0 else -1
    require((n,q)==expected,'active 5 controller trap')

# Rational counterpart: exact letter replay, retained denominator; reclosure differs.
z=F(19,101)
for a,expected in zip([1,1,5],[F(79,101),F(169,101),F(19,101)]):
    z=(3*z+1)/2**a
    require(z==expected and z.denominator==101,'rational cycle')
    require(z.numerator%2==1 and z.denominator%2==1,'exact rational valuation')
closures=[]
for word in ([1,1,5,1],[1,3,3,1],[3,1,3,1]):
    j,d,c=word_data(3,word)
    closures.append(F(c,2**d-3**j))
require(closures==[F(37,35),F(221,175),F(11,7)],'reclosure not transition')

# Exact full-support model: countdown j has w_j=2^-j, mean2 and tail2^-K.
for K in range(1,31):
    cost=sum(F(j,2**j) for j in range(1,K+1))+F(K,2**K)
    require(cost==2-F(1,2**(K-1)) and cost<=2,'truncated atomic mean')
    for j in range(1,31):
        w=F(1,2**j); tail=F(1,2**K)
        require(not(tail<w) or j<=K,'atomic target tail')
        require(j<=2//w,'budget extraction integer bound')
# All halt, but finite mean not automatic: length2^j and weight2^-j give each term1.
require(sum(F(2**j,2**j) for j in range(1,41))==40,'infinite-mean partial sums')
# Strict rational decrease does not force a floor decrease.
require(F(3,2)>F(4,3) and int(F(3,2))==int(F(4,3)),'rounding hostile')
for denominator in range(1,12):
    eps=F(1,denominator)
    for numerator in range(1,50):
        x=F(numerator,7)
        y=x-eps
        if y>=0:
            require(x//eps-y//eps>=1,'uniform-gap rounding')

# Bounded-reset compiler: arbitrary finite budget profile, every legal local edge
# and every decreasing-source reset, independently computed prefix sums.
profile = [0] + [(7*n*n+3*n)%19 for n in range(1,31)]
def offset(n):
    return sum(profile[m]+1 for m in range(1,n))
reset_edges=0
for n in range(1,31):
    for r in range(1,profile[n]+1):
        for s in range(r):
            require(offset(n)+s < offset(n)+r, 'compiler internal rank')
    for m in range(1,n):
        for refreshed in range(profile[m]+1):
            require(offset(m)+refreshed < offset(n), 'compiler paid reset')
            reset_edges+=1
require(controller_edge(3,1)==(5,0,1),'unpaid first-return hostile')
require(2*(3-5)==-4,'unpaid original-source ledger')
# Actual controller entry budgets also instantiate the finite compiler exactly.
actual_budget=[0]+[rank(2*j-1,1) for j in range(1,11)]
actual_offset=[sum(actual_budget[m]+1 for m in range(1,j)) for j in range(1,11)]
for j in range(2,11):
    for h in range(1,j):
        require(actual_offset[h-1]+actual_budget[h]<actual_offset[j-1],
                'controller budget reset compiler conditional on source decrease')

# Closed logarithmic entry-rank bound and explicit paid-reset offset.
require(2*27**20<28**20,'exact twenty-contraction inequality')
for n in range(1,1002,2):
    L=n.bit_length()
    p=potential(n,1)
    require(p<=3*2**(L-1),'bit-length initial potential envelope')
    require(F(27,28)**(20*(L-1)+1)*p<3,'geometric budget witness')
    require(rank(n,1)<=60*L-56,'logarithmic native return bound')
for m in range(1,1001):
    f=lambda n:60*n*n.bit_length()
    require(f(m+1)-f(m)>=60*m.bit_length(),'closed offset adjacent gap')
for source in range(3,202,2):
    for child in range(1,source,2):
        require(f(source)-f(child)>60*child.bit_length()-56,
                'closed offset dominates refreshed rank envelope')
require(rank(1,1)==4 and rank(13,1)==163<=184,'displayed quantitative cases')

# Malformed domains must reject even under python -O.
rejections=0
for f,args in [(step,(0,)),(step,(2,)),(step,(True,)),(step,(3,4)),
               (rank,(3,-1)),(budget,(F(2),)),(budget,(3,)),
               (word_data,(3,[0])),(word_data,(3,[True]))]:
    try:
        f(*args)
    except ValueError:
        rejections+=1
    else:
        raise RuntimeError('malformed input accepted')
require(rejections==9,'all malformed controls')
print(json.dumps({'status':'pass','checks':checks,'controller_edges':edge_count,'compiled_reset_edges':reset_edges,
 'controller_complete_runs':404,'max_observed_return_steps':max_steps,
 'affine_native_steps':affine_count,'malformed_rejections':rejections,
 'scope':'finite exact controls; infinite statements depend on proofs; Collatz OPEN'},sort_keys=True,indent=2))
