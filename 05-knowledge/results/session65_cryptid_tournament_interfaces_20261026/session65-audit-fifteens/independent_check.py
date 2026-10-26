#!/usr/bin/env python3
"""Independent exact audit: stdlib only, no imported research implementation."""
from collections import Counter
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import hashlib
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def cycle_lengths(p):
    unseen = set(range(len(p)))
    sizes = []
    while unseen:
        v = min(unseen)
        w, k = v, 0
        while w in unseen:
            unseen.remove(w)
            w, k = p[w], k + 1
        require(w == v, 'vertex cycle closure')
        sizes.append(k)
    return tuple(sorted(sizes))


def pair_orbits(p):
    unseen = set(combinations(range(len(p)), 2))
    out = []
    while unseen:
        start = min(unseen)
        ordered = start
        length = 0
        while tuple(sorted(ordered)) in unseen:
            unseen.remove(tuple(sorted(ordered)))
            ordered = (p[ordered[0]], p[ordered[1]])
            length += 1
        require(tuple(sorted(ordered)) == start, 'pair closure')
        out.append((length, 1 if ordered == start else -1))
    return out


def convolve(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def polynomial(p):
    value = [1]
    for length, sign in pair_orbits(p):
        value = convolve(value, [comb(length, k)*(1+sign*(-1)**k) for k in range(length+1)])
    return Counter({k:v for k,v in enumerate(value) if v})


def brute_distribution(p):
    """All directed edge bitstrings; compare each with its pushforward bitstring."""
    pairs = list(combinations(range(len(p)), 2))
    indices = {e:i for i,e in enumerate(pairs)}
    bit_images = []
    flip = 0
    for u,v in pairs:
        a,b = p[u],p[v]
        q = 1 << indices[tuple(sorted((a,b)))]
        bit_images.append(q)
        if a > b:
            flip |= q
    # Linear part of permutation on all bitstrings, evaluated by least-bit recursion.
    pushes = [0] * (1 << len(pairs))
    out = Counter()
    for mask in range(len(pushes)):
        if mask:
            bit = mask & -mask
            pushes[mask] = pushes[mask ^ bit] ^ bit_images[bit.bit_length()-1]
        out[(mask ^ pushes[mask] ^ flip).bit_count()] += 1
    return out


def structure(p):
    sizes = cycle_lengths(p)
    orbits = pair_orbits(p)
    negatives = [L for L,e in orbits if e == -1]
    require(sorted(negatives) == sorted(s//2 for s in sizes if s%2 == 0), 'negative orbit antipodes')
    poly = polynomial(p)
    minimum = sum(s%2 == 0 for s in sizes)
    require(min(poly) == minimum, 'minimum even-cycle count')
    require(all(k%2 == minimum%2 for k in poly), 'global parity')
    minimizers = 2**len(orbits)
    for L in negatives:
        minimizers *= L
    require(poly[minimum] == minimizers, 'exact minimizer count')
    require(sum(poly.values()) == 2**comb(len(p),2), 'all orientations accounted for')
    if all(s <= 2 for s in sizes):
        f, b = sizes.count(1), sizes.count(2)
        m = b*(f+b-1)
        predicted = Counter({b+2*k: 2**(comb(f,2)+b+m)*comb(m,k) for k in range(m+1)})
        require(poly == predicted, 'involution specialization')
    return sizes, poly


def tournament_tests():
    exhaustive = []
    tested = 0
    for n in range(6):
        for p in permutations(range(n)):
            sizes, pred = structure(p)
            got = brute_distribution(p)
            require(pred == got, 'brute full distribution')
            tested += sum(got.values())
        exhaustive.append({'n':n, 'permutations':factorial(n), 'tournaments_per_permutation':2**comb(n,2)})
    representatives = {}
    for p in permutations(range(6)):
        representatives.setdefault(cycle_lengths(p), p)
    type_rows = []
    for sizes,p in sorted(representatives.items()):
        _,pred = structure(p)
        got = brute_distribution(p)
        require(pred == got, 'n6 type full distribution')
        tested += sum(got.values())
        type_rows.append({'cycle_type':sizes,'distribution':dict(sorted(got.items()))})
    structurally_checked = 0
    for n in (6,7):
        for p in permutations(range(n)):
            structure(p)
            structurally_checked += 1
    return {'all_permutations_n_le_5':exhaustive,'all_types_n6':type_rows,'all_orientation_permutation_pairs_checked':tested,'additional_permutation_structure_checks':structurally_checked}


# Independent Hasse coefficient route over F_37[t]/(t^2-2), no Legendre enumeration.
P=37
ZERO=(0,0)
ONE=(1,0)
def add(x,y): return ((x[0]+y[0])%P,(x[1]+y[1])%P)
def neg(x): return ((-x[0])%P,(-x[1])%P)
def sub(x,y): return add(x,neg(y))
def mul(x,y): return ((x[0]*y[0]+2*x[1]*y[1])%P,(x[0]*y[1]+x[1]*y[0])%P)
def fpow(x,n):
    result=ONE
    while n:
        if n%2: result=mul(result,x)
        x=mul(x,x)
        n//=2
    return result
def inv(x):
    require(x != ZERO,'invert nonzero')
    return fpow(x,P*P-2)
def scale(c,x): return (c*x[0]%P,c*x[1]%P)
def curve(j):
    if j == ZERO: return ZERO,ONE
    if j == (1728%P,0): return ONE,ZERO
    k=mul(j,inv(sub(j,(1728%P,0))))
    return scale(-3,k),scale(-2,k)
def hasse(A,B):
    h=ZERO
    for i in range(19):
        q=36-3*i
        r=18-i-q
        if q>=0 and r>=0:
            coefficient=factorial(18)//(factorial(i)*factorial(q)*factorial(r))
            h=add(h,scale(coefficient,mul(fpow(A,q),fpow(B,r))))
    return h

def arithmetic_tests():
    require(pow(2,18,37)==36,'quadratic extension irreducible')
    roots=[]
    for a in range(P):
        for b in range(P):
            j=(a,b)
            A,B=curve(j)
            den=add(scale(4,fpow(A,3)),scale(27,fpow(B,2)))
            require(den!=ZERO,'smooth representative')
            j_check=mul(scale(1728*4,fpow(A,3)),inv(den))
            require(j_check==j,'j formula')
            if hasse(A,B)==ZERO:
                roots.append(j)
    require(sorted(roots)==[(3,10),(3,27),(8,0)],'supersingular Hasse zero set')
    require([fpow(j,37) for j in roots] == [(3,27),(3,10),(8,0)],'Frobenius action')
    prime_field_roots=[]
    for j in range(37):
        A,B=curve((j,0))
        points=1
        for x in range(37):
            rhs=(x**3+A[0]*x+B[0])%37
            leg=pow(rhs,18,37)
            points+=1+(0 if rhs==0 else 1 if leg==1 else -1)
        if points==38:
            prime_field_roots.append(j)
    require(prime_field_roots==[8],'independent trace zero count')
    require(brute_distribution((0,2,1))==Counter({1:4,3:4}),'37 defect distribution')
    return {'p':37,'extension':'t^2=2','hasse_coefficient_all_j_checked':37**2,'supersingular_j_coefficients':roots,'prime_field_trace_zero_j':prime_field_roots,'defect_distribution':{1:4,3:4}}


def ternary(n):
    if n==0: return '0'
    s=''
    while n:
        n,d=divmod(n,3)
        s=str(d)+s
    return s

def ternary_tests():
    omissions={str(d):[] for d in range(3)}
    for n in range(1001):
        s=ternary(2**n)
        for d in omissions:
            if d not in s:
                omissions[d].append(n)
    require(omissions['0']==[0,1,2,3,4,15],'Sloane exceptions <=1000')
    require(omissions['2']==[0,2,8],'Erdos exceptions <=1000')
    require(omissions['1']==[1,3,9],'shifted digit-one exceptions <=1000')
    require(ternary(2**15)=='1122221122','exact 2^15 ternary')
    require(ternary(2**6)=='2101','exact 64 ternary')
    require(ternary(2**15).count('2')==6,'persistence intermediate')
    tail_rows=[]
    for k in range(1,11):
        modulus=3**k
        period=2*3**(k-1)
        powers=[pow(2,n,modulus) for n in range(period)]
        require(len(set(powers))==period,'primitive root')
        tails=[ternary(a).zfill(k) for a in powers]
        zero_free=sum('0' not in s for s in tails)
        two_free=sum('2' not in s for s in tails)
        require(zero_free==2**k and two_free==2**(k-1),'sieve tail counts')
        branch=3**(k-1)-1
        require(branch%2==0,'allones even parity')
        require(pow(2,branch,modulus)==(modulus-1)//2,'allones branch')
        require(pow(2,period,3*modulus)==1+modulus,'first lift')
        for n,a in enumerate(powers):
            lifts=[pow(2,n+j*period,3*modulus) for j in range(3)]
            require(len(set(x//modulus for x in lifts))==3,'all digit lifts')
            base_digit=lifts[0]//modulus
            require([x//modulus for x in lifts]==[(base_digit+j*(a%3))%3 for j in range(3)],'carry formula')
        tail_rows.append({'k':k,'period':period,'zero_free':zero_free,'two_free':two_free,'allones_least_exponent':branch})
    for k in range(1,101):
        m=3**k
        n=3**(k-1)-1
        require(pow(2,n,m)==(m-1)//2,'allones depth100')
        require((3**k-1)-n == 2*3**(k-1),'branch compatibility')
    return {'exponent_bound':1000,'digit_omission_exponents':omissions,'all_exponent_residue_lifts_depths':tail_rows,'allones_branch_tested_to_depth':100,'persistence_chain_ternary':['2'*15,ternary(2**15),ternary(64),'0']}


def main():
    out={'status':'PASS: independent finite-exact audit, not conjecture proof','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'tournaments':tournament_tests(),'arithmetic':arithmetic_tests(),'ternary':ternary_tests()}
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
