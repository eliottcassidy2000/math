"""Exact divisor-domain and marked Pythagorean-fibre controls.

Classical/inherited classifications are not novelty claims. No primality
or convergence conclusion is inferred from a prime-support quotient.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, lcm, prod
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(n, minimum=1):
    need(type(n) is int and n >= minimum, 'exact integer in domain')
    return n


def factor(n):
    integer(n)
    answer = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                n //= p
                a += 1
            answer.append((p,a))
        p = 3 if p == 2 else p+2
    if n > 1:
        answer.append((n,1))
    return tuple(answer)


def prime(p):
    integer(p,2)
    need(factor(p) == ((p,1),), 'prime label required')
    return p


def divisors(n):
    integer(n)
    ds = [1]
    for p,a in factor(n):
        ds = [d*p**j for d in ds for j in range(a+1)]
    return tuple(sorted(ds))


def profile_counts(exponents, k=2):
    need(type(exponents) is tuple and len(exponents)>0, 'nonempty literal profile')
    for a in exponents:
        integer(a)
    integer(k)
    r = len(exponents)
    f = prod(a+1 for a in exponents)-2
    s = prod(min(a+1,k) for a in exponents)-1-int(all(a<k for a in exponents))
    u = r-int(exponents == (1,))
    return f,s,u


def counts(n, k=2):
    integer(n,2)
    return profile_counts(tuple(a for _,a in factor(n)),k)


def shape(p,q,r,a):
    prime(p),prime(q),prime(r),integer(a,2)
    need(len({p,q,r}) == 3, 'three distinct prime labels')
    return p**a*q*r


def predecessor_data(n):
    """L is lcm(p-1), not the Carmichael function for arbitrary prime powers."""
    integer(n,2)
    fs = factor(n)
    bound = lcm(*(p-1 for p,_ in fs))
    squarefree = all(a == 1 for _,a in fs)
    return {'source':n,'factorization':fs,'predecessor_factorization':factor(n-1),
            'prime_predecessors':tuple(p-1 for p,_ in fs),
            'prime_predecessor_lcm':bound,
            'slack':str(Fraction(n-1,bound)),
            'korselt_composite':squarefree and len(fs)>=2 and (n-1)%bound == 0}


def coprime_splits(n):
    integer(n,2)
    powers = tuple(p**a for p,a in factor(n))
    pairs = []
    for choices in product((0,1),repeat=len(powers)):
        u = prod(x for x,chosen in zip(powers,choices) if chosen)
        v = n//u
        if u<v:
            pairs.append((u,v))
    return tuple(sorted(pairs))


def euclid_pair(m,n):
    integer(m,2),integer(n)
    need(m>n and gcd(m,n)==1 and (m-n)%2==1, 'primitive positive Euclid pair')
    return (m*m-n*n,2*m*n,m*m+n*n)


def branch(pair, letter):
    need(type(pair) is tuple and len(pair)==2, 'marked parameter pair')
    m,n = pair
    euclid_pair(m,n)
    need(type(letter) is str and letter in ('U','A','D'), 'Berggren branch label')
    return {'U':(2*m-n,m),'A':(2*m+n,m),'D':(m+2*n,n)}[letter]


def ancestry(m,n):
    euclid_pair(m,n)
    reversed_word = []
    while (m,n)!=(2,1):
        old = m+n
        if m<2*n:
            letter,(m,n) = 'U',(n,2*n-m)
        elif m<3*n:
            letter,(m,n) = 'A',(n,m-2*n)
        else:
            letter,(m,n) = 'D',(m-2*n,n)
        euclid_pair(m,n)
        need(m+n<old,'strict parameter descent')
        reversed_word.append(letter)
    return ''.join(reversed(reversed_word))


def decode_branch(word):
    need(type(word) is str,'literal branch word')
    pair = (2,1)
    for letter in word:
        pair = branch(pair,letter)
    return pair


@dataclass(frozen=True)
class Triangle:
    odd_leg: int
    even_leg: int
    hypotenuse: int
    parameters: tuple
    branch_word: str


def triangle(m,n):
    a,b,c = euclid_pair(m,n)
    return Triangle(a,b,c,(m,n),ancestry(m,n))


def audit_triangle(receipt, role, source):
    need(type(receipt) is Triangle,'literal marked triangle')
    integer(source)
    need(type(role) is str and role in ('odd_leg','even_leg','hypotenuse'),'source role')
    for x in (receipt.odd_leg,receipt.even_leg,receipt.hypotenuse):
        integer(x)
    need(type(receipt.parameters) is tuple and len(receipt.parameters)==2,'Euclid sidecar')
    need(type(receipt.branch_word) is str,'literal ancestry word')
    expected = triangle(*receipt.parameters)
    need(receipt==expected and decode_branch(receipt.branch_word)==receipt.parameters,'exact primitive geometry and ancestry')
    need(getattr(receipt,role)==source,'same source in the declared role')
    return True


def odd_leg_fibre(n):
    integer(n,3)
    need(n%2==1,'odd leg')
    return tuple(triangle((u+v)//2,(v-u)//2) for u,v in coprime_splits(n))


def even_leg_fibre(n):
    integer(n,2)
    need(n%2==0,'even leg')
    if n%4:
        return ()
    return tuple(triangle(v,u) for u,v in coprime_splits(n//2))


def norm_pairs(n):
    integer(n)
    answer=[]
    for a in range(isqrt(n)+1):
        b=isqrt(n-a*a)
        if a<=b and a*a+b*b==n:
            answer.append((a,b))
    return tuple(answer)


def hypotenuse_fibre(n):
    integer(n,3)
    return tuple(triangle(b,a) for a,b in norm_pairs(n)
                 if a>0 and gcd(a,b)==1 and (b-a)%2==1)


def main():
    checks=0
    def check(ok,msg):
        nonlocal checks
        checks+=1
        if not ok:
            raise RuntimeError(msg)
    def reject(fn,*args):
        try:
            fn(*args)
        except ValueError:
            check(True,'hostile rejected')
        else:
            check(False,'invalid source accepted')

    for n in range(2,1501):
        proper=tuple(d for d in divisors(n) if 1<d<n)
        for k in (2,4):
            f,s,u=counts(n,k)
            direct_s=sum(all(a<k for _,a in factor(d)) for d in proper)
            direct_u=sum(factor(d)==((d,1),) for d in proper)
            check((f,s,u)==(len(proper),direct_s,direct_u),'independent divisor predicates')
        fs=factor(n)
        profile=sorted(a for _,a in fs)
        f,s,u=counts(n)
        check((f==s+u)==(profile in ([1],[3],[1,1,2])),'inherited squarefree balance')
        for k in (2,3,4,5):
            f,s,u=counts(n,k)
            expected=profile in ([1],[k+1],[1,1,k]) or (k>=3 and profile==[2,k])
            check((f==s+u)==expected,'inherited k-free balance')

    profile_controls=0
    for r in range(1,5):
        for exponents in product(range(1,7),repeat=r):
            f,s,u=profile_counts(exponents,4)
            profile=sorted(exponents)
            check((f==s+u)==(profile in ([1],[5],[2,4],[1,1,4])),'four-free profile universe')
            profile_controls+=1
    for p,q,r in ((2,3,5),(3,5,13),(5,2,13),(3,11,17)):
        for a in range(2,13):
            n=shape(p,q,r,a)
            f,s,u=counts(n)
            check((f,s,u)==(4*a+2,7,3) and f-s-u==4*(a-2),'marked exponent transport')
            check(counts(n,a)==(4*a+2,4*a-1,3),'correctly changed freeness threshold')
    check(shape(2,3,5,4)==240 and counts(240)==(18,7,3),'minimal p4qr hostile')
    check(counts(60)==(10,7,3) and counts(240,4)==(18,15,3),'positive old and changed-domain cases')

    pred={n:predecessor_data(n) for n in (561,1105,1729,85,65)}
    check(pred[561]['prime_predecessor_lcm']==80 and pred[561]['slack']=='7','561 predecessor slack')
    check(pred[561]['predecessor_factorization']==((2,4),(5,1),(7,1)),'actual 2,5,7 carrier')
    check(all(pred[n]['korselt_composite'] for n in (561,1105,1729)),'three certified Korselt inputs')
    check(not pred[85]['korselt_composite'] and 84%16!=0,'radical inclusion loses valuation')
    for n in (561,1105,1729):
        f,s,u=counts(n)
        check((f,s,u)==(6,6,3) and f-s-u==-3,'squarefree Carmichael is not divisor balance')
        check(all(pow(a,n,n)==a for a in range(n)),'literal all-residue Fermat control')
    check(1729==1**3+12**3==9**3+10**3,'independent taxicab identity')

    for n in range(3,802,2):
        fibre=odd_leg_fibre(n)
        check(len(fibre)==2**(len(factor(n))-1),'complete odd-leg factor partitions')
        check(len({t.parameters for t in fibre})==len(fibre),'distinct marked geometries')
        for t in fibre:
            check(audit_triangle(t,'odd_leg',n),'source-owned odd-leg receipt')
            u,v=isqrt(t.hypotenuse-t.even_leg),isqrt(t.hypotenuse+t.even_leg)
            check(u*v==n and gcd(u,v)==1 and u*u==t.hypotenuse-t.even_leg
                  and v*v==t.hypotenuse+t.even_leg,'lossless coprime split inverse')
    for m in range(2,71):
        for n in range(1,m):
            if gcd(m,n)==1 and (m-n)%2:
                t=triangle(m,n)
                check(t in odd_leg_fibre(t.odd_leg),'independent Euclid converse')
                check(t in even_leg_fibre(t.even_leg),'independent even-leg converse')
                check(t in hypotenuse_fibre(t.hypotenuse),'independent norm converse')

    odd561=odd_leg_fibre(561)
    even560=even_leg_fibre(560)
    hyp1105=hypotenuse_fibre(1105)
    check(len(odd561)==len(even560)==len(hyp1105)==4,'three distinct four-element fibres')
    check(hypotenuse_fibre(561)==hypotenuse_fibre(1729)==(),'inert-prime hypotenuse obstruction')
    check(norm_pairs(561)==norm_pairs(1729)==(),'squarefree norm obstruction')
    check(any((t.odd_leg,t.even_leg,t.hypotenuse)==(47,1104,1105) for t in hyp1105),'positive 1105 bridge')
    check(even_leg_fibre(70)==(),'radical quotient destroys even-leg admissibility')
    check(norm_pairs(561*561)==((0,561),(264,495)),'nonprimitive hypotenuse boundary')
    check(gcd(264,495)==33 and (264//33,495//33,561//33)==(8,15,17),'content sidecar')
    check(all(gcd(a,b)%3==0 for a,b in norm_pairs(shape(3,5,13,2))),'squarecontent p2')
    check(all(gcd(a,b)%9==0 for a,b in norm_pairs(shape(3,5,13,4))),'squarecontent p4')

    reject(counts,True)
    reject(counts,1)
    reject(shape,2,2,5,4)
    reject(shape,4,3,5,4)
    reject(shape,2,3,5,4.0)
    reject(odd_leg_fibre,560)
    reject(odd_leg_fibre,1)
    reject(even_leg_fibre,561)
    reject(branch,(2,1),'B')
    reject(audit_triangle,replace(odd561[0],odd_leg=561.0),'odd_leg',561)
    reject(audit_triangle,odd561[0],'hypotenuse',561)
    reject(audit_triangle,replace(odd561[0],parameters=(True,2)),'odd_leg',561)
    reject(audit_triangle,replace(odd561[0],branch_word=''),'odd_leg',561)

    def compact(t):
        runs=[]
        for ch in t.branch_word:
            if runs and runs[-1][0]==ch:
                runs[-1][1]+=1
            else:
                runs.append([ch,1])
        return {'triple':[t.odd_leg,t.even_leg,t.hypotenuse],
                'euclid':t.parameters,'ancestry_runs':runs}
    report={'checks':checks,'status':'Exact source-role interfaces; inherited divisor and Pythagorean classifications, no new prime-density or Collatz claim',
            'direct_integer_range':[2,1500],'four_free_profile_controls':profile_controls,
            'p4qr_hostile':{'N':240,'F_S_U':counts(240),'F_S4_U':counts(240,4)},
            'predecessors':pred,'odd_leg_561':[compact(t) for t in odd561],
            'even_leg_560':[compact(t) for t in even560],
            'hypotenuse_1105':[compact(t) for t in hyp1105]}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
