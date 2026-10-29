#!/usr/bin/env python3
"""Independent adversarial controls: no producer imports, no assert statements."""
import hashlib
from collections import Counter
import itertools
import json
from pathlib import Path
from math import gcd, isqrt


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]


def walk(x,t,m):
    for _ in range(t):
        x=(x*x-x+1)%m
    return x


def jet(t,m):
    # Coefficients of f^t(2+z) mod z^3, computed without expanding f^t.
    a,b,c=2,1,0
    for _ in range(t):
        a,b,c=(a*a-a+1)%m,((2*a-1)*b)%m,((2*a-1)*c+b*b)%m
    return a,b,c


def direct_order(a,m):
    x=a%m
    for d in range(1,m+1):
        if x==1:
            return d
        x=x*a%m
    raise RuntimeError('not a unit')


def fq(a,p):
    a%=p*p
    need(gcd(a,p)==1,'quotient unit domain')
    return (pow(a,p-1,p*p)-1)//p


def endpoint(n,word,sign):
    for a in word:
        z=3*n+sign
        if z%(2**a) or (z//2**a)%2!=1:
            return None
        n=z//2**a
    return n


def main():
    p=19
    need(jet(4,p**5)[0] == 1807,'fourfold exact return')
    need(jet(4,p**5)[1] == 16575,'fourfold derivative')
    need(jet(12,p**4)[0] == 2+17*p**3,'twelvefold defect')
    need(jet(12,p)[1]==1,'twelvefold tangent identity')
    need(fq(2,p)!=0,'19 not binary Wieferich')
    periods=[]
    for k in range(1,7):
        m=p**k
        L=4 if k<=2 else 12*p**(k-3)
        need(walk(2,L,m)==2,'return at predicted period')
        for q in (2,3,19):
            if L%q==0:
                need(walk(2,L//q,m)!=2,'prime-divisor minimality')
        periods.append(L)
    # Independently expose each mod-19^3 fiber at the next three levels.
    base={walk(2,t,19**3) for t in range(12)}
    # A degree < q polynomial vanishes at a primitive q=19^t root iff
    # coefficients are equal on every coset modulo q/19. Test imbalance,
    # independently of producer polynomial-remainder reduction.
    fourier_count=0
    for b in range(19**3):
        q=19**3//gcd(b,19**3)
        if q==1:
            need(len(base)==12,'constant Fourier value')
        else:
            unit=b//gcd(b,19**3)
            coeff=Counter((unit*c)%q for c in base)
            r=min(coeff)%(q//19)
            profile=[coeff[r+j*(q//19)] for j in range(19)]
            need(len(set(profile))>1,'nonconstant cyclotomic coefficient coset')
        fourier_count+=1
    fibers=[]
    for k in range(3,7):
        counts={c:0 for c in base}; x=2
        L=12*19**(k-3)
        for _ in range(L):
            need(x%(19**3) in counts,'support projects to base')
            counts[x%(19**3)]+=1
            x=(x*x-x+1)%(19**k)
        need(set(counts.values())=={19**(k-3)},'complete base fibers')
        fibers.append([k,len(counts),min(counts.values())])
    # Direct multiplication, not order-from-Euler/factor stripping.
    rows=[]
    for p in (1093,3511):
        need(p in primes(p),'prime proof by trial division')
        for a in (2,3):
            d=direct_order(a,p)
            ds=direct_order(a,p*p)
            v=0
            z=pow(a,d)-1
            while z%p==0:
                v+=1; z//=p
            rows.append([p,a,d,ds,v,fq(a,p)])
    need(rows==[[1093,2,364,364,2,0],[1093,3,7,7651,1,312],
                [3511,2,1755,1755,2,0],[3511,3,270,947970,1,7]],'inherited clock rows')
    subgroup_rows=[]
    for p,expected in ((1093,364),(3511,3510)):
        # Closure under generators by breadth-first search, not pairwise powers.
        reached={1}; queue=[1]
        for x in queue:
            for a in (2,3):
                y=x*a%p
                if y not in reached:
                    reached.add(y); queue.append(y)
        need(len(reached)==expected,'joint subgroup by graph closure')
        d=direct_order(3,p)
        for k in range(1,5):
            z=pow(3,d*p**(k-1),p**(k+1))
            need((z-1)%p**k==0 and z!=1,'ternary full lifting kernel')
        subgroup_rows.append([p,len(reached),p*len(reached)])
    clock_count=affine_count=0
    for p in (5,7,11,13,19,1093,3511):
        m=p*p
        for J,D in itertools.product(range(-8,9),repeat=2):
            alpha=pow(3,J,m)*pow(2,-D,m)%m
            need(fq(alpha,p)==(J*fq(3,p)-D*fq(2,p))%p,'signed clocks')
            if alpha%p==1:
                need((alpha==1)==((J*fq(3,p)-D*fq(2,p))%p==0),'conditional slope lift iff')
            clock_count+=1
    # Every lift of every affine identity at these primes; direct all-point test.
    for p in (3,5,7):
        m=p*p
        for u,v in itertools.product(range(p),repeat=2):
            a,b=1+p*u,p*v
            fixed=sum((a*x+b)%m==x for x in range(m))
            need(fixed==(m if u==v==0 else (p if u else 0)),'full fixed-point count')
            for x in range(m):
                z=x
                for t in range(1,p+1):
                    z=(a*z+b)%m
                    need(z==(x+t*p*(u*x+v))%m,'iterate expansion')
                need(z==x,'period divides p')
            affine_count+=1
    word_count=source_count=lift_count=0
    for sign in (-1,1):
        for J in range(1,4):
            for word in itertools.product(range(1,5),repeat=J):
                D=sum(word); M=2**(D+1); A=3**J
                # Carry by evaluating the rational affine numerator at zero.
                C=0; den=1
                for a in word:
                    C=3*C+den; den*=2**a
                r=((den-sign*C)*pow(A,-1,M))%M
                accepted=[]
                # Test ALL positive odd representatives, not just predicted r.
                for n in range(1,M,2):
                    y=endpoint(n,word,sign)
                    need((y is not None)==(n==r),'exact cylinder both directions')
                    if y is not None:
                        accepted.append(n)
                        need(den*y==A*n+sign*C,'native endpoint')
                    source_count+=1
                need(accepted==[r],'unique exact source')
                for m in (5,7,25,35):
                    points=[]; prefixes=[set() for _ in range(J+1)]
                    for t in range(-m,0):
                        n=r+M*t; x=n
                        prefixes[0].add(x%m)
                        for i,a in enumerate(word):
                            x=(3*x+sign)//2**a
                            prefixes[i+1].add(x%m)
                        points.append((n%m,x%m))
                    need(all(len(s)==m for s in prefixes),'all-prefix CRT permutations signed parameters')
                    need(len(set(points))==m,'joint graph size')
                    need(all((den*y-A*n-sign*C)%m==0 for n,y in points),'joint graph constraint')
                for p in (5,7,11,13,19):
                    m=p*p
                    alpha=A*pow(den,-1,m)%m
                    beta=sign*C*pow(den,-1,m)%m
                    if alpha%p==1 and beta%p==0:
                        actual=(alpha==1 and beta==0)
                        predicted=((J*fq(3,p)-D*fq(2,p))%p==0 and C%m==0)
                        need(actual==predicted,'full word identity lift iff')
                        lift_count+=1
                word_count+=1
    # Hostile domains: p=3 destroys endpoint uniformity; p=2 destroys source uniformity.
    need(len({(1+6*t)%3 for t in range(3)})==1,'ternary nonunit hostile')
    need(len({(1+8*t)%2 for t in range(2)})==1,'dyadic nonunit hostile')
    root=Path(__file__).resolve().parents[1]
    inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for d in ('orders','recursions') for p in sorted((root/d).glob('*')) if p.suffix in ('.md','.py')}
    return dict(status='PASS',fourier_nonzero_checks=fourier_count,periods19=periods,fiber_counts=fibers,joint_subgroups=subgroup_rows,clock_rows=rows,signed_clock_checks=clock_count,
                full_affine_checks=affine_count,exact_signed_words=word_count,
                exhaustive_odd_sources=source_count,actual_word_identity_lift_checks=lift_count,
                source_sha256=inputs,scope='Independent finite controls; proof audit in AUDIT.md. No producer imports.')

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
