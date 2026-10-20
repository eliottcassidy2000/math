#!/usr/bin/env python3
"""Independent exact controls for branch integrality; no third-party packages."""
from fractions import Fraction
from collections import Counter
from math import factorial,gcd
import json,argparse

def require(condition, detail=None):
    if not condition:
        raise AssertionError(detail)

def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for j in range(minimum,n+1):
        for tail in partitions(n-j,j):
            yield (j,)+tail

def odd_partitions(n, minimum=1):
    if n == 0:
        yield ()
    for j in range(minimum,n+1,2):
        for tail in odd_partitions(n-j,j):
            yield (j,)+tail

def v2(a):
    require(a>0)
    return (a & -a).bit_length()-1

def zpart(mu):
    z=1
    for j,m in Counter(mu).items():
        z*=j**m*factorial(m)
    return z

def o2(mu):
    return sum(j//2 for j in mu)+sum(gcd(a,b) for i,a in enumerate(mu) for b in mu[i+1:])

def level(mu):
    return len({v2(j) for j in mu})==1

def fixed_formula(mu):
    if not level(mu):
        return 0
    exponent=o2(mu)-len(mu)+(1 if mu[0]%2 else 0)
    require(exponent>=0)
    return 1<<exponent

def rank_basis(vectors):
    basis={}
    for x in vectors:
        while x:
            pivot=x.bit_length()-1
            if pivot in basis:
                x^=basis[pivot]
            else:
                basis[pivot]=x
                break
    return basis

def fixed_linear(mu):
    n=sum(mu)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    lookup={ij:k for k,ij in enumerate(pairs)}
    perm=[]
    start=0
    for length in mu:
        perm.extend(list(range(start+1,start+length))+[start])
        start+=length
    columns=[]
    twist=0
    for k,(i,j) in enumerate(pairs):
        a,b=perm[i],perm[j]
        target=lookup[tuple(sorted((a,b)))]
        columns.append((1<<k)^(1<<target))
        if a>b:
            twist|=1<<target
    for vertex in range(n-1):
        columns.append(sum(1<<k for k,(i,j) in enumerate(pairs) if vertex in (i,j)))
    basis=rank_basis(columns)
    while twist:
        pivot=twist.bit_length()-1
        if pivot not in basis:
            return 0
        twist^=basis[pivot]
    return 1<<(len(pairs)-len(basis))

def run(max_n,linear_n):
    records=[]
    failures=[]
    odd_terms=0
    linear_types=0
    # Frozen inherited positive and hostile controls, n = 1,...,16.
    known_total=[1,1,1,2,2,6,12,79,792,19576]
    known_lev={2:Fraction(1,2),4:1,6:2,8:15,10:280,12:23464,14:6819360,16:6981258440}
    for n in range(1,max_n+1):
        odd=Fraction()
        minimum=None
        for mu in odd_partitions(n):
            odd_terms+=1
            k=len(mu)
            e=o2(mu)-k+1
            val=e-v2(zpart(mu))
            minimum=val if minimum is None else min(minimum,val)
            if val<0:
                failures.append([n,list(mu),val])
            require(e >= (n-k)//2+(k-1)*(k-2)//2)
            require(v2(zpart(mu)) <= v2(factorial(k)))
            odd+=Fraction(1<<e,zpart(mu))
        lev=Fraction()
        # All-even level partitions correspond to odd partitions of n/2^a.
        a=1
        while (1<<a)<=n:
            scale=1<<a
            if n%scale==0:
                for odd_mu in odd_partitions(n//scale):
                    mu=tuple(scale*x for x in odd_mu)
                    lev+=Fraction(fixed_formula(mu),zpart(mu))
            a+=1
        if n==2:
            require(odd==lev==Fraction(1,2))
        else:
            require(odd.denominator==lev.denominator==1)
        if n<=len(known_total):
            require(odd+lev==known_total[n-1])
        if n in known_lev:
            require(lev==known_lev[n])
        if n<=linear_n:
            total_linear=Fraction()
            odd_linear=Fraction()
            for mu in partitions(n):
                linear_types+=1
                f=fixed_linear(mu)
                require(f==fixed_formula(mu),(mu,f,fixed_formula(mu)))
                total_linear+=Fraction(f,zpart(mu))
                if all(j%2 for j in mu):
                    odd_linear+=Fraction(f,zpart(mu))
            require(odd==odd_linear and odd+lev==total_linear)
        records.append({'n':n,'N_odd':str(odd),'N_lev':str(lev),'minimum_term_v2':minimum})
    require(failures==[[2,[1,1],-1]],failures)
    return {'status':'ALL CHECKS PASSED','max_n':max_n,'linear_n':linear_n,
            'odd_partition_summands':odd_terms,'independent_linear_cycle_types':linear_types,
            'negative_term_v2_exceptions':failures,'records':records}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-n',type=int,default=64)
    ap.add_argument('--linear-n',type=int,default=10)
    ap.add_argument('--out',default='branch_checks.json')
    args=ap.parse_args()
    result=run(args.max_n,args.linear_n)
    with open(args.out,'w') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(result['status'])
    print('exact odd-branch terms:',result['odd_partition_summands'])
    print('independent affine-linear cycle types:',result['independent_linear_cycle_types'])
    print('sole nonintegral 2-adic term:',result['negative_term_v2_exceptions'])
    for row in result['records'][:16]:
        print(row['n'],row['N_odd'],row['N_lev'],row['minimum_term_v2'])
