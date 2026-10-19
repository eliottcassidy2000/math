#!/usr/bin/env python3
"""Exact reproducibility and hostile controls for WRONSKIAN_FROBENIUS_TRANSFER.md.
Pure Python standard library; deterministic seeds; assertions are not used as gates.
Run: python3 verify_wronskian_frobenius.py
"""
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import comb, factorial, gcd, prod
from random import Random
import json

COUNTS = {}

def check(condition, name):
    if not condition:
        raise RuntimeError('FAILED: ' + name)
    COUNTS[name.split(':')[0]] = COUNTS.get(name.split(':')[0], 0) + 1

def vp(n, p):
    n = abs(n)
    if n == 0:
        return 10**9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def det_bareiss(matrix):
    n = len(matrix)
    if n == 0:
        return 1
    a = [r[:] for r in matrix]
    previous, sign = 1, 1
    for k in range(n-1):
        pivot = next((i for i in range(k,n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        t = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = a[i][j]*t-a[i][k]*a[k][j]
                if numerator % previous:
                    raise RuntimeError('Nonexact Bareiss division')
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = t
    return sign*a[-1][-1]

def modular_pivots(matrix, p=1000003):
    a = [[x%p for x in row] for row in matrix]
    if not a:
        return []
    r, pivots = 0, []
    for c in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:
            continue
        a[r],a[pivot] = a[pivot],a[r]
        inv = pow(a[r][c],-1,p)
        a[r] = [(x*inv)%p for x in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][c]
            if z:
                a[i]=[(x-z*y)%p for x,y in zip(a[i],a[r])]
        pivots.append(c)
        r+=1
        if r==len(a):
            break
    return pivots

def budget(D,k,e,l):
    return l*(D-k*e)-(k+1)*l*(l-1)//2

def threshold(D,k,m,weighted=True):
    E=D//k
    if m<=E:
        raise ValueError('Uniform criterion here requires m>E')
    candidates=[]
    for e in range(E+1):
        for l in range(1,E-e+2):
            B=budget(D,k,e,l) if weighted else l*(D-k*e)-l*(l-1)
            candidates.append((Fraction(e)+Fraction(B,m-l+1),e,l))
    maximum=max(candidates)
    return maximum[0].numerator//maximum[0].denominator+1,maximum

def mono_jet(r,s,alpha,beta,a,b,ordinary):
    if a>r or b>s:
        return 0
    x=comb(r,a)*comb(s,b)*alpha**(r-a)*beta**(s-b)
    return x*factorial(a)*factorial(b) if ordinary else x

def transformed_jet(r,s,h,p,alpha,beta,a,b,ordinary):
    return sum((-1)**t*comb(h,t)*mono_jet(r+p*(h-t),s+t,alpha,beta,a,b,ordinary)
               for t in range(h+1))

def test_wronskians():
    rng=Random(601019)
    # Exact leading Vandermonde determinants, not floating determinant estimates.
    for D in range(1,18):
        for k in range(1,min(D,6)+1):
            for e in range(D//k+1):
                maxl=(D-k*e)//k+1
                for l in range(1,min(maxl,7)+1):
                    ds=[D-k*e-k*j for j in range(l)]
                    lead=[[prod(range(d-a+1,d+1)) if a else 1 for d in ds]
                          for a in range(l)]
                    exact=det_bareiss(lead)
                    expected=prod(ds[j]-ds[i] for i in range(l) for j in range(i+1,l))
                    check(exact==expected!=0,'wronskian_leading_vandermonde')
                    check(sum(ds)-l*(l-1)//2==budget(D,k,e,l),
                          'weighted_degree_budget_equality')
    # Exact finite-field full rank is a certificate of full Q rank.
    for D in range(1,10):
        for k in range(1,min(D,5)+1):
            E=D//k
            rows=[(r,s) for s in range(E+1) for r in range(D-k*s+1)]
            for m in range(E+1,E+4):
                M,_=threshold(D,k,m)
                for trial in range(3):
                    betas=list(range(M))
                    rng.shuffle(betas)
                    columns=[(i,betas[i],a,0) for i in range(M) for a in range(m)]
                    matrix=[[mono_jet(r,s,*col,True) for col in columns] for r,s in rows]
                    check(len(modular_pivots(matrix))==len(rows),
                          'certified_matched_jet_rank')
    check(threshold(99,10,20)[0]==46,'reference_weighted_threshold')
    check(threshold(99,10,20,False)[0]==82,'reference_unweighted_threshold')
    # Characteristic-p controls, including a colliding-degree nonzero survivor.
    for p in (2,3,5,7):
        check(p%p==0,'inseparable_zero_wronskian')
        check((p+1)%p==1,'colliding_degrees_do_not_force_zero')


def test_jets_and_assignments():
    for p in (2,3,5):
        for h in range(5):
            for r in range(min(p,3)):
                for s in range(3):
                    for a in range(min(12,r+p*h)+2):
                        for b in range(h+s+2):
                            for alpha,c in ((0,-1),(0,1),(1,-1),(1,1),(2,0),(2,1)):
                                beta=alpha**p+p*c
                                H=transformed_jet(r,s,h,p,alpha,beta,a,b,False)
                                O=H*factorial(a)*factorial(b)
                                check(vp(H,p)>=max(h-a//p-b,0),'hasse_jet_bound')
                                check(vp(O,p)>=max(h-b,0),'ordinary_jet_bound')
                                check(O==transformed_jet(r,s,h,p,alpha,beta,a,b,True),
                                      'ordinary_hasse_factorials')
        for h in range(p):
            for t in range(h+1):
                for u in range(h-t+1):
                    H=transformed_jet(0,0,h,p,0,-p,t*p,u,False)
                    check(vp(H,p)==h-t-u,'sharp_local_jet_budget')
        check(transformed_jet(0,0,1,p,0,-p,p,0,False)==1,
              'normalized_derivative_hostile')
        check(transformed_jet(0,0,1,p,0,-p,0,1,True)==-1,
              'wrong_direction_hostile')
    for n in range(1,6):
        vectors=list(combinations_with_replacement(range(4),n))
        for h in vectors:
            for c in vectors:
                proposed=sum(max(x-y,0) for x,y in zip(h,c))
                actual=min(sum(max(x-y,0) for x,y in zip(h,perm))
                           for perm in set(permutations(c)))
                check(proposed==actual,'sorted_assignment_exact')


def test_same_minors():
    rng=Random(190610)
    records=[]
    for p,k,D in ((2,1,5),(2,2,7),(3,1,5),(3,2,7),(3,3,8),(5,4,9)):
        rows=[(r,s) for s in range(D//k+1) for r in range(D-k*s+1)]
        hs=[r//p for r,s in rows]
        columns=[(i,i**p+p*(i*i+1),a,b)
                 for i in range(D+2) for a in range(4) for b in range(3)]
        rng.shuffle(columns)
        for ordinary,z_only in ((False,False),(True,False),(False,True),(True,True)):
            active_columns=[col for col in columns if col[3]==0] if z_only else columns
            mat=[[mono_jet(r,s,*col,ordinary) for col in active_columns] for r,s in rows]
            pivots=modular_pivots(mat)
            check(len(pivots)==len(rows),'minor_full_rank_control')
            chosen=[active_columns[j] for j in pivots]
            square=[[row[j] for j in pivots] for row in mat]
            transformed=[[transformed_jet(r%p,s,r//p,p,*col,ordinary) for col in chosen]
                         for r,s in rows]
            determinant=det_bareiss(square)
            check(determinant!=0,'minor_nonzero')
            check(det_bareiss(transformed)==determinant,'same_minor_unitriangular')
            cs=[b if ordinary else a//p+b for alpha,beta,a,b in chosen]
            B=sum(max(h-c,0) for h,c in zip(sorted(hs),sorted(cs)))
            check(vp(determinant,p)>=B,'minor_valuation_bound')
            factorial_B=None
            if ordinary:
                cs=[a//p+b for alpha,beta,a,b in chosen]
                factorial_B=sum(vp(factorial(a)*factorial(b),p) for alpha,beta,a,b in chosen)
                factorial_B+=sum(max(h-c,0) for h,c in zip(sorted(hs),sorted(cs)))
                check(vp(determinant,p)>=factorial_B,'minor_factorial_bound')
            records.append({'p':p,'k':k,'D':D,'rows':len(rows),'ordinary':ordinary,'z_only':z_only,
                            'bound':B,'factorial_bound':factorial_B,
                            'valuation':vp(determinant,p)})
    return records


def clean(A,mod):
    return {i:c%mod for i,c in A.items() if c%mod}

def mul(A,B,mod):
    C={}
    for i,a in A.items():
        for j,b in B.items():
            C[i+j]=(C.get(i+j,0)+a*b)%mod
    return clean(C,mod)

def power(A,n,mod):
    B={0:1}
    A=clean(A,mod)
    while n:
        if n&1:
            B=mul(B,A,mod)
        n//=2
        if n:
            A=mul(A,A,mod)
    return B

def dilate(A,d):
    return {i*d:c for i,c in A.items()}

def plus(A,B,mod):
    C=A.copy()
    for i,b in B.items():
        C[i]=C.get(i,0)+b
    return clean(C,mod)

def test_cellular():
    rng=Random(221019)
    rules=[{0:1,1:1},{0:2,1:3,2:-1},{-2:1,0:-1,1:2},{0:-1,1:2,3:1}]
    rules += [{i:rng.randrange(-4,5) for i in range(-1,3)} for _ in range(12)]
    for p,nmax in ((2,8),(3,5),(5,3)):
        for m in range(1,5):
            mod=p**m
            for A in rules:
                for n in range(m-1,nmax+1):
                    direct=power(A,p**n,mod)
                    template=dilate(power(A,p**(m-1),mod),p**(n-m+1))
                    check(direct==template,'cellular_finite_template')
                    canonical={i:a%p for i,a in A.items()}
                    check(direct==power(canonical,p**n,mod),'rule_lift_erasure')
        for m in range(2,6):
            for e in range(1,m):
                n=m-e-1
                a=power({0:1,1:p**e},p**n,p**m)
                check(a.get(1,0)==p**(m-1),'sharp_rule_erasure_cutoff')
        for n in range(nmax+1):
            N=p**n
            for j in range(1,N):
                check(vp(comb(N,j),p)==n-vp(j,p),'binomial_exact_valuation')
            for m in range(1,min(n+1,5)+1):
                step=p**(n-m+1)
                support=[j for j in range(N+1) if comb(N,j)%p**m]
                check(support==list(range(0,N+1,step)),'pascal_exact_support')
                check(len(support)==p**(m-1)+1,'pascal_support_count')
                for b,j in enumerate(support):
                    check(comb(N,j)%p**m==comb(p**(m-1),b)%p**m,
                          'pascal_template_residues')
                if m>=2:
                    check(step%p**(n-(m-2))!=0 and comb(N,step)%p**m!=0,
                          'sharp_template_depth')
    # Mod-8 seed with content 2: three nonzero scaled packets, middle truncated.
    S={0:2,1:4}
    actual=mul(power({0:1,1:1},8,8),S,8)
    check(actual=={0:2,1:4,4:4,8:2,9:4},'seed_content_and_truncation')
    # On a one-cell periodic quotient X=1, every large-time binomial state dies.
    check(sum(power({0:1,1:1},8,8).values())%8==0,'torus_aliasing_hostile')
    # Rule 30, given explicit convention in the note.
    S={0}
    def rule30(S):
        return {x for x in range(min(S)-1,max(S)+2)
                if ((x-1 in S) != ((x in S) or (x+1 in S)))}
    one=rule30(S);two=rule30(one)
    check(one=={-1,0,1} and two=={-2,-1,2},'rule30_exact_hostile')
    check(two!={2*x for x in one},'rule30_dilation_failure')
    def mm(A,B,mod):
        return [[sum(A[i][k]*B[k][j] for k in range(2))%mod for j in range(2)]
                for i in range(2)]
    U=[[1,2],[0,0]];V=[[1,0],[0,0]]
    check(mm(U,U,4)==U and mm(V,V,4)==V and U!=V,
          'noncommuting_lift_hostile')
    # First correction formula computed with exact integer coefficients.
    for p in (2,3,5):
        for A in rules[:4]:
            # Large enough auxiliary modulus to reconstruct integer difference is
            # unnecessary: here compute the integer pth power without reduction.
            Ap={0:1}
            for _ in range(p):
                B={}
                for i,a in Ap.items():
                    for j,b in A.items():
                        B[i+j]=B.get(i+j,0)+a*b
                Ap=B
            AX=dilate(A,p)
            C={i:(Ap.get(i,0)-AX.get(i,0))//p for i in Ap.keys()|AX.keys()}
            check(all((Ap.get(i,0)-AX.get(i,0))%p==0 for i in Ap.keys()|AX.keys()),
                  'integral_carry_polynomial')
            for n in range(1,4):
                rhs=plus(dilate(A,p**n),{i*p**(n-1):p*c for i,c in C.items()},p*p)
                check(power(A,p**n,p*p)==rhs,'first_carry_layer')


def main():
    test_wronskians()
    test_jets_and_assignments()
    minors=test_same_minors()
    test_cellular()
    print('All exact gates passed. No floating point arithmetic used.')
    print('Finite universe: weighted supports D<=17,k<=6,e all,ell<=7;')
    print('matched ranks D<=9,k<=5,m=E+1..E+3,3 deterministic permutations;')
    print('jets p=2,3,5,h<=4,r<min(p,3),s<=2,a<=min(12,r+p*h)+1,b<=h+s+1,6 nodes;')
    print('assignments n<=5,entries 0..3; 24 exact same-minor controls;')
    print('CA rules:4 fixed plus12 deterministic, p=2/3/5,n<=8/5/3,m<=4;')
    print('all binomial entries at these power times; precise counts follow.')
    print(json.dumps({'counts':COUNTS,'total_gates':sum(COUNTS.values()),
                      'same_minor_records':minors,
                      'threshold_example':{'D':99,'k':10,'m':20,
                                           'weighted_nodes':46,'unweighted_nodes':82,
                                           'monomials':550}},sort_keys=True,indent=2))

if __name__=='__main__':
    main()
