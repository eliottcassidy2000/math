#!/usr/bin/env python3
"""Exact complete finite-universe controls for the scalar template classification.
Run beside verify_wronskian_frobenius.py: python3 verify_template_depth.py
No third-party imports, floats, or optimization-sensitive assert gates.
"""
from itertools import product
from math import gcd, comb
import json
from verify_wronskian_frobenius import check, COUNTS, vp, power, mul, dilate, clean


def integer_power(A,n):
    B={0:1}
    for _ in range(n):
        C={}
        for i,a in B.items():
            for j,b in A.items():
                C[i+j]=C.get(i+j,0)+a*b
        B={i:c for i,c in C.items() if c}
    return B

def carry(A,p):
    AP=integer_power(A,p)
    AX=dilate(A,p)
    numerators={i:AP.get(i,0)-AX.get(i,0) for i in AP.keys()|AX.keys()}
    check(all(c%p==0 for c in numerators.values()),'carry_integrality')
    return {i:c//p for i,c in numerators.items() if c}

def compatible(A,p,m):
    return power(A,p,p**m)==clean(dilate(A,p),p**m)

def residue_classes(A,p):
    return {i%p for i,a in A.items() if a%p}

def test_carry_obstruction():
    universe=[]
    for p,D in ((2,5),(3,4),(5,3),(7,2)):
        count=0
        for coeffs in product(range(p),repeat=D+1):
            A={i:a for i,a in enumerate(coeffs) if a}
            C=carry(A,p)
            offgrid=any(i%p and c%p for i,c in C.items())
            expected=len(residue_classes(A,p))<=1
            check((not offgrid)==expected,'intrinsic_first_obstruction_iff')
            # Independently check the exact integral derivative identity modulo p.
            fp={i-1:(i*a)%p for i,a in A.items() if i*a%p}
            lhs=clean({i-1:i*c for i,c in C.items()},p)
            term1=mul(power(A,p-1,p),fp,p)
            term2={i+p-1:-c for i,c in power(fp,p,p).items()}
            rhs=term1.copy()
            for i,c in term2.items():
                rhs[i]=rhs.get(i,0)+c
            check(lhs==clean(rhs,p),'carry_derivative_identity')
            # A deliberately changed lift; transverse class must be unchanged.
            B={-1:1,0:2,2:3}
            changed=A.copy()
            for i,b in B.items():
                changed[i]=changed.get(i,0)+p*b
            C1=carry(changed,p)
            desired=C.copy()
            for i,b in B.items():
                desired[p*i]=desired.get(p*i,0)-b
            check(clean(C1,p)==clean(desired,p),'change_of_lift_law')
            if expected:
                correction={i//p:c%p for i,c in C.items() if c%p}
                repaired=A.copy()
                for i,b in correction.items():
                    repaired[i]=repaired.get(i,0)+p*b
                check(compatible(repaired,p,2),'constructed_compatible_first_lift')
            count+=1
        universe.append({'p':p,'degree_bound':D,'all_mod_p_rules':count})
    return universe


def test_granularity():
    universe=[]
    for p,D in ((2,6),(3,4),(5,3)):
        count=0
        for coeffs in product(range(p),repeat=D+1):
            A={i:a for i,a in enumerate(coeffs) if a}
            S=sorted(A)
            for m in range(1,4):
                for n in (m-1,m):
                    K=power(A,p**n,p**m)
                    if not S:
                        check(not K,'zero_kernel_classification')
                    elif len(S)==1:
                        check(set(K)=={S[0]*p**n},'monomial_kernel_classification')
                    else:
                        d=S[0]
                        s=min(vp(i-d,p) for i in S[1:])
                        actual=min(vp(i-min(K),p) for i in K if i!=min(K))
                        check(actual==n-m+1+s,'complete_kernel_granularity')
                        scale=p**s
                        B={(i-d)//scale:a for i,a in A.items()}
                        j=min(i for i in B if i%p)
                        N=p**(m-1)
                        exact=integer_power(B,N).get(j,0)
                        expected=N*B[j]*B[0]**(N-1)
                        check(exact==expected,'minimal_exponent_exact_coefficient')
            count+=1
        universe.append({'p':p,'degree_bound':D,'all_mod_p_rules':count,
                         'm':'1..3','n':'m-1 and m'})
    return universe


def expected_patterns(p,m,D):
    stride=p**(m-1)
    return {f for f in product(range(p),repeat=D+1)
            if len({i%stride for i,a in enumerate(f) if a})<=1}

def test_all_lifts():
    records=[]
    # Every coefficient vector in each stated modulus/degree box is visited.
    for p,m,D in ((2,2,3),(2,3,3),(2,3,4),(2,4,2),(3,2,3),(3,3,2),(5,2,2)):
        mod=p**m
        found=set(); compatible_count=0
        for coeffs in product(range(mod),repeat=D+1):
            A={i:a for i,a in enumerate(coeffs) if a}
            ok=compatible(A,p,m)
            predicted=len({i%p**(m-1) for i,a in A.items() if a%p})<=1
            check(not ok or predicted,'exhaustive_lift_necessity')
            if ok:
                found.add(tuple(a%p for a in coeffs))
                compatible_count+=1
        expected=expected_patterns(p,m,D)
        check(found==expected,'exhaustive_lift_existence_iff')
        # Reconstruct one explicit degree-preserving lift for every admissible f.
        for f in expected:
            S=[i for i,a in enumerate(f) if a]
            if not S:
                A={}
            else:
                d=min(S);N=p**(m-1)
                B={(i-d)//N:f[i] for i in S}
                A={i+d:a for i,a in integer_power(B,N).items()}
            check(clean(A,p)=={i:a for i,a in enumerate(f) if a},'explicit_lift_reduces_correctly')
            check(compatible(A,p,m),'explicit_all_precision_lift')
            check(not A or min(A)>=0 and max(A)<=D,'lift_preserves_degree_interval')
        records.append({'p':p,'m':m,'D':D,'all_lifts':mod**(D+1),
                        'compatible_lifts':compatible_count,
                        'admissible_mod_p_patterns':len(found),
                        'all_mod_p_patterns':p**(D+1)})
    return records


def test_finite_logarithm():
    # Full root/derivative census through 199 plus the two supplied Wieferich controls.
    from math import isqrt
    primes=[p for p in range(3,200) if all(p%d for d in range(2,isqrt(p)+1))]
    roots={}
    for p in primes:
        coeffs=[0]+[(comb(p,j)//p)%p for j in range(1,p)]
        check(all(coeffs[j]==(-1)**(j-1)*pow(j,-1,p)%p for j in range(1,p)),
              'finite_logarithm_coefficients')
        H=lambda r:sum(c*pow(r,j,p) for j,c in enumerate(coeffs))%p
        D=lambda r:sum(j*c*pow(r,j-1,p) for j,c in enumerate(coeffs) if j)%p
        DD=lambda r:sum(j*(j-1)*c*pow(r,j-2,p) for j,c in enumerate(coeffs) if j>=2)%p
        check(H(0)==0 and H(p-1)==0 and D(0)==1 and D(p-1)==p-1,
              'finite_logarithm_simple_endpoints')
        found=[]
        for r in range(1,p-1):
            check(D(r)==0,'finite_logarithm_critical_points')
            check(DD(r)==pow(r*(r+1),-1,p),'finite_logarithm_nonzero_second_derivative')
            if H(r)==0:
                check(D(r)==0 and DD(r)!=0,'finite_logarithm_exact_double_roots')
                found.append(r)
        roots[str(p)]=found
        for r in range(1,p-1):
            check(H((-1-r)%p)==H(r),'finite_logarithm_reflection')
            check((pow(r,p,p)*H(pow(r,-1,p))-H(r))%p==0,
                  'finite_logarithm_inversion')
        if p>3:
            # Work exactly in Z[omega], omega^2+omega+1=0.
            def wm(u,v):
                a,b=u;c,d=v
                return (a*c-b*d,a*d+b*c-b*d)
            def wp(u,n):
                v=(1,0)
                for _ in range(n):
                    v=wm(v,u)
                return v
            w=(0,1);onew=(1,1)
            a=wp(onew,p);b=wp(w,p)
            check((a[0]-1-b[0],a[1]-b[1])==(0,0),
                  'cyclotomic_exact_factor')
            a=wp(onew,p-1);b=wp(w,p-1)
            derivative=(a[0]-b[0],a[1]-b[1])
            if p%3==1:
                a=wp(onew,p-2);b=wp(w,p-2)
                second=((p-1)*(a[0]-b[0]),(p-1)*(a[1]-b[1]))
                check(derivative==(0,0) and second==(p-1,0),
                      'cyclotomic_exact_double_factor')
            else:
                check(derivative!= (0,0),'cyclotomic_exact_simple_factor')
        q=(pow(2,p-1,p*p)-1)//p
        check(H(1)==2*q%p,'finite_logarithm_fermat_specialization')
    for p in (1093,3511):
        check(pow(2,p-1,p*p)==1,'wieferich_known_positive_controls')
        # At r=1 exact identities prove multiplicity 2 without storing all coefficients.
        H1=((pow(2,p,p*p)-2)//p)%p
        check(H1==0 and pow(2,p-1,p)-1==0 and pow(2,-1,p)!=0,
              'wieferich_double_root_controls')
    return roots


def main():
    first=test_carry_obstruction()
    grid=test_granularity()
    lifts=test_all_lifts()
    roots=test_finite_logarithm()
    print('All exact scalar-template and lift-classification gates passed.')
    print(json.dumps({'counts':COUNTS,'total_gates':sum(COUNTS.values()),
                      'first_obstruction_universes':first,'granularity_universes':grid,
                      'exhaustive_lift_universes':lifts,
                      'finite_logarithm_nontrivial_roots_below_200':roots},
                     sort_keys=True,indent=2))

if __name__=='__main__':
    main()
