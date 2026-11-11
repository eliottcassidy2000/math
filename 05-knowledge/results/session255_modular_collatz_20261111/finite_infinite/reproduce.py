#!/usr/bin/env python3
"""Exact new-claim controls. Python standard library; SymPy for determinants."""
import itertools, json, math
from fractions import Fraction
from pathlib import Path
import sympy as S


def affine(word, sigma=1):
    a, B = 1, 0
    for i, b in enumerate(word):
        m = 3 if b else 1
        a, B = m*a, m*B + sigma*b*2**i
    return a, B


def step(n, sigma=1):
    return (3*n+sigma)//2 if n%2 else n//2


def residue_step(r, b, q, sigma=1):
    return ((3*r+sigma) if b else r)*pow(2,-1,q)%q


def matrix(q, u, sigma=1):
    A=S.zeros(q)
    for r in range(q):
        for b in (0,1):
            A[r,residue_step(r,b,q,sigma)] += u if b else 1
    return A


def run():
    checks={}
    # Unique dyadic cylinders, direct signed integer execution, affine receipts.
    count=0
    for sigma in (-1,1):
        for k in range(1,10):
            seen={}
            for n in range(2**k):
                x=n; w=[]
                for _ in range(k):
                    w.append(x%2); x=step(x,sigma)
                a,B=affine(w,sigma)
                assert 2**k*x == a*n+B
                assert tuple(w) not in seen
                seen[tuple(w)]=n
                assert n == (-B*pow(a,-1,2**k))%2**k
                count+=1
            assert len(seen)==2**k
    checks['signed_dyadic_cylinders']=count
    # Direct finite-state path enumeration independent of gcd fixed-point formula.
    cases=0; snapshots={}
    for sigma in (-1,1):
        for depth in range(1,5):
            q=3**depth
            for k in range(1,10):
                coeff=[0]*(k+1)
                for w in itertools.product((0,1),repeat=k):
                    s=sum(w); fixed=0
                    for r in range(q):
                        y=r
                        for b in w: y=residue_step(y,b,q,sigma)
                        fixed += (y==r)
                    a,B=affine(w,sigma); g=math.gcd(2**k-a,q)
                    assert fixed==(g if B%g==0 else 0)
                    coeff[s]+=fixed
                expected=[math.gcd(2**k-1,q)]+[math.comb(k,s) for s in range(1,k+1)]
                assert coeff==expected
                cases+=1
                if sigma==1 and k in (2,6,9): snapshots[f'q={q},k={k}']=coeff
    checks['signed_prime_power_trace_polynomials']=cases
    z,u=S.symbols('z u'); determinants={}
    for sigma in (-1,1):
        for depth in (1,2):
            q=3**depth
            target=(1-(1+u)*z)*S.prod(1-z**(2*3**j) for j in range(depth))
            got=(S.eye(q)-z*matrix(q,u,sigma)).det(method='domain-ge')
            assert S.expand(got-target)==0
            determinants[f'q={q},sigma={sigma}']=str(S.factor(got))
    checks['symbolic_determinants']=determinants
    # Mixed-prime extension: ternary inflation changes only the even-word trace.
    mixed=0
    for sigma in (-1,1):
        for m,depth in ((5,1),(5,2),(7,1)):
            q=3**depth*m
            for k in range(1,8):
                difference=[0]*(k+1)
                for w in itertools.product((0,1),repeat=k):
                    counts=[]
                    for modulus in (q,m):
                        fixed=0
                        for r in range(modulus):
                            y=r
                            for b in w:y=residue_step(y,b,modulus,sigma)
                            fixed+=(y==r)
                        counts.append(fixed)
                    difference[sum(w)]+=counts[0]-counts[1]
                assert difference==[math.gcd(2**k-1,m)*(math.gcd(2**k-1,3**depth)-1)]+[0]*k
                mixed+=1
    checks['mixed_prime_trace_polynomials']=mixed
    for sigma in (-1,1):
        for value in (0,1,2):
            d15=(S.eye(15)-z*matrix(15,value,sigma)).det(method='domain-ge')
            d5=(S.eye(5)-z*matrix(5,value,sigma)).det(method='domain-ge')
            assert S.expand(d15-d5*(1-z**2)*(1-z**4)**2)==0
    checks['mixed_prime_determinant_specializations']=6
    # Arbitrary-driver extension: the no-consecutive-ones two-state DFA.
    D0=S.Matrix([[1,0],[1,0]]); D1=S.Matrix([[0,1],[0,0]])
    driver_controls={}
    for sigma in (-1,1):
        E=matrix(3,0,sigma); O=matrix(3,1,sigma)-E
        M=S.kronecker_product(E,D0)+u*S.kronecker_product(O,D1)
        expected=(1-z-u*z**2)*(1-z**2)
        got=(S.eye(6)-z*M).det(method='domain-ge')
        assert S.expand(got-expected)==0
        assert S.expand((S.eye(2)-z*(D0+u*D1)).det()-(1-z-u*z**2))==0
        beta=S.Matrix([1,1,0,0,0,0]) # accept either driver state at residue zero
        resolvent=(S.eye(6)-z*M).inv()*beta
        series={}
        for residue in (0,1):
            value=S.cancel(resolvent[2*residue])
            series[str(residue)]=str(value)
            # Independent direct recurrence for length-zero through length-six.
            row=S.zeros(1,6); row[2*residue]=1
            coefficients=[]
            for k in range(7):
                coefficients.append(S.expand((row*beta)[0])); row=row*M
            expansion=S.series(value,z,0,7).removeO().expand()
            assert all(S.expand(expansion.coeff(z,k)-coefficients[k])==0 for k in range(7))
        assert series['0']!=series['1']
        driver_controls[str(sigma)]={'determinant':str(S.factor(got)),'accepted_series_by_start_residue':series}
    checks['finite_driver_no_consecutive_ones']=driver_controls


    # CRT synchronized branch versus independent choices: false mixed transitions.
    def crt(a,b): return next(r for r in range(15) if r%3==a and r%5==b)
    sync={residue_step(0,b,15) for b in (0,1)}
    independent={crt(residue_step(0,b,3),residue_step(0,c,5)) for b,c in itertools.product((0,1),repeat=2)}
    assert sync=={0,8} and independent-sync=={3,5}
    checks['crt_hostile']={'synchronized':sorted(sync),'independent':sorted(independent),'spurious':sorted(independent-sync)}
    # Multiplication modulo 8 preserves low bits; division loses one bit.
    checks['times_three_mod8']=[3*r%8 for r in range(8)]
    assert 1%8==9%8 and step(1)%8!=step(9)%8
    checks['mod8_division_hostile']={'inputs':[1,9],'outputs':[step(1),step(9)]}
    # A periodic rational point survives every ternary quotient but is not integral.
    w=(1,0,0); a,B=affine(w); x=Fraction(B,2**len(w)-a)
    assert x==Fraction(1,5)
    for depth in range(1,9):
        q=3**depth; r=pow(5,-1,q); y=r
        for b in w:y=residue_step(y,b,q)
        assert y==r
    assert all((8-3)*r%5 != B%5 for r in range(5))
    checks['integer_realization_hostile']={'word':'100','rational_root':str(x),'all_ternary_levels':'accepted by theorem','mod5':'rejected'}
    checks['trace_snapshots']=snapshots
    out=Path(__file__).with_name('EXACT_RESULTS.json')
    out.write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':run()
