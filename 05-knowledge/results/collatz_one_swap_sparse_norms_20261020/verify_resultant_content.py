#!/usr/bin/env python3
"""Exact symbolic certificates and independent positional one-swap tests.
No repository mathematical-code imports. Run from this directory.
"""
import argparse
import json
import math
from pathlib import Path
import sympy as sp

X, Z = sp.symbols('X Z')

def check(ok, message):
    if not ok:
        raise ArithmeticError(message)

def ceildiv(a, b):
    return -(-a // b)

def strip6(n):
    n = abs(int(n))
    check(n != 0, 'zero resultant')
    for ell in (2, 3):
        while n % ell == 0:
            n //= ell
    return n

def mechanical(q, p):
    return [((i + 1) * p) // q - (i * p) // q for i in range(q)]

def carry_positional(word):
    ones = [i for i, bit in enumerate(word) if bit]
    return sum(2**i * 3**(len(ones)-1-j) for j, i in enumerate(ones))

def carry_composition(word):
    # Affine numerator recurrence after appending the chronological next bit.
    B = 0
    for i, bit in enumerate(word):
        B = (3 if bit else 1) * B + bit * 2**i
    return B

def polynomial_coefficients(p, r, t, epsilon):
    eta = [ceildiv((t + m)*r, p) - ceildiv(t*r, p) for m in range(p)]
    h = max(0, -eta[1]-epsilon)
    coeff = [3**(p-1-m) * 2**(h+eta[m]+(epsilon if m == 1 else 0)) for m in range(p)]
    check(all(isinstance(c, int) for c in coeff), 'nonintegral polynomial')
    return coeff, h

def polynomial_value(coeff, x, modulus=None):
    b = 0
    for c in reversed(coeff):
        b = b*x+c
        if modulus:
            b %= modulus
    return b

def reduced_trinomial(p, a, epsilon):
    if epsilon == 1:
        return Z**(a+1)-Z**a+1
    return Z**(p-a)-Z+1

def norm_data(p, a, epsilon):
    Q = reduced_trinomial(p, a, epsilon)
    R = int(sp.resultant(Z**p-2, Q, Z))
    return R, strip6(R)

def gcd_period(coeff, p, r, M):
    if M == 1:
        return 1, [1]
    T = int(sp.n_order(2, M))
    values = []
    for k in range(T):
        x = pow(2, k, M)
        d = (pow(2, r, M)*pow(x, p, M)-pow(3, p, M)) % M
        values.append(math.gcd(M, math.gcd(d, polynomial_value(coeff, x, M))))
    minimal = next(d for d in sp.divisors(T) if all(values[i] == values[i % d] for i in range(T)))
    return int(minimal), values[:minimal]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--catalogue', type=Path, default=None)
    args = parser.parse_args()
    norms = {}
    universal = []
    for p in range(2, 13):
        for a in range(p):
            for epsilon in (-1, 1):
                R, M = norm_data(p, a, epsilon)
                norms[p, a, epsilon] = M
                check(M <= 3**p-2, f'norm bound {(p,a,epsilon)}')
                universal.append(dict(p=p, a=a, epsilon=epsilon, resultant=R, M=M))
    # Compare the entire polynomial construction with the universal norm, not
    # merely integer samples. This also checks residue-independent type labels.
    symbolic_checks = 0
    for p in range(2, 11):
        for r in range(1, p):
            if math.gcd(p, r) != 1:
                continue
            for t in range(p):
                a = (-(t+1)*r) % p
                for epsilon in (-1, 1):
                    coeff, h = polynomial_coefficients(p, r, t, epsilon)
                    P = sum(c*X**j for j, c in enumerate(coeff))
                    R = sp.resultant(2**r*X**p-3**p, P, X)
                    check(strip6(R) == norms[p, a, epsilon], f'resultant mismatch {(p,r,t,epsilon)}')
                    symbolic_checks += 1
    direct_checks = crossing_checks = k1_checks = balanced_returns = 0
    direct_max = None
    for p in range(2, 13):
        for q in range(p+1, 181):
            if math.gcd(p,q) != 1:
                continue
            k, r = divmod(q,p)
            c = mechanical(q,p)
            # Odd indices start at the last 1 of the zero-intercept word.
            s = [(j*k+ceildiv(j*r,p)-1) % q for j in range(p)]
            check(set(s) == {i for i,b in enumerate(c) if b}, 'odd indexing')
            expected_occurrences = sum(c[i] != c[(i+1)%q] for i in range(q))
            got = 0
            rotations = {tuple(c[v:]+c[:v]) for v in range(q)}
            for j in range(p):
                t = (j-1) % p
                a = (-j*r) % p
                for epsilon in (-1,1):
                    target = (s[j]+epsilon) % q
                    if c[target]:
                        continue
                    got += 1
                    w = c.copy()
                    w[s[j]], w[target] = 0,1
                    anchor = s[t]
                    rotated = w[anchor:]+w[:anchor]
                    coeff,h = polynomial_coefficients(p,r,t,epsilon)
                    P = polynomial_value(coeff, 2**k)
                    B = carry_positional(w)
                    Br = carry_positional(rotated)
                    D = 2**q-3**p
                    g = math.gcd(abs(D),B)
                    check(B == carry_composition(w), 'carry implementations')
                    check(P == 2**h*Br, f'polynomial value {(p,q,j,epsilon)}')
                    check(g == math.gcd(abs(D),P), 'rotation/scaling invariance')
                    M = norms[p,a,epsilon]
                    check(M % g == 0, f'norm divisibility {(p,q,j,epsilon,g,M)}')
                    check(g == math.gcd(M,math.gcd(abs(D),P)), 'exact capped gcd')
                    if a == 0 and epsilon == -1:
                        check(g == math.gcd(3**p-2, 2**k*3**(r-1)-1), 'minus endpoint formula')
                    if a == p-1 and epsilon == 1:
                        check(g == math.gcd(3**p-2**(p-1), 2**(k+r)-3**(r+1)), 'plus endpoint formula')
                    if D > 0 and g == D:
                        check(2**(q-1) < 3**p < 2**q, 'critical strip implication')
                    if tuple(w) in rotations:
                        balanced_returns += 1
                        check(g == 1, 'balanced return not free')
                    crossing_checks += int({s[j],target} == {0,q-1})
                    k1_checks += int(k == 1)
                    direct_checks += 1
                    if direct_max is None or g > direct_max['g']:
                        direct_max = dict(p=p,q=q,j=j,epsilon=epsilon,g=g)
            check(got == expected_occurrences, 'missed or duplicate directed swap')
    classes = []
    period_checks = 0
    for p in (3,4,5):
        for r in range(1,p):
            if math.gcd(p,r) != 1:
                continue
            for t in range(p):
                a = (-(t+1)*r) % p
                for epsilon in (-1,1):
                    coeff,h = polynomial_coefficients(p,r,t,epsilon)
                    M=norms[p,a,epsilon]
                    T,values=gcd_period(coeff,p,r,M)
                    hits={str(g):[k for k,v in enumerate(values) if v == g] for g in sorted(set(values)) if g>1}
                    classes.append(dict(p=p,r=r,t=t,a=a,epsilon=epsilon,M=M,minimal_period=T,nonunit_classes=hits))
                    for k in range(2,2+2*T):
                        g=math.gcd(abs(2**(p*k+r)-3**p),polynomial_value(coeff,2**k))
                        check(g == values[k%T], 'period replay')
                        period_checks += 1
    sharp = []
    for p in range(2,9):
        M=3**p-2
        k=int(sp.n_order(2,M))
        if k<2:
            k*=2
        r=1
        coeff,h=polynomial_coefficients(p,r,p-1,-1)
        g=math.gcd(2**(p*k+r)-3**p,polynomial_value(coeff,2**k))
        check(g == M, 'sharpness replay')
        sharp.append(dict(p=p,k=k,q=p*k+1,g=g))
    catalogue = dict(universal_norms=universal,complete_periodic_classes=classes)
    if args.catalogue:
        args.catalogue.write_text(json.dumps(catalogue,indent=2,sort_keys=True)+'\n')
    receipt = dict(status='PASS',symbolic_resultant_comparisons=symbolic_checks,
        symbolic_universe='2<=p<=10; every coprime 1<=r<p; all p roots and both directions',
        universal_norm_count=len(universal),universal_norm_universe='2<=p<=12; all a and both directions',
        direct_swap_checks=direct_checks,direct_universe='2<=p<=12; p<q<=180; gcd(p,q)=1; every unequal cyclic adjacent swap occurrence',
        crossing_cut_checks=crossing_checks,k_equals_1_checks=k1_checks,balanced_return_checks=balanced_returns,
        maximum_content_in_direct_universe=direct_max,period_classes=len(classes),period_replay_checks=period_checks,
        sharpness_examples=sharp)
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
