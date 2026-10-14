"""Exact rational-run-bank audit. No imports from other repository modules.

Writes nothing; all checks survive python -O. Scope is printed as JSON.
"""
from fractions import Fraction
import json

CHECKS = 0

def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)

def valuations(n, count):
    if type(n) is not int or n <= 0 or n % 2 != 1:
        raise ValueError('positive odd integer required')
    if type(count) is not int or count < 0:
        raise ValueError('nonnegative exact count required')
    word = []
    for _ in range(count):
        if n == 1:
            break
        z = 3*n+1
        a = 0
        while z % 2 == 0:
            z //= 2
            a += 1
        word.append(a)
        n = z
    return tuple(word), n

def block(k, m):
    P, Q = 3**(k+m), 2**(k+2*m)
    B = Q+P-2**(k+1)*3**m
    u = pow(3**k, -1, 4**m)
    n = 2**(k+1)*u-1
    return P, Q, B, n

def schedule_certificate(a, b):
    check(3**(a+b) < 2**(b+2*a), 'strict slope contraction')
    seeds = [(r, (a*r+b-1)//b) for r in range(1,b+1)]
    L = max(Fraction(3**(k+m), 2**(k+2*m)) for k,m in seeds)
    Z = (1+L)/(1-L)
    K = 1
    while 2**(K+1)-1 < Z:
        K += 1
    short = []
    for k in range(1,K):
        m = (a*k+b-1)//b
        P,Q,B,n = block(k,m)
        check((Q-P)*n > B, 'finite native threshold gate')
        short.append([k,m,n,str(Fraction(B,Q-P))])
    density = sum((Fraction(1,2**(k+2*m)) for k,m in seeds),Fraction())
    density /= 1-Fraction(1,2**(b+2*a))
    for k in range(1,129):
        m = (a*k+b-1)//b
        P,Q,B,n0 = block(k,m)
        for t in range(8):
            n = n0+2*Q*t
            w,h = valuations(n,k+m)
            check(w == (1,)*k+(2,)*m, 'literal exact valuation word')
            check(Q*h == P*n+B and 1<h<n, 'native payment and carrier')
        r = (k-1)%b+1
        j = (k-r)//b
        mr = (a*r+b-1)//b
        check(m == mr+j*a, 'schedule period')
        check(Fraction(P,Q) <= L, 'global multiplier bound')
    # Full native residue universes, independent parity replay.
    census = []
    for k in range(1,5):
        m=(a*k+b-1)//b
        P,Q,B,n0=block(k,m)
        found=[]
        for n in range(1,2*Q,2):
            w,h=valuations(n,k+m)
            if w == (1,)*k+(2,)*m:
                found.append(n)
        check(found == [n0], 'complete exact native cell')
        census.append([k,m,2*Q,n0])
    # Exact Mersenne rejection for each tested exponent; symbolic proof in report.
    for E in range(2,129):
        k=E-1
        m=(a*k+b-1)//b
        check((3**k-1) % (4**m) != 0, 'Mersenne complement')
    return {'a':a,'b':b,'density':str(density),'max_multiplier':str(L),
            'uniform_threshold':str(Z),'tail_k_from':K,
            'short_certificate':short,'complete_cell_census':census,
            'literal_family_k':[1,128],'members_per_k':8}

if __name__=='__main__':
    banks=[schedule_certificate(2,1),schedule_certificate(3,2),schedule_certificate(10,7)]
    check(Fraction(banks[1]['density'])==Fraction(3,85),'closed new bank density')
    check(valuations(967,5)==((1,1,2,2,2),919),'strict new source witness')
    check(valuations(967,6)[0] != (1,1,2,2,2,2),'old bank excludes witness')
    check(valuations(1,10)==((),1),'root receives no padding')
    print(json.dumps({'status':'PROVED scoped rational-run banks; universal coverage and child grounding OPEN',
        'rational_run_banks':banks,'strict_extension_witness':{'source':967,'endpoint':919,'word':[1,1,2,2,2]},
        'density_convention':'Relative to positive odd integers; divide by 2 for all positive integers',
        'checks':CHECKS},indent=2,sort_keys=True))
