#!/usr/bin/env python3
"""Standalone exact audit. No repository imports; no external side effects."""
from fractions import Fraction
from itertools import product
from collections import Counter
import json
import sys

COUNTS = Counter()
def check(ok, label):
    COUNTS[label] += 1
    if not ok:
        raise RuntimeError('FAILED: ' + label)

def require_word(s, nonempty=False):
    if not isinstance(s, str) or (nonempty and not s) or any(c not in '01' for c in s):
        raise ValueError('expected a binary word' + (' of positive length' if nonempty else ''))

def tm(n):
    return ''.join(str(i.bit_count() % 2) for i in range(n))

def subst(s):
    return ''.join('01' if c == '0' else '10' for c in s)

def morph(s, U, V):
    return ''.join(U if c == '0' else V for c in s)

def iterate(A, n):
    w = ''
    for _ in range(n):
        e = A % 2
        w += str(e)
        A = (3*A + e)//2
    return w, A

def carry(w):
    C = 0
    power = 1
    for e in w:
        C = 3*C + int(e)*power
        power *= 2
    return C

def v2(n):
    if n <= 0:
        raise ValueError('positive valuation argument required')
    return (n & -n).bit_length()-1

def certificate(A, prefix, U, V):
    if type(A) is not int or A <= 0:
        raise ValueError('positive ordinary integer source required')
    require_word(prefix)
    require_word(U, True)
    require_word(V, True)
    actual_prefix, B = iterate(A, len(prefix))
    for i, (e, f) in enumerate(zip(actual_prefix, prefix)):
        if e != f:
            return {'first_mismatch': i, 'horizon': len(prefix), 'kind': 'prefix'}
    a = len(prefix)
    L = len(U)+len(V)
    n = 1
    k = L
    while 2**(a+5*k) < 3**(a+3*k)*(A+1):
        n += 1
        k *= 2
    seed = '011'
    for _ in range(n):
        seed = subst(seed)
    period = morph(seed, U, V)
    p = len(period)
    C = carry(period)
    D = 3**p-2**p
    h = v2(D*B+C)
    if p != 3*k or h >= 5*k:
        raise RuntimeError('the proved finite bound failed')
    return {'first_mismatch': a+h, 'horizon': a+5*k, 'kind': 'periodic-shadow',
            'substitution_depth': n, 'k': k, 'period_length': p,
            'period': period, 'tail_match': h, 'shadow_denominator_unreduced': D,
            'valuation_argument': D*B+C}

def target_prefix(prefix, U, V, N):
    if N <= len(prefix):
        return prefix[:N]
    count = (N-len(prefix)+min(len(U),len(V))-1)//min(len(U),len(V))
    return (prefix+morph(tm(count), U,V))[:N]

def verify_certificate(A, prefix, U, V):
    out = certificate(A,prefix,U,V)
    H = out['horizon']
    actual,_=iterate(A,H)
    target=target_prefix(prefix,U,V,H)
    first=next((i for i,(x,y) in enumerate(zip(actual,target)) if x!=y),None)
    check(first==out['first_mismatch'] and first<H, 'source-bound exact first mismatch')
    return out

def greedy(n):
    z=Fraction(1)
    w=''
    for _ in range(n):
        z*=Fraction(3,2)
        e=z.numerator//z.denominator
        z-=e
        w+=str(e)
    return w

def follower(w,d):
    q=0
    resets=0
    for c in w:
        if c==d[q]:
            q+=1
        elif c=='0' and d[q]=='1':
            q=0
            resets+=1
        else:
            return False,q,resets
    return True,q,resets

def safe_slack(w):
    m=Fraction(1)
    for c in w:
        m=3*m/2-int(c)
        if m<=0:
            return False
        m=min(m,1)
    return True

def polynomial(w,z):
    s=Fraction(0)
    for c in reversed(w):
        s=s*z+int(c)
    return s

def main():
    check(tm(8192)==subst(tm(4096)), 'independent Thue-Morse definition')
    for n in range(10):
        a,b='01101','011'
        for _ in range(n): a,b=subst(a),subst(b)
        check(a==(b*2)[:len(a)], 'substituted overlap')
        check(a==tm(len(a)), 'actual fixed-point overlap')
    words=[''.join(x) for l in range(1,5) for x in product('01',repeat=l)]
    z=Fraction(2,3)
    for U in words:
        for V in words:
            check(morph(tm(128),U,V)==morph(tm(64),U+V,V+U), 'equal-block recoding')
            delta=polynomial(V+U,z)-polynomial(U+V,z)
            check(bool(delta)==(U+V!=V+U), 'noncommuting specialization guard')
            for n in range(1,5):
                seed='01101'; pseed='011'
                for _ in range(n):seed,pseed=subst(seed),subst(pseed)
                t,p=morph(seed,U,V),morph(pseed,U,V)
                k=2**(n-1)*(len(U)+len(V))
                check(len(t)==5*k and len(p)==3*k and t==(p*2)[:len(t)], 'unequal-length clock')
    pairs=[('0','1'),('000','100'),('0','100'),('100','101000000'),
           ('01','010'),('00','000'),('1','11'),('1010','0001')]
    for A in range(1,1025):
        for U,V in pairs:
            verify_certificate(A,'',U,V)
    for A in range(1,257):
        for prefix in ('0','1','100','00100','10101'):
            for U,V in pairs[:4]:verify_certificate(A,prefix,U,V)
    d=greedy(2000)
    loops=[d[:l-1]+'0' for l in range(1,18) if d[l-1]=='1']
    for U in loops:
        for V in loops:
            w=morph(tm(64),U,V)
            ok,q,resets=follower(w,d)
            check(ok and q==0 and resets>=64 and safe_slack(w), 'safe reset-bank control')
    runway=[]
    U,V='000','100'
    for N in (8,16,32,64,128,256):
        w=target_prefix('',U,V,N)
        C=carry(w)
        A=(-C*pow(3**N,-1,2**N))%(2**N)
        check(A>0 and iterate(A,N)[0]==w and safe_slack(w), 'moving positive-source cylinder')
        out=verify_certificate(A,'',U,V)
        check(out['first_mismatch']>=N, 'moving-source boundary preserved')
        runway.append({'N':N,'A':str(A),'bit_length':A.bit_length(),
                       'first_mismatch':out['first_mismatch'],'horizon':out['horizon']})
    for bad in (True,False,0,-1,1.0,'1'):
        try: certificate(bad,'','0','1')
        except ValueError: check(True,'invalid source rejected')
        else: check(False,'invalid source rejected')
    for args in (('', '1'),('0',''),('2','1'),('0',None)):
        try: certificate(1,'',*args)
        except ValueError: check(True,'invalid morphism rejected')
        else: check(False,'invalid morphism rejected')
    examples=[]
    for A in (1,13,124173,3532045,10**30+1):
        out=verify_certificate(A,'','000','100')
        examples.append({k:out[k] for k in ('first_mismatch','horizon','k','period_length')})
        examples[-1]['A']=str(A)
    payload={'status':'PASS','checks':sum(COUNTS.values()),'check_groups':dict(sorted(COUNTS.items())),
             'source_universe':'A=1..1024, eight morphisms; A=1..256, five prefixes, four morphisms',
             'morphism_universe':'all 900 ordered nonempty binary image pairs with each length <=4',
             'safe_reset_loop_lengths':list(map(len,loops)), 'examples':examples, 'moving_source_controls':runway,
             'scope':'Finite exact audits only. All-parameter source exclusion is proved in THUE_MORSE_FINITE_SOURCE_EXCLUSION.md. These computations do not certify transcendence.'}
    print(json.dumps(payload,indent=2,sort_keys=True))

if __name__=='__main__':main()
