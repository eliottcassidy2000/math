#!/usr/bin/env python3
"""Exact inverse aliquot enumeration by increasing prime factors; no source scan."""
import json, math
from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
from sympy import primerange, isprime, factorint, divisor_sigma

def need(b, msg):
    if not b: raise RuntimeError(msg)

def fiber(t):
    if t < 2: raise ValueError('finite-target domain t>=2')
    primes=list(primerange(2,t+1))
    out=set(); calls=0
    def visit(m, sig, start):
        nonlocal calls
        calls+=1
        sm=sig-m
        # Last prime to exponent one: linear equation, with strict ordering.
        if sm:
            num=t-sig
            if num>0 and num%sm==0:
                q=num//sm
                if (start==0 or q>primes[start-1]) and isprime(q): out.add(m*q)
        for i in range(start,len(primes)):
            p=primes[i]
            # No p or larger can be a last factor of exponent >=2, or
            # a non-last factor, once both monotone lower bounds exceed t.
            terminal2=sm*p*p+sig*(p+1)
            continuation1=sm*p*p+2*sig*p+sig
            if min(terminal2,continuation1)>t: break
            pp=p; ps=1+p; exponent=1
            while True:
                n=m*pp; sn=sig*ps; value=sn-n
                if value>t: break
                if exponent>=2 and value==t: out.add(n)
                if i+1<len(primes):
                    q=primes[i+1]
                    if value*q+sn<=t: visit(n,sn,i+1)
                pp*=p; ps+=pp; exponent+=1
    visit(1,1,0)
    return sorted(out),calls

def sieve(N):
    s=[0]*(N+1)
    for d in range(1,N//2+1):
        for n in range(2*d,N+1,d): s[n]+=d
    return s

def main():
    s=sieve(299**2)
    for t in range(2,301):
        got,_=fiber(t)
        expected=[n for n in range(1,(t-1)**2+1) if s[n]==t]
        need(got==expected, f'fiber mismatch {t}: {got}, {expected}')
    need(fiber(2)[0]==[] and fiber(5)[0]==[], 'untouchable controls')
    need(fiber(12)[0]==[121], 'prime-square exception')
    need(fiber(13)[0]==[27,35], 'nonsemiprime hostile')
    need(int(divisor_sigma(4))-4==3, 'repeated-prime hostile')
    for n in range(1,10001):
        v=int(divisor_sigma(n))
        square_or_twice=(math.isqrt(n)**2==n or (n%2==0 and math.isqrt(n//2)**2==n//2))
        need((v%2==1)==square_or_twice,'square parity')
    for a in range(1,31):
        need(int(divisor_sigma(2**a))-2**a==2**a-1,'Mersenne family')
    results={}
    for t in [1001,1000001,1001000]:
        f,c=fiber(t)
        need(all(int(divisor_sigma(n))-n==t for n in f),'reencoding')
        pairs=[(int(p),t-1-int(p)) for p in primerange(2,(t-1)//2+1) if p<t-1-p and isprime(t-1-p)]
        semis=sorted(p*q for p,q in pairs)
        need(set(semis)<=set(f),'Goldbach missing')
        for p,q in pairs:
            need((t-1)**2-4*p*q==(q-p)**2,'gap-square identity')
            need(int(divisor_sigma(p*q))-p*q==t,'distinct semiprime identity')
        results[str(t)]={'factorization':str(factorint(t)), 's_target':int(divisor_sigma(t))-t,'fiber':f,'count':len(f),'semiprime_count':len(semis),'non_semiprime_sources':sorted(set(f)-set(semis)),'calls':c}
        print(t,len(f),len(semis),c,flush=True)
    with (DATA_DIR / 'results.json').open('w') as h: json.dump(results,h,indent=2,sort_keys=True)
    print('PASS: complete fibers agree with independent divisor sieve for targets 2..300')
if __name__=='__main__': main()
