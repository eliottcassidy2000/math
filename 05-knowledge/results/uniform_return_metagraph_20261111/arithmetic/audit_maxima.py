"""Independent arithmetic routes for every entry of the million-target prefix.

Even sources: additive proper-divisor sieve, not factorization.
Odd squares: trial factorization of roots, not SPF.
Odd extrema: prime-list predecessor search and squared rational bound.
"""
from array import array
from bisect import bisect_left
from hashlib import sha256
from math import isqrt
from pathlib import Path
import argparse
import json
import maxima


def gate(ok, message):
    if not ok:
        raise RuntimeError(message)


def independent_primes(bound):
    sieve=[True]*(bound+1)
    sieve[0]=sieve[1]=False
    for p in range(2,isqrt(bound)+1):
        if sieve[p]:
            for k in range(p*p,bound+1,p):
                sieve[k]=False
    return sieve,[p for p in range(2,bound+1) if sieve[p]]


def square_aliquot(root):
    a=root
    sigma=1
    p=3
    while p*p<=a:
        exponent=0
        while a%p==0:
            a//=p
            exponent+=1
        if exponent:
            sigma*=sum(p**j for j in range(2*exponent+1))
        p+=2
    if a>1:
        sigma*=1+a+a*a
    return sigma-root*root


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--limit',type=int,default=1_000_000)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    X=args.limit
    got,_,_=maxima.enumerate_maxima(X)
    prime,primes=independent_primes(X)
    expected=array('Q',[0])*(X+1)
    s=array('Q',[0])*(2*X+1)
    for d in range(1,X+1):
        for m in range(2*d,2*X+1,d):
            s[m]+=d
    for m in range(4,2*X+1,2):
        n=int(s[m])
        if 2<=n<=X and n%2==0:
            expected[n]=max(expected[n],m)
    for p in primes:
        if p+1<=X:
            expected[p+1]=p*p
    root=3
    odd_square_count=0
    while root**3<=X*X:
        if not prime[root]:
            n=square_aliquot(root)
            odd_square_count+=1
            if 2<=n<=X:
                gate(n%2==0,'square parity')
                expected[n]=max(expected[n],root*root)
        root+=2
    fallback=[]
    for n in range(5,X+1,2):
        best=0
        index=bisect_left(primes,(n-1)//2)-1
        while index>=0:
            p=primes[index]
            q=n-1-p
            if q<=X and prime[q]:
                best=p*q
                break
            index-=1
        # Exact m > r^3, independent of producer's integer-square-root cap.
        if best and (2*best+3*n-2)**2>n*n*(4*n-3):
            expected[n]=best
        else:
            fallback.append(n)
            gate((n-1)**2<=2*X,'audit fallback source universe exceeded')
            expected[n]=max((m for m in range(4,(n-1)**2+1) if s[m]==n),default=0)
    gate(got[2:]==expected[2:],'million-entry disagreement')
    receipt={
        'status':'PASS', 'targets':[2,X], 'entry_comparisons':X-1,
        'complete_additive_divisor_sieve_bound':2*X,
        'trial_factored_composite_odd_roots':odd_square_count,
        'global_bound_fallback_targets':fallback,
        'table_sha256':maxima.normalized_hash(expected),
        'scope':'Independent arithmetic implementations; shared written source-type proof, not independent discovery.'}
    encoded=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(encoded,encoding='utf-8')
    print(encoded)


if __name__=='__main__':
    main()
