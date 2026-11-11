"""Certified A135244 prefix by extremal source types, without Goldbach assumptions.

Run: python -B maxima.py --limit 1000000 --output maxima_receipt.json
The fallback is complete for every input; observed cheap fallback is finite evidence.
"""
from array import array
from collections import Counter
from hashlib import sha256
from math import isqrt
from pathlib import Path
import argparse
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def prime_sieve(bound):
    a = bytearray(b'\1') * (bound + 1)
    a[:2] = b'\0\0'
    for p in range(2, isqrt(bound) + 1):
        if a[p]:
            a[p*p:bound+1:p] = b'\0' * ((bound-p*p)//p + 1)
    return a


def spf_sieve(bound):
    a = array('I', range(bound + 1))
    for p in range(2, isqrt(bound) + 1):
        if a[p] == p:
            for j in range(p*p, bound+1, p):
                if a[j] == j:
                    a[j] = p
    return a


def factor_sigma(n, spf, power=1):
    result = 1
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        result *= (p ** (power*e+1)-1)//(p-1)
    return result


def integer_cuberoot(n):
    lo, hi = 0, 1 << ((n.bit_length()+2)//3)
    while lo < hi:
        mid = (lo+hi+1)//2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid-1
    return lo


def cubic_cap(n):
    """floor(r^3) where r^2+r+1=n: all sources with Omega>=3 lie below it."""
    return (isqrt(n*n*(4*n-3))-3*n+2)//2


def dense_divisor_sums(bound):
    """Independent additive formulation, used only for the small complete audit."""
    s = array('Q', [0]) * (bound+1)
    for d in range(1, bound//2+1):
        for m in range(2*d, bound+1, d):
            s[m] += d
    return s


def enumerate_maxima(limit):
    require(limit >= 3, 'limit >=3 required')
    prime = prime_sieve(2*limit)
    spf = spf_sieve(2*limit)
    answer = array('Q', [0]) * (limit+1)
    kind = bytearray(limit+1)
    # 1 prime square, 2 closest distinct semiprimes, 3 even source,
    # 4 composite odd square, 5 exhaustive residual, 0 empty.
    for n in range(2, limit+1):
        if prime[n-1]:
            answer[n] = (n-1)**2
            kind[n] = 1
    even_sources = 0
    for m in range(6, 2*limit-5, 2):
        n = factor_sigma(m, spf)-m
        if 2 <= n <= limit and n % 2 == 0:
            even_sources += 1
            if m > answer[n]:
                answer[n], kind[n] = m, 3
    root_bound = integer_cuberoot(limit*limit)
    odd_composite_roots = 0
    square_hits = []
    for a in range(9, root_bound+1, 2):
        if prime[a]:
            continue
        odd_composite_roots += 1
        n = factor_sigma(a, spf, 2)-a*a
        if 2 <= n <= limit:
            require(n % 2 == 0, 'odd-square parity')
            if a*a > answer[n]:
                answer[n], kind[n] = a*a, 4
                square_hits.append((n,a*a))
    fallback = []
    odd_pair_scan_steps = 0
    for n in range(5, limit+1, 2):
        best = 0
        for p in range((n-2)//2, 1, -1):
            odd_pair_scan_steps += 1
            q = n-1-p
            if prime[p] and prime[q]:
                best = p*q
                break
        if best > cubic_cap(n):
            answer[n], kind[n] = best, 2
        else:
            fallback.append((n,best,cubic_cap(n)))
            answer[n], kind[n] = best, (2 if best else 0)
    # Only these residual targets need higher-factor sources. Every source of
    # Omega>=3 is at most cubic_cap(n); Omega2 has already been exhausted.
    if fallback:
        b = max(t[2] for t in fallback)
        fallback_spf = spf if b <= 2*limit else spf_sieve(b)
        wanted = {n for n,_,_ in fallback}
        for m in range(4,b+1):
            n = factor_sigma(m,fallback_spf)-m
            if n in wanted and m > answer[n]:
                answer[n], kind[n] = m, 5
    return answer, kind, {
        'limit':limit, 'even_source_bound':2*limit-6,
        'odd_composite_square_root_bound':root_bound,
        'even_source_hits':even_sources,
        'odd_composite_roots_tested':odd_composite_roots,
        'odd_pair_scan_steps':odd_pair_scan_steps,
        'fallback_targets':fallback,
        'types':dict(sorted(Counter(int(x) for x in kind[2:]).items())),
        'type_4_first20':[(n,int(answer[n])) for n in range(2,limit+1) if kind[n]==4][:20],
    }


def audit_prefix(answer, limit=500):
    bound=(limit-1)**2
    s=dense_divisor_sums(bound)
    truth=[0]*(limit+1)
    omega=bytearray(bound+1)
    spf=spf_sieve(bound)
    layer_checks=0
    equality=[]
    for m in range(2,bound+1):
        omega[m]=omega[m//spf[m]]+1
        n=int(s[m])
        if 2<=n<=limit:
            truth[n]=m
        if omega[m]>=3:
            require(m<=cubic_cap(n), 'Omega3 cap failure')
            layer_checks+=1
        if omega[m]>=2:
            # Compare all k<=Omega without floating roots, using the weakest
            # suffix-product consequence m^(k-1)<=(s(m)-1)^k.
            for k in range(2,omega[m]+1):
                require(m**(k-1)<=(n-1)**k,'factor layer bound')
        if omega[m]==3 and cubic_cap(n)==m:
            equality.append((m,n))
    require(list(answer[2:limit+1])==truth[2:], 'full source-sieve audit')
    if limit>=15:
        require(answer[5]==0 and answer[7]==8 and answer[15]==33, 'small residual controls')
    if bound>=9:
        require(s[9]==4 and 1+3+3!=s[9], 'repeated-prime hostile')
    if limit>=12:
        require(answer[12]==121 and 121>2*12, 'even-target linear-bound hostile')
    return {'targets':[2,limit], 'complete_source_bound':bound,
            'Omega3_sources_checked':layer_checks,'integer_cap_equalities_first20':equality[:20]}


def normalized_hash(answer):
    h=sha256()
    for n in range(2,len(answer)):
        h.update(f'{n} {answer[n]}\n'.encode())
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--limit',type=int,default=1_000_000)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--table',type=Path)
    args=ap.parse_args()
    values,kinds,receipt=enumerate_maxima(args.limit)
    receipt['independent_audit']=audit_prefix(values,min(500,args.limit))
    receipt['status']='FINITE-EXACT certified prefix; no global Goldbach or maximum-shape claim'
    receipt['table_sha256']=normalized_hash(values)
    receipt['first59']=[int(v) for v in values[2:61]]
    receipt['last10']=[(n,int(values[n])) for n in range(max(2,args.limit-9),args.limit+1)]
    receipt['type_key']={'0':'empty','1':'prime square','2':'closest distinct semiprimes',
                         '3':'even source','4':'composite odd square','5':'residual exhaustive'}
    if args.table:
        with args.table.open('w',encoding='utf-8',newline='\n') as f:
            for n in range(2,len(values)):
                f.write(f'{n} {values[n]}\n')
    encoded=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(encoded,encoding='utf-8')
    print(encoded)


if __name__=='__main__':
    main()
